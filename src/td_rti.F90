!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.

#include "global.h"

module td_rti
  use mix
  use td_exp

  implicit none

  type td_rti_type
    integer           :: method         ! which evolution method to use
    type(td_exp_type) :: te             ! how to apply the propagator (e^{-i H \Delta t})


    real(r8), pointer :: v_old(:, :, :) ! storage of the KS potential of previous iterations
  end type td_rti_type

  integer, parameter :: OLD_REVERSAL         = 1, &
                        REVERSAL             = 2, &
                        APP_REVERSAL         = 3, &
                        EXPONENTIAL_MIDPOINT = 4  
  
contains
  subroutine td_rti_init(m, st, tr)
    type(mesh_type),   intent(in)    :: m
    type(states_type), intent(in)    :: st
    type(td_rti_type), intent(inout) :: tr

    call oct_parse_int("TDEvolutionMethod", REVERSAL, tr%method)
    select case(tr%method)
    case(OLD_REVERSAL);         message(1) = 'Info: Evolution method:  Old-Style.'
    case(REVERSAL);             message(1) = 'Info: Evolution method:  Enforced Time-Reversal Symmetry'
    case(APP_REVERSAL);         message(1) = 'Info: Evolution method:  Approx.Enforced Time-Reversal Symmetry' 
    case(EXPONENTIAL_MIDPOINT); message(1) = 'Info: Evolution method:  Exponential Midpoint Rule.'
    case default
      write(message(1), '(a,i6,a)') "Input: '", tr%method, "' is not a valid TDEvolutionMethod"
      message(2) = '(1 <= TDEvolutionMethod <= 4)'
      call write_fatal(2)
    end select
    call write_info(1)
    
    allocate(tr%v_old(m%np, st%nspin, 3)) ! allocate memory to store the old KS potentials
    call td_exp_init(m, tr%te)            ! initialize propagator

  end subroutine td_rti_init

  subroutine td_rti_end(tr)
    type(td_rti_type), intent(inout) :: tr

    ASSERT(associated(tr%v_old)) ! sanity check

    deallocate(tr%v_old)         ! clean ols KS potentials
    nullify(tr%v_old)

    call td_exp_end(tr%te)       ! clean propagator method
  end subroutine td_rti_end

  subroutine td_rti_run_zero_iter(h, tr)
    type(hamiltonian_type), intent(in) :: h
    type(td_rti_type), intent(inout) :: tr
    tr%v_old(:, :, 2) = h%vhxc(:, :)
    tr%v_old(:, :, 3) = tr%v_old(:, :, 2)
  end subroutine td_rti_run_zero_iter

  subroutine td_rti_dt(h, m, st, sys, tr, t, dt)
    type(hamiltonian_type), intent(inout) :: h
    type(mesh_type), intent(in) :: m
    type(states_type), intent(inout) :: st
    type(system_type), intent(in) :: sys
    type(td_rti_type), intent(inout) :: tr
    real(r8), intent(in) :: t, dt
    
    call push_sub('td_rti')
    
    tr%v_old(:, :, 3) = tr%v_old(:, :, 2)
    tr%v_old(:, :, 2) = tr%v_old(:, :, 1)
    tr%v_old(:, :, 1) = h%vhxc(:, :)

    select case(tr%method)
    case(OLD_REVERSAL)
      call td_rti1
    case(REVERSAL)
      call td_rti2
    case(APP_REVERSAL)
      if(t<3*dt) then
        call td_rti2
      else
        call td_rti3
      endif
    case(EXPONENTIAL_MIDPOINT)
      if(t<3*dt) then
        call td_rti2
      else
        call td_rti4
      endif
    end select

    call pop_sub()
  contains

    ! Warning: this subroutine should only be used with LDA/GGA functionals
    subroutine td_rti1
      integer is, ik, ist
      real(r8), allocatable :: aux(:,:)
      complex(r8), allocatable :: zpsi1(:,:,:,:)
    
      call push_sub('td_rti1')
      
      allocate(aux(m%np, st%nspin))
      
      call xpolate_pot(dt/2._r8, m%np, st%nspin, &
           tr%v_old(:, :, 3), tr%v_old(:, :, 2), tr%v_old(:, :, 1), aux)
      
      h%vhxc = aux
      allocate(zpsi1(m%np, st%dim, st%st_start:st%st_end, st%nik))
      zpsi1 = st%zpsi
      do ik = 1, st%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, sys, h, st%zpsi(:,:, ist, ik), ik, dt, t-dt)
        end do
      end do
      st%zpsi = zpsi1
      deallocate(zpsi1)
      
      call zcalcdens(st, m%np, aux, .true.)
      st%rho = (st%rho +  aux) / 2.0_r8
      deallocate(aux)
      
      call zhamiltonian_setup(h, m, st, sys)
      
      do ik = 1, st%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, sys, h, st%zpsi(:,:, ist, ik), ik, dt, t-dt)
        end do
      end do
      
      call pop_sub()
    end subroutine td_rti1
    
    subroutine td_rti2
      real(r8), allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
      complex(r8), allocatable :: zpsi1(:,:,:,:)
      integer is, ik, ist
      
      call push_sub('td_rti2')
      
      if(.not.h%ip_app) then
        allocate(zpsi1(m%np, st%dim, st%st_start:st%st_end, st%nik))
        zpsi1 = st%zpsi ! store zpsi
        
        allocate(vhxc_t1(m%np, st%nspin), vhxc_t2(m%np, st%nspin))
        vhxc_t1 = h%vhxc
        
        ! propagate dt with H(t-dt)
        do ik = 1, st%nik
          do ist = st%st_start, st%st_end
            call td_exp_dt(tr%te, sys, h, st%zpsi(:,:, ist, ik), ik, dt, t-dt)
          end do
        end do
        
        call zcalcdens(st, m%np, st%rho, .true.)
        call zhamiltonian_setup(h, m, st, sys)
        
        st%zpsi = zpsi1
        deallocate(zpsi1)
        
        vhxc_t2 = h%vhxc
        h%vhxc = vhxc_t1
      end if
      
      ! propagate dt/2 with H(t-dt)
      do ik = 1, st%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, sys, h, st%zpsi(:,:, ist, ik), ik, dt/2._r8, t-dt)
        end do
      end do
      
      ! propagate dt/2 with H(t)
      h%vhxc = vhxc_t2
      do ik = 1, st%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, sys, h, st%zpsi(:,:, ist, ik), ik, dt/2._r8, t)
        end do
      end do
      
      if(.not.h%ip_app) deallocate(vhxc_t1, vhxc_t2)
      
      call pop_sub()
    end subroutine td_rti2
    
    subroutine td_rti3
      integer is, ik, ist
      real(r8), allocatable :: aux(:,:)

      call push_sub('td_rti3')
      
      if(.not.h%ip_app) then
        allocate(aux(m%np, st%nspin))
        call xpolate_pot(dt, m%np, st%nspin, &
             tr%v_old(:, :, 3), tr%v_old(:, :, 2), tr%v_old(:, :, 1), aux)
      end if
      
      do ik = 1, st%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, sys, h, st%zpsi(:,:, ist, ik), ik, dt/2._r8, t-dt)
        end do
      end do
      
      if(.not.h%ip_app) h%vhxc = aux
      do ik = 1, st%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, sys, h, st%zpsi(:,:, ist, ik), ik, dt/2._r8, t)
        end do
      end do
      
      if(.not.h%ip_app) deallocate(aux)
      
      call pop_sub()
    end subroutine td_rti3
    
    subroutine td_rti4
      integer :: ist, ik
      
      call push_sub('td_rti4')
      
      if(.not.h%ip_app) then
        call xpolate_pot(dt/2._r8, m%np, st%nspin, &
             tr%v_old(:, :, 3), tr%v_old(:, :, 2), tr%v_old(:, :, 1), h%vhxc)
      end if
      
      do ik = 1, st%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, sys, h, st%zpsi(:,:, ist, ik), ik, dt, t - dt/2._r8)
        end do
      end do
      
      call pop_sub()
    end subroutine td_rti4
    
    subroutine xpolate_pot(t, np, dim, pot2, pot1, pot0, pot)
      real(r8), intent(in)  :: t
      integer, intent(in)   :: np, dim
      real(r8), intent(in)  :: pot0(np, dim), pot1(np, dim), pot2(np, dim)
      real(r8), intent(out) :: pot(np, dim)
    
      !pot = pot0 + t/dt       * (3._r8/2._r8*pot0 + 1._r8/2._r8*pot2 - 2.0_r8*pot1) + &
      !             t**2/dt**2 * (1._r8/2._r8*pot0 + 1._r8/2._r8*pot2 -        pot1)
      call dcopy(np*dim,                                                    pot0(1, 1), 1, pot(1, 1), 1)
      call daxpy(np*dim, (t/dt)*(3._r8/2._r8) + (t**2/dt**2)*(1._r8/2._r8), pot0(1, 1), 1, pot(1, 1), 1)
      call daxpy(np*dim, (t/dt)*(-2._r8)      + (t**2/dt**2)*(-1._r8),      pot1(1, 1), 1, pot(1, 1), 1)
      call daxpy(np*dim, (t/dt)*(1._r8/2._r8) + (t**2/dt**2)*(1._r8/2._r8), pot2(1, 1), 1, pot(1, 1), 1)
      
    end subroutine xpolate_pot
    
  end subroutine td_rti_dt

end module td_rti

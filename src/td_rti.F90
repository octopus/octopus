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
    real(r8), pointer :: vmagnus(:, :, :) ! auxiliary function to store the Magnus potentials.
  end type td_rti_type

  integer, parameter :: REVERSAL             = 2, &
                        APP_REVERSAL         = 3, &
                        EXPONENTIAL_MIDPOINT = 4, &
                        MAGNUS               = 5 

  real(r8), parameter, private :: scf_threshold = 1.0e-3_r8
  
contains
  subroutine td_rti_init(m, st, tr)
    type(mesh_type),   intent(in)    :: m
    type(states_type), intent(in)    :: st
    type(td_rti_type), intent(inout) :: tr

    call oct_parse_int("TDEvolutionMethod", REVERSAL, tr%method)
    select case(tr%method)
    case(REVERSAL);             message(1) = 'Info: Evolution method:  Enforced Time-Reversal Symmetry'
    case(APP_REVERSAL);         message(1) = 'Info: Evolution method:  Approx.Enforced Time-Reversal Symmetry' 
    case(EXPONENTIAL_MIDPOINT); message(1) = 'Info: Evolution method:  Exponential Midpoint Rule.'
    case(MAGNUS);               message(1) = 'Info: Evolution method:  Magnus expansion.'
       allocate(tr%vmagnus(m%np, st%nspin, 2))
    case default
      write(message(1), '(a,i6,a)') "Input: '", tr%method, "' is not a valid TDEvolutionMethod"
      message(2) = '(2 <= TDEvolutionMethod <= 4)'
      call write_fatal(2)
    end select
    call write_info(1)
    
    allocate(tr%v_old(m%np, st%nspin, 0:3)) ! allocate memory to store the old KS potentials
    call td_exp_init(m, tr%te)            ! initialize propagator

  end subroutine td_rti_init

  subroutine td_rti_end(tr)
    type(td_rti_type), intent(inout) :: tr

    ASSERT(associated(tr%v_old)) ! sanity check
    deallocate(tr%v_old)         ! clean ols KS potentials
    nullify(tr%v_old)

    if(tr%method == MAGNUS) then
      ASSERT(associated(tr%vmagnus))
      deallocate(tr%vmagnus); nullify(tr%vmagnus)
    endif

    call td_exp_end(tr%te)       ! clean propagator method
  end subroutine td_rti_end

  subroutine td_rti_run_zero_iter(h, tr)
    type(hamiltonian_type), intent(in) :: h
    type(td_rti_type), intent(inout) :: tr
    tr%v_old(:, :, 2) = h%vhxc(:, :)
    tr%v_old(:, :, 3) = tr%v_old(:, :, 2)
    tr%v_old(:, :, 1) = h%vhxc(:, :)
  end subroutine td_rti_run_zero_iter

  subroutine td_rti_dt(h, m, st, sys, tr, t, dt)
    type(hamiltonian_type), intent(inout) :: h
    type(mesh_type), intent(in) :: m
    type(states_type), intent(inout) :: st
    type(system_type), intent(in) :: sys
    type(td_rti_type), intent(inout) :: tr
    real(r8), intent(in) :: t, dt

    integer :: is
    logical :: self_consistent
    complex(r8), allocatable :: zpsi1(:, :, :, :)
    
    call push_sub('td_rti')

    self_consistent = .false.
    if(t<3*dt .and. (.not.h%ip_app)) then
       self_consistent = .true.
       allocate(zpsi1(m%np, st%dim, st%st_start:st%st_end, st%nik))
       zpsi1 = st%zpsi
    endif

    tr%v_old(:, :, 3) = tr%v_old(:, :, 2)
    tr%v_old(:, :, 2) = tr%v_old(:, :, 1)
    tr%v_old(:, :, 1) = h%vhxc(:, :)
    call dextrapolate(2, m%np*st%nspin, tr%v_old(:, :, 1:3), tr%v_old(:, :, 0), dt, dt)

    select case(tr%method)
    case(REVERSAL);             call td_rti2
    case(APP_REVERSAL);         call td_rti3
    case(EXPONENTIAL_MIDPOINT); call td_rti4
    case(MAGNUS);               call td_rti5
    end select

    if(self_consistent) then
      do
        tr%v_old(:, :, 3) = tr%v_old(:, :, 0)

        call zcalcdens(st, m%np, st%rho, .true.)
        call zh_calc_vhxc(h, m, st, sys)
        tr%v_old(:, :, 0) = h%vhxc
        h%vhxc = tr%v_old(:, :, 1)

        if( maxval((/ (dmf_nrm2(m, tr%v_old(:, is, 3)-tr%v_old(:, is, 0)),is=1,st%nspin) /)) < scf_threshold) exit

        st%zpsi = zpsi1
        select case(tr%method)
         case(REVERSAL);             call td_rti2
         case(APP_REVERSAL);         call td_rti3
         case(EXPONENTIAL_MIDPOINT); call td_rti4
         case(MAGNUS);               call td_rti5
        end select
      enddo
      deallocate(zpsi1)
    endif

    call pop_sub()
  contains

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
        call zh_calc_vhxc(h, m, st, sys)
        
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
      call push_sub('td_rti3')
      
      do ik = 1, st%nik
        do ist = st%st_start, st%st_end
           call td_exp_dt(tr%te, sys, h, st%zpsi(:,:, ist, ik), ik, dt/2._r8, t-dt)
        end do
      end do
      
      h%vhxc = tr%v_old(:, :, 0)
      do ik = 1, st%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, sys, h, st%zpsi(:,:, ist, ik), ik, dt/2._r8, t)
        end do
      end do
      
      call pop_sub()
    end subroutine td_rti3
    
    subroutine td_rti4
      integer :: ist, ik
      call push_sub('td_rti4')

      if(.not.h%ip_app) then
        call dextrapolate(2, m%np*st%nspin, tr%v_old(:, :, 0:2), h%vhxc, dt, -dt/M_TWO)
      end if
      
      do ik = 1, st%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, sys, h, st%zpsi(:,:, ist, ik), ik, dt, t - dt/2._r8)
        end do
      end do
      
      call pop_sub()
    end subroutine td_rti4

    subroutine td_rti5
      integer :: j, is, ist, k
      real(r8) :: time(2), x(3), f(3)
      real(r8), allocatable :: vaux(:, :, :)

      call push_sub('td_rti5')

      allocate(vaux(m%np, st%nspin, 2))

      time(1) = (M_HALF-sqrt(M_THREE)/6.0_r8)*dt
      time(2) = (M_HALF+sqrt(M_THREE)/6.0_r8)*dt

      if(.not.h%ip_app) then
        do j = 1, 2
           call dextrapolate(2, m%np*st%nspin, tr%v_old(:, :, 0:2), vaux(:, :, j), dt, time(j)-dt)
        enddo
      else
        vaux = M_ZERO
      endif

      do j = 1, 2
        if(h%ep%no_lasers > 0) then
          select case(h%gauge)
          case(1) ! length gauge
            call laser_field(h%ep%no_lasers, h%ep%lasers, t-dt+time(j), f)
            do k = 1, m%np
               call mesh_xyz(m, k, x)
               do is = 1, st%spin_channels
                  vaux(k, is, j) = vaux(k, is, j) + sum(x*f)
               enddo
            end do
          case(2) ! velocity gauge
            write(message(1),'(a)') 'Inconsistency in td_rti5'
            call write_fatal(1)
          end select
        endif
      enddo

      tr%vmagnus(:, :, 2)  = M_HALF*(vaux(:, :, 1) + vaux(:, :, 2))
      tr%vmagnus(:, :, 1) = (sqrt(M_THREE)/12.0_r8)*dt*(vaux(:, :, 2) - vaux(:, :, 1))

      do k = 1, st%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, sys, h, st%zpsi(:,:, ist, k), k, dt, M_ZERO, &
                         vmagnus = tr%vmagnus)
        end do
      end do

      deallocate(vaux)
      call pop_sub()
    end subroutine td_rti5
    
  end subroutine td_rti_dt

end module td_rti

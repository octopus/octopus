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

subroutine td_rti(h, m, st, sys, td, t)
  type(hamiltonian_type), intent(inout) :: h
  type(mesh_type), intent(in) :: m
  type(states_type), intent(inout) :: st
  type(system_type), intent(in) :: sys
  type(td_type), intent(inout) :: td
  real(r8), intent(in) :: t

  integer :: is  
  sub_name = 'td_rti'; call push_sub()

!!$  call dcopy(st%nspin*m%np, td%v_old(1, 1, 2), 1, td%v_old(1, 1, 3), 1)
!!$  call dcopy(st%nspin*m%np, td%v_old(1, 1, 1), 1, td%v_old(1, 1, 2), 1)
  td%v_old(:, :, 3) = td%v_old(:, :, 2)
  td%v_old(:, :, 2) = td%v_old(:, :, 1)
  select case(st%ispin)
  case(UNPOLARIZED, SPIN_POLARIZED)
    do is = 1, st%nspin
       td%v_old(:, is, 1) = h%vhartree(:) + h%vxc(:, is)
    enddo
  case(SPINORS)
     td%v_old(:, 1, 1) = h%vhartree(:) + h%vxc(:, 1)
     td%v_old(:, 2, 1) = h%vhartree(:) + h%vxc(:, 2)
     td%v_old(:, 3, 1) = h%vxc(:, 3)
     td%v_old(:, 4, 1) = h%vxc(:, 4)
  end select
  
  select case(td%evolution_method)
  case(OLD_REVERSAL)
    call td_rti1
  case(REVERSAL)
    call td_rti2
  case(APP_REVERSAL)
    if(t<3*td%dt) then
      call td_rti2
    else
      call td_rti3
    endif
  case(EXPONENTIAL_MIDPOINT)
    if(t<3*td%dt) then
      call td_rti2
    else
      call td_rti4
    endif
  case(SIMPLE_EXP)
    call td_rti0
  end select

  call pop_sub(); return
contains

  subroutine td_rti0
    integer is, ik, ist

    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(h, sys, td, ik, st%zpsi(:,:, ist, ik), td%dt, t-td%dt/M_TWO)
      end do
    end do
  end subroutine td_rti0

  ! Warning: this subroutine should only be used with LDA/GGA functionals
  subroutine td_rti1
    integer is, ik, ist
    real(r8), allocatable :: aux(:,:)
    complex(r8), allocatable :: zpsi1(:,:,:,:)
    
    sub_name = 'td_rti1'; call push_sub()
    
    allocate(aux(m%np, st%nspin))
    
    call xpolate_pot(td%dt/2._r8, td%dt, m%np, st%nspin, &
         td%v_old(:, :, 3), td%v_old(:, :, 2), td%v_old(:, :, 1), aux)
    
    h%VHartree = 0._r8; h%Vxc = aux
    allocate(zpsi1(0:m%np, st%dim, st%st_start:st%st_end, st%nik))
    zpsi1 = st%zpsi
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(h, sys, td, ik, st%zpsi(:,:, ist, ik), td%dt, t-td%dt)
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
        call td_dtexp(h, sys, td, ik, st%zpsi(:,:, ist, ik), td%dt, t-td%dt)
      end do
    end do
    
    call pop_sub(); return
  end subroutine td_rti1

  subroutine td_rti2
    real(r8), allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
    complex(r8), allocatable :: zpsi1(:,:,:,:)
    integer is, ik, ist

    sub_name = 'td_rti2'; call push_sub()

    if(.not.h%ip_app) then
      allocate(zpsi1(0:m%np, st%dim, st%st_start:st%st_end, st%nik))
      zpsi1 = st%zpsi ! store zpsi
    
      allocate(vhxc_t1(m%np, st%nspin))
      do is = 1, min(2,st%nspin)
        Vhxc_t1(1:m%np, is) = h%Vhartree(1:m%np) + h%Vxc(1:m%np, is)
      end do
      if(st%nspin == 4) vhxc_t1(:, 3:4) = h%vxc(:, 3:4)
      
      ! propagate dt with H(t-dt)
      do ik = 1, st%nik
        do ist = st%st_start, st%st_end
          call td_dtexp(h, sys, td, ik, st%zpsi(:,:, ist, ik), td%dt, t-td%dt)
        end do
      end do
    
      call zcalcdens(st, m%np, st%rho, .true.)
      call zhamiltonian_setup(h, m, st, sys)
    
      st%zpsi = zpsi1
      deallocate(zpsi1)
    
      ! store Vhxc at t
      allocate(vhxc_t2(m%np, st%nspin))
      do is = 1, min(2,st%nspin)
        Vhxc_t2(1:m%np, is) = h%Vhartree(1:m%np) + h%Vxc(1:m%np, is)
      end do
      if(st%nspin == 4) vhxc_t2(:, 3:4) = h%vxc(:, 3:4)

      h%Vhartree = 0._r8
      h%Vxc = Vhxc_t1
    end if

    ! propagate dt/2 with H(t-dt)
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(h, sys, td, ik, st%zpsi(:,:, ist, ik), td%dt/2._r8, t-td%dt)
      end do
    end do
    
    if(.not.h%ip_app) then
      deallocate(vhxc_t1)
    
      h%Vxc = Vhxc_t2
    end if

    ! propagate dt/2 with H(t)
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(h, sys, td, ik, st%zpsi(:,:, ist, ik), td%dt/2._r8, t)
      end do
    end do
    
    if(.not.h%ip_app) deallocate(vhxc_t2)

    call pop_sub(); return
  end subroutine td_rti2

  subroutine td_rti3
    integer is, ik, ist
    real(r8), allocatable :: aux(:,:)

    sub_name = 'td_rti3'; call push_sub()

    if(.not.h%ip_app) then
      allocate(aux(m%np, st%nspin))

      call xpolate_pot(td%dt, td%dt, m%np, st%nspin, &
           td%v_old(:, :, 3), td%v_old(:, :, 2), td%v_old(:, :, 1), aux)
      
      ! propagate dt/2 with H(t-dt)
      h%Vhartree = 0._r8
      h%vxc = td%v_old(:, :, 1)
    end if

    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(h, sys, td, ik, st%zpsi(:,:, ist, ik), td%dt/2._r8, t-td%dt)
      end do
    end do
    
    if(.not.h%ip_app) h%vxc = aux
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(h, sys, td, ik, st%zpsi(:,:, ist, ik), td%dt/2._r8, t)
      end do
    end do
 
    if(.not.h%ip_app) deallocate(aux)

    call pop_sub(); return
  end subroutine td_rti3

  subroutine td_rti4
    integer :: ist, ik

    sub_name = 'td_rti4'; call push_sub()

    call xpolate_pot(td%dt/2._r8, td%dt, m%np, st%nspin, &
         td%v_old(:, :, 3), td%v_old(:, :, 2), td%v_old(:, :, 1), h%vxc)
    
    h%vhartree = 0._r8
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(h, sys, td, ik, st%zpsi(:,:, ist, ik), td%dt, t - td%dt/2._r8)
      end do
    end do

    call pop_sub(); return
  end subroutine td_rti4

  subroutine xpolate_pot(t, dt, np, dim, pot2, pot1, pot0, pot)
    real(r8), intent(in)  :: t, dt
    integer, intent(in)   :: np, dim
    real(r8), intent(in)  :: pot0(np, dim), pot1(np, dim), pot2(np, dim)
    real(r8), intent(out) :: pot(np, dim)

    !pot = pot0 + t/dt       * (3._r8/2._r8*pot0 + 1._r8/2._r8*pot2 - 2.0_r8*pot1) + &
    !             t**2/dt**2 * (1._r8/2._r8*pot0 + 1._r8/2._r8*pot2 -        pot1)
    call dcopy(np*dim,                                                    pot0, 1, pot, 1)
    call daxpy(np*dim, (t/dt)*(3._r8/2._r8) + (t**2/dt**2)*(1._r8/2._r8), pot0, 1, pot, 1)
    call daxpy(np*dim, (t/dt)*(-2._r8)      + (t**2/dt**2)*(-1._r8),      pot1, 1, pot, 1)
    call daxpy(np*dim, (t/dt)*(1._r8/2._r8) + (t**2/dt**2)*(1._r8/2._r8), pot2, 1, pot, 1)
    
  end subroutine xpolate_pot

end subroutine td_rti

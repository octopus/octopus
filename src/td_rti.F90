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

subroutine td_rti(sys, h, td, t)
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h
  type(td_type), intent(inout) :: td
  real(r8), intent(in) :: t
  
  sub_name = 'td_rti'; call push_sub()
  
  select case(td%evolution_method)
  case(1)
    call td_rti1(sys%m, sys%st)
  case(2)
    call td_rti2(sys%m, sys%st)
#ifdef THREE_D
  case(3) ! split operator
    call td_rti3(sys, h, td, t)
#endif
  case(4)
    call td_rti4(sys%m, sys%st)
  end select

  call pop_sub()
  return

contains

  subroutine td_dtexp(ik, zpsi, timestep, t)
    integer, intent(in) :: ik
    complex(r8), intent(inout) :: zpsi(0:sys%m%np, sys%st%dim)
    real(r8), intent(in) :: timestep, t

    integer, parameter :: order = 4

    complex(r8) :: zfact
    complex(r8), allocatable :: zpsi1(:,:), hzpsi1(:,:), grad(:,:)
    real(r8) :: x(3), f(3)
    integer i, k, idim

    sub_name = 'td_dtexp'; call push_sub()

    allocate(zpsi1(0:sys%m%np, sys%st%dim), hzpsi1(sys%m%np, sys%st%dim))

    zfact = 1._r8
    zpsi1 = zpsi
    do i = 1, order ! forth order method
      zfact = zfact*(-M_zI*timestep)/i
      call zHpsi(h, sys, ik, zpsi1, hzpsi1)
      
      ! apply lasers
      if(td%no_lasers > 0) then
        select case(td%gauge)
        case(1) ! length gauge
          call laser_field(td%no_lasers, td%lasers, t, f)
          
          do k = 1, sys%m%np
            call mesh_xyz(sys%m, k, x)
            hzpsi1(k,:) = hzpsi1(k,:) + sum(x*f) * zpsi1(k,:)
          end do
          
        case(2) ! velocity gauge
          call laser_vector_field(td%no_lasers, td%lasers, t, f)
          allocate(grad(3, sys%m%np))
          do idim = 1, sys%st%dim
            call zmesh_derivatives(sys%m, zpsi1(:, idim), grad=grad)
            do k = 1, sys%m%np
              hzpsi1(k, idim) = hzpsi1(k, idim) - M_zI * sum(f(:)*grad(:, k)) + &
                   sum(f**2)/2._r8 * zpsi1(k, idim)
            end do
          end do
          deallocate(grad)
        end select
      end if
      
      ! absorbing potential
      if(td%ab .eq. 1) then
        do idim = 1, sys%st%dim
          hzpsi1(:, idim) = hzpsi1(:, idim) + M_zI*td%ab_pot(:)*zpsi1(1:, idim)
        end do
      end if
      
      zpsi(1:,:) = zpsi(1:,:) + zfact*hzpsi1(:,:)
      
      if(i .ne. order) zpsi1(1:,:) = hzpsi1(:,:)
    end do
    
    deallocate(zpsi1, hzpsi1)
    call pop_sub(); return
  end subroutine td_dtexp
  
  ! Warning: this subroutine should only be used with LDA/GGA functionals
  subroutine td_rti1(m, st)
    type(mesh_type), intent(IN) :: m
    type(states_type), intent(inout) :: st

    integer is, ik, ist
    real(r8), allocatable :: aux(:,:)
    complex(r8), allocatable :: zpsi1(:,:,:,:)
    
    sub_name = 'td_rti1'; call push_sub()

    allocate(aux(m%np, st%nspin))
    
    do is = 1, sys%st%nspin
      aux(:, is) = 1.875_r8*(h%VHartree(:) + h%Vxc(:, is)) &
           -1.25_r8*td%v_old1(:, is) + 0.375_r8*td%v_old2(:, is)
    end do
    
    td%v_old2 = td%v_old1
    do is = 1, st%nspin
      td%v_old1(:, is) = h%VHartree(:) + h%Vxc(:, is)
    end do
    
    h%VHartree = 0._r8; h%Vxc = aux
    allocate(zpsi1(0:m%np, st%dim, st%st_start:st%st_end, st%nik))
    zpsi1 = st%zpsi
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(ik, st%zpsi(:,:, ist, ik), td%dt, t-td%dt)
      end do
    end do
    st%zpsi = zpsi1
    deallocate(zpsi1)
    
    call zcalcdens(st, m%np, aux, .true.)
    st%rho = (st%rho +  aux) / 2.0_r8
    deallocate(aux)
    
    call zhamiltonian_setup(h, sys)
    
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(ik, sys%st%zpsi(:,:, ist, ik), td%dt, t-td%dt)
      end do
    end do
    
    call pop_sub(); return
  end subroutine td_rti1

  subroutine td_rti2(m, st)
    type(mesh_type), intent(inout) :: m
    type(states_type), intent(inout) :: st
    
    real(r8), allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
    complex(r8), allocatable :: zpsi1(:,:,:,:)
    integer is, ik, ist
#if defined(HAVE_MPI) && defined(MPI_TD)
    real(r8), allocatable :: reduce_rho(:,:) ! temporary to do MPI_reduce
    integer :: ierr
#endif

    sub_name = 'td_rti2'; call push_sub()

    allocate(zpsi1(0:m%np, st%dim, st%st_start:st%st_end, st%nik))
    zpsi1 = st%zpsi ! store zpsi
    
    allocate(vhxc_t1(m%np, st%nspin))
    do is = 1, st%nspin ! store Vhxc
      Vhxc_t1(:, is) = h%VHartree(:) + h%Vxc(:, is)
    end do
    
    ! propagate dt with H(t-dt)
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(ik, st%zpsi(:,:, ist, ik), td%dt, t-td%dt)
      end do
    end do
    
    call zcalcdens(st, m%np, st%rho, .true.)
    call zhamiltonian_setup(h, sys)
    
    st%zpsi = zpsi1
    deallocate(zpsi1)
    
    ! store Vhxc at t
    allocate(vhxc_t2(m%np, st%nspin))
    do is = 1, st%nspin
      Vhxc_t2(:, is) = h%VHartree(:) + h%Vxc(:, is)
    end do
    
    ! propagate dt/2 with H(t-dt)
    h%Vhartree = 0._r8
    h%Vxc = Vhxc_t1
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(ik, st%zpsi(:,:, ist, ik), td%dt/2._r8, t-td%dt)
      end do
    end do
    
    deallocate(vhxc_t1)
    
    ! propagate dt/2 with H(t)
    h%Vxc = Vhxc_t2
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(ik, st%zpsi(:,:, ist, ik), td%dt/2._r8, t)
      end do
    end do
    
    deallocate(vhxc_t2)

    call pop_sub(); return
  end subroutine td_rti2

  subroutine td_rti4(m, st)
    type(mesh_type), intent(inout) :: m
    type(states_type), intent(inout) :: st
    
    integer is, ik, ist
#if defined(HAVE_MPI) && defined(MPI_TD)
    real(r8), allocatable :: reduce_rho(:,:) ! temporary to do MPI_reduce
    integer :: ierr
#endif

    real(r8), allocatable :: aux(:,:)

    sub_name = 'td_rti2'; call push_sub()

    allocate(aux(m%np, st%nspin))
    
    do is = 1, sys%st%nspin
      aux(:, is) = 1.875_r8*(h%VHartree(:) + h%Vxc(:, is)) &
           -1.25_r8*td%v_old1(:, is) + 0.375_r8*td%v_old2(:, is)
    end do
    
    td%v_old2 = td%v_old1
    do is = 1, st%nspin
      td%v_old1(:, is) = h%VHartree(:) + h%Vxc(:, is)
    end do
    
    ! propagate dt/2 with H(t-dt)
    h%Vhartree = 0._r8
    h%vxc = td%v_old1
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(ik, st%zpsi(:,:, ist, ik), td%dt/2._r8, t-td%dt)
      end do
    end do
    
    h%vxc = aux
    do ik = 1, st%nik
      do ist = st%st_start, st%st_end
        call td_dtexp(ik, st%zpsi(:,:, ist, ik), td%dt/2._r8, t)
      end do
    end do
    
    call pop_sub(); return
  end subroutine td_rti4

end subroutine td_rti

#ifdef THREE_D
#include "td_rti3.F90"
#endif

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

#include "config_F90.h"

module lcao
  use global
  use liboct
  use spline
  use mesh
  use system
  use hamiltonian
  use states
  use mix

  implicit none

  private
  public :: lcao_dens, lcao_init, lcao_wf, lcao_end

type lcao_type
  integer           :: mode
  integer           :: dim
  R_TYPE , pointer  :: psis(:, :, :, :)
  R_TYPE , pointer  :: hamilt    (:, :, :), &
                       k_plus_psv(:, :, :), &
                       s         (:, :, :)
  logical, pointer  :: atoml(:,:)
end type

type(lcao_type) :: lcao_data

integer, parameter :: MEM_INTENSIVE = 0, &
                      CPU_INTENSIVE = 1

contains

!builds a density which is the sum of the atomic densities
subroutine lcao_dens(sys, nspin, rho)
  type(system_type), intent(IN) :: sys
  integer, intent(in) :: nspin
  real(r8), intent(out) :: rho(sys%m%np, nspin)
  
  integer :: ia, is
  real(r8) :: r
  type(specie_type), pointer :: s
  type(atom_type),   pointer :: a

  sub_name = 'lcao_dens'; call push_sub()
  
  rho = 0._r8
  do ia = 1, sys%natoms
    a => sys%atom(ia) ! shortcuts
    s => a%spec
    
    select case(s%label(1:5))
    case('jelli', 'point')
      call from_jellium(sys%m, rho(:, 1))
    case('usdef')
      call from_userdef(sys%m, rho(:, 1))
    case default
      call from_pseudopotential(sys%m)
    end select
  end do

  ! we now renormalize the density (necessary if we have a charged system)
  r = dmesh_integrate(sys%m, rho(:, 1))
  write(message(1),'(a,f13.6)')'Info: Unnormalized total charge = ', r
  call write_info(1)
  r = sys%st%qtot/r
  do is = 1, nspin
    rho(:, is) = r*rho(:, is)
  end do
  r = dmesh_integrate(sys%m, rho(:, 1))
  write(message(1),'(a,f13.6)')'Info: Renormalized total charge = ', r
  call write_info(1)

  call pop_sub()
contains
  subroutine from_jellium(m, rho)
    type(mesh_type), intent(in) :: m
    real(r8), intent(inout) :: rho(m%np)

    integer :: i, in_points
    real(r8) :: r

    do i = 1, m%np
      call mesh_r(m, i, r, a=a%x)
      if(r <= s%jradius) then
        in_points = in_points + 1
      end if
    end do
    
    if(in_points > 0) then
      do i = 1, m%np
        call mesh_r(m, i, r, a=a%x)
        if(r <= s%jradius) then
          rho(i) = rho(i) + real(s%Z_val, r8)/(in_points*m%vol_pp)
        end if
      end do
    end if
  end subroutine from_jellium
  
  subroutine from_userdef(m, rho)
    type(mesh_type), intent(in) :: m
    real(r8), intent(inout) :: rho(m%np)

    integer :: i

    do i = 1, m%np
      call mesh_r(m, i, r, a=a%x)
      rho(i) = rho(i) + real(s%Z_val, r8)/(m%np*m%vol_pp)
    end do
  end subroutine from_userdef

  subroutine from_pseudopotential(m)
    type(mesh_type), intent(in) :: m

    integer :: i, l, is
    real(r8) :: r
    R_TYPE :: psi1, psi2

    do i = 1, m%np
      call mesh_r(m, i, r, a=a%x)
#if defined(THREE_D)
      do l = 0 , s%ps%L_max_occ
        if(r >= r_small) then
          select case(sys%st%spin_channels)
          case(1)
             psi1 = splint(s%ps%Ur(l, 1), r)
             rho(i, 1) = rho(i, 1) + s%ps%occ(l, 1)*psi1*psi1*(r**(2*l))/(4*M_PI)
          case(2)
             psi1 = splint(s%ps%ur(l, 1), r)
             psi2 = splint(s%ps%ur(l, 2), r)
             rho(i, 1) = rho(i, 1) + (s%ps%occ(l, 1)*psi1*psi1 + s%ps%occ(l, 2)*psi2*psi2) * &
                                     (r**(2*l))/(4*M_PI)
             rho(i, sys%st%nspin) = rho(i, sys%st%nspin) + &
                                    (s%ps%occ(l, 1)*psi1*psi1 - s%ps%occ(l, 2)*psi2*psi2) * &
                                    (r**(2*l))/(4*M_PI)
          end select          
        end if
      end do
#elif defined(ONE_D)
      rho(i, 1) = rho(i, 1) + s%z_val*exp(-r**2)/sqrt(m_pi)
#endif
    end do
    
  end subroutine from_pseudopotential
end subroutine lcao_dens

subroutine lcao_init(sys, h)
  type(system_type), intent(IN)      :: sys
  type(hamiltonian_type), intent(IN) :: h

  integer :: norbs, i, ispin, a, ik, n1, i1, l1, lm, lm1, d1, n2, i2, l2, lm2, d2
  integer, parameter :: orbs_local = 2

  R_TYPE, allocatable :: psi1(:,:), psi2(:,:), hpsi(:,:)
  real(r8) :: s

  real(r8) :: uVpsi

  sub_name = 'lcao_init'; call push_sub

  ! Counting
  allocate(lcao_data%atoml(sys%natoms, 0:3))
  lcao_data%atoml = .true.
  norbs = 0
  atoms_loop: do i1 = 1, sys%natoms
      l_loop: do l1 = 0, sys%atom(i1)%spec%ps%L_max_occ
           if(sum(sys%atom(i1)%spec%ps%occ(l1, :)).ne.0.0_r8) then
              norbs = norbs + (2*l1+1)
           else
              lcao_data%atoml(i1, l1) = .false.
           endif
      end do l_loop
  enddo atoms_loop

  select case(sys%st%ispin)
   case(1);
   case(2); ! No need to multiply by two, since each spin-channel goes to a k-subspace.
   case(3); norbs = norbs * 2
  end select

  lcao_data%dim = norbs
  if(norbs < sys%st%nst) then
    write(message(1), '(a)') 'Internal bug: LCAO basis dimension, norbs, is smaller than'
    write(message(2), '(a)') 'number of required states.'
    call write_fatal(2)
  endif
  write(message(1), '(a,i6)') 'Info: LCAO basis dimension: ', lcao_data%dim
  write(message(2), '(a)')    '      (not considering spin or k-points)'
  call write_info(2)

  allocate(lcao_data%psis(0:sys%m%np, sys%st%dim, norbs, sys%st%nik))
  lcao_data%psis = 0._r8
  do ik = 1, sys%st%nik
     n1 = 1
     do i1 = 1, sys%natoms
        do l1 = 0, sys%atom(i1)%spec%ps%L_max_occ
           if(.not. lcao_data%atoml(i1, l1)) cycle
           do lm1 = -l1, l1
              do d1 = 1, sys%st%dim
                 ispin = states_spin_channel(sys%st%ispin, ik, d1)
                 call get_wf(sys, i1, l1, lm1, ispin, lcao_data%psis(:, d1, n1, ik))
                 n1 = n1 + 1
              end do
           end do
        end do
     end do
  end do
!!$  end if

  ! Allocation of variables
  allocate(lcao_data%hamilt    (sys%st%nik, norbs, norbs), &
           lcao_data%s         (sys%st%nik, norbs, norbs), &
           lcao_data%k_plus_psv(sys%st%nik, norbs, norbs))

  call pop_sub()
end subroutine lcao_init

subroutine lcao_end
  sub_name = 'lcao_end'; call push_sub()

  if(lcao_data%mode == MEM_INTENSIVE) deallocate(lcao_data%psis)
  deallocate(lcao_data%hamilt, lcao_data%s, lcao_data%atoml)

  call pop_sub()
end subroutine lcao_end

subroutine lcao_wf(sys, h)
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(in) :: h
  
  integer, parameter :: orbs_local = 2

  integer :: a, idim, i, ispin, lm, ik, n1, n2, i1, i2, l1, l2, lm1, lm2, d1, d2
  integer :: norbs, mode
  R_TYPE, allocatable :: hpsi(:,:)
  R_TYPE, allocatable :: s(:,:)
  real(r8) :: uvpsi

  ! variables for dsyev (LAPACK)
  integer :: lwork, info
  R_TYPE, allocatable :: work(:)
  real(r8), allocatable :: rwork(:), w(:)

  sub_name = 'lcao_wf'; call push_sub()

  norbs = lcao_data%dim
  mode  = lcao_data%mode

  ! Hamiltonian and overlap matrices.
  allocate(hpsi(sys%m%np, sys%st%dim))
  do ik = 1, sys%st%nik
    do n1 = 1, lcao_data%dim
       call R_FUNC(Hpsi)(h, sys%m, sys%st, sys, ik, lcao_data%psis(:, :, n1, ik), hpsi(:, :))
       do n2 = 1, lcao_data%dim
              lcao_data%hamilt(ik, n1, n2) = &
                        R_FUNC(states_dotp)(sys%m, sys%st%dim, &
                                            hpsi, lcao_data%psis(1:, : ,n2, ik))
              lcao_data%s(ik, n1, n2) = &
                        R_FUNC(states_dotp)(sys%m, sys%st%dim, &
                                            lcao_data%psis(1:, :, n1, ik), lcao_data%psis(1:, : ,n2, ik))
       enddo
    enddo
  enddo

  do ik = 1, sys%st%nik

    lwork = 3*norbs - 1
    allocate(work(lwork), w(norbs), rwork(lwork), s(norbs, norbs))
    s(1:norbs, 1:norbs) = lcao_data%s(ik, 1:norbs, 1:norbs)
#ifdef COMPLEX_WFNS
    call zhegv (1, 'v', 'u', norbs, lcao_data%hamilt(ik, :, :), norbs, s, norbs, w, work, lwork, rwork, info)
#else
    call dsygv (1, 'v', 'u', norbs, lcao_data%hamilt(ik, :, :), norbs, s, norbs, w, work, lwork, info)
#endif
    if(info.ne.0) then
      write(message(1),'(a,i5)') 'LAPACK "zhegv/dsygv" returned error code ', info
      call write_fatal(1)
    endif
    sys%st%eigenval(1:sys%st%nst, ik) = w(1:sys%st%nst)
    deallocate(work, w, s, rwork)
    sys%st%R_FUNC(psi)(:,:,:, ik) = R_TOTYPE(0.0_r8)
    do n1 = 1, lcao_data%dim
       do n2 = 1, sys%st%nst
          do d1 = 1, sys%st%dim
            sys%st%R_FUNC(psi)(:, d1, n2, ik) = sys%st%R_FUNC(psi)(:, d1, n2, ik) + &
              lcao_data%hamilt(ik, n1, n2) * lcao_data%psis(:, d1, n1, ik)
          enddo
       enddo
    enddo
  end do

  deallocate(hpsi)
  call pop_sub(); return
end subroutine lcao_wf

subroutine get_wf(sys, i, l, lm, ispin, psi)
  type(system_type), intent(IN) :: sys
  integer, intent(in)   :: i, l, lm, ispin
  R_TYPE, intent(out) :: psi(0:sys%m%np)
    
  integer :: j, d2
  real(r8) :: x(3), a(3), r, p, ylm, g(3)
  type(spline_type), pointer :: s

  sub_name = 'get_wf'; call push_sub()
    
  a = sys%atom(i)%x
  psi(0) = 0.0_r8
  if(sys%atom(i)%spec%local) then
    ! add a couple of harmonic oscilator functions
  else
    s => sys%atom(i)%spec%ps%ur(l, ispin)
    do j = 1, sys%m%np
      call mesh_r(sys%m, j, r, x=x, a=a)
      p = splint(s, r)
      ylm = oct_ylm(x(1), x(2), x(3), l, lm)
      if(r > 0._r8) then
        psi(j) = p * ylm * r**l
      end if
    end do
  end if

  call pop_sub(); return    
end subroutine get_wf

end module lcao

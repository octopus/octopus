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

module lcao
  use hamiltonian

  implicit none

  private
  public :: lcao_dens, lcao_init, lcao_wf, lcao_end

type lcao_type
  !integer           :: mode
  integer           :: state ! 0 => non-initialized;
                             ! 1 => initialized (k, s and v1 matrices filled)
  integer           :: dim
  R_TYPE , pointer  :: psis(:, :, :, :)
  ! hamilt stores the Hamiltonian in the LCAO subspace;
  ! s is the overlap matrix;
  ! k is the kinetic + spin orbit operator matrix;
  ! v is the potential.
  R_TYPE , pointer  :: hamilt    (:, :, :), &
                       s         (:, :, :), &
                       k         (:, :, :), &
                       v         (:, :, :)
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

  call push_sub('lcao_dens')
  
  rho = M_ZERO
  do ia = 1, sys%natoms
    a => sys%atom(ia) ! shortcuts
    s => a%spec

    if(conf%dim==1.or.conf%dim==2) then
      call from_userdef(sys%m)
    else
      select case(s%label(1:5))
      case('usdef')
        call from_userdef(sys%m)
      case('jelli', 'point')
        call from_jellium(sys%m)
      case default
        call from_pseudopotential(sys%m)
      end select
    end if
  end do

  ! we now renormalize the density (necessary if we have a charged system)
  r = M_ZERO
  do is = 1, sys%st%spin_channels
     r = r + dmf_integrate(sys%m, rho(:, is))
  end do
  write(message(1),'(a,f13.6)')'Info: Unnormalized total charge = ', r
  call write_info(1)
  r = sys%st%qtot/r
  rho = r*rho
  r = M_ZERO
  do is = 1, sys%st%spin_channels
     r = r + dmf_integrate(sys%m, rho(:, is))
  end do
  write(message(1),'(a,f13.6)')'Info: Renormalized total charge = ', r
  call write_info(1)

  call pop_sub()
contains
  subroutine from_userdef(m)
    type(mesh_type), intent(in) :: m

    integer :: i

    do i = 1, m%np
      call mesh_r(m, i, r, a=a%x)
      rho(i, 1) = rho(i, 1) + real(s%Z_val, r8)/(m%np*m%vol_pp)
    end do
  end subroutine from_userdef

  subroutine from_jellium(m)
    type(mesh_type), intent(in) :: m

    integer :: i, in_points
    real(r8) :: r

    in_points = 0
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
          rho(i, 1) = rho(i, 1) + real(s%Z_val, r8)/(in_points*m%vol_pp)
        end if
      end do
    end if
  end subroutine from_jellium
  
  subroutine from_pseudopotential(m)
    type(mesh_type), intent(in) :: m
    
    integer :: i, n, is
    real(r8) :: r
    R_TYPE :: psi1, psi2
    do i = 1, m%np
      call mesh_r(m, i, r, a=a%x)
      do n = 1, s%ps%conf%p
        if(r >= r_small) then
          select case(sys%st%spin_channels)
          case(1)
            psi1 = splint(s%ps%Ur(n, 1), r)
            rho(i, 1) = rho(i, 1) + s%ps%conf%occ(n, 1)*psi1*psi1 /(4*M_PI)  
          case(2)
            ! We will build a spin-unpolarized density, also for spin-polarized calculations.
            psi1 = splint(s%ps%Ur(n, 1), r)
            rho(i, 1) = rho(i, 1) + M_HALF*sum(s%ps%conf%occ(n, :))*psi1*psi1 / (M_FOUR*M_PI)
            rho(i, 2) = rho(i, 1)
          end select
        end if
      end do
    end do

  end subroutine from_pseudopotential
end subroutine lcao_dens

subroutine lcao_init(sys, h)
  type(system_type), intent(IN)      :: sys
  type(hamiltonian_type), intent(IN) :: h

  integer :: norbs, i, ispin, a, ik, n1, i1, l, l1, lm, lm1, d1, n2, i2, l2, lm2, d2
  integer, parameter :: orbs_local = 2

  R_TYPE, allocatable :: psi1(:,:), psi2(:,:), hpsi(:,:)
  real(r8) :: s

  real(r8) :: uVpsi

  if(conf%dim.ne.3) return
  if(lcao_data%state == 1) return

  call push_sub('lcao_init')

  ! Counting
  allocate(lcao_data%atoml(sys%natoms, 6))
  lcao_data%atoml = .true.
  norbs = 0
  atoms_loop: do i1 = 1, sys%natoms
    l_loop: do l1 = 1, sys%atom(i1)%spec%ps%conf%p
      l = sys%atom(i1)%spec%ps%conf%l(l1)
      if(sum(sys%atom(i1)%spec%ps%conf%occ(l1, :)).ne.M_ZERO) then
        norbs = norbs + (2*l+1)
      else
        lcao_data%atoml(i1, l1) = .false.
      endif
    end do l_loop
  end do atoms_loop

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
  call write_info(1)

  allocate(lcao_data%psis(sys%m%np, sys%st%dim, norbs, sys%st%nik))
  lcao_data%psis = 0._r8
  do ik = 1, sys%st%nik
     n1 = 1
     do i1 = 1, sys%natoms
        do l1 = 1, sys%atom(i1)%spec%ps%conf%p
           l = sys%atom(i1)%spec%ps%conf%l(l1)
           if(.not. lcao_data%atoml(i1, l1)) cycle
           do lm1 = -l, l
              do d1 = 1, sys%st%dim
                 ispin = states_spin_channel(sys%st%ispin, ik, d1)
                 call get_wf(sys, i1, l1, lm1, ispin, lcao_data%psis(:, d1, n1, ik))
                 n1 = n1 + 1
              end do
           end do
        end do
     end do
  end do

  ! Allocation of variables
  allocate(lcao_data%hamilt (norbs, norbs, sys%st%nik), &
           lcao_data%s      (norbs, norbs, sys%st%nik), &
           lcao_data%k      (norbs, norbs, sys%st%nik), &
           lcao_data%v      (norbs, norbs, sys%st%nik))

  ! Overlap and kinetic+so matrices.
  allocate(hpsi(sys%m%np, sys%st%dim))
  do ik = 1, sys%st%nik
    do n1 = 1, lcao_data%dim
      call X(kinetic) (h, ik, sys%m, sys%st, lcao_data%psis(:, :, n1, ik), hpsi(:, :))
      ! Relativistic corrections...
      select case(h%reltype)
      case(NOREL)
#if defined(COMPLEX_WFNS) && defined(R_TCOMPLEX)
      case(SPIN_ORBIT)
        call zso (h, sys%m, sys%natoms, sys%atom, sys%st%dim, ik, lcao_data%psis(:, :, n1, ik), hpsi(:, :))
#endif
      case default
        message(1) = 'Error: Internal.'
      call write_fatal(1)
      end select 
 
      do n2 = n1, lcao_data%dim
        lcao_data%k(n1, n2, ik) = X(states_dotp)(sys%m, sys%st%dim, &
             hpsi, lcao_data%psis(1:, : ,n2, ik))
        lcao_data%s(n1, n2, ik) = X(states_dotp)(sys%m, sys%st%dim, &
             lcao_data%psis(:, :, n1, ik), lcao_data%psis(:, : ,n2, ik))
      end do
      
    end do
  end do
  deallocate(hpsi)

  lcao_data%state = 1
  call pop_sub(); return
end subroutine lcao_init

subroutine lcao_end
  call push_sub('lcao_end')

  if(associated(lcao_data%hamilt)) then
    deallocate(lcao_data%hamilt, lcao_data%s, lcao_data%k, lcao_data%v)
  endif
  if(associated(lcao_data%psis)) deallocate(lcao_data%psis)

  lcao_data%state = 0  
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

  if(conf%dim.ne.3) return
  call push_sub('lcao_wf')

  norbs = lcao_data%dim

  ! Hamiltonian and overlap matrices.
  allocate(hpsi(sys%m%np, sys%st%dim))
  do ik = 1, sys%st%nik
    do n1 = 1, lcao_data%dim
      hpsi = M_ZERO
      call X(vlpsi) (h, sys%m, sys%st, ik, lcao_data%psis(:, :, n1, ik), hpsi(:, :))
      if(sys%nlpp) call X(vnlpsi) (ik, sys%m, sys%st, sys, lcao_data%psis(:, :, n1, ik), hpsi(:, :))
      do n2 = n1, lcao_data%dim
        lcao_data%v(n1, n2, ik) = X(states_dotp)(sys%m, sys%st%dim, &
                                                      hpsi, lcao_data%psis(1:, : ,n2, ik))
        lcao_data%hamilt(n1, n2, ik) = lcao_data%k(n1, n2, ik) + lcao_data%v(n1 , n2, ik)
      end do
    end do
  end do
  
  do ik = 1, sys%st%nik

    lwork = 5*norbs
    allocate(work(lwork), w(norbs), rwork(max(1,3*norbs-2)), s(norbs, norbs))
    s(1:norbs, 1:norbs) = lcao_data%s(1:norbs, 1:norbs, ik)
#ifdef COMPLEX_WFNS
    call zhegv (1, 'v', 'u', norbs, lcao_data%hamilt(1, 1, ik), norbs, s(1, 1), norbs, w(1), work(1), lwork, rwork(1), info)
#else
    call dsygv (1, 'v', 'u', norbs, lcao_data%hamilt(1, 1, ik), norbs, s(1, 1), norbs, w(1), work(1), lwork, info)
#endif
    if(info.ne.0) then
      write(message(1),'(a,i5)') 'LAPACK "zhegv/dsygv" returned error code ', info
      call write_fatal(1)
    endif
    sys%st%eigenval(1:sys%st%nst, ik) = w(1:sys%st%nst)
    deallocate(work, w, s, rwork)
    sys%st%X(psi)(:,:,:, ik) = R_TOTYPE(0.0_r8)

    ! Change of base
    call X(gemm)('N', 'N', sys%m%np*sys%st%dim, sys%st%nst, lcao_data%dim, &
                 R_TOTYPE(M_ONE),                                         &
                 lcao_data%psis(1, 1, 1, ik), sys%m%np*sys%st%dim,        &
                 lcao_data%hamilt(1, 1, ik), norbs,                       &
                 R_TOTYPE(M_ZERO),                                        &
                 sys%st%X(psi)(1, 1, 1, ik), sys%m%np*sys%st%dim)
 
   end do

  deallocate(hpsi)
  call pop_sub()
end subroutine lcao_wf

subroutine get_wf(sys, i, l, lm, ispin, psi)
  type(system_type), intent(IN) :: sys
  integer, intent(in)   :: i, l, lm, ispin
  R_TYPE, intent(out) :: psi(sys%m%np)
    
  integer :: j, d2, ll
  real(r8) :: x(3), a(3), r, p, ylm, g(3)
  type(spline_type), pointer :: s

  call push_sub('get_wf')
    
  a = sys%atom(i)%x
  if(sys%atom(i)%spec%local) then
    ! add a couple of harmonic oscilator functions
  else
    s => sys%atom(i)%spec%ps%Ur(l, ispin)

    ll = sys%atom(i)%spec%ps%conf%l(l)
    do j = 1, sys%m%np
      call mesh_r(sys%m, j, r, x=x, a=a)
      p = splint(s, r)
      ylm = oct_ylm(x(1), x(2), x(3), ll, lm)
      psi(j) = p * ylm
    end do
  end if
 
  call pop_sub()
end subroutine get_wf

end module lcao

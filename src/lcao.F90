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
  use global
  use lib_oct
  use lib_oct_gsl_spline
  use lib_basic_alg
  use lib_adv_alg
  use mesh
  use specie
  use atom
  use geometry
  use hamiltonian

  implicit none

  private
  public :: lcao_init, lcao_wf, lcao_end

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

subroutine lcao_init(m, st, geo, h)
  type(mesh_type),        intent(IN) :: m
  type(states_type),      intent(IN) :: st
  type(geometry_type),    intent(IN) :: geo
  type(hamiltonian_type), intent(IN) :: h

  integer :: norbs, ispin, ik, n1, i1, l, l1, lm1, d1, n2
  integer, parameter :: orbs_local = 2

  R_TYPE, allocatable :: hpsi(:,:)

  if(conf%dim.ne.3) return
  if(lcao_data%state == 1) return

  call push_sub('lcao_init')

  ! Counting
  allocate(lcao_data%atoml(geo%natoms, 6))
  lcao_data%atoml = .true.
  norbs = 0
  atoms_loop: do i1 = 1, geo%natoms
    l_loop: do l1 = 1, geo%atom(i1)%spec%ps%conf%p
      l = geo%atom(i1)%spec%ps%conf%l(l1)
      if(sum(geo%atom(i1)%spec%ps%conf%occ(l1, :)).ne.M_ZERO) then
        norbs = norbs + (2*l+1)
      else
        lcao_data%atoml(i1, l1) = .false.
      endif
    end do l_loop
  end do atoms_loop

  select case(st%d%ispin)
  case(UNPOLARIZED);
  case(SPIN_POLARIZED); ! No need to multiply by two, since each spin-channel goes to a k-subspace.
  case(SPINORS); norbs = norbs * 2
  end select

  lcao_data%dim = norbs
  if(norbs < st%nst) then
    write(message(1), '(a)') 'Internal bug: LCAO basis dimension, norbs, is smaller than'
    write(message(2), '(a)') 'number of required states.'
    call write_fatal(2)
  endif
  write(message(1), '(a,i6)') 'Info: LCAO basis dimension: ', lcao_data%dim
  call write_info(1)

  allocate(lcao_data%psis(m%np, st%dim, norbs, st%nik))
  lcao_data%psis = M_ZERO
  do ik = 1, st%nik
     n1 = 1
     do i1 = 1, geo%natoms
        do l1 = 1, geo%atom(i1)%spec%ps%conf%p
           l = geo%atom(i1)%spec%ps%conf%l(l1)
           if(.not. lcao_data%atoml(i1, l1)) cycle
           do lm1 = -l, l
              do d1 = 1, st%dim
                 ispin = states_spin_channel(st%d%ispin, ik, d1)
                 call atom_get_wf(m, geo%atom(i1), l1, lm1, ispin, lcao_data%psis(:, d1, n1, ik))
                 n1 = n1 + 1
              end do
           end do
        end do
     end do
  end do

  ! Allocation of variables
  allocate(lcao_data%hamilt (norbs, norbs, st%nik), &
           lcao_data%s      (norbs, norbs, st%nik), &
           lcao_data%k      (norbs, norbs, st%nik), &
           lcao_data%v      (norbs, norbs, st%nik))

  ! Overlap and kinetic+so matrices.
  allocate(hpsi(m%np, st%dim))
  do ik = 1, st%nik
    do n1 = 1, lcao_data%dim
      call X(kinetic) (h, m, lcao_data%psis(:, :, n1, ik), hpsi(:, :), ik)
      ! Relativistic corrections...
      select case(h%reltype)
      case(NOREL)
#if defined(COMPLEX_WFNS) && defined(R_TCOMPLEX)
      case(SPIN_ORBIT)
        call zso (h, m, lcao_data%psis(:, :, n1, ik), hpsi(:, :), st%dim, ik)
#endif
      case default
        message(1) = 'Error: Internal.'
      call write_fatal(1)
      end select
 
      do n2 = n1, lcao_data%dim
        lcao_data%k(n1, n2, ik) = X(states_dotp)(m, st%dim, hpsi, lcao_data%psis(1:, : ,n2, ik))
        lcao_data%s(n1, n2, ik) = X(states_dotp)(m, st%dim, lcao_data%psis(:, :, n1, ik), &
                                                 lcao_data%psis(:, : ,n2, ik))
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

subroutine lcao_wf(m, st, h)
  type(mesh_type),        intent(IN)    :: m
  type(states_type),      intent(inout) :: st
  type(hamiltonian_type), intent(IN)    :: h

  integer, parameter :: orbs_local = 2

  integer :: np, dim, nst, ik, n1, n2
  integer :: norbs
  R_TYPE, allocatable :: hpsi(:,:)
  FLOAT, allocatable :: ev(:)

  if(conf%dim.ne.3) return
  call push_sub('lcao_wf')

  norbs = lcao_data%dim
  np = m%np
  dim = st%dim
  nst = st%nst

  ! Hamiltonian and overlap matrices.
  allocate(hpsi(np, dim))
  do ik = 1, st%nik
    do n1 = 1, lcao_data%dim
      hpsi = M_ZERO
      call X(vlpsi) (h, m, lcao_data%psis(:, :, n1, ik), hpsi(:, :), ik)
      if (h%ep%nvnl > 0) call X(vnlpsi) (h, m, lcao_data%psis(:, :, n1, ik), hpsi(:, :), ik)
      do n2 = n1, lcao_data%dim
        lcao_data%v(n1, n2, ik) = X(states_dotp)(m, dim, hpsi, lcao_data%psis(1:, : ,n2, ik))
        lcao_data%hamilt(n1, n2, ik) = lcao_data%k(n1, n2, ik) + lcao_data%v(n1 , n2, ik)
      end do
    end do
  end do
  
  do ik = 1, st%nik
    allocate(ev(norbs))
    call lalg_geneigensolve(norbs, lcao_data%hamilt(1:norbs, 1:norbs, ik), &
         lcao_data%s(1:norbs, 1:norbs, ik), ev)

    st%eigenval(1:nst, ik) = ev(1:nst)
    deallocate(ev)

    st%X(psi)(:,:,:, ik) = R_TOTYPE(M_ZERO)

    ! Change of base
    call X(lalg_gemm)(np*dim, nst, norbs, R_TOTYPE(M_ONE), lcao_data%psis(1:np, 1:dim, 1:norbs, ik), &
                      lcao_data%hamilt(1:norbs, 1:nst, ik), &
                      R_TOTYPE(M_ZERO),st%X(psi)(1:np, 1:dim, 1:nst, ik))

   end do

  deallocate(hpsi)
  call pop_sub()
end subroutine lcao_wf

end module lcao

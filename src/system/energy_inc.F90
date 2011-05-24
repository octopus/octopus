!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!!
!! $Id$


! ---------------------------------------------------------
! calculates the eigenvalues of the real orbitals
subroutine X(calculate_eigenvalues)(hm, der, st, time)
  type(hamiltonian_t), intent(inout) :: hm
  type(derivatives_t), intent(inout) :: der
  type(states_t),      intent(inout) :: st
  FLOAT,   optional,   intent(in)    :: time

  R_TYPE, allocatable :: hpsi(:, :, :)
  R_TYPE, allocatable :: eigen(:)
  integer :: ik, minst, maxst, ib
  type(batch_t) :: hpsib
  type(profile_t), save :: prof

  PUSH_SUB(X(calculate_eigenvalues))

  if(hm%theory_level == CLASSICAL) then
    st%eigenval = M_ZERO
    POP_SUB(X(calculate_eigenvalues))
    return
  end if

  call profiling_in(prof, "EIGENVALUE_CALC")

  if(in_debug_mode) then
    write(message(1), '(a)') 'Debug: Calculating eigenvalues.'
    call messages_info(1)
  end if

  ! FIXME: for TD open boundaries this is wrong. But the GS case like above
  ! is also wrong.
  ! The correct way to calculate the eigenvalue here is:
  !      / Psi_L | H_LL H_LC 0    | Psi_L \
  ! e = <  Psi_C | H_CL H_CC H_CR | Psi_C  >
  !      \ Psi_R | 0    H_RC H_RR | Psi_R /
  ! But I am not sure how to calculate this right now.

  SAFE_ALLOCATE(eigen(st%st_start:st%st_end))

  do ik = st%d%kpt%start, st%d%kpt%end
    do ib = st%block_start, st%block_end

      minst = st%block_range(ib, 1)
      maxst = st%block_range(ib, 2)

      call batch_init(hpsib, st%d%dim, st%block_size(ib))
      call X(batch_new)(hpsib, minst, maxst, der%mesh%np)
      call X(hamiltonian_apply_batch)(hm, der, st%psib(ib, ik), hpsib, ik, time)
      call X(mesh_batch_dotp_vector)(der%mesh, st%psib(ib, ik), hpsib, eigen(minst:maxst))
      call X(batch_delete)(hpsib)
      call batch_end(hpsib)

    end do
    
    st%eigenval(st%st_start:st%st_end, ik) = eigen(st%st_start:st%st_end)

  end do

  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(eigen)

  call profiling_out(prof)
  POP_SUB(X(calculate_eigenvalues))
end subroutine X(calculate_eigenvalues)

! ---------------------------------------------------------
FLOAT function X(electronic_kinetic_energy)(hm, gr, st) result(t0)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st

  integer :: ik, ist
  R_TYPE, allocatable :: tpsi(:, :), psi(:, :)
  FLOAT, allocatable  :: t(:, :)

  PUSH_SUB(X(electronic_kinetic_energy))

  SAFE_ALLOCATE(tpsi(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(t(st%st_start:st%st_end, 1:st%d%nik))

  t = M_ZERO

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      tpsi = R_TOTYPE(M_ZERO)
      call X(hamiltonian_apply)(hm, gr%der, st%X(psi)(:, :, ist, ik), tpsi, ist, ik, terms = TERM_KINETIC)
      t(ist, ik) = X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), tpsi)
    end do
  end do
  
  t0 = states_eigenvalues_sum(st, t)
  
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(tpsi)
  SAFE_DEALLOCATE_A(t)
  POP_SUB(X(electronic_kinetic_energy))
end function X(electronic_kinetic_energy)

! ---------------------------------------------------------
FLOAT function X(electronic_external_energy)(hm, gr, st) result(v)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st

  integer :: ik, ist
  R_TYPE, allocatable :: vpsi(:, :)
  FLOAT, allocatable :: t(:, :)

  PUSH_SUB(X(electronic_external_energy))

  SAFE_ALLOCATE(vpsi(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(t(st%st_start:st%st_end, 1:st%d%nik))
  t = M_ZERO

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      call X(hamiltonian_apply)(hm, gr%der, st%X(psi)(:, :, ist, ik), vpsi, ist, ik, &
        terms = TERM_NON_LOCAL_POTENTIAL + TERM_LOCAL_EXTERNAL)
      t(ist, ik) = X(mf_dotp) (gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), vpsi)
    end do
  end do

  v = states_eigenvalues_sum(st, t)

  SAFE_DEALLOCATE_A(vpsi)
  SAFE_DEALLOCATE_A(t)
  POP_SUB(X(electronic_external_energy))
end function X(electronic_external_energy)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

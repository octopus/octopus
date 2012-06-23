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
!> calculates the eigenvalues of the real orbitals
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

      minst = states_block_min(st, ib)
      maxst = states_block_max(st, ib)

      call batch_copy(st%psib(ib, ik), hpsib, reference = .false.)

      if(hamiltonian_apply_packed(hm, der%mesh)) then
        call batch_pack(st%psib(ib, ik))
        call batch_pack(hpsib, copy = .false.)
      end if

      call X(hamiltonian_apply_batch)(hm, der, st%psib(ib, ik), hpsib, ik, time)
      call X(mesh_batch_dotp_vector)(der%mesh, st%psib(ib, ik), hpsib, eigen(minst:maxst))

      if(hamiltonian_apply_packed(hm, der%mesh)) call batch_unpack(st%psib(ib, ik), copy = .false.)
      
      call batch_end(hpsib, copy = .false.)

    end do
    
    st%eigenval(st%st_start:st%st_end, ik) = eigen(st%st_start:st%st_end)

  end do

  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(eigen)

  call profiling_out(prof)
  POP_SUB(X(calculate_eigenvalues))
end subroutine X(calculate_eigenvalues)

! ---------------------------------------------------------
R_TYPE function X(energy_calc_electronic)(hm, der, st, terms, cproduct) result(energy)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(inout) :: der
  type(states_t),      intent(inout) :: st
  integer,             intent(in)    :: terms
  logical, optional,   intent(in)    :: cproduct

  integer :: ik, ist, ib, minst, maxst
  type(batch_t) :: hpsib
  R_TYPE, allocatable  :: tt(:, :)
  logical :: cproduct_
  
  cproduct_ = optional_default(cproduct, .false.)
 
  PUSH_SUB(X(energy_calc_electronic))

  SAFE_ALLOCATE(tt(st%st_start:st%st_end, 1:st%d%nik))

  tt = M_ZERO

  do ik = st%d%kpt%start, st%d%kpt%end
    do ib = st%block_start, st%block_end
      minst = states_block_min(st, ib)
      maxst = states_block_max(st, ib)

      call batch_copy(st%psib(ib, ik), hpsib, reference = .false.)

      call X(hamiltonian_apply_batch)(hm, der, st%psib(ib, ik), hpsib, ik, terms = terms)
      call X(mesh_batch_dotp_vector)(der%mesh, st%psib(ib, ik), hpsib, tt(minst:maxst, ik), cproduct = cproduct_)

      call batch_end(hpsib, copy = .false.)

    end do
  end do
  
#ifdef R_TCOMPLEX
  energy = zstates_eigenvalues_sum(st, tt)
#else  
  energy = states_eigenvalues_sum(st, real(tt, REAL_PRECISION))
#endif  
  
  SAFE_DEALLOCATE_A(tt)
  POP_SUB(X(energy_calc_electronic))
end function X(energy_calc_electronic)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$


! ---------------------------------------------------------
!> calculates the eigenvalues of the orbitals
subroutine X(calculate_eigenvalues)(hm, der, st, time)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(inout) :: der
  type(states_t),      intent(inout) :: st
  FLOAT,   optional,   intent(in)    :: time

  R_TYPE, allocatable :: eigen(:, :)
  integer :: ik
  logical :: cmplxscl

  PUSH_SUB(X(calculate_eigenvalues))
  
  cmplxscl = hm%cmplxscl%space

  if(hm%theory_level == CLASSICAL) then
    st%eigenval = M_ZERO
    POP_SUB(X(calculate_eigenvalues))
    return
  end if

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

  SAFE_ALLOCATE(eigen(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
  call X(calculate_expectation_values)(hm, der, st, eigen, time = time)

  st%eigenval(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end) = &
    real(eigen(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end), REAL_PRECISION)
#ifdef R_TCOMPLEX    
  if(cmplxscl) st%zeigenval%Im(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end) = &
    aimag(eigen(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
#endif

  SAFE_DEALLOCATE_A(eigen)

  POP_SUB(X(calculate_eigenvalues))
end subroutine X(calculate_eigenvalues)

subroutine X(calculate_expectation_values)(hm, der, st, eigen, time, terms)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(inout) :: der
  type(states_t),      intent(inout) :: st
  R_TYPE,              intent(out)   :: eigen(st%st_start:, st%d%kpt%start:) !< (:st%st_end, :st%d%kpt%end)
  FLOAT,   optional,   intent(in)    :: time
  integer, optional,   intent(in)    :: terms

  integer :: ik, minst, maxst, ib
  type(batch_t) :: hpsib
  type(profile_t), save :: prof
  logical :: cmplxscl

  PUSH_SUB(X(calculate_expectation_values))
  
  call profiling_in(prof, "EIGENVALUE_CALC")

  cmplxscl = hm%cmplxscl%space

  do ik = st%d%kpt%start, st%d%kpt%end
    do ib = st%group%block_start, st%group%block_end

      minst = states_block_min(st, ib)
      maxst = states_block_max(st, ib)

      call batch_copy(st%group%psib(ib, ik), hpsib, reference = .false.)

      if(hamiltonian_apply_packed(hm, der%mesh)) then
        call batch_pack(st%group%psib(ib, ik))
        if(st%have_left_states) call batch_pack(st%psibL(ib, ik))
        call batch_pack(hpsib, copy = .false.)
      end if

      call X(hamiltonian_apply_batch)(hm, der, st%group%psib(ib, ik), hpsib, ik, time = time, terms = terms)
      if(st%have_left_states) then
        call X(mesh_batch_dotp_vector)(der%mesh, st%psibL(ib, ik), hpsib, eigen(minst:maxst, ik), cproduct = cmplxscl)
      else
        call X(mesh_batch_dotp_vector)(der%mesh, st%group%psib(ib, ik), hpsib, eigen(minst:maxst, ik), cproduct = cmplxscl)        
      end if
      if(hamiltonian_apply_packed(hm, der%mesh)) then
        call batch_unpack(st%group%psib(ib, ik), copy = .false.)
        if(st%have_left_states) call batch_unpack(st%psibL(ib, ik), copy = .false.)
      end if

      call batch_end(hpsib, copy = .false.)

    end do
  end do

  call profiling_out(prof)
  POP_SUB(X(calculate_expectation_values))
end subroutine X(calculate_expectation_values)

! ---------------------------------------------------------
R_TYPE function X(energy_calc_electronic)(hm, der, st, terms) result(energy)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(inout) :: der
  type(states_t),      intent(inout) :: st
  integer,             intent(in)    :: terms

  integer :: ik, ib, minst, maxst
  R_TYPE, allocatable  :: tt(:, :)
 
  PUSH_SUB(X(energy_calc_electronic))

  SAFE_ALLOCATE(tt(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))

  call X(calculate_expectation_values)(hm, der, st, tt, terms = terms)

  if(hm%cmplxscl%space) then
#ifdef R_TCOMPLEX
    energy = zstates_eigenvalues_sum(st, tt)
#else
    message(1) = "Internal error in energy_calc_electronic, real states but complex scaling."
    call messages_fatal(1)
#endif
  else
    energy = states_eigenvalues_sum(st, real(tt, REAL_PRECISION))
  endif

  SAFE_DEALLOCATE_A(tt)
  POP_SUB(X(energy_calc_electronic))
end function X(energy_calc_electronic)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

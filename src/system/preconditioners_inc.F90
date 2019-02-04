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

! --------------------------------------------------------- 
subroutine X(preconditioner_apply)(pre, gr, hm, ik, a, b, omega)
  type(preconditioner_t), intent(in)    :: pre
  type(grid_t), target,   intent(in)    :: gr
  type(hamiltonian_t),    intent(in)    :: hm
  integer,                intent(in)    :: ik
  R_TYPE,                 intent(inout) :: a(:,:)
  R_TYPE,                 intent(inout) :: b(:,:)
  R_TYPE,       optional, intent(in)    :: omega
  
  integer :: idim
  R_TYPE  :: omega_
  type(profile_t), save :: preconditioner_prof
  
  call profiling_in(preconditioner_prof, "PRECONDITIONER")
  PUSH_SUB(X(preconditioner_apply))

  omega_ = M_ZERO
  if(present(omega)) omega_ = omega

  select case(pre%which)
  case(PRE_NONE)
    do idim = 1, hm%d%dim
      call lalg_copy(gr%mesh%np, a(:,idim), b(:,idim))
    end do

  case(PRE_FILTER)
    do idim = 1, hm%d%dim
      call X(derivatives_perform)(pre%op, gr%der, a(:, idim), b(:, idim))
    end do
  end select

  POP_SUB(X(preconditioner_apply))
  call profiling_out(preconditioner_prof)
end subroutine X(preconditioner_apply)

! ----------------------------------------

subroutine X(preconditioner_apply_batch)(pre, gr, hm, ik, aa, bb, omega)
  type(preconditioner_t), intent(in)    :: pre
  type(grid_t),           intent(in)    :: gr
  type(hamiltonian_t),    intent(in)    :: hm
  integer,                intent(in)    :: ik
  type(batch_t),          intent(inout) :: aa
  type(batch_t),          intent(inout) :: bb
  R_TYPE,       optional, intent(in)    :: omega(:)

  integer :: ii
  type(profile_t), save :: prof
  R_TYPE, allocatable :: psia(:, :), psib(:, :)

  PUSH_SUB(X(preconditioner_apply_batch))
  call profiling_in(prof, 'PRECONDITIONER_BATCH')

  if(pre%which == PRE_FILTER) then
    call X(derivatives_batch_perform)(pre%op, gr%der, aa, bb)
  else if(pre%which == PRE_NONE) then
    call batch_copy_data(gr%der%mesh%np, aa, bb)
  end if

  call profiling_out(prof)
  POP_SUB(X(preconditioner_apply_batch))
end subroutine X(preconditioner_apply_batch)
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

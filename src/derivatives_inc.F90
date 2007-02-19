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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

! This module calculates the derivatives (gradients, laplacians, etc.) 
! of a function. Note that the function whose derivative is to be calculated
! *has* to be defined (1:m%np_part), while the (1:m%np) values of the derivative
! are calculated. This was made to simplify the parallel mode, and has to be
! followed by the rest of the code.

! ---------------------------------------------------------
! The transpose of the Laplacian.
subroutine X(derivatives_laplt)(der, f, lapl, have_ghost_)
  type(der_discr_t), intent(in)    :: der
  R_TYPE,            intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,            intent(out)   :: lapl(:)  ! lapl(m%np)
  logical, optional, intent(in)    :: have_ghost_

  logical :: have_ghost

  call push_sub('derivatives_inc.Xderivatives_laplt')

  ASSERT(ubound(f,    DIM=1) == der%m%np_part)
  ASSERT(ubound(lapl, DIM=1) >= der%m%np)

  have_ghost = .false.
  if(present(have_ghost_)) have_ghost = have_ghost_

  if(.not.have_ghost.and.der%zero_bc) then
    call X(zero_bc)(der%m, f)
  end if

  call X(nl_operator_operate) (der%laplt, f, lapl)

  call pop_sub()
end subroutine X(derivatives_laplt)


! ---------------------------------------------------------
subroutine X(derivatives_lapl)(der, f, lapl, have_ghost_, ghost_update)
  type(der_discr_t), target, intent(in)    :: der
  R_TYPE, target,            intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,                    intent(out)   :: lapl(:)  ! lapl(m%np)
  logical, optional,         intent(in)    :: have_ghost_
  logical, optional,         intent(in)    :: ghost_update

  logical             :: have_ghost

  call push_sub('derivatives_inc.Xderivatives_lapl')

  ASSERT(ubound(f,    DIM=1) == der%m%np_part)
  ASSERT(ubound(lapl, DIM=1) >= der%m%np)

  have_ghost = .false.
  if(present(have_ghost_)) have_ghost = have_ghost_

  ! If the derivatives are defined with non-constant weights, then we do not
  ! need extra points.
  ! Otherwise, set the boundary points to zero (if we have zero boundary conditions).
  if(.not.have_ghost.and.der%zero_bc) then
    call X(zero_bc)(der%m, f)
  end if

  call X(nl_operator_operate) (der%lapl, f, lapl, ghost_update=ghost_update)

  call pop_sub()
end subroutine X(derivatives_lapl)


! ---------------------------------------------------------
subroutine X(derivatives_grad)(der, f, grad, ghost_update)
  type(der_discr_t), intent(in)    :: der
  R_TYPE, target,    intent(inout) :: f(:)       ! f(m%np_part)
  R_TYPE,            intent(out)   :: grad(:,:)  ! grad(m%np, m%sb%dim)
  logical, optional, intent(in)    :: ghost_update

  integer :: i

  call push_sub('derivatives_inc.Xderivatives_grad')
  
  ASSERT(ubound(f,    DIM=1) == der%m%np_part)
  ASSERT(ubound(grad, DIM=1) >= der%m%np)
  ASSERT(ubound(grad, DIM=2) >= der%m%sb%dim)

  if(der%zero_bc) then
    call X(zero_bc)(der%m, f)
  end if

  grad(:,:) = R_TOTYPE(M_ZERO)
  do i = 1, der%m%sb%dim
    call X(nl_operator_operate) (der%grad(i), f, grad(:,i), ghost_update=ghost_update)
  end do

  call pop_sub()
end subroutine X(derivatives_grad)


! ---------------------------------------------------------
subroutine X(derivatives_div)(der, f, div, ghost_update)
  type(der_discr_t), intent(in)    :: der
  R_TYPE,            intent(inout) :: f(:,:)   ! f(m%np_part, m%sb%dim)
  R_TYPE,            intent(out)   :: div(:)   ! div(m%np)
  logical, optional, intent(in)    :: ghost_update

  R_TYPE, allocatable :: tmp(:)
  integer             :: i

  call push_sub('derivatives_inc.Xderivatives_div')

  ASSERT(ubound(f,   DIM=1) == der%m%np_part)
  ASSERT(ubound(div, DIM=1) >= der%m%np)

  if(der%zero_bc) then
    do i = 1, der%m%sb%dim
      call X(zero_bc)(der%m, f(:, i))
    end do
  end if

  ALLOCATE(tmp(der%m%np), der%m%np)

  div(:) = R_TOTYPE(M_ZERO)
  do i = 1, der%m%sb%dim
    call X(nl_operator_operate) (der%grad(i), f(:,i), tmp, ghost_update=ghost_update)
    div(:) = div(:) + tmp(:)
  end do

  deallocate(tmp)

  call pop_sub()
end subroutine X(derivatives_div)


! ---------------------------------------------------------
subroutine X(derivatives_curl)(der, f, curl, ghost_update)
  type(der_discr_t), intent(in)    :: der
  R_TYPE,            intent(inout) :: f(:,:)    ! f(m%np_part, der%m%sb%dim)
  R_TYPE,            intent(out)   :: curl(:,:) ! curl(m%np, der%m%sb%dim)
  logical, optional, intent(in)    :: ghost_update

  R_TYPE, allocatable :: tmp(:)
  integer             :: i

  call push_sub('derivatives_inc.Xderivatives_div')

  ASSERT(der%m%sb%dim == 3)
  ASSERT(ubound(f,    DIM=1) == der%m%np_part)
  ASSERT(ubound(curl, DIM=1) >= der%m%np)
  ASSERT(ubound(curl, DIM=2) == der%m%sb%dim)

  if(der%zero_bc) then
    do i = 1, der%m%sb%dim
      call X(zero_bc)(der%m, f(:, i))
    end do
  end if

  ALLOCATE(tmp(der%m%np), der%m%np)

  curl(:,:) = R_TOTYPE(M_ZERO)

  call X(nl_operator_operate) (der%grad(3), f(:,1), tmp, ghost_update=ghost_update)
  curl(:,2) = curl(:,2) + tmp(:)
  call X(nl_operator_operate) (der%grad(2), f(:,1), tmp, ghost_update=ghost_update)
  curl(:,3) = curl(:,3) - tmp(:)

  call X(nl_operator_operate) (der%grad(3), f(:,2), tmp, ghost_update=ghost_update)
  curl(:,1) = curl(:,1) - tmp(:)
  call X(nl_operator_operate) (der%grad(1), f(:,2), tmp, ghost_update=ghost_update)
  curl(:,3) = curl(:,3) + tmp(:)

  call X(nl_operator_operate) (der%grad(2), f(:,3), tmp, ghost_update=ghost_update)
  curl(:,1) = curl(:,1) + tmp(:)
  call X(nl_operator_operate) (der%grad(1), f(:,3), tmp, ghost_update=ghost_update)
  curl(:,2) = curl(:,2) - tmp(:)

  deallocate(tmp)

  call pop_sub()
end subroutine X(derivatives_curl)


! ---------------------------------------------------------
! Set all boundary points in f to zero to implement zero
! boundary conditions for the derivatives.
subroutine X(zero_bc)(m, f)
  type(mesh_t), intent(in)    :: m
  R_TYPE,       intent(inout) :: f(:)

  integer :: bndry_start, bndry_end
  integer :: p

  call push_sub('derivatives_inc.Xzero_bc')

  p = m%vp%partno

  ! The boundary points are at different locations depending on the presence
  ! of ghost points due to domain parallelization.
  if(m%parallel_in_domains) then
    bndry_start = m%vp%np_local(p) + m%vp%np_ghost(p) + 1
    bndry_end   = m%vp%np_local(p) + m%vp%np_ghost(p) + m%vp%np_bndry(p)
  else
    bndry_start = m%np+1
    bndry_end   = m%np_part
  end if
  f(bndry_start:bndry_end) = R_TOTYPE(M_ZERO)

  call pop_sub()
end subroutine X(zero_bc)

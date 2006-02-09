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
!!
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
  R_TYPE,               intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,               intent(out)   :: lapl(:)  ! lapl(m%np)
  logical, optional,    intent(in)    :: have_ghost_

  logical :: have_ghost

  ASSERT(ubound(f,    DIM=1) == der%m%np_part)
  ASSERT(ubound(lapl, DIM=1) >= der%m%np)

  have_ghost = .false.
  if(present(have_ghost_)) have_ghost = have_ghost_

  if(.not.have_ghost.and.der%zero_bc) then
    f(der%m%np+1:der%m%np_part) = R_TOTYPE(M_ZERO)
  end if

  call X(nl_operator_operate) (der%laplt, f, lapl)

end subroutine X(derivatives_laplt)


! ---------------------------------------------------------
subroutine X(derivatives_lapl)(der, f, lapl, have_ghost_)
  type(der_discr_t), intent(in)    :: der
  R_TYPE, target,       intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,               intent(out)   :: lapl(:)  ! lapl(m%np)
  logical, optional,    intent(in)    :: have_ghost_

  logical :: have_ghost

  call push_sub('derivatives_inc.Xderivatives_lapl')

  ASSERT(ubound(f,    DIM=1) == der%m%np_part)
  ASSERT(ubound(lapl, DIM=1) >= der%m%np)

  have_ghost = .false.
  if(present(have_ghost_)) have_ghost = have_ghost_

  ! If the derivatives are defined with non-constant weights, then we do not
  ! need extra points.
  if( (.not.(have_ghost)).and.der%zero_bc) then
    f(der%m%np+1:der%m%np_part) = R_TOTYPE(M_ZERO)
  end if

  call X(nl_operator_operate) (der%lapl, f, lapl)

  call pop_sub()

end subroutine X(derivatives_lapl)


! ---------------------------------------------------------
subroutine X(derivatives_grad)(der, f, grad)
  type(der_discr_t), intent(in)    :: der
  R_TYPE, target,       intent(inout) :: f(:)       ! f(m%np_part)
  R_TYPE,               intent(out)   :: grad(:,:)  ! grad(m%np, calc_dim)

  integer :: i

  ASSERT(ubound(f,    DIM=1) == der%m%np_part)
  ASSERT(ubound(grad, DIM=1) >= der%m%np)
  ASSERT(ubound(grad, DIM=2) >= calc_dim)

  if(der%zero_bc) then
    f(der%m%np+1:der%m%np_part) = R_TOTYPE(M_ZERO)
  end if

  grad(:,:) = R_TOTYPE(M_ZERO)
  do i = 1, calc_dim
    call X(nl_operator_operate) (der%grad(i), f, grad(:,i))
  end do

end subroutine X(derivatives_grad)


! ---------------------------------------------------------
subroutine X(derivatives_div)(der, f, div)
  type(der_discr_t), intent(in)    :: der
  R_TYPE,               intent(inout) :: f(:,:)   ! f(m%np_part, calc_dim)
  R_TYPE,               intent(out)   :: div(:)   ! div(m%np)

  R_TYPE, allocatable :: tmp(:)
  integer :: i

  ASSERT(ubound(f,   DIM=1) == der%m%np_part)
  ASSERT(ubound(div, DIM=1) >= der%m%np)

  if(der%zero_bc) then
    f(der%m%np+1:der%m%np_part,:) = R_TOTYPE(M_ZERO)
  end if
  ALLOCATE(tmp(der%m%np), der%m%np)

  div(:) = R_TOTYPE(M_ZERO)
  do i = 1, calc_dim
    call X(nl_operator_operate) (der%grad(i), f(:,i), tmp)
     div(:) = div(:) + tmp(:)
  end do

  deallocate(tmp)
end subroutine X(derivatives_div)


! ---------------------------------------------------------
subroutine X(derivatives_curl)(der, f, curl)
  type(der_discr_t), intent(in)    :: der
  R_TYPE,               intent(inout) :: f(:,:)    ! f(m%np_part, calc_dim)
  R_TYPE,               intent(out)   :: curl(:,:) ! curl(m%np, calc_dim)

  R_TYPE, allocatable :: tmp(:)

  ASSERT(calc_dim == 3)
  ASSERT(ubound(f,    DIM=1) == der%m%np_part)
  ASSERT(ubound(curl, DIM=1) >= der%m%np)
  ASSERT(ubound(curl, DIM=2) == calc_dim)

  if(der%zero_bc) then
    f(der%m%np+1:der%m%np_part,:) = R_TOTYPE(M_ZERO)
  end if
  ALLOCATE(tmp(der%m%np), der%m%np)

  curl(:,:) = R_TOTYPE(M_ZERO)

  call X(nl_operator_operate) (der%grad(3), f(:,1), tmp)
  curl(:,2) = curl(:,2) + tmp(:)
  call X(nl_operator_operate) (der%grad(2), f(:,1), tmp)
  curl(:,3) = curl(:,3) - tmp(:)

  call X(nl_operator_operate) (der%grad(3), f(:,2), tmp)
  curl(:,1) = curl(:,1) - tmp(:)
  call X(nl_operator_operate) (der%grad(1), f(:,2), tmp)
  curl(:,3) = curl(:,3) + tmp(:)

  call X(nl_operator_operate) (der%grad(2), f(:,3), tmp)
  curl(:,1) = curl(:,1) + tmp(:)
  call X(nl_operator_operate) (der%grad(1), f(:,3), tmp)
  curl(:,2) = curl(:,2) - tmp(:)

  deallocate(tmp)

end subroutine X(derivatives_curl)

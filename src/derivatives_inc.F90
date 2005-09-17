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

! ---------------------------------------------------------
! The transpose of the Laplacian.
! WARNING: We should avoid the unnecessary use of zero-ed ghost points,
! also for the case of the Laplacian itself.
subroutine X(derivatives_laplt)(der, f, lapl, have_ghost_)
  type(der_discr_type), intent(in)  :: der
  R_TYPE, target,       intent(in)  :: f(:)     ! f(m%np)
  R_TYPE,               intent(out) :: lapl(:)  ! lapl(m%np)
  logical, optional,    intent(in)  :: have_ghost_

  R_TYPE, pointer :: fp(:)
  logical :: have_ghost

  have_ghost = .false.
  if(present(have_ghost_)) have_ghost = have_ghost_

  if(.not.have_ghost.and.der%zero_bc) then
     allocate(fp(der%m%np_part))
     fp(1:der%m%np)              = f(1:der%m%np)
     fp(der%m%np+1:der%m%np_part) = R_TOTYPE(M_ZERO)
  else
     fp => f
  end if

  call X(nl_operator_operate) (der%laplt, fp, lapl)

  if(.not.have_ghost.and.der%zero_bc) deallocate(fp)
end subroutine X(derivatives_laplt)

! ---------------------------------------------------------
subroutine X(derivatives_lapl)(der, f, lapl, have_ghost_)
  type(der_discr_type), intent(in)    :: der
  R_TYPE, target,       intent(inout) :: f(:)     ! f(m%np)
  R_TYPE,               intent(out)   :: lapl(:)  ! lapl(m%np)
  logical, optional,    intent(in)    :: have_ghost_

  R_TYPE, pointer:: fp(:)
  logical :: have_ghost

  have_ghost = .false.
  if(present(have_ghost_)) have_ghost = have_ghost_

  ! If the derivatives are defined with non-constant weights, then we do not
  ! need extra points.
  if( (.not.(have_ghost))   .and.   der%zero_bc    .and.  der%lapl%const_w  ) then
     allocate(fp(der%m%np_part))
     fp(1:der%m%np)              = f(1:der%m%np)
     fp(der%m%np+1:der%m%np_part) = R_TOTYPE(M_ZERO)
     call X(nl_operator_operate) (der%lapl, fp, lapl)
     deallocate(fp)
  else
     call X(nl_operator_operate) (der%lapl, f, lapl)
  endif

end subroutine X(derivatives_lapl)



! ---------------------------------------------------------
subroutine X(derivatives_grad)(sb, der, f, grad)
  type(simul_box_type), intent(in)  :: sb
  type(der_discr_type), intent(in)  :: der
  R_TYPE, target,       intent(in)  :: f(:)       ! f(m%np)
  R_TYPE,               intent(out) :: grad(:,:)  ! grad(m%np, sb%dim)

  R_TYPE, pointer :: fp(:)
  integer :: i

  if(der%zero_bc) then
     allocate(fp(der%m%np_part))
     fp(1:der%m%np)              = f(1:der%m%np)
     fp(der%m%np+1:der%m%np_part) = R_TOTYPE(M_ZERO)
  else
     fp => f
  end if

  do i = 1, sb%dim
     call X(nl_operator_operate) (der%grad(i), fp, grad(:,i))
  end do

  if(der%zero_bc) deallocate(fp)
end subroutine X(derivatives_grad)


! ---------------------------------------------------------
subroutine X(derivatives_div)(sb, der, f, div)
  type(simul_box_type), intent(in)  :: sb
  type(der_discr_type), intent(in)  :: der
  R_TYPE, target,       intent(in)  :: f(:,:)   ! f(m%np, sb%dim)
  R_TYPE,               intent(out) :: div(:)   ! div(m%np)

  R_TYPE, pointer     :: fp(:)
  R_TYPE, allocatable :: tmp(:)
  integer :: i

  if(der%zero_bc) then
     allocate(fp(der%m%np_part))
     fp(der%m%np+1:der%m%np_part) = R_TOTYPE(M_ZERO)
  end if
  allocate(tmp(der%m%np))

  div(:) = R_TOTYPE(M_ZERO)
  do i = 1, sb%dim
     if(der%zero_bc) then
        fp(1:der%m%np) = f(1:der%m%np,i)
     else
        fp => f(:,i)
     end if

     call X(nl_operator_operate) (der%grad(i), fp, tmp)
     div(:) = div(:) + tmp(:)
  end do

  deallocate(tmp)
  if(der%zero_bc) deallocate(fp)
end subroutine X(derivatives_div)


! ---------------------------------------------------------
subroutine X(derivatives_curl)(sb, der, f, curl)
  type(simul_box_type), intent(in)  :: sb
  type(der_discr_type), intent(in)  :: der
  R_TYPE, target,       intent(in)  :: f(:,:)    ! f(m%np, sb%dim)
  R_TYPE,               intent(out) :: curl(:,:) ! curl(m%np, sb%dim)

  R_TYPE, pointer     :: fp(:)
  R_TYPE, allocatable :: tmp(:)

  ASSERT(sb%dim == 3)

  if(der%zero_bc) then
     allocate(fp(der%m%np_part))
     fp(der%m%np+1:der%m%np_part) = R_TOTYPE(M_ZERO)
  end if
  allocate(tmp(der%m%np))

  curl(:,:) = R_TOTYPE(M_ZERO)

  call get_f(1)
  call X(nl_operator_operate) (der%grad(3), fp, tmp)
  curl(:,2) = curl(:,2) + tmp(:)
  call X(nl_operator_operate) (der%grad(2), fp, tmp)
  curl(:,3) = curl(:,3) - tmp(:)

  call get_f(2)
  call X(nl_operator_operate) (der%grad(3), fp, tmp)
  curl(:,1) = curl(:,1) - tmp(:)
  call X(nl_operator_operate) (der%grad(1), fp, tmp)
  curl(:,3) = curl(:,3) + tmp(:)

  call get_f(3)
  call X(nl_operator_operate) (der%grad(2), fp, tmp)
  curl(:,1) = curl(:,1) + tmp(:)
  call X(nl_operator_operate) (der%grad(1), fp, tmp)
  curl(:,2) = curl(:,2) - tmp(:)

  deallocate(tmp)
  if(der%zero_bc) deallocate(fp)

contains
  subroutine get_f(ii)
    integer :: ii
    if(der%zero_bc) then
       fp(1:der%m%np) = f(1:der%m%np,ii)
    else
       fp => f(:,ii)
    end if
  end subroutine get_f
end subroutine X(derivatives_curl)

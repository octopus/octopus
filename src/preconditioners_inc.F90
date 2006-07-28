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
subroutine X(apply_precond_smoothing)(this, a, b)
  type(preconditioner_smoothing_t), intent(inout) :: this
  R_TYPE, intent(inout) :: a(:)
  R_TYPE, intent(inout) :: b(:)
  
  call X(nl_operator_operate) (this%op, a(:), b(:))
  
end subroutine X(apply_precond_smoothing)


! --------------------------------------------------------- 
subroutine X(apply_precond_smoothing_wfs)(this, a, b, n)
  type(preconditioner_smoothing_t), intent(inout) :: this
  R_TYPE, intent(inout) :: a(:,:)
  R_TYPE, intent(inout) :: b(:,:)
  integer, intent(in)   :: n

  integer :: i

  do i=1, n
    call X(nl_operator_operate) (this%op, a(:,i), b(:,i))
  end  do
  
end subroutine X(apply_precond_smoothing_wfs)

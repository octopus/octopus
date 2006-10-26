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

! --------------------------------------------------------- 
subroutine X(preconditioner_apply)(this, a, b)
  type(preconditioner_t), intent(in)    :: this
  R_TYPE,                 intent(inout) :: a(:)
  R_TYPE,                 intent(out)   :: b(:)
  
  select case(this%which)
  case(PRECONDITIONER_NONE)
    b(:) = a(:)
  case(PRECONDITIONER_SMOOTHING)
    call X(nl_operator_operate) (this%op, a(:), b(:))
  end select
  
end subroutine X(preconditioner_apply)

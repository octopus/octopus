!! Copyright (C) 2009 X. Andrade
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
!! $Id: symm_op_inc.F90 3479 2007-11-09 08:36:10Z xavier $

! -------------------------------------------------------------------------------
pure function X(symm_op_apply)(this, aa) result(bb)
  type(symm_op_t),  intent(in)  :: this
  R_TYPE,            intent(in)  :: aa(:)
  R_TYPE                         :: bb(1:3)
  
  bb(1:3) = matmul(aa(1:3), dble(this%rotation(1:3, 1:3))) + this%translation(1:3)
  
end function X(symm_op_apply)

! -------------------------------------------------------------------------------
pure function X(symm_op_apply_inv)(this, aa) result(bb)
  type(symm_op_t),  intent(in)  :: this
  R_TYPE,            intent(in)  :: aa(:)
  R_TYPE                         :: bb(1:3)
  
  bb(1:3) = matmul(dble(this%rotation(1:3, 1:3)), aa(1:3))
  
end function X(symm_op_apply_inv)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

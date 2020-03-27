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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

! -------------------------------------------------------------------------------
pure function X(symm_op_apply_red)(this, aa) result(bb)
  type(symm_op_t),  intent(in)  :: this
  R_TYPE,           intent(in)  :: aa(:) !< (this%dim)
  R_TYPE                        :: bb(1:this%dim)
  
  bb(1:this%dim) = matmul(TOFLOAT(this%rot_red(1:this%dim, 1:this%dim)),&
                       aa(1:this%dim)) + this%trans_red(1:this%dim)
 
end function X(symm_op_apply_red)

! -------------------------------------------------------------------------------
pure function X(symm_op_apply_inv_red)(this, aa) result(bb)
  type(symm_op_t),  intent(in)  :: this
  R_TYPE,           intent(in)  :: aa(:) !< (this%dim)
  R_TYPE                        :: bb(1:this%dim)
  
  bb(1:this%dim) = aa(1:this%dim) - this%trans_red(1:this%dim)
  bb(1:this%dim) = matmul(TOFLOAT(this%rot_red_inv(1:this%dim, 1:this%dim)), &
                        bb(1:this%dim))

end function X(symm_op_apply_inv_red)

! -------------------------------------------------------------------------------
pure function X(symm_op_apply_transpose_red)(this, aa) result(bb)
  type(symm_op_t),  intent(in)  :: this
  R_TYPE,           intent(in)  :: aa(:) !< (this%dim)
  R_TYPE                        :: bb(1:this%dim)

  bb(1:this%dim) = matmul(aa(1:this%dim), TOFLOAT(this%rot_red(1:this%dim, 1:this%dim)))

end function X(symm_op_apply_transpose_red)


! -------------------------------------------------------------------------------
pure function X(symm_op_apply_cart)(this, aa) result(bb)
  type(symm_op_t),  intent(in)  :: this
  R_TYPE,           intent(in)  :: aa(:) !< (this%dim)
  R_TYPE                        :: bb(1:this%dim)

  bb(1:this%dim) = matmul(this%rot_cart(1:this%dim, 1:this%dim), aa(1:this%dim)) &
                            + this%trans_cart(1:this%dim)

end function X(symm_op_apply_cart)

! -------------------------------------------------------------------------------
pure function X(symm_op_apply_inv_cart)(this, aa) result(bb)
  type(symm_op_t),  intent(in)  :: this
  R_TYPE,           intent(in)  :: aa(:) !< (this%dim)
  R_TYPE                        :: bb(1:this%dim)

  bb(1:this%dim) = aa(1:this%dim) - this%trans_cart(1:this%dim)
  bb(1:this%dim) = matmul(bb(1:this%dim), this%rot_cart(1:this%dim, 1:this%dim))

end function X(symm_op_apply_inv_cart)

! -------------------------------------------------------------------------------
logical pure function X(symm_op_invariant_cart)(this, aa, prec) result(invariant)
  type(symm_op_t),  intent(in)  :: this
  R_TYPE,           intent(in)  :: aa(:)
  FLOAT,            intent(in)  :: prec

  R_TYPE :: cc(1:3)

  cc = symm_op_apply_cart(this, aa)
  invariant = all(abs(cc(1:this%dim) - aa(1:this%dim)) < prec)

end function X(symm_op_invariant_cart)

! -------------------------------------------------------------------------------
logical pure function X(symm_op_invariant_red)(this, aa, prec) result(invariant)
  type(symm_op_t),  intent(in)  :: this
  R_TYPE,           intent(in)  :: aa(:)
  FLOAT,            intent(in)  :: prec

  R_TYPE :: cc(1:3)

  cc = symm_op_apply_red(this, aa)
  invariant = all(abs(cc(1:this%dim) - aa(1:this%dim)) < prec)

end function X(symm_op_invariant_red)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

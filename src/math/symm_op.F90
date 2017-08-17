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

#include "global.h"

module symm_op_oct_m
  use global_oct_m
  use messages_oct_m

  implicit none

  private
  
  public ::                            &
       symm_op_t,                      &
       symm_op_init,                   &
       symm_op_copy,                   &
       symm_op_end,                    &
       symm_op_apply_red,              &
       symm_op_apply_inv_red,          &
       symm_op_invariant_red,          &
       symm_op_has_translation,        &
       symm_op_rotation_matrix_red,    &
       symm_op_translation_vector_red, &
       symm_op_apply_cart,             &
       symm_op_apply_inv_cart,         &
       symm_op_invariant_cart,         &
       symm_op_rotation_matrix_cart,   &
       symm_op_translation_vector_cart,&
       symm_op_is_identity

  type symm_op_t
    private
    integer  :: rot_red(1:3, 1:3)
    integer  :: rot_red_inv(1:3, 1:3)
    FLOAT    :: rot_cart(1:3, 1:3)
    FLOAT    :: trans_red(1:3)
    FLOAT    :: trans_cart(1:3)
  end type symm_op_t
  
  interface symm_op_apply_red
    module procedure isymm_op_apply_red, dsymm_op_apply_red, zsymm_op_apply_red
  end interface symm_op_apply_red
 
  interface symm_op_apply_cart
    module procedure dsymm_op_apply_cart, zsymm_op_apply_cart
  end interface symm_op_apply_cart

  interface symm_op_apply_inv_red
    module procedure isymm_op_apply_inv_red, dsymm_op_apply_inv_red, zsymm_op_apply_inv_red
  end interface symm_op_apply_inv_red

  interface symm_op_apply_inv_cart
    module procedure dsymm_op_apply_inv_cart, zsymm_op_apply_inv_cart
  end interface symm_op_apply_inv_cart

contains
 
  !TODO: We should handle the low-dimensional case 
  subroutine symm_op_init(this, rot, rlattice, klattice, trans)
    type(symm_op_t),     intent(out) :: this
    integer,             intent(in)  :: rot(:, :)
    FLOAT,               intent(in)  :: rlattice(:, :), klattice(:, :)
    FLOAT, optional,     intent(in)  :: trans(:)

    PUSH_SUB(symm_op_init)

    !What comes out of spglib is the rotation matrix in fractional coordinates
    !This is not a rotation matrix!
    this%rot_red(1:3, 1:3) = rot(1:3, 1:3)
    !The inverse operation is only given by its inverse, not the transpose  
    !By Contrustion the reduced matrix also has a determinant of +1 
    ! Calculate the inverse of the matrix
    this%rot_red_inv(1,1) = +(rot(2,2)*rot(3,3) - rot(2,3)*rot(3,2))
    this%rot_red_inv(2,1) = -(rot(2,1)*rot(3,3) - rot(2,3)*rot(3,1))
    this%rot_red_inv(3,1) = +(rot(2,1)*rot(3,2) - rot(2,2)*rot(3,1))
    this%rot_red_inv(1,2) = -(rot(1,2)*rot(3,3) - rot(1,3)*rot(3,2))
    this%rot_red_inv(2,2) = +(rot(1,1)*rot(3,3) - rot(1,3)*rot(3,1))
    this%rot_red_inv(3,2) = -(rot(1,1)*rot(3,2) - rot(1,2)*rot(3,1))
    this%rot_red_inv(1,3) = +(rot(1,2)*rot(2,3) - rot(1,3)*rot(2,2))
    this%rot_red_inv(2,3) = -(rot(1,1)*rot(2,3) - rot(1,3)*rot(2,1))
    this%rot_red_inv(3,3) = +(rot(1,1)*rot(2,2) - rot(1,2)*rot(2,1))

    !This is a proper rotation matrix
    !R_{cart} = A*R_{red}*A^{-1}
    !Where A is the matrix containing the lattice vectors as column
    this%rot_cart(1:3, 1:3)  = matmul(rlattice(1:3, 1:3), &
                    matmul(rot(1:3, 1:3),transpose(klattice(1:3, 1:3))/ (M_TWO * M_PI)))

    if(present(trans)) then
      this%trans_red(1:3) = trans(1:3)
      this%trans_cart(1:3) = matmul(rlattice,trans(1:3))
    else
      this%trans_red(1:3) = M_ZERO
      this%trans_cart(1:3) = M_ZERO
    end if

    if(sum(abs(matmul(this%rot_cart,transpose(this%rot_cart)) &
              -reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3, 3/))))>CNST(1.0e-6)) then
      message(1) = "Internal error: This matrix is not a rotation matrix"
      write(message(2),'(3(3f7.3,2x))') this%rot_cart
      call messages_fatal(2) 
    end if
    POP_SUB(symm_op_init)
  end subroutine symm_op_init
  
  ! -------------------------------------------------------------------------------
  subroutine symm_op_copy(inp, outp)
    type(symm_op_t),     intent(in) :: inp
    type(symm_op_t),     intent(out) :: outp

    PUSH_SUB(symm_op_copy)

    outp%rot_red(1:3, 1:3) =  inp%rot_red(1:3, 1:3)
    outp%rot_red_inv(1:3, 1:3) =  inp%rot_red_inv(1:3, 1:3)
    outp%trans_red(1:3) =  inp%trans_red(1:3)
    outp%rot_cart(1:3, 1:3) =  inp%rot_cart(1:3, 1:3)
    outp%trans_cart(1:3) =  inp%trans_cart(1:3)

    POP_SUB(symm_op_copy)
  end subroutine symm_op_copy
  
  ! -------------------------------------------------------------------------------
  subroutine symm_op_end(this)
    type(symm_op_t),  intent(inout) :: this

    PUSH_SUB(symm_op_end)

    !nothing to do for the moment

    POP_SUB(symm_op_end)
  end subroutine symm_op_end
  
  ! -------------------------------------------------------------------------------
  logical pure function symm_op_invariant_cart(this, aa, prec) result(invariant)
    type(symm_op_t),  intent(in)  :: this
    FLOAT,            intent(in)  :: aa(:)
    FLOAT,            intent(in)  :: prec

    FLOAT :: cc(1:3)

    cc = symm_op_apply_cart(this, aa)
    invariant = all(abs(cc(1:3) - aa(1:3)) < prec)

  end function symm_op_invariant_cart

  ! -------------------------------------------------------------------------------
 logical pure function symm_op_invariant_red(this, aa, prec) result(invariant)
    type(symm_op_t),  intent(in)  :: this
    FLOAT,            intent(in)  :: aa(:)
    FLOAT,            intent(in)  :: prec

    FLOAT :: cc(1:3)

    cc = symm_op_apply_red(this, aa)
    invariant = all(abs(cc(1:3) - aa(1:3)) < prec)

  end function symm_op_invariant_red

  ! -------------------------------------------------------------------------------
  logical pure function symm_op_has_translation(this, prec) result(has)
    type(symm_op_t),  intent(in)  :: this
    FLOAT,            intent(in)  :: prec

    has = any(abs(this%trans_red(1:3)) >= prec)

  end function symm_op_has_translation

  ! -------------------------------------------------------------------------------

  function symm_op_rotation_matrix_red(this) result(matrix)
    type(symm_op_t),  intent(in)  :: this
    integer                       :: matrix(1:3, 1:3)

    matrix(1:3, 1:3) = this%rot_red(1:3, 1:3)

  end function symm_op_rotation_matrix_red

  ! -------------------------------------------------------------------------------

  function symm_op_rotation_matrix_cart(this) result(matrix)
    type(symm_op_t),  intent(in)  :: this
    integer                       :: matrix(1:3, 1:3)

    matrix(1:3, 1:3) = this%rot_cart(1:3, 1:3)

  end function symm_op_rotation_matrix_cart


  ! -------------------------------------------------------------------------------
 
  function symm_op_translation_vector_red(this) result(vector)
    type(symm_op_t),  intent(in)  :: this
    FLOAT                         :: vector(1:3)

    vector(1:3) = this%trans_red(1:3)

  end function symm_op_translation_vector_red

  ! -------------------------------------------------------------------------------

  function symm_op_translation_vector_cart(this) result(vector)
    type(symm_op_t),  intent(in)  :: this
    FLOAT                         :: vector(1:3)

    vector(1:3) = this%trans_cart(1:3)

  end function symm_op_translation_vector_cart

  ! -------------------------------------------------------------------------------

  logical pure function symm_op_is_identity(this) result(is_identity)
    type(symm_op_t),  intent(in)  :: this

    is_identity = .true.
    is_identity = is_identity .and. all(abs(this%trans_red) < CNST(1.0e-5))
    is_identity = is_identity .and. all(abs(this%rot_red(:, 1) - (/ CNST(1.0), CNST(0.0), CNST(0.0)/)) < CNST(1.0e-5))
    is_identity = is_identity .and. all(abs(this%rot_red(:, 2) - (/ CNST(0.0), CNST(1.0), CNST(0.0)/)) < CNST(1.0e-5))
    is_identity = is_identity .and. all(abs(this%rot_red(:, 3) - (/ CNST(0.0), CNST(0.0), CNST(1.0)/)) < CNST(1.0e-5))

  end function symm_op_is_identity

  ! -------------------------------------------------------------------------------
  
  pure function isymm_op_apply_red(this, aa) result(bb)
    type(symm_op_t),  intent(in)  :: this
    integer,          intent(in)  :: aa(:) !< (3)
    integer                       :: bb(1:3)

    bb(1:3) = nint(dsymm_op_apply_red(this, real(aa, REAL_PRECISION)))
    
  end function isymm_op_apply_red

  ! -------------------------------------------------------------------------------

  function isymm_op_apply_inv_red(this, aa) result(bb)
    type(symm_op_t),  intent(in)  :: this
    integer,          intent(in)  :: aa(:) !< (3)
    integer                       :: bb(1:3)

    bb(1:3) = nint(dsymm_op_apply_inv_red(this, real(aa, REAL_PRECISION)))
    
  end function isymm_op_apply_inv_red

#include "undef.F90"
#include "real.F90"
#include "symm_op_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "symm_op_inc.F90"

end module symm_op_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

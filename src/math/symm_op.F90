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
  use lattice_vectors_oct_m
  use messages_oct_m

  implicit none

  private
  
  public ::                            &
       symm_op_t,                      &
       symm_op_init,                   &
       symm_op_copy,                   &
       symm_op_apply_red,              &
       symm_op_apply_inv_red,          &
       symm_op_apply_transpose_red,    &
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
    integer  :: dim
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

  interface symm_op_apply_transpose_red
    module procedure dsymm_op_apply_transpose_red, zsymm_op_apply_transpose_red
  end interface symm_op_apply_transpose_red

  interface symm_op_apply_inv_cart
    module procedure dsymm_op_apply_inv_cart, zsymm_op_apply_inv_cart
  end interface symm_op_apply_inv_cart

  interface symm_op_invariant_cart
    module procedure dsymm_op_invariant_cart, zsymm_op_invariant_cart
  end interface symm_op_invariant_cart

  interface symm_op_invariant_red
    module procedure dsymm_op_invariant_red, zsymm_op_invariant_red
  end interface symm_op_invariant_red

contains
 
  subroutine symm_op_init(this, rot, latt, dim, trans)
    type(symm_op_t),         intent(out) :: this
    integer,                 intent(in)  :: rot(:, :)
    type(lattice_vectors_t), intent(in)  :: latt
    integer,                 intent(in)  :: dim
    FLOAT, optional,         intent(in)  :: trans(:)

    integer :: idim

    PUSH_SUB(symm_op_init)

    this%dim = dim

    !What comes out of spglib is the rotation matrix in fractional coordinates
    !This is not a rotation matrix!
    this%rot_red = 0
    this%rot_red(1:dim, 1:dim) = rot(1:dim, 1:dim)
    do idim = dim+1,3
      this%rot_red(idim,idim) = 1
    end do

    !The inverse operation is only given by its inverse, not the transpose  
    !By Contrustion the reduced matrix also has a determinant of +1 
    ! Calculate the inverse of the matrix
    this%rot_red_inv(1,1) = +(this%rot_red(2,2)*this%rot_red(3,3) - this%rot_red(2,3)*this%rot_red(3,2))
    this%rot_red_inv(2,1) = -(this%rot_red(2,1)*this%rot_red(3,3) - this%rot_red(2,3)*this%rot_red(3,1))
    this%rot_red_inv(3,1) = +(this%rot_red(2,1)*this%rot_red(3,2) - this%rot_red(2,2)*this%rot_red(3,1))
    this%rot_red_inv(1,2) = -(this%rot_red(1,2)*this%rot_red(3,3) - this%rot_red(1,3)*this%rot_red(3,2))
    this%rot_red_inv(2,2) = +(this%rot_red(1,1)*this%rot_red(3,3) - this%rot_red(1,3)*this%rot_red(3,1))
    this%rot_red_inv(3,2) = -(this%rot_red(1,1)*this%rot_red(3,2) - this%rot_red(1,2)*this%rot_red(3,1))
    this%rot_red_inv(1,3) = +(this%rot_red(1,2)*this%rot_red(2,3) - this%rot_red(1,3)*this%rot_red(2,2))
    this%rot_red_inv(2,3) = -(this%rot_red(1,1)*this%rot_red(2,3) - this%rot_red(1,3)*this%rot_red(2,1))
    this%rot_red_inv(3,3) = +(this%rot_red(1,1)*this%rot_red(2,2) - this%rot_red(1,2)*this%rot_red(2,1))

    !This is a proper rotation matrix
    !R_{cart} = A*R_{red}*A^{-1}
    !Where A is the matrix containing the lattice vectors as column
    this%rot_cart = M_ZERO
    this%rot_cart(1:dim, 1:dim)  = matmul(latt%rlattice(1:dim, 1:dim), &
                    matmul(rot(1:dim, 1:dim), latt%rlattice_inverse(1:dim, 1:dim)))
    do idim = dim+1,3
      this%rot_cart(idim,idim) = M_ONE
    end do

    this%trans_red(1:3) = M_ZERO
    this%trans_cart(1:3) = M_ZERO
    if(present(trans)) then
      this%trans_red(1:dim) = trans(1:dim)
      this%trans_cart(1:dim) = latt%red_to_cart(trans(1:dim))
    end if

    !We check that the rotation matrix in cartesian space is indeed a rotation matrix
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
 
    outp%dim = inp%dim
    outp%rot_red(1:3, 1:3) =  inp%rot_red(1:3, 1:3)
    outp%rot_red_inv(1:3, 1:3) =  inp%rot_red_inv(1:3, 1:3)
    outp%trans_red(1:3) =  inp%trans_red(1:3)
    outp%rot_cart(1:3, 1:3) =  inp%rot_cart(1:3, 1:3)
    outp%trans_cart(1:3) =  inp%trans_cart(1:3)

    POP_SUB(symm_op_copy)
  end subroutine symm_op_copy

  ! -------------------------------------------------------------------------------
  logical pure function symm_op_has_translation(this, prec) result(has)
    type(symm_op_t),  intent(in)  :: this
    FLOAT,            intent(in)  :: prec

    has = any(abs(this%trans_red(1:this%dim)) >= prec)

  end function symm_op_has_translation

  ! -------------------------------------------------------------------------------

  function symm_op_rotation_matrix_red(this) result(matrix)
    type(symm_op_t),  intent(in)  :: this
    integer                       :: matrix(1:this%dim, 1:this%dim)

    matrix(1:this%dim, 1:this%dim) = this%rot_red(1:this%dim, 1:this%dim)

  end function symm_op_rotation_matrix_red

  ! -------------------------------------------------------------------------------

  function symm_op_rotation_matrix_cart(this) result(matrix)
    type(symm_op_t),  intent(in)  :: this
    FLOAT                         :: matrix(1:this%dim, 1:this%dim)

    matrix(1:this%dim, 1:this%dim) = this%rot_cart(1:this%dim, 1:this%dim)

  end function symm_op_rotation_matrix_cart


  ! -------------------------------------------------------------------------------
 
  function symm_op_translation_vector_red(this) result(vector)
    type(symm_op_t),  intent(in)  :: this
    FLOAT                         :: vector(1:this%dim)

    vector(1:this%dim) = this%trans_red(1:this%dim)

  end function symm_op_translation_vector_red

  ! -------------------------------------------------------------------------------

  function symm_op_translation_vector_cart(this) result(vector)
    type(symm_op_t),  intent(in)  :: this
    FLOAT                         :: vector(1:this%dim)

    vector(1:this%dim) = this%trans_cart(1:this%dim)

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
    integer                       :: bb(1:this%dim)

    bb(1:this%dim) = nint(dsymm_op_apply_red(this, TOFLOAT(aa)))
    
  end function isymm_op_apply_red

  ! -------------------------------------------------------------------------------

  function isymm_op_apply_inv_red(this, aa) result(bb)
    type(symm_op_t),  intent(in)  :: this
    integer,          intent(in)  :: aa(:) !< (3)
    integer                       :: bb(1:this%dim)

    bb(1:this%dim) = nint(dsymm_op_apply_inv_red(this, TOFLOAT(aa)))
    
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

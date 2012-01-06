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
!! $Id: symm_op.F90 3479 2007-11-09 08:36:10Z xavier $

#include "global.h"

module symm_op_m
  use global_m
  use messages_m
  use profiling_m

  implicit none

  private
  
  public ::                        &
       symm_op_init,               &
       symm_op_init_from_dump,     &
       symm_op_dump,               &
       symm_op_copy,               &
       symm_op_end,                &
       symm_op_apply,              &
       symm_op_apply_inv,          &
       symm_op_invariant,          &
       symm_op_has_translation,    &
       symm_op_rotation_matrix,    &
       symm_op_translation_vector, &
       symm_op_t

  type symm_op_t
    private
    integer :: rotation(1:3, 1:3)
    FLOAT   :: translation(1:3)
  end type symm_op_t
  
  interface symm_op_apply
    module procedure isymm_op_apply, dsymm_op_apply, zsymm_op_apply
  end interface

  interface symm_op_apply_inv
    module procedure isymm_op_apply_inv, dsymm_op_apply_inv, zsymm_op_apply_inv
  end interface

contains
  
  subroutine symm_op_init(this, rot, trans)
    type(symm_op_t),     intent(out) :: this
    integer,             intent(in)  :: rot(:, :)
    FLOAT, optional,     intent(in)  :: trans(:)

    PUSH_SUB(symm_op_init)

    this%rotation(1:3, 1:3) = rot(1:3, 1:3)

    if(present(trans)) then
      this%translation(1:3) = trans(1:3)
    else
      this%translation(1:3) = M_ZERO
    end if

    POP_SUB(symm_op_init)
  end subroutine symm_op_init
  
  ! -------------------------------------------------------------------------------
  subroutine symm_op_init_from_dump(this, iunit)
    type(symm_op_t), intent(inout) :: this
    integer,         intent(in)    :: iunit
    !
    integer :: gb
    !
    PUSH_SUB(symm_op_init_from_dump)
    read(iunit) gb
    ASSERT(gb==GUARD_BITS)
    read(iunit) this%rotation
    read(iunit) this%translation
    read(iunit) gb
    ASSERT(gb==GUARD_BITS)
    POP_SUB(symm_op_init_from_dump)
    return
  end subroutine symm_op_init_from_dump

  ! -------------------------------------------------------------------------------
  subroutine symm_op_dump(this, iunit)
    type(symm_op_t), intent(in) :: this
    integer,         intent(in) :: iunit
    !
    PUSH_SUB(symm_op_dump)
    write(iunit) GUARD_BITS
    write(iunit) this%rotation
    write(iunit) this%translation
    write(iunit) GUARD_BITS
    POP_SUB(symm_op_dump)
    return
  end subroutine symm_op_dump

  ! -------------------------------------------------------------------------------
  subroutine symm_op_copy(inp, outp)
    type(symm_op_t),     intent(in) :: inp
    type(symm_op_t),     intent(out) :: outp

    PUSH_SUB(symm_op_copy)

    outp%rotation(1:3, 1:3) =  inp%rotation(1:3, 1:3)
    outp%translation(1:3) =  inp%translation(1:3)

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
  logical pure function symm_op_invariant(this, aa, prec) result(invariant)
    type(symm_op_t),  intent(in)  :: this
    FLOAT,            intent(in)  :: aa(:)
    FLOAT,            intent(in)  :: prec

    FLOAT :: cc(1:3)

    cc = symm_op_apply(this, aa)
    invariant = all(abs(cc(1:3) - aa(1:3)) < prec)

  end function symm_op_invariant

  ! -------------------------------------------------------------------------------
  logical pure function symm_op_has_translation(this, prec) result(has)
    type(symm_op_t),  intent(in)  :: this
    FLOAT,            intent(in)  :: prec

    has = any(abs(this%translation(1:3)) >= prec)

  end function symm_op_has_translation

  ! -------------------------------------------------------------------------------

  function symm_op_rotation_matrix(this) result(matrix)
    type(symm_op_t),  intent(in)  :: this
    FLOAT                         :: matrix(1:3, 1:3)

    matrix(1:3, 1:3) = this%rotation(1:3, 1:3)

  end function symm_op_rotation_matrix

  ! -------------------------------------------------------------------------------
 
  function symm_op_translation_vector(this) result(vector)
    type(symm_op_t),  intent(in)  :: this
    FLOAT                         :: vector(1:3)

    vector(1:3) = this%translation(1:3)

  end function symm_op_translation_vector

#include "undef.F90"
#include "integer.F90"
#include "symm_op_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "symm_op_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "symm_op_inc.F90"

end module symm_op_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

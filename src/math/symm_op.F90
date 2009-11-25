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
  
  public ::                     &
       symm_op_init,            &
       symm_op_end,             &
       symm_op_apply,           &
       symm_op_apply_inv,       &
       symm_op_invariant,       &
       symm_op_t

  type symm_op_t
    private
    integer :: rotation(1:3, 1:3)
    FLOAT   :: translation(1:3)
  end type symm_op_t
  
contains
  
  subroutine symm_op_init(this, rot, trans)
    type(symm_op_t),     intent(out) :: this
    integer,             intent(in)  :: rot(:, :)
    FLOAT,               intent(in)  :: trans(:)

    this%rotation(1:3, 1:3) = rot(1:3, 1:3)
    this%translation(1:3) = trans(1:3)
  end subroutine symm_op_init
  
  ! -------------------------------------------------------------------------------

  subroutine symm_op_end(this)
    type(symm_op_t),  intent(inout) :: this

    !nothing to do for the moment
  end subroutine symm_op_end

  ! -------------------------------------------------------------------------------
  pure function symm_op_apply(this, aa) result(bb)
    type(symm_op_t),  intent(in)  :: this
    FLOAT,            intent(in)  :: aa(:)
    FLOAT                         :: bb(1:3)

    bb(1:3) = matmul(aa(1:3), dble(this%rotation(1:3, 1:3))) + this%translation(1:3)
      
  end function symm_op_apply

  ! -------------------------------------------------------------------------------
  pure function symm_op_apply_inv(this, aa) result(bb)
    type(symm_op_t),  intent(in)  :: this
    FLOAT,            intent(in)  :: aa(:)
    FLOAT                         :: bb(1:3)

    bb(1:3) = matmul(dble(this%rotation(1:3, 1:3)), aa(1:3))
      
  end function symm_op_apply_inv
  
  ! -------------------------------------------------------------------------------
  logical pure function symm_op_invariant(this, aa, prec) result(invariant)
    type(symm_op_t),  intent(in)  :: this
    FLOAT,            intent(in)  :: aa(:)
    FLOAT,            intent(in)  :: prec

    FLOAT :: cc(1:3)

    cc = symm_op_apply(this, aa)
    invariant = all(abs(cc(1:3) - aa(1:3)) < prec)

  end function symm_op_invariant

end module symm_op_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

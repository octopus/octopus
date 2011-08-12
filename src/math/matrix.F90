!! Copyright (C) 2011 X. Andrade
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
!! $Id: matrix.F90 7610 2011-03-25 03:41:01Z xavier $

#include "global.h"

module matrix_m
  use global_m
  use messages_m
  use profiling_m
  use types_m

  implicit none

  private
  public ::                        &
    matrix_t,                      &
    matrix_init,                   &
    matrix_end

  type matrix_t
    private
    integer        :: dim(1:2)
    FLOAT, pointer :: dmat(:, :)
    FLOAT, pointer :: zmat(:, :)
  end type matrix_t

contains

  ! ---------------------------------------------------------
  subroutine matrix_init(this, dim1, dim2, type)
    type(matrix_t),             intent(out) :: this
    integer,                    intent(in)  :: dim1
    integer,                    intent(in)  :: dim2
    type(type_t),               intent(in)  :: type

    PUSH_SUB(matrix_init)

    this%dim(1:2) = (/dim1, dim2/)

    nullify(this%dmat)
    nullify(this%zmat)

    ASSERT(type == TYPE_FLOAT .or. type == TYPE_CMPLX)

    if(type == TYPE_FLOAT) then
      SAFE_ALLOCATE(this%dmat(dim1, dim2))
    else
      SAFE_ALLOCATE(this%zmat(dim1, dim2))
    end if

    POP_SUB(matrix_init)
  end subroutine matrix_init

  ! ---------------------------------------------------------

  subroutine matrix_end(this)
    type(matrix_t),             intent(inout) :: this
    
    PUSH_SUB(matrix_end)

    SAFE_DEALLOCATE_P(this%dmat)
    SAFE_DEALLOCATE_P(this%zmat)

    POP_SUB(matrix_end)
  end subroutine matrix_end

  ! ---------------------------------------------------------

#include "undef.F90"
#include "real.F90"

#include "matrix_inc.F90"

#include "undef.F90"
#include "complex.F90"

#include "matrix_inc.F90"

end module matrix_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

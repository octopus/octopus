!! Copyright (C) 2011-2015 X. Andrade
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
!! $Id: matrix.F90 7610 2011-03-25 03:41:01Z xavier $

#include "global.h"

module matrix_oct_m
  use global_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m
  use types_oct_m

  implicit none

  private
  public ::                        &
    matrix_t,                      &
    matrix_init,                   &
    matrix_end,                    &
    matrix_set_zero,               &
    matrix_set_block,              &
    matrix_get_block,              &
    matrix_print,                  &
    matrix_type
  
  type matrix_t
    integer            :: dim(1:2)
    FLOAT, allocatable :: dmat(:, :)
    CMPLX, allocatable :: zmat(:, :)
    type(type_t)       :: type
    type(mpi_grp_t)    :: mpi_grp
  end type matrix_t

  interface matrix_init
    module procedure matrix_init_empty, dmatrix_init_data, zmatrix_init_data
  end interface matrix_init

  interface matrix_set_block
    module procedure dmatrix_set_block, zmatrix_set_block
  end interface matrix_set_block

  interface matrix_get_block
    module procedure dmatrix_get_block, zmatrix_get_block
  end interface matrix_get_block
  
contains

  ! ---------------------------------------------------------
  subroutine matrix_init_empty(this, dim1, dim2, type, mpi_grp)
    type(matrix_t),             intent(out) :: this    !< the object to be initialized
    integer,                    intent(in)  :: dim1    !< the first dimension of the matrix
    integer,                    intent(in)  :: dim2    !< the second dimension of the matrix
    type(type_t),               intent(in)  :: type    !< the type of the elements of the matrix TYPE_FLOAT or TYPE_CMPLX
    type(mpi_grp_t),            intent(in)  :: mpi_grp !< the group of processors that shares this matrix

    PUSH_SUB(matrix_init_empty)

    this%dim(1:2) = (/dim1, dim2/)

    this%type = type

    this%mpi_grp = mpi_grp
    
    ASSERT(type == TYPE_FLOAT .or. type == TYPE_CMPLX)

    if(type == TYPE_FLOAT) then
      SAFE_ALLOCATE(this%dmat(1:dim1, 1:dim2))
    else
      SAFE_ALLOCATE(this%zmat(1:dim1, 1:dim2))
    end if

    POP_SUB(matrix_init_empty)
  end subroutine matrix_init_empty

  ! ---------------------------------------------------------

  subroutine matrix_end(this)
    type(matrix_t),             intent(inout) :: this !< the object to be destroyed
    
    PUSH_SUB(matrix_end)

    SAFE_DEALLOCATE_A(this%dmat)
    SAFE_DEALLOCATE_A(this%zmat)

    POP_SUB(matrix_end)
  end subroutine matrix_end

  ! ---------------------------------------------------------

  subroutine matrix_set_zero(this)
    type(matrix_t), intent(inout) :: this

    PUSH_SUB(matrix_set_zero)

    if(this%type == TYPE_FLOAT) then
      this%dmat = CNST(0.0)
    else
      this%zmat = CNST(0.0)
    end if
    
    POP_SUB(matrix_set_zero)
  end subroutine matrix_set_zero

  ! ---------------------------------------------------------

  subroutine matrix_print(this)
    type(matrix_t), intent(in) :: this

    integer :: ii
    
    PUSH_SUB(matrix_print)

    do ii = 1, this%dim(1)
      if(this%type == TYPE_FLOAT) then
        print*, this%dmat(ii, 1:this%dim(2))
      else
        print*, this%zmat(ii, 1:this%dim(2))
      end if
    end do
    
    POP_SUB(matrix_print)
  end subroutine matrix_print

  ! ---------------------------------------------------------

  type(type_t) function matrix_type(this)
    type(matrix_t), intent(in) :: this

    matrix_type = this%type
  end function matrix_type

  
#include "undef.F90"
#include "real.F90"

#include "matrix_inc.F90"

#include "undef.F90"
#include "complex.F90"

#include "matrix_inc.F90"

end module matrix_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

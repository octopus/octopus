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

 ! ---------------------------------------------------------
subroutine X(matrix_init_data)(this, dim1, dim2, data, mpi_grp)
  type(matrix_t),             intent(out) :: this
  integer,                    intent(in)  :: dim1
  integer,                    intent(in)  :: dim2
  R_TYPE,                     intent(in)  :: data(:, :)
  type(mpi_grp_t),            intent(in)  :: mpi_grp !< the group of processors that shares this matrix

  PUSH_SUB(X(matrix_init_data))

  ASSERT(all(ubound(data) == (/dim1, dim2/)))

  this%dim(1:2) = (/dim1, dim2/)

  this%type = R_TYPE_VAL

  this%mpi_grp = mpi_grp

  ASSERT(this%type == TYPE_FLOAT .or. this%type == TYPE_CMPLX)

  SAFE_ALLOCATE(this%X(mat)(1:dim1, 1:dim2))

  this%X(mat)(1:dim1, 1:dim2) = data(1:dim1, 1:dim2)

  POP_SUB(X(matrix_init_data))
end subroutine X(matrix_init_data)

 ! ---------------------------------------------------------

subroutine X(matrix_set_block)(this, min1, max1, min2, max2, data)
  type(matrix_t),             intent(inout) :: this
  integer,                    intent(in)    :: min1
  integer,                    intent(in)    :: max1
  integer,                    intent(in)    :: min2
  integer,                    intent(in)    :: max2
  R_TYPE,                     intent(in)    :: data(:, :)

  PUSH_SUB(X(matrix_init_data))

  !    print*, min1, max1, this%dim(1)
  !    print*, min2, max2, this%dim(2)

  ASSERT(this%type == R_TYPE_VAL)
  ASSERT(min1 <= max1)
  ASSERT(min2 <= max2)
  ASSERT(0 < min1 .and. min1 <= this%dim(1))
  ASSERT(0 < max1 .and. max1 <= this%dim(1))
  ASSERT(0 < min2 .and. min2 <= this%dim(2))
  ASSERT(0 < max2 .and. max2 <= this%dim(2))

  this%X(mat)(min1:max1, min2:max2) = data(1:max1 - min1 + 1, 1:max2 - min2 + 1)

  POP_SUB(X(matrix_init_data))
end subroutine X(matrix_set_block)

 ! ---------------------------------------------------------

subroutine X(matrix_get_block)(this, min1, max1, min2, max2, data)
  type(matrix_t),             intent(in)    :: this
  integer,                    intent(in)    :: min1
  integer,                    intent(in)    :: max1
  integer,                    intent(in)    :: min2
  integer,                    intent(in)    :: max2
  R_TYPE,                     intent(inout) :: data(:, :)

  PUSH_SUB(X(matrix_init_data))

  !    print*, min1, max1, this%dim(1)
  !    print*, min2, max2, this%dim(2)

  ASSERT(this%type == R_TYPE_VAL)
  ASSERT(min1 <= max1)
  ASSERT(min2 <= max2)
  ASSERT(0 < min1 .and. min1 <= this%dim(1))
  ASSERT(0 < max1 .and. max1 <= this%dim(1))
  ASSERT(0 < min2 .and. min2 <= this%dim(2))
  ASSERT(0 < max2 .and. max2 <= this%dim(2))

  data(1:max1 - min1 + 1, 1:max2 - min2 + 1) = this%X(mat)(min1:max1, min2:max2)

  POP_SUB(X(matrix_init_data))
end subroutine X(matrix_get_block)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

!! Copyright (C) 2021 S. Ohlmann
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

module merge_sorted_oct_m
  use global_oct_m
  use heap_oct_m
  use messages_oct_m
  implicit none
  private
  public :: merge_sorted_arrays

contains
  ! Merge a number of sorted arrays
  ! The sorted arrays are linearly ordered in array and their sizes are given in sizes.
  ! The merged array is returned and optionally also an index map that can be used
  ! to also merge another data array.
  !
  ! This routine uses a minheap to to a k-way merge.
  subroutine merge_sorted_arrays(array, sizes, merged, index_map)
    integer,           intent(in)    :: array(:)
    integer,           intent(in)    :: sizes(:)
    integer,           intent(inout) :: merged(:)
    integer, optional, intent(inout) :: index_map(:)

    type(heap_t) :: heap

    PUSH_SUB(merge_sorted_arrays)

    call heap%init(array, sizes)
    call heap%merge(merged, index_map)
    call heap%end()

    POP_SUB(merge_sorted_arrays)
  end subroutine
end module merge_sorted_oct_m

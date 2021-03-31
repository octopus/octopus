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
  use profiling_oct_m
  implicit none
  private
  public :: merge_sorted_arrays

  interface merge_sorted_arrays
    module procedure merge_sorted_arrays_4
    module procedure merge_sorted_arrays_8
  end interface


contains
  ! Merge a number of sorted arrays
  ! The sorted arrays are linearly ordered in array and their sizes are given in sizes.
  ! The merged array is returned and optionally also an index map that can be used
  ! to also merge another data array.
  !
  ! This routine uses a minheap to to a k-way merge.
  ! There are two versions, one for integer(4) and one for integer(8).
  subroutine merge_sorted_arrays_8(array, sizes, merged, index_map)
    integer(8),           intent(in)    :: array(:)
    integer(8),           intent(in)    :: sizes(:)
    integer(8),           intent(inout) :: merged(:)
    integer(8), optional, intent(inout) :: index_map(:)

    type(heap_t) :: heap

    PUSH_SUB(merge_sorted_arrays_8)

    call heap%init(array, sizes)
    call heap%merge(merged, index_map)
    call heap%end()

    POP_SUB(merge_sorted_arrays_8)
  end subroutine merge_sorted_arrays_8

  ! do type conversions before and after the usage of the heap
  subroutine merge_sorted_arrays_4(array, sizes, merged, index_map)
    integer,           intent(in)    :: array(:)
    integer,           intent(in)    :: sizes(:)
    integer,           intent(inout) :: merged(:)
    integer, optional, intent(inout) :: index_map(:)

    type(heap_t) :: heap
    integer(8), allocatable :: larray(:), lsizes(:), lmerged(:), lindex_map(:)

    PUSH_SUB(merge_sorted_arrays_4)

    SAFE_ALLOCATE(larray(1:ubound(array, dim=1)))
    SAFE_ALLOCATE(lsizes(1:ubound(sizes, dim=1)))
    SAFE_ALLOCATE(lmerged(1:ubound(merged, dim=1)))

    call heap%init(larray, lsizes)
    if (present(index_map)) then
      SAFE_ALLOCATE(lindex_map(1:ubound(index_map, dim=1)))
      call heap%merge(lmerged, lindex_map)
      index_map(:) = int(lindex_map(:))
      SAFE_DEALLOCATE_A(lindex_map)
    else
      call heap%merge(lmerged)
    end if
    call heap%end()

    merged(:) = int(lmerged(:), 4)

    SAFE_DEALLOCATE_A(larray)
    SAFE_DEALLOCATE_A(lsizes)
    SAFE_DEALLOCATE_A(lmerged)

    POP_SUB(merge_sorted_arrays_4)
  end subroutine merge_sorted_arrays_4
end module merge_sorted_oct_m

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

module heap_oct_m
  use global_oct_m
  use messages_oct_m
  implicit none
  private
  public :: heap_t

  ! A min-heap type for implementing a k-way merge
  ! The min-heap is built on the first indices of the sorted arrays.
  ! It allows to always extract the minimum element of the sorted arrays
  ! and the heap structure is then maintained by swapping the indices and
  ! the corresponding sizes.
  ! To merge the sorted arrays, the minimum of the heap is extracted and
  ! replaced by the next element of that array; then the min-heap property
  ! is restored by down-sifting (minheapify).
  ! The time complexity of this merge is O(n log k) where n is the full size
  ! of the arrays and k is the number of arrays.
  type, private :: heap_t
    integer(8), pointer :: a(:)            !< a list of sorted arrays
    integer(8) :: length                   !< the number of sorted arrays -> this is the length of the heap
    integer(8), allocatable :: indices(:)  !< the starting indices of the sorted arrays
    integer(8), allocatable :: sizes(:)    !< the sizes of the sorted arrays
  contains
    procedure :: init => heap_init
    procedure :: end => heap_end
    procedure :: merge => heap_merge
  end type heap_t

contains
  ! Initialize the heap
  subroutine heap_init(heap, list, sizes)
    class(heap_t),      intent(inout) :: heap      !< the heap type
    integer(8), target, intent(in)    :: list(:)   !< the list of numbers to be merged
    integer(8),         intent(in)    :: sizes(:)  !< the sizes of the arrays
    integer(8) :: i

    PUSH_SUB(heap_init)
    ! we keep a pointer to access the data; we do not need to modify it
    heap%a => list

    ! get the number of arrays to merge
    heap%length = ubound(sizes, dim=1)
    allocate(heap%sizes(heap%length))
    allocate(heap%indices(heap%length))
    do i = 1, heap%length
      heap%sizes(i) = sizes(i)
    end do
    ! compute the starting indices
    heap%indices(1) = 1
    do i = 2, heap%length
      heap%indices(i) = heap%indices(i-1) + heap%sizes(i-1)
    end do
    ! sanity check
    if (sum(heap%sizes) /= ubound(list, dim=1)) then
      message(1) = "Error! Mismatch in size of array when initializing heap!"
      call messages_fatal(1)
    end if

    ! build the heap
    call build_minheap(heap)
    POP_SUB(heap_init)
  end subroutine

  subroutine heap_end(heap)
    class(heap_t), intent(inout) :: heap

    PUSH_SUB(heap_end)
    nullify(heap%a)
    heap%length = 0
    deallocate(heap%sizes)
    deallocate(heap%indices)
    POP_SUB(heap_end)
  end subroutine

  subroutine heap_merge(heap, merged, index_map)
    class(heap_t),        intent(inout) :: heap
    integer(8),           intent(inout) :: merged(:)
    integer(8), optional, intent(inout) :: index_map(:)  !< index map for sorting another data array

    integer(8) :: i

    PUSH_SUB(heap_merge)
    if (sum(heap%sizes) /= ubound(merged, dim=1)) then
      message(1) = "Error! Mismatch in size of array when doing k-way merge!"
      call messages_fatal(1)
    end if

    do i = 1, sum(heap%sizes)
      ! extract the minimum element of the heap
      merged(i) = heap%a(heap%indices(1))
      if (present(index_map)) then
        index_map(i) = heap%indices(1)
      end if
      ! replace by the next element of that array and reduce its size
      heap%indices(1) = heap%indices(1) + 1
      heap%sizes(1) = heap%sizes(1) - 1
      if (heap%sizes(1) == 0) then
        ! if this array is already empty, swap it with the last one and
        ! reduce the size of the heap by one
        call swap(heap, 1_8, heap%length)
        heap%length = heap%length - 1
      end if
      ! down-sift (restore min-heap property)
      call minheapify(heap, 1_8)
    end do
    POP_SUB(heap_merge)
  end subroutine

  ! build the minheap by downsifting all parent nodes
  subroutine build_minheap(heap)
    class(heap_t), intent(inout) :: heap

    integer(8) :: i

    ! this goes bottom up for all parent nodes
    do i = heap%length/2, 1, -1
      call minheapify(heap, i)
    end do
  end subroutine

  ! downsifting from index i (restore minheap property)
  recursive subroutine minheapify(heap, i)
    class(heap_t),    intent(inout) :: heap
    integer(8),       intent(in)    :: i

    integer(8) :: left, right, smallest

    left = 2*i
    right = 2*i + 1
    smallest = i

    ! get the smallest value of both children (if they are still in the heap)
    if (left <= heap%length) then
      if (is_smaller(heap, left, smallest)) then
        smallest = left
      end if
    end if
    if (right <= heap%length) then
      if (is_smaller(heap, right, smallest)) then
        smallest = right
      end if
    end if
    ! we only swap and downsift if one of the children is smaller
    if (smallest /= i) then
      call swap(heap, i, smallest)
      call minheapify(heap, smallest)
    end if
  end subroutine

  ! swap two of the sorted arrays in the heap
  subroutine swap(heap, i, j)
    class(heap_t),    intent(inout) :: heap
    integer(8),       intent(in)    :: i
    integer(8),       intent(in)    :: j

    ! swap index and associated size
    call swap_int(heap%indices, i, j)
    call swap_int(heap%sizes, i, j)
  end subroutine

  ! swap two integers
  subroutine swap_int(a, i, j)
    integer(8), intent(inout) :: a(:)
    integer(8), intent(in)    :: i
    integer(8), intent(in)    :: j

    integer(8) :: tmp

    tmp = a(i)
    a(i) = a(j)
    a(j) = tmp
  end subroutine

  ! compare the current value of two arrays in the heap
  logical function is_smaller(heap, i, j)
    class(heap_t),    intent(inout) :: heap
    integer(8),       intent(in)    :: i
    integer(8),       intent(in)    :: j

    is_smaller = heap%a(heap%indices(i)) < heap%a(heap%indices(j))
  end function
end module heap_oct_m

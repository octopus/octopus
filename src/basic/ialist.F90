!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: io.F90 3613 2007-11-29 16:47:41Z xavier $

#include "global.h"

! This module implements a simple associative list for integer keys and
! values. It is used by the separate changing hash-table implementation.

module ialist_m
  use global_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::        &
    ialist_t,      &
    iacons_t,      &
    ialist_init,   &
    ialist_end,    &
    ialist_drop,   &
    ialist_delete, &
    ialist_insert, &
    ialist_lookup

  type ialist_t
    type(iacons_t), pointer :: head
    integer                 :: length
  end type ialist_t

  type iacons_t
    integer                 :: key
    integer                 :: val
    type(iacons_t), pointer :: next
  end type iacons_t

contains

  ! ---------------------------------------------------------
  ! Set up an empty associative list.
  subroutine ialist_init(l)
    type(ialist_t), intent(out) :: l

    PUSH_SUB(ialist_init)

    l%length = 0
    nullify(l%head)

    POP_SUB(ialist_init)
  end subroutine ialist_init


  ! ---------------------------------------------------------
  ! Drop the head of the list.
  subroutine ialist_drop(l)
    type(ialist_t), intent(inout) :: l

    type(iacons_t), pointer :: old_head

    PUSH_SUB(ialist_drop)

    if(l%length.gt.0) then
      old_head => l%head
      l%head => l%head%next
      SAFE_DEALLOCATE_P(old_head)
      l%length = l%length - 1
    end if

    POP_SUB(ialist_drop)
  end subroutine ialist_drop


  ! ---------------------------------------------------------
  ! Delete the pair with the given key from the list.
  subroutine ialist_delete(key, l)
    integer,        intent(in)    :: key
    type(ialist_t), intent(inout) :: l

    integer                 :: i
    type(iacons_t), pointer :: ptr

    PUSH_SUB(ialist_delete)

    ! Deletions are only possible from nonempty lists.
    if (l%length.ge.1) then
      ! We take the head as special case.
      if(l%head%key.eq.key) then
        call ialist_drop(l)
      else
        ! Process tail of list, if head had a different key.
        if(l%length.gt.1) then
          ptr => l%head
          do i = 1, l%length - 2
            if(ptr%next%key.eq.key) then
              SAFE_DEALLOCATE_P(ptr%next)
              ptr%next => ptr%next%next
              l%length = l%length - 1
              exit
            else
              ptr => ptr%next
            end if
          end do
        end if
      end if
    end if

    POP_SUB(ialist_delete)
  end subroutine ialist_delete


  ! ---------------------------------------------------------
  ! Insert a (key, val) pair in the list. If key is already present,
  ! its value is updated.
  subroutine ialist_insert(key, val, l)
    integer,        intent(in) :: key
    integer,        intent(in) :: val
    type(ialist_t), intent(inout) :: l
  
    integer                 :: i
    logical                 :: found
    type(iacons_t), pointer :: ptr

    ! List is empty.
    if(l%length.eq.0) then
      SAFE_ALLOCATE(ptr)
      ptr%key  =  key
      ptr%val  =  val
      l%head   => ptr
      l%length =  1
    ! List is not empty.
    else
      ! Look for key in list.
      ptr => l%head
      found = .false.
      do i = 1, l%length
        ! If found, replace key`s value.
        if(ptr%key.eq.key) then
          ptr%val = val
          found = .true.
          exit
        end if
        ptr => ptr%next
      end do
      ! If not found, prepend a new cons to the list.
      if(.not. found) then
        SAFE_ALLOCATE(ptr)
        ptr%key  =  key
        ptr%val  =  val
        ptr%next => l%head
        l%head   => ptr
        l%length =  l%length + 1
      end if
    end if
  end subroutine ialist_insert


  ! ---------------------------------------------------------
  ! Get the value of key. If found is present and .false. the return
  ! value of ialist_lookup is meaningless (i.e. undefined). For this reason,
  ! always pass found if you do not know, for different reasons, that key
  ! is member of the list.
  integer function ialist_lookup(key, l, found)
    integer,           intent(in) :: key
    type(ialist_t),    intent(in) :: l
    logical, optional, intent(out) :: found

    integer                 :: i
    type(iacons_t), pointer :: ptr

    ! we define a default return value to avoid problems with checks
    ! for uninitialized values
    ialist_lookup = -1

    ptr => l%head
    do i = 1, l%length
      if(ptr%key.eq.key) then
        ialist_lookup = ptr%val
        if(present(found)) then
          found = .true.
        end if
        return
      else
        ptr => ptr%next
      end if
    end do
    if(present(found)) then
      found = .false.
    end if
  end function ialist_lookup


  ! ---------------------------------------------------------
  ! Deallocate all cons of the list.
  subroutine ialist_end(l)
    type(ialist_t), intent(inout) :: l

    integer :: i

    PUSH_SUB(ialist_end)

    do i = 1, l%length
      call ialist_drop(l)
    end do
    l%length = 0
    nullify(l%head)

    POP_SUB(ialist_end)
  end subroutine ialist_end
end module ialist_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

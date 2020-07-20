!! Copyright (C) 2019 M. Oliveira
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

module linked_list_oct_m
  use global_oct_m
  use list_node_oct_m
  use messages_oct_m
  implicit none

  private
  public :: linked_list_t,  &
            list_iterator_t


  !---------------------------------------------------------------------------
  !> This class implements a linked list of unlimited polymorphic values. This
  !> allows the storage of any type of data. Iterating over the list is done
  !> using the associated list_counter_t. There are two ways of iterating. The
  !> first is by using a "do while" construct:
  !>
  !>  call iter%start(list)
  !>  do while (iter%has_next())
  !>     value => iter%get_next()
  !>     ...
  !>  end do
  !>
  !> The second method is with a simple "do":
  !>
  !>  call iter%start(list)
  !>  do
  !>     if (.not. iter%has_next()) exit
  !>     value => iter%get_next()
  !>     ...
  !>  end do
  !>
  type :: linked_list_t
    private
    class(list_node_t), pointer :: first_node => null()
    class(list_node_t), pointer :: last_node => null()
  contains
    procedure :: add => linked_list_add_node
    final     :: linked_list_finalize
  end type linked_list_t

  type :: list_iterator_t
    private
    class(list_node_t), pointer :: next_node => null()
  contains
    procedure :: start    => list_iterator_start
    procedure :: has_next => list_iterator_has_next
    procedure :: get_next => list_iterator_get_next
  end type list_iterator_t

contains

  ! ---------------------------------------------------------
  subroutine linked_list_add_node(this, value)
    class(linked_list_t) :: this
    class(*),             target        :: value

    class(list_node_t), pointer :: new_node

    PUSH_SUB(linked_list_add_node)

    if (.not. associated(this%first_node)) then
      this%first_node => list_node(value, this%first_node)
      this%last_node => this%first_node
    else
      new_node => list_node(value, this%last_node%next())
      call this%last_node%set_next(new_node)
      this%last_node => new_node
    end if

    POP_SUB(linked_list_add_node)
  end subroutine linked_list_add_node

  ! ---------------------------------------------------------
  subroutine linked_list_finalize(this)
    type(linked_list_t), intent(inout) :: this

    class(list_node_t), pointer :: current, next

    PUSH_SUB(linked_list_finalize)

    current => this%first_node
    do while (associated(current))
      next => current%next()
      deallocate(current)
      current => next
    end do

    POP_SUB(linked_list_finalize)
  end subroutine linked_list_finalize

  ! ---------------------------------------------------------
  subroutine list_iterator_start(this, list)
    class(list_iterator_t), intent(inout)      :: this
    class(linked_list_t),   intent(in), target :: list

    PUSH_SUB(list_iterator_start)

    this%next_node => list%first_node

    POP_SUB(list_iterator_start)
  end subroutine list_iterator_start

  ! ---------------------------------------------------------
  logical function list_iterator_has_next(this)
    class(list_iterator_t), intent(in) :: this

    PUSH_SUB(list_iterator_has_next)

    list_iterator_has_next = associated(this%next_node)

    POP_SUB(list_iterator_has_next)
  end function list_iterator_has_next

  ! ---------------------------------------------------------
  function list_iterator_get_next(this) result(value)
    class(list_iterator_t), intent(inout) :: this
    class(*),              pointer        :: value

    PUSH_SUB(list_iterator_get_next)

    value => this%next_node%get()
    this%next_node => this%next_node%next()

    POP_SUB(list_iterator_get_next)
  end function list_iterator_get_next

end module linked_list_oct_m

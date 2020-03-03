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
  public :: linked_list_t, &
            list_counter_t


  !---------------------------------------------------------------------------
  !> This class implements a linked list of unlimited polymorphic values. This
  !> allows the storeage of any type of data. Iterating over the list is done
  !> using the associated list_counter_t. There are two ways of iterating. The
  !> first is by using a "do while" construct:
  !>
  !>  counter = list%counter_start()
  !>  do while (counter%iterate())
  !>     value => list%get(counter)
  !>     ...
  !>  end do
  !>
  !> The second method is with a simple "do":
  !>
  !>  counter = list%counter_start()
  !>  do
  !>     call counter%increment()
  !>     value => list%get(counter)
  !>     ...
  !>     if (counter%has_next()) exit
  !>  end do
  !>
  type :: linked_list_t
    private
    class(list_node_t), pointer :: first_node => null()
    class(list_node_t), pointer :: last_node => null()
  contains
    procedure :: add => linked_list_add_node
    procedure :: start_counter => linked_list_start_counter
    procedure, nopass :: get => linked_list_get
    final     :: linked_list_finalize
  end type linked_list_t

  type :: list_counter_t
    private
    class(list_node_t), pointer :: current_node => null()
    class(list_node_t), pointer :: first_node => null()
  contains
    procedure :: iterate => counter_iterate
    procedure :: increment => counter_increment
    procedure :: has_next => counter_has_next
  end type list_counter_t

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
  function linked_list_get(counter) result(value)
    type(list_counter_t), intent(in) :: counter
    class(*),             pointer    :: value

    PUSH_SUB(linked_list_get)

    value => counter%current_node%get()

    POP_SUB(linked_list_get)
  end function linked_list_get

  ! ---------------------------------------------------------
  type(list_counter_t) function linked_list_start_counter(this) result(counter)
    class(linked_list_t), intent(in), target :: this

    counter%first_node => this%first_node
    nullify(counter%current_node)

  end function linked_list_start_counter

  ! ---------------------------------------------------------
  logical function counter_iterate(this)
    class(list_counter_t), intent(inout) :: this

    PUSH_SUB(counter_iterate)

    ! Get the next node in the list
    call this%increment()

    ! Are we done?
    counter_iterate = this%has_next()

    POP_SUB(counter_iterate)
  end function counter_iterate

  ! ---------------------------------------------------------
  subroutine counter_increment(this)
    class(list_counter_t), intent(inout) :: this

    PUSH_SUB(counter_increment)

    if (associated(this%current_node)) then
      this%current_node => this%current_node%next()
    else
      this%current_node => this%first_node
    end if

    POP_SUB(counter_increment)
  end subroutine counter_increment

  ! ---------------------------------------------------------
  logical function counter_has_next(this)
    class(list_counter_t), intent(in) :: this

    PUSH_SUB(counter_has_next)

    counter_has_next = associated(this%current_node)

    POP_SUB(counter_has_next)
  end function counter_has_next

end module linked_list_oct_m

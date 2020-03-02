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
  public :: linked_list_t

  type :: linked_list_t
    private
    class(list_node_t), pointer :: first_node => null()
    class(list_node_t), pointer :: last_node => null()
    class(list_node_t), pointer :: current_node => null()
  contains
    procedure :: add => linked_list_add_node
    procedure :: next => linked_list_next
    procedure :: current => linked_list_current
    procedure :: rewind => linked_list_rewind
    procedure :: has_more_values => linked_list_has_more_values
    procedure :: iterate => linked_list_iterate
    final     :: linked_list_finalize
  end type linked_list_t

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
  function linked_list_current(this)
    class(linked_list_t), intent(in) :: this
    class(*),             pointer    :: linked_list_current

    PUSH_SUB(linked_list_current)

    linked_list_current => this%current_node%get()

    POP_SUB(linked_list_current)
  end function linked_list_current

  ! ---------------------------------------------------------
  subroutine linked_list_next(this)
    class(linked_list_t), intent(inout) :: this

    PUSH_SUB(linked_list_next)

    this%current_node => this%current_node%next()

    POP_SUB(linked_list_next)
  end subroutine linked_list_next

  ! ---------------------------------------------------------
  logical function linked_list_has_more_values(this)
    class(linked_list_t), intent(in) :: this

    PUSH_SUB(linked_list_has_more_values)

    linked_list_has_more_values = associated(this%current_node)

    POP_SUB(linked_list_has_more_values)
  end function linked_list_has_more_values

  ! ---------------------------------------------------------
  subroutine linked_list_rewind(this)
    class(linked_list_t), intent(inout) :: this

    PUSH_SUB(linked_list_rewind)

    this%current_node => this%first_node

    POP_SUB(linked_list_rewind)
  end subroutine linked_list_rewind

  ! ---------------------------------------------------------
  logical function linked_list_iterate(this, iteration_counter, value)
    class(linked_list_t), intent(in)        :: this
    type(list_node_t),    pointer           :: iteration_counter
    class(*),             pointer, optional :: value

    PUSH_SUB(linked_list_iterate)

    ! Get the next node in the list
    if (associated(iteration_counter)) then
      iteration_counter => iteration_counter%next()
    else
      iteration_counter => this%first_node
    end if

    ! Are we done?
    linked_list_iterate = associated(iteration_counter)

    ! Get a pointer to the value stored in the list
    if (present(value) .and. linked_list_iterate) then
      value => iteration_counter%get()
    end if

    POP_SUB(linked_list_iterate)
  end function linked_list_iterate

  ! ---------------------------------------------------------
  subroutine linked_list_finalize(this)
    type(linked_list_t), intent(inout) :: this

    class(list_node_t), pointer :: next

    PUSH_SUB(linked_list_finalize)

    call this%rewind()
    do while (associated(this%current_node))
      next => this%current_node%next()
      deallocate(this%current_node)
      this%current_node => next
    end do

    POP_SUB(linked_list_finalize)
  end subroutine linked_list_finalize

end module linked_list_oct_m

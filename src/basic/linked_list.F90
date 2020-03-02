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
    class(list_node_t), pointer, public :: first_node => null()
    class(list_node_t), pointer :: last_node => null()
  contains
    procedure :: add => linked_list_add_node
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

end module linked_list_oct_m

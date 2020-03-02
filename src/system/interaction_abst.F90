!! Copyright (C) 2020 M. Oliveira
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

module interaction_abst_oct_m
  use global_oct_m
  use linked_list_oct_m
  use list_node_oct_m
  use messages_oct_m
  implicit none

  private
  public ::               &
    interaction_abst_t,   &
    interaction_list_t

  !> The only purpose of the following class is to act as a surrogate and
  !> avoid circular dependencies between the interactions and the systems.
  type, abstract :: interaction_abst_t
    private
  end type interaction_abst_t

  !> Wrapper class for linked lists to guarantee that only interactions are
  !> stored in the list. Note that we cannot simply make this class an extension
  !> of linked_list_t, as in that case we would no be able to change the
  !> interface of the methods, which is exactly what we want here.
  type interaction_list_t
    private
    type(linked_list_t) :: list
  contains
    procedure :: add => interaction_list_add_node
    procedure :: iterate => interaction_list_iterate
  end type interaction_list_t

contains

  ! ---------------------------------------------------------
  subroutine interaction_list_add_node(this, value)
    class(interaction_list_t)         :: this
    class(interaction_abst_t), target :: value

    PUSH_SUB(interaction_list_add_node)

    call this%list%add(value)

    POP_SUB(interaction_list_add_node)
  end subroutine interaction_list_add_node

  ! ---------------------------------------------------------
  logical function interaction_list_iterate(this, iteration_counter, value)
    class(interaction_list_t), intent(in)        :: this
    type(list_node_t),         pointer           :: iteration_counter
    class(interaction_abst_t), pointer, optional :: value

    class(*), pointer :: ptr

    PUSH_SUB(interaction_list_iterate)

    interaction_list_iterate = this%list%iterate(iteration_counter, ptr)

    if (present(value) .and. interaction_list_iterate) then
      select type (ptr)
      class is (interaction_abst_t)
        value => ptr
      class default
        ASSERT(.false.)
      end select
    end if

    POP_SUB(interaction_list_iterate)
  end function interaction_list_iterate

end module interaction_abst_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

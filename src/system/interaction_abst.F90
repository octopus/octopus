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

  type, extends(linked_list_t) :: interaction_list_t
    private
  contains
    procedure :: add => interaction_list_add_node
    procedure :: get_interaction => interaction_list_get
  end type interaction_list_t

contains

  ! ---------------------------------------------------------
  subroutine interaction_list_add_node(this, value)
    class(interaction_list_t)        :: this
    class(*),                 target :: value

    PUSH_SUB(interaction_list_add_node)

    select type (value)
    class is (interaction_abst_t)
      call this%linked_list_t%add(value)
    class default
      ASSERT(.false.)
    end select

    POP_SUB(interaction_list_add_node)
  end subroutine interaction_list_add_node

  ! ---------------------------------------------------------
  function interaction_list_get(this, counter) result(value)
    class(interaction_list_t), intent(in) :: this
    type(list_counter_t),      intent(in) :: counter
    class(interaction_abst_t), pointer    :: value

    class(*), pointer :: ptr

    PUSH_SUB(interaction_list_get)

    ptr => this%get(counter)
    select type (ptr)
    class is (interaction_abst_t)
      value => ptr
    class default
      ASSERT(.false.)
    end select

    POP_SUB(interaction_list_get)
  end function interaction_list_get

end module interaction_abst_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

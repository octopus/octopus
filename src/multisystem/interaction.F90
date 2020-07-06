!! Copyright (C) 2020 M. Oliveira, Heiko Appel
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

module interaction_oct_m
  use clock_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  implicit none

  private
  public ::                &
    interaction_t,         &
    interaction_end,       &
    interaction_list_t,    &
    interaction_iterator_t

  type, abstract :: interaction_t
    private
    !> The interaction requires access to some quantities from a system to be evaluated.
    integer,              public :: n_system_quantities  !< Number of quantities needed from the system
    integer, allocatable, public :: system_quantities(:) !< Identifiers of the quantities needed from the system
    type(clock_t), public :: clock !< Clock storing the time at which the interaction was last updated.
    character(len=:), public, allocatable :: label
  contains
    procedure :: init_clock => interaction_init_clock
    procedure(interaction_update),    deferred :: update
    procedure(interaction_calculate), deferred :: calculate
  end type interaction_t

  abstract interface
    logical function interaction_update(this, namespace, requested_time)
      import interaction_t
      import clock_t
      import namespace_t
      class(interaction_t),      intent(inout) :: this
      type(namespace_t),         intent(in)    :: namespace
      class(clock_t),            intent(in)    :: requested_time
    end function interaction_update

    subroutine interaction_calculate(this, namespace)
      import interaction_t
      import namespace_t
      class(interaction_t),      intent(inout) :: this
      type(namespace_t),         intent(in)    :: namespace
    end subroutine interaction_calculate
  end interface

  !> These classes extend the list and list iterator to make an interaction list
  type, extends(linked_list_t) :: interaction_list_t
    private
  contains
    procedure :: add => interaction_list_add_node
  end type interaction_list_t

  type, extends(linked_list_iterator_t) :: interaction_iterator_t
    private
  contains
    procedure :: get_next => interaction_iterator_get_next
  end type interaction_iterator_t

contains

  ! ---------------------------------------------------------
  subroutine interaction_init_clock(this, label, dt, algo_dt)
    class(interaction_t),      intent(inout) :: this
    character(len=*),          intent(in)    :: label
    FLOAT,                     intent(in)    :: dt
    FLOAT,                     intent(in)    :: algo_dt

    PUSH_SUB(interaction_init_clock)

    this%clock = clock_t(label, dt, algo_dt, initial_tick=-1)

    POP_SUB(interaction_init_clock)
  end subroutine interaction_init_clock

  ! ---------------------------------------------------------
  subroutine interaction_end(this)
    class(interaction_t), intent(inout) :: this

    PUSH_SUB(interaction_end)

    SAFE_DEALLOCATE_A(this%system_quantities)
    if (allocated(this%label)) then
      deallocate(this%label)
    end if

    POP_SUB(interaction_end)
  end subroutine interaction_end

  ! ---------------------------------------------------------
  subroutine interaction_list_add_node(this, interaction)
    class(interaction_list_t)         :: this
    class(interaction_t),      target :: interaction

    PUSH_SUB(interaction_list_add_node)

    call this%add_ptr(interaction)

    POP_SUB(interaction_list_add_node)
  end subroutine interaction_list_add_node

  ! ---------------------------------------------------------
  function interaction_iterator_get_next(this) result(interaction)
    class(interaction_iterator_t), intent(inout) :: this
    class(interaction_t),          pointer       :: interaction

    PUSH_SUB(interaction_iterator_get_next)

    select type (ptr => this%get_next_ptr())
    class is (interaction_t)
      interaction => ptr
    class default
      ASSERT(.false.)
    end select

    POP_SUB(interaction_iterator_get_next)
  end function interaction_iterator_get_next

end module interaction_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

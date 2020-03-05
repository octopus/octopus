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

#include "global.h"

module interaction_abst_oct_m
  use clock_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  implicit none

  private
  public ::               &
    interaction_abst_t,   &
    interaction_abst_end, &
    interaction_iterator_t

  !> The only purpose of the following class is to act as a surrogate and
  !> avoid circular dependencies between the interactions and the systems.
  type, abstract :: interaction_abst_t
    private

    integer, public :: n_system_observables
    integer, public :: n_partner_observables
    integer, allocatable, public :: system_observables(:)
    integer, allocatable, public :: partner_observables(:)

    type(clock_t), public :: clock
  contains
    procedure :: init_clock => interaction_init_clock
    procedure(interaction_update), deferred :: update
  end type interaction_abst_t


  !> This class extends the list iterator and adds one method to get the
  !> interaction as a pointer of type class(interaction_abst_t).
  type, extends(list_iterator_t) :: interaction_iterator_t
    private
  contains
    procedure :: get_next_interaction => interaction_iterator_get_next
  end type interaction_iterator_t

  abstract interface
    logical function interaction_update(this, clock)
      import interaction_abst_t
      import clock_t
      class(interaction_abst_t), intent(inout) :: this
      class(clock_t),            intent(in)    :: clock
    end function interaction_update
  end interface

contains

  ! ---------------------------------------------------------
  subroutine interaction_init_clock(this, namespace, dt, algo_dt)
    class(interaction_abst_t), intent(inout) :: this
    type(namespace_t),         intent(in)    :: namespace
    FLOAT,                     intent(in)    :: dt
    FLOAT,                     intent(in)    :: algo_dt

    PUSH_SUB(interaction_init_clock)

    this%clock = clock_t(namespace%get(), dt, algo_dt)

    POP_SUB(interaction_init_clock)
  end subroutine interaction_init_clock

  ! ---------------------------------------------------------
  subroutine interaction_abst_end(this)
    class(interaction_abst_t), intent(inout) :: this

    PUSH_SUB(interaction_abst_end)

    SAFE_DEALLOCATE_A(this%system_observables)
    SAFE_DEALLOCATE_A(this%partner_observables)

    POP_SUB(interaction_abst_end)

  end subroutine interaction_abst_end

  ! ---------------------------------------------------------
  function interaction_iterator_get_next(this) result(value)
    class(interaction_iterator_t), intent(inout) :: this
    class(interaction_abst_t),     pointer       :: value

    class(*), pointer :: ptr

    PUSH_SUB(interaction_iterator_get_next)

    ptr => this%get_next()
    select type (ptr)
    class is (interaction_abst_t)
      value => ptr
    class default
      ASSERT(.false.)
    end select

    POP_SUB(interaction_iterator_get_next)
  end function interaction_iterator_get_next

end module interaction_abst_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

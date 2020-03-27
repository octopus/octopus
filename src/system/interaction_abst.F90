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
  use clock_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use profiling_oct_m
  implicit none

  private
  public ::               &
    interaction_abst_t,   &
    interaction_abst_end, &
    interaction_iterator_t

  !> An interaction is a unidirectional relationship beween two systems. One of the
  !! systems owns the interaction and feels it`s effects. The other system is
  !! refered to as the interaction partner.
  type, abstract :: interaction_abst_t
    private
    !> The interaction requires access to some quantities to be evaluated, both
    !> from the system and from the partner.
    integer,              public :: n_system_quantities  !< Number of quantities needed from the system
    integer, allocatable, public :: system_quantities(:) !< Identifiers of the quantities needed from the system

    integer,              public :: n_partner_quantities !< Number of quantities needed from the partner
    integer, allocatable, public :: partner_quantities(:)!< Identifiers of the quantities needed from the parner

    type(clock_t), public :: clock !< Clock storing the time at which the interaction was last updated.
  contains
    procedure :: init_clock => interaction_init_clock
    procedure(interaction_update), deferred :: update
  end type interaction_abst_t

  !> This class extends the list iterator and adds one method to get the
  !! interaction as a pointer of type class(interaction_abst_t).
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
  subroutine interaction_init_clock(this, label, dt, algo_dt)
    class(interaction_abst_t), intent(inout) :: this
    character(len=*),          intent(in)    :: label
    FLOAT,                     intent(in)    :: dt
    FLOAT,                     intent(in)    :: algo_dt

    PUSH_SUB(interaction_init_clock)

    this%clock = clock_t(label, dt, algo_dt)

    POP_SUB(interaction_init_clock)
  end subroutine interaction_init_clock

  ! ---------------------------------------------------------
  subroutine interaction_abst_end(this)
    class(interaction_abst_t), intent(inout) :: this

    PUSH_SUB(interaction_abst_end)

    SAFE_DEALLOCATE_A(this%system_quantities)
    SAFE_DEALLOCATE_A(this%partner_quantities)

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

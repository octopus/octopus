!! Copyright (C) 2019 N. Tancogne-Dejean
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

module system_abst_oct_m
  use global_oct_m
  use interaction_abst_oct_m
  use messages_oct_m
  use namespace_oct_m
  use linked_list_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m
  use clock_oct_m

  implicit none

  private
  public ::               &
    system_abst_t

  integer, public, parameter ::        &
    TOTAL_CURRENT                = 1,  &
    FORCE                        = 2

  type, abstract :: system_abst_t
    private
    type(namespace_t),   public :: namespace

    class(propagator_abst_t), pointer, public :: prop
    type(clock_t),          public :: clock

  contains
    procedure :: dt_operation =>  system_dt_operation
    procedure :: set_propagator => system_set_propagator
    procedure :: init_clock => system_init_clock
    procedure(system_add_interaction_partner),       deferred :: add_interaction_partner
    procedure(system_has_interaction),               deferred :: has_interaction
    procedure(system_do_td_op),                      deferred :: do_td_operation
    procedure(system_update_interaction_as_partner), deferred :: update_interaction_as_partner
    procedure(system_update_interactions),           deferred :: update_interactions
    procedure(system_write_td_info),                 deferred :: write_td_info
  end type system_abst_t

  abstract interface
    subroutine system_add_interaction_partner(this, partner)
      import system_abst_t
      class(system_abst_t),     intent(inout) :: this
      class(system_abst_t),     intent(in)    :: partner
    end subroutine system_add_interaction_partner

    logical function system_has_interaction(this, interaction)
      import system_abst_t
      import interaction_abst_t
      class(system_abst_t),      intent(in) :: this
      class(interaction_abst_t), intent(in) :: interaction
    end function system_has_interaction

    subroutine system_do_td_op(this, operation)
      import system_abst_t
      class(system_abst_t), intent(inout) :: this
      integer             , intent(in)    :: operation
    end subroutine system_do_td_op

    subroutine system_update_interaction_as_partner(this, interaction)
      import system_abst_t
      import interaction_abst_t
      class(system_abst_t),      intent(in)    :: this
      class(interaction_abst_t), intent(inout) :: interaction
    end subroutine system_update_interaction_as_partner

    subroutine system_update_interactions(this)
      import system_abst_t
      class(system_abst_t),      intent(inout) :: this
    end subroutine system_update_interactions

    subroutine system_write_td_info(this)
      import system_abst_t
      class(system_abst_t), intent(in) :: this
    end subroutine system_write_td_info

  end interface

contains

  ! ---------------------------------------------------------
  subroutine system_dt_operation(this)
    class(system_abst_t),     intent(inout) :: this

    integer :: tdop

    PUSH_SUB(system_dt_operation)

    tdop = this%prop%get_td_operation()
    select case(tdop)
    case(FINISHED)
      if (debug%info) then
        message(1) = "Debug: Propagation step finished for " + trim(this%namespace%get())
        call messages_info(1)
      end if
      call this%prop%finished()
      call this%clock%increment()
      !DO OUTPUT HERE AND BROADCAST NEEDED QUANTITIES
      !ONLY IF WE ARE NOT YET FINISHED

    case(UPDATE_INTERACTIONS)
      if (debug%info) then
        message(1) = "Debug: Propagation step - Updating interactions for " + trim(this%namespace%get())
        call messages_info(1)
      end if

      call this%update_interactions()

      call this%prop%next()

    case default
      call this%do_td_operation(tdop)
    end select

    POP_SUB(system_dt_operation)
  end subroutine system_dt_operation

  ! ---------------------------------------------------------
  subroutine system_set_propagator(this, propagator)
    class(system_abst_t),             intent(inout) :: this
    class(propagator_abst_t), target, intent(in)    :: propagator

    PUSH_SUB(system_set_propagator)

    this%prop => propagator

    POP_SUB(system_set_propagator)
  end subroutine system_set_propagator

  ! ---------------------------------------------------------
  subroutine system_init_clock(this, dt, smallest_algo_dt)
    class(system_abst_t), intent(inout) :: this
    FLOAT,                intent(in)    :: dt, smallest_algo_dt

    PUSH_SUB(system_init_clock)

    this%clock = clock_t(this%namespace, dt, smallest_algo_dt)

    POP_SUB(system_init_clock)
  end subroutine system_init_clock

end module system_abst_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

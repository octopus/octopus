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
  use clock_oct_m
  use global_oct_m
  use interaction_abst_oct_m
  use messages_oct_m
  use namespace_oct_m
  use linked_list_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m
  use quantity_oct_m
  use space_oct_m
  implicit none

  private
  public ::               &
    system_abst_t

  type, abstract :: system_abst_t
    private
    type(namespace_t),   public :: namespace
    type(space_t), public :: space

    class(propagator_abst_t), pointer, public :: prop
    type(clock_t), public :: clock

    integer :: accumulated_loop_ticks

    type(linked_list_t), public :: interactions !< List we all the interactions of this system with other systems


    type(quantity_t), public :: quantities(MAX_QUANTITIES) !< Array of all possible quantities.
                                                           !< The elements of the array are accessed using the
                                                           !< quantity`s identifiers.
  contains
    procedure :: dt_operation =>  system_dt_operation
    procedure :: set_propagator => system_set_propagator
    procedure :: init_clocks => system_init_clocks
    procedure :: reset_clocks => system_reset_clocks
    procedure :: update_exposed_quantities => system_update_exposed_quantities
    procedure(system_add_interaction_partner),        deferred :: add_interaction_partner
    procedure(system_has_interaction),                deferred :: has_interaction
    procedure(system_do_td_op),                       deferred :: do_td_operation
    procedure(system_write_td_info),                  deferred :: write_td_info
    procedure(system_is_tolerance_reached),           deferred :: is_tolerance_reached
    procedure(system_store_current_status),           deferred :: store_current_status
    procedure(system_update_quantity),                deferred :: update_quantity
    procedure(system_update_exposed_quantity),        deferred :: update_exposed_quantity
    procedure(system_set_pointers_to_interaction),    deferred :: set_pointers_to_interaction
  end type system_abst_t

  abstract interface
    ! ---------------------------------------------------------
    subroutine system_add_interaction_partner(this, partner)
      import system_abst_t
      class(system_abst_t), target,    intent(inout) :: this
      class(system_abst_t),            intent(in)    :: partner
    end subroutine system_add_interaction_partner

    ! ---------------------------------------------------------
    logical function system_has_interaction(this, interaction)
      import system_abst_t
      import interaction_abst_t
      class(system_abst_t),      intent(in) :: this
      class(interaction_abst_t), intent(in) :: interaction
    end function system_has_interaction

    ! ---------------------------------------------------------
    subroutine system_do_td_op(this, operation)
      import system_abst_t
      class(system_abst_t), intent(inout) :: this
      integer             , intent(in)    :: operation
    end subroutine system_do_td_op

    ! ---------------------------------------------------------
    subroutine system_write_td_info(this)
      import system_abst_t
      class(system_abst_t), intent(in) :: this
    end subroutine system_write_td_info

    ! ---------------------------------------------------------
    logical function system_is_tolerance_reached(this, tol)
      import system_abst_t
      class(system_abst_t), intent(in) :: this
      FLOAT,                intent(in) :: tol
    end function system_is_tolerance_reached

    ! ---------------------------------------------------------
    subroutine system_store_current_status(this)
      import system_abst_t
      class(system_abst_t), intent(inout) :: this
    end subroutine system_store_current_status

    ! ---------------------------------------------------------
    logical function system_update_quantity(this, iq, clock)
      import system_abst_t
      import clock_t
      class(system_abst_t),      intent(inout) :: this
      integer,                   intent(in)    :: iq
      class(clock_t),            intent(in)    :: clock
    end function system_update_quantity

    ! ---------------------------------------------------------
    logical function system_update_exposed_quantity(this, iq, clock)
      import system_abst_t
      import clock_t
      class(system_abst_t),      intent(inout) :: this
      integer,                   intent(in)    :: iq
      class(clock_t),            intent(in)    :: clock
    end function system_update_exposed_quantity

    ! ---------------------------------------------------------
    subroutine system_set_pointers_to_interaction(this, inter)
      import system_abst_t
      import interaction_abst_t
      class(system_abst_t),   target,   intent(in)    :: this
      class(interaction_abst_t),        intent(inout) :: inter
    end subroutine system_set_pointers_to_interaction

  end interface

contains

  ! ---------------------------------------------------------
  subroutine system_dt_operation(this)
    class(system_abst_t),     intent(inout) :: this

    integer :: tdop, iq
    logical :: all_updated, obs_updated
    class(interaction_abst_t), pointer :: interaction
    type(interaction_iterator_t) :: iter

    PUSH_SUB(system_dt_operation)

    tdop = this%prop%get_td_operation()

    if (debug%info) then
      write(message(1), '(a,a,1X,a)') "Debug: ", trim(propagator_step_debug_message(tdop)), trim(this%namespace%get())
      call messages_info(1)
    end if

    select case (tdop)
    case (FINISHED)
      if (.not. this%prop%step_is_done()) then
        call this%clock%increment()
      end if
      call this%prop%finished()

    case (UPDATE_INTERACTIONS)
      !We increment by one algorithmic step
      call this%prop%clock%increment()

      !Loop over all interactions
      all_updated = .true.
      call iter%start(this%interactions)
      do while (iter%has_next())
        interaction => iter%get_next_interaction()

        ! Update the system quantities that will be needed for computing the interaction
        obs_updated = .true.
        do iq = 1, interaction%n_system_quantities
          obs_updated = this%update_quantity(interaction%system_quantities(iq), this%prop%clock) .and. obs_updated
        end do

        if (obs_updated) then
          ! Try to update the interaction
          all_updated = interaction%update(this%prop%clock) .and. all_updated
        else
          ! We are not able to update the interaction
          all_updated = .false.
          if (debug%info) then
            write(message(1), '(a,a)') "Debug: -- Internal quantities are not up-to-date, so could not update interactions for ", &
              trim(this%namespace%get())
            call messages_info(1)
          end if
        end if
      end do

      !Move to next propagator step if all interactions have been
      !updated. Otherwise try again later.
      if(all_updated) then
        this%accumulated_loop_ticks = this%accumulated_loop_ticks + 1
        call this%prop%next()
      else
        call this%prop%clock%decrement()
      end if

    case (START_SCF_LOOP)
      ASSERT(this%prop%predictor_corrector)
 
      call this%prop%save_scf_start()
      this%accumulated_loop_ticks = 0

      if (debug%info) then
        write(message(1), '(a,i3,a)') "Debug: -- SCF iter ", this%prop%scf_count, " for " + trim(this%namespace%get())
        call messages_info(1)
      end if

    case (END_SCF_LOOP)
      !Here we first check if we did the maximum number of steps.
      !Otherwise, we need check the tolerance 
      if(this%prop%scf_count == this%prop%max_scf_count) then
        if (debug%info) then
          message(1) = "Debug: -- Max SCF Iter reached for " + trim(this%namespace%get())
          call messages_info(1)
        end if
        call this%prop%next()
      else
        !We reset the pointer to the begining of the scf loop
        if(this%is_tolerance_reached(this%prop%scf_tol)) then
          if (debug%info) then
            message(1) = "Debug: -- SCF tolerance reached for " + trim(this%namespace%get())
            call messages_info(1)
          end if
          call this%prop%next()
        else
          !We rewind the instruction stack
          call this%prop%rewind_scf_loop()

          !We reset the clocks
          call this%reset_clocks(this%accumulated_loop_ticks)
          this%accumulated_loop_ticks = 0
          if (debug%info) then
            write(message(1), '(a,i3,a)') "Debug: -- SCF iter ", this%prop%scf_count, " for " + trim(this%namespace%get())
           call messages_info(1)
         end if 
        end if
      end if

    case (STORE_CURRENT_STATUS)
      call this%store_current_status()
      call this%prop%next()

    case default
      call this%do_td_operation(tdop)
      call this%prop%next()
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
  subroutine system_init_clocks(this, dt, smallest_algo_dt)
    class(system_abst_t), intent(inout) :: this
    FLOAT,                intent(in)    :: dt, smallest_algo_dt

    integer :: iq
    type(interaction_iterator_t) :: iter
    class(interaction_abst_t), pointer :: interaction

    PUSH_SUB(system_init_clocks)

    ! Internal clock
    this%clock = clock_t(this%namespace%get(), dt, smallest_algo_dt)

    ! Propagator clock
    this%prop%clock = clock_t(this%namespace%get(), dt/this%prop%algo_steps, smallest_algo_dt)

    ! Interaction clocks
    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next_interaction()
      call interaction%init_clock(this%namespace%get(), dt, smallest_algo_dt)
    end do

    ! Internal quantities clocks
    where (this%quantities%internal)
      this%quantities%clock = clock_t(this%namespace%get(), dt, smallest_algo_dt)
    end where

    POP_SUB(system_init_clocks)
  end subroutine system_init_clocks

  ! ---------------------------------------------------------
  subroutine system_reset_clocks(this, accumulated_ticks)
    class(system_abst_t),      intent(inout) :: this
    integer,                   intent(in)    :: accumulated_ticks

    integer :: it, iq
    type(interaction_iterator_t) :: iter
    class(interaction_abst_t), pointer :: interaction

    PUSH_SUB(system_reset_clocks)

    do it = 1, accumulated_ticks
      ! Propagator clock
      call this%prop%clock%decrement()

      ! Interaction clocks
      call iter%start(this%interactions)
      do while (iter%has_next())
        interaction => iter%get_next_interaction()
        call interaction%clock%decrement()
      end do

      ! Internal quantities clocks
      do iq = 1, MAX_QUANTITIES
        if (this%quantities(iq)%internal) call this%quantities(iq)%clock%decrement()
      end do
    end do

    POP_SUB(system_reset_clocks)
  end subroutine system_reset_clocks

  ! ---------------------------------------------------------
  logical function system_update_exposed_quantities(this, clock, n_quantities, quantities) result(updated)
    class(system_abst_t),      intent(inout) :: this
    type(clock_t),             intent(in)    :: clock
    integer,                   intent(in)    :: n_quantities
    integer,                   intent(in)    :: quantities(:)

    integer :: iq

    if (this%clock < clock .and. this%clock%is_earlier_with_step(clock)) then
      !This is not the best moment to update the quantities
      updated = .false.
    else
      !This is the best moment to update the quantities
      updated = .true.
      do iq = 1, n_quantities
        updated = this%update_exposed_quantity(quantities(iq), clock) .and. updated
      end do
    end if

  end function system_update_exposed_quantities

end module system_abst_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

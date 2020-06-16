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
  use interaction_partner_oct_m
  use interaction_with_partner_oct_m
  use messages_oct_m
  use namespace_oct_m
  use linked_list_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m
  use propagator_beeman_oct_m
  use propagator_exp_mid_oct_m
  use propagator_verlet_oct_m
  use quantity_oct_m
  use space_oct_m
  use varinfo_oct_m
  implicit none

  private
  public ::               &
    system_abst_t,        &
    system_list_t,        &
    system_iterator_t

  type, extends(interaction_partner_t), abstract :: system_abst_t
    private
    type(space_t), public :: space

    class(propagator_abst_t), pointer, public :: prop

    integer :: accumulated_loop_ticks

    integer :: interaction_quantity_timing  !< parameter to determine if interactions
      !< should use the quantities at the exact time or if retardation is allowed

    type(interaction_list_t), public :: interactions !< List with all the interactions of this system
  contains
    procedure :: dt_operation =>  system_dt_operation
    procedure :: init_clocks => system_init_clocks
    procedure :: reset_clocks => system_reset_clocks
    procedure :: update_exposed_quantities => system_update_exposed_quantities
    procedure :: init_propagator => system_init_propagator
    procedure :: update_interactions => system_update_interactions
    procedure :: propagation_start => system_propagation_start
    procedure :: propagation_finish => system_propagation_finish
    procedure :: has_reached_final_propagation_time => system_has_reached_final_propagation_time
    procedure :: propagation_step_finish => system_propagation_step_finish
    procedure :: propagation_step_is_done => system_propagation_step_is_done
    procedure(system_add_interaction_partner),        deferred :: add_interaction_partner
    procedure(system_has_interaction),                deferred :: has_interaction
    procedure(system_initial_conditions),             deferred :: initial_conditions
    procedure(system_do_td_op),                       deferred :: do_td_operation
    procedure(system_iteration_info),                 deferred :: iteration_info
    procedure(system_is_tolerance_reached),           deferred :: is_tolerance_reached
    procedure(system_store_current_status),           deferred :: store_current_status
    procedure(system_update_quantity),                deferred :: update_quantity
    procedure(system_update_interactions_start),      deferred :: update_interactions_start
    procedure(system_update_interactions_finish),     deferred :: update_interactions_finish
    procedure(system_output_start),                   deferred :: output_start
    procedure(system_output_write),                   deferred :: output_write
    procedure(system_output_finish),                  deferred :: output_finish
  end type system_abst_t

  abstract interface
    ! ---------------------------------------------------------
    subroutine system_add_interaction_partner(this, partner)
      import system_abst_t
      class(system_abst_t), target,    intent(inout) :: this
      class(system_abst_t),            intent(inout) :: partner
    end subroutine system_add_interaction_partner

    ! ---------------------------------------------------------
    logical function system_has_interaction(this, interaction)
      import system_abst_t
      import interaction_abst_t
      class(system_abst_t),      intent(in) :: this
      class(interaction_abst_t), intent(in) :: interaction
    end function system_has_interaction

    ! ---------------------------------------------------------
    subroutine system_initial_conditions(this, from_scratch)
      import system_abst_t
      class(system_abst_t), intent(inout) :: this
      logical,              intent(in)    :: from_scratch
    end subroutine system_initial_conditions

    ! ---------------------------------------------------------
    subroutine system_do_td_op(this, operation)
      import system_abst_t
      class(system_abst_t), intent(inout) :: this
      integer             , intent(in)    :: operation
    end subroutine system_do_td_op

    ! ---------------------------------------------------------
    subroutine system_iteration_info(this)
      import system_abst_t
      class(system_abst_t), intent(in) :: this
    end subroutine system_iteration_info

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
    subroutine system_update_quantity(this, iq, requested_time)
      import system_abst_t
      import clock_t
      class(system_abst_t),      intent(inout) :: this
      integer,                   intent(in)    :: iq
      class(clock_t),            intent(in)    :: requested_time
    end subroutine system_update_quantity

    ! ---------------------------------------------------------
    subroutine system_update_interactions_start(this)
      import system_abst_t
      class(system_abst_t), intent(inout) :: this
    end subroutine system_update_interactions_start

    ! ---------------------------------------------------------
    subroutine system_update_interactions_finish(this)
      import system_abst_t
      class(system_abst_t), intent(inout) :: this
    end subroutine system_update_interactions_finish

    ! ---------------------------------------------------------
    subroutine system_output_start(this)
      import system_abst_t
      class(system_abst_t), intent(inout) :: this
    end subroutine system_output_start

    ! ---------------------------------------------------------
    subroutine system_output_write(this, iter)
      import system_abst_t
      class(system_abst_t), intent(inout) :: this
      integer,              intent(in)    :: iter
    end subroutine system_output_write

    ! ---------------------------------------------------------
    subroutine system_output_finish(this)
      import system_abst_t
      class(system_abst_t), intent(inout) :: this
    end subroutine system_output_finish
  end interface

  !> These class extends the list and list iterator to create a system list
  type, extends(linked_list_t) :: system_list_t
    private
  contains
    procedure :: add => system_list_add_node
  end type system_list_t
  
  type, extends(linked_list_iterator_t) :: system_iterator_t
    private
  contains
    procedure :: get_next => system_iterator_get_next
  end type system_iterator_t

contains

  ! ---------------------------------------------------------
  subroutine system_dt_operation(this)
    class(system_abst_t),     intent(inout) :: this

    integer :: tdop
    logical :: all_updated

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
      ! We increment by one algorithmic step
      call this%prop%clock%increment()

      ! Try to update all the interactions
      all_updated = this%update_interactions(this%prop%clock)

      ! Move to next propagator step if all interactions have been
      ! updated. Otherwise try again later.
      if(all_updated) then
        this%accumulated_loop_ticks = this%accumulated_loop_ticks + 1
        call this%prop%next()
      else
        call this%prop%clock%decrement()
      end if

    case (START_SCF_LOOP)
      ASSERT(this%prop%predictor_corrector)

      call this%prop%save_scf_start()
      this%prop%inside_scf = .true.
      this%accumulated_loop_ticks = 0

      if (debug%info) then
        write(message(1), '(a,i3,a)') "Debug: -- SCF iter ", this%prop%scf_count, " for " + trim(this%namespace%get())
        call messages_info(1)
      end if

    case (END_SCF_LOOP)
      ! Here we first check if we did the maximum number of steps.
      ! Otherwise, we need check the tolerance
      if(this%prop%scf_count == this%prop%max_scf_count) then
        if (debug%info) then
          message(1) = "Debug: -- Max SCF Iter reached for " + trim(this%namespace%get())
          call messages_info(1)
        end if
        this%prop%inside_scf = .false.
        call this%prop%next()
      else
        ! We reset the pointer to the begining of the scf loop
        if(this%is_tolerance_reached(this%prop%scf_tol)) then
          if (debug%info) then
            message(1) = "Debug: -- SCF tolerance reached for " + trim(this%namespace%get())
            call messages_info(1)
          end if
          this%prop%inside_scf = .false.
          call this%prop%next()
        else
          ! We rewind the instruction stack
          call this%prop%rewind_scf_loop()

          ! We reset the clocks
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
  subroutine system_init_clocks(this, smallest_algo_dt)
    class(system_abst_t), intent(inout) :: this
    FLOAT,                intent(in)    :: smallest_algo_dt

    type(interaction_iterator_t) :: iter
    class(interaction_abst_t), pointer :: interaction

    PUSH_SUB(system_init_clocks)

    ! Internal clock
    this%clock = clock_t(this%namespace%get(), this%prop%dt, smallest_algo_dt)

    ! Propagator clock
    this%prop%clock = clock_t(this%namespace%get(), this%prop%dt/this%prop%algo_steps, smallest_algo_dt)

    ! Interaction clocks
    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next()
      call interaction%init_clock(this%namespace%get(), this%prop%dt, smallest_algo_dt)
    end do

    ! Required quantities clocks
    where (this%quantities%required)
      this%quantities%clock = clock_t(this%namespace%get(), this%prop%dt/this%prop%algo_steps, smallest_algo_dt)
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
        interaction => iter%get_next()
        call interaction%clock%decrement()
      end do

      ! Internal quantities clocks
      do iq = 1, MAX_QUANTITIES
        if (this%quantities(iq)%required) call this%quantities(iq)%clock%decrement()
      end do
    end do

    POP_SUB(system_reset_clocks)
  end subroutine system_reset_clocks

  ! ---------------------------------------------------------
  ! this function is called as partner from the interaction
  logical function system_update_exposed_quantities(partner, requested_time, interaction) result(allowed_to_update)
    class(system_abst_t),      intent(inout) :: partner
    type(clock_t),             intent(in)    :: requested_time
    class(interaction_abst_t), intent(inout) :: interaction

    logical :: ahead_in_time, behind_in_time, right_in_time, need_to_copy
    integer :: iq, q_id

    PUSH_SUB(system_update_exposed_quantities)

    if (debug%info) then
      write(message(1), '(a,a)') "Debug: ----- updating exposed quantities for partner ", trim(partner%namespace%get())
      call messages_info(1)
    end if

    select type (interaction)
    class is (interaction_with_partner_t)

      if(partner%prop%inside_scf) then
        ! we are inside an SCF cycle and therefore are not allowed to expose any quantities.
        allowed_to_update = .false.
      else
        allowed_to_update = .true.
        need_to_copy = .true.
        do iq = 1, interaction%n_partner_quantities
          ! Get the requested quantity ID
          q_id = interaction%partner_quantities(iq)

          ! All needed quantities must have been marked as required. If not, then fix your code!
          ASSERT(partner%quantities(q_id)%required)

          ! First update the exposed quantities that are not protected
          if (.not.partner%quantities(q_id)%protected) then
            if (partner%quantities(q_id)%clock%is_later_with_step(requested_time)) then
              ! We can update because the partner will reach this time in its current timestep
              ! This is not a protected quantity, so we update it
              call partner%update_exposed_quantity(q_id, requested_time)
            end if
          end if

          ! Now compare the times
          ahead_in_time = partner%quantities(q_id)%clock > requested_time
          right_in_time = partner%quantities(q_id)%clock == requested_time

          if(partner%interaction_quantity_timing == OPTION__INTERACTIONQUANTITYTIMING__TIMING_EXACT) then
            ! only allow interaction at exactly the same time
            allowed_to_update = allowed_to_update .and. right_in_time
            need_to_copy = allowed_to_update
          else if(partner%interaction_quantity_timing == OPTION__INTERACTIONQUANTITYTIMING__TIMING_RETARDED) then
            ! allow retarded interaction
            allowed_to_update = allowed_to_update .and. &
              (right_in_time .or. ahead_in_time)
            need_to_copy = need_to_copy .and. .not. ahead_in_time
          else
            call messages_not_implemented("Method for interaction quantity timing")
          end if

          ! Debug stuff
          if (debug%info) then
            write(message(1), '(a,i3)') "Debug: ------ updating exposed quantities ", q_id
            write(message(2), '(a,i3,a,i3)') "Debug: ------ requested time is ", requested_time%get_tick(), &
              " and partner time is ", partner%quantities(q_id)%clock%get_tick()
            call messages_info(2)
          end if

        end do

        ! If the quantities have been updated, we copy them to the interaction
        if (allowed_to_update .and. need_to_copy) then
          call partner%copy_quantities_to_interaction(interaction)
        end if
      end if

    class default
      message(1) = "A system can only expose quantities to an interaction as a partner."
      call messages_fatal(1, namespace=partner%namespace)
    end select


    POP_SUB(system_update_exposed_quantities)
  end function system_update_exposed_quantities

  ! ---------------------------------------------------------
  logical function system_update_interactions(this, requested_time) result(all_updated)
    class(system_abst_t),      intent(inout) :: this
    type(clock_t),             intent(in)    :: requested_time !< Requested time for the update

    integer :: iq, q_id
    class(interaction_abst_t), pointer :: interaction
    type(interaction_iterator_t) :: iter

    PUSH_SUB(system_update_interactions)

    ! Some systems might need to perform some specific operations before the update
    call this%update_interactions_start()

    !Loop over all interactions
    all_updated = .true.
    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next()

      if (.not. interaction%clock == requested_time) then
        ! Update the system quantities that will be needed for computing the interaction
        do iq = 1, interaction%n_system_quantities
          ! Get requested quantity ID
          q_id = interaction%system_quantities(iq)

          ! All needed quantities must have been marked as required. If not, then fix your code!
          ASSERT(this%quantities(q_id)%required)

          ! We do not need to update the protected quantities, the propagator takes care of that
          if (this%quantities(q_id)%protected) cycle

          if (.not. this%quantities(q_id)%clock == requested_time) then
            ! The requested quantity is not at the requested time, so we try to update it

            ! Sanity check: it should never happen that the quantity is in advance
            ! with respect to the requested time.
            if (this%quantities(q_id)%clock > requested_time) then
              message(1) = "The quantity clock is in advance compared to the requested time."
              call messages_fatal(1, namespace=this%namespace)
            end if

            call this%update_quantity(q_id, requested_time)
          end if

        end do

        ! We can now try to update the interaction
        all_updated = interaction%update(this%namespace, requested_time) .and. all_updated
      end if
    end do

    ! Some systems might need to perform some specific operations after the update
    call this%update_interactions_finish()

    POP_SUB(system_update_interactions)
  end function system_update_interactions

  ! ---------------------------------------------------------
  subroutine system_init_propagator(this, smallest_algo_dt)
    class(system_abst_t),      intent(inout) :: this
    FLOAT,                     intent(inout) :: smallest_algo_dt

    integer :: prop

    PUSH_SUB(system_init_propagator)

    call messages_experimental('Multisystem propagator framework')

    !%Variable TDSystemPropagator
    !%Type integer
    !%Default verlet
    !%Section Time-Dependent::Propagation
    !%Description
    !% A variable to set the propagator in the multisystem framework.
    !% This is a temporary solution, and should be replaced by the
    !% TDPropagator variable.
    !%Option verlet 1
    !% (Experimental) Verlet propagator.
    !%Option beeman 2
    !% (Experimental) Beeman propagator without predictor-corrector.
    !%Option beeman_scf 3
    !% (Experimental) Beeman propagator with predictor-corrector scheme.
    !%Option exp_mid 4
    !% (Experimental) Exponential midpoint propagator without predictor-corrector.
    !%Option exp_mid_scf 5
    !% (Experimental) Exponential midpoint propagator with predictor-corrector scheme.
    !%End
    call parse_variable(this%namespace, 'TDSystemPropagator', PROP_VERLET, prop)
    if(.not.varinfo_valid_option('TDSystemPropagator', prop)) call messages_input_error(this%namespace, 'TDSystemPropagator')
    call messages_print_var_option(stdout, 'TDSystemPropagator', prop)

    select case(prop)
    case(PROP_VERLET)
      this%prop => propagator_verlet_t(this%namespace)
    case(PROP_BEEMAN)
      this%prop => propagator_beeman_t(this%namespace, .false.)
    case(PROP_BEEMAN_SCF)
      this%prop => propagator_beeman_t(this%namespace, .true.)
    case(PROP_EXPMID)
      this%prop => propagator_exp_mid_t(this%namespace, .false.)
    case(PROP_EXPMID_SCF)
      this%prop => propagator_exp_mid_t(this%namespace, .true.)
    end select

    call this%prop%rewind()

    ! Check if this propagators dt is smaller then the current smallest dt.
    ! If so, replace the current smallest dt by the one from this propagator.
    smallest_algo_dt = min(smallest_algo_dt, this%prop%dt/this%prop%algo_steps)

    !%Variable InteractionQuantityTiming
    !%Type integer
    !%Default exact
    !%Section Time-Dependent::Propagation
    !%Description
    !% A parameter to determine if interactions should use the quantities
    !% at the exact time or if retardation is allowed.
    !%Option timing_exact 1
    !% Only allow interactions at exactly the same times
    !%Option timing_retarded 2
    !% Allow retarded interactions
    !%End
    call parse_variable(this%namespace, 'InteractionQuantityTiming', &
      OPTION__INTERACTIONQUANTITYTIMING__TIMING_EXACT, &
      this%interaction_quantity_timing)
    if(.not.varinfo_valid_option('InteractionQuantityTiming', this%interaction_quantity_timing)) then
      call messages_input_error(this%namespace, 'InteractionQuantityTiming')
    end if
    call messages_print_var_option(stdout, 'InteractionQuantityTiming', &
      this%interaction_quantity_timing)

    POP_SUB(system_init_propagator)
  end subroutine system_init_propagator

  ! ---------------------------------------------------------
  subroutine system_propagation_start(this)
    class(system_abst_t),      intent(inout) :: this

    logical :: all_updated

    PUSH_SUB(system_propagation_start)

    ! Update interactions at initial time
    all_updated = this%update_interactions(this%clock)
    if (.not. all_updated) then
      message(1) = "Unable to update interactions when initializing the propagation."
      call messages_fatal(1, namespace=this%namespace)
    end if

    ! System-specific and propagator-specific initialization step
    call this%do_td_operation(this%prop%start_step)

    ! Start output
    call this%output_start()

    ! Write information for first iteration
    call this%iteration_info()

    POP_SUB(system_propagation_start)
  end subroutine system_propagation_start

  ! ---------------------------------------------------------
  subroutine system_propagation_finish(this)
    class(system_abst_t),      intent(inout) :: this

    PUSH_SUB(system_propagation_finish)

    ! Finish output
    call this%output_finish()

    ! System-specific and propagator-specific finalization step
    call this%do_td_operation(this%prop%final_step)

    POP_SUB(system_propagation_finish)
  end subroutine system_propagation_finish

  ! ---------------------------------------------------------
  logical function system_has_reached_final_propagation_time(this)
    class(system_abst_t),      intent(inout) :: this

    PUSH_SUB(system_has_reached_final_propagation_time)

    ! Fixme: should be changed to final propagation time
    system_has_reached_final_propagation_time = (this%clock%get_sim_time() >= this%prop%max_td_steps*this%prop%dt)

    POP_SUB(system_has_reached_final_propagation_time)
  end function system_has_reached_final_propagation_time

  ! ---------------------------------------------------------
  subroutine system_propagation_step_finish(this, iteration)
    class(system_abst_t),      intent(inout) :: this
    integer,                   intent(in)    :: iteration

    PUSH_SUB(system_propagation_step_finish)

    ! Print information about the current iteration and write output
    call this%output_write(iteration)
    call this%iteration_info()

    ! Reset propagator for next step
    call this%prop%rewind()

    POP_SUB(system_propagation_step_finish)
  end subroutine system_propagation_step_finish

  ! ---------------------------------------------------------
  logical function system_propagation_step_is_done(this)
    class(system_abst_t),      intent(inout) :: this

    PUSH_SUB(system_propagation_step_is_done)

    system_propagation_step_is_done = this%prop%step_is_done()

    POP_SUB(system_propagation_step_is_done)
  end function system_propagation_step_is_done

  ! ---------------------------------------------------------
  subroutine system_list_add_node(this, system)
    class(system_list_t)         :: this
    class(system_abst_t), target :: system

    PUSH_SUB(system_list_add_node)

    call this%add_ptr(system)

    POP_SUB(system_list_add_node)
  end subroutine system_list_add_node

  ! ---------------------------------------------------------
  function system_iterator_get_next(this) result(system)
    class(system_iterator_t), intent(inout) :: this
    class(system_abst_t),     pointer       :: system

    PUSH_SUB(system_iterator_get_next)

    select type (ptr => this%get_next_ptr())
    class is (system_abst_t)
      system => ptr
    class default
      ASSERT(.false.)
    end select

    POP_SUB(system_iterator_get_next)
  end function system_iterator_get_next

end module system_abst_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

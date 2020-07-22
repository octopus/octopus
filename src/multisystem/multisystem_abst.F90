!! Copyright (C) 2019-2020 M. Oliveira
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

module multisystem_abst_oct_m
  use clock_oct_m
  use global_oct_m
  use ghost_interaction_oct_m
  use interaction_oct_m
  use interaction_with_partner_oct_m
  use io_oct_m
  use loct_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use system_oct_m
  use system_replica_oct_m
  use system_factory_abst_oct_m
  implicit none

  private
  public ::               &
    multisystem_abst_t,   &
    multisystem_abst_init

  type, extends(system_t), abstract :: multisystem_abst_t
    type(system_list_t) :: list
  contains
    procedure :: dt_operation =>  multisystem_abst_dt_operation
    procedure :: init_clocks => multisystem_abst_init_clocks
    procedure :: reset_clocks => multisystem_abst_reset_clocks
    procedure :: init_propagator => multisystem_abst_init_propagator
    procedure :: propagation_start => multisystem_abst_propagation_start
    procedure :: propagation_finish => multisystem_abst_propagation_finish
    procedure :: has_reached_final_propagation_time => multisystem_abst_has_reached_final_propagation_time
    procedure :: propagation_step_finish => multisystem_abst_propagation_step_finish
    procedure :: propagation_step_is_done => multisystem_abst_propagation_step_is_done
    procedure :: init_all_interactions => multisystem_abst_init_all_interactions
    procedure :: init_interaction => multisystem_abst_init_interaction
    procedure :: write_interaction_graph => multisystem_abst_write_interaction_graph
    procedure :: initial_conditions => multisystem_abst_initial_conditions
    procedure :: do_td_operation => multisystem_abst_do_td_operation
    procedure :: iteration_info => multisystem_abst_iteration_info
    procedure :: output_start => multisystem_abst_output_start
    procedure :: output_write => multisystem_abst_output_write
    procedure :: output_finish => multisystem_abst_output_finish
    procedure :: is_tolerance_reached => multisystem_abst_is_tolerance_reached
    procedure :: store_current_status => multisystem_abst_store_current_status
    procedure :: update_quantity => multisystem_abst_update_quantity
    procedure :: update_exposed_quantity => multisystem_abst_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => multisystem_abst_copy_quantities_to_interaction
  end type multisystem_abst_t

contains


  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_abst_init(system, namespace, factory, system_replica)
    class(multisystem_abst_t),         intent(inout) :: system
    type(namespace_t),            intent(in) :: namespace
    class(system_factory_abst_t), intent(in) :: factory
    type(system_replica_t),       intent(inout) :: system_replica

    integer :: isys, system_type
    character(len=128) :: system_name, replica_name
    type(block_t) :: blk
    class(system_t), pointer :: sys
    integer :: system_replicas_default, jj

    PUSH_SUB(multisystem_abst_init)

    system%namespace = namespace
    system%system_replica = system_replica

    !%Variable SystemReplicas
    !%Type integer
    !%Section Time-Dependent::Propagation
    !%Description
    !% Number of system replicas
    !%End

    system_replicas_default = 0
    call parse_variable(system%namespace, 'SystemReplicas', system_replicas_default, system%system_replica%n_replicas)
    write(message(1), '(a,a,a,i6)') 'Namespace: ', trim(system%namespace%get()), ' SystemReplicas:', &
            system%system_replica%n_replicas
    call messages_info(1)


    !%Variable Systems
    !%Type block
    !%Section System
    !%Description
    !% List of systems that will be treated in the calculation.
    !% The first column should be a string containing the system name.
    !% The second column should be the system type. See below for a list of
    !% available system types.
    !%Option electronic 1
    !% An electronic system. (NOT IMPLEMENTED)
    !%Option maxwell 2
    !% A maxwell system.
    !%Option classical_particle 3
    !% A classical particle. Used for testing purposes only.
    !%Option charged_particle 4
    !% A charged classical particle.
    !%Option multisystem 5
    !% A system containing other systems.
    !%End
    if (parse_block(system%namespace, 'Systems', blk) == 0) then

      do isys = 1, parse_block_n(blk)
        ! Parse system name and type
        call parse_block_string(blk, isys - 1, 0, system_name)
        if (len_trim(system_name) == 0) then
          call messages_input_error(system%namespace, 'Systems', 'All systems must have a name.')
        end if
        call parse_block_integer(blk, isys - 1, 1, system_type)

        ! Create folder to store system files.
        ! Needs to be done before creating the system as this in turn might create subfolders.
        call io_mkdir(system_name, namespace=system%namespace)

        ! Create system
        sys => factory%create(system%namespace, system_name, system_type, system%system_replica)
        if (.not. associated(sys)) then
          call messages_input_error(system%namespace, 'Systems', 'Unknown system type.')
        end if

        ! Add system to list of systems
        call system%list%add(sys)

        ! Add replicas to the list of systems
        do jj = 1, system%system_replica%n_replicas
           write(replica_name,'(a,a,i8.8)') trim(system_name), '-', jj
           call io_mkdir(replica_name, namespace=system%namespace)
           system%system_replica%is_replica = .true.
           sys => factory%create(system%namespace, replica_name, system_type, system%system_replica)

           if (.not. associated(sys)) then
             call messages_input_error(system%namespace, 'Systems', 'Unknown system type.')
           end if
           call system%list%add(sys)
        end do

      end do
      call parse_block_end(blk)
    else
      message(1) = "Input error while reading block Systems."
      call messages_fatal(1, namespace=system%namespace)
    end if

    POP_SUB(multisystem_abst_init)
  end subroutine multisystem_abst_init

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_abst_dt_operation(this)
    class(multisystem_abst_t),     intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_dt_operation)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%dt_operation()
    end do

    POP_SUB(multisystem_abst_dt_operation)
  end subroutine multisystem_abst_dt_operation

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_abst_init_clocks(this, smallest_algo_dt)
    class(multisystem_abst_t), intent(inout) :: this
    FLOAT,                intent(in)    :: smallest_algo_dt

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_init_clocks)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%init_clocks(smallest_algo_dt)
    end do

    POP_SUB(multisystem_abst_init_clocks)
  end subroutine multisystem_abst_init_clocks

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_abst_reset_clocks(this, accumulated_ticks)
    class(multisystem_abst_t),      intent(inout) :: this
    integer,                   intent(in)    :: accumulated_ticks

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_reset_clocks)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%reset_clocks(accumulated_ticks)
    end do

    POP_SUB(multisystem_abst_reset_clocks)
  end subroutine multisystem_abst_reset_clocks

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_abst_init_propagator(this, smallest_algo_dt)
    class(multisystem_abst_t),      intent(inout) :: this
    FLOAT,                     intent(inout) :: smallest_algo_dt

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_init_propagator)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%init_propagator(smallest_algo_dt)
    end do

    POP_SUB(multisystem_abst_init_propagator)
  end subroutine multisystem_abst_init_propagator

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_abst_propagation_start(this)
    class(multisystem_abst_t),      intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_propagation_start)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%propagation_start()
    end do

    POP_SUB(multisystem_abst_propagation_start)
  end subroutine multisystem_abst_propagation_start

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_abst_propagation_finish(this)
    class(multisystem_abst_t),      intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_propagation_finish)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%propagation_finish()
    end do

    POP_SUB(multisystem_abst_propagation_finish)
  end subroutine multisystem_abst_propagation_finish

  ! ---------------------------------------------------------------------------------------
  recursive logical function multisystem_abst_has_reached_final_propagation_time(this)
    class(multisystem_abst_t),      intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_has_reached_final_propagation_time)

    multisystem_abst_has_reached_final_propagation_time = .true.
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      multisystem_abst_has_reached_final_propagation_time = multisystem_abst_has_reached_final_propagation_time .and. &
        system%has_reached_final_propagation_time()
    end do

    POP_SUB(multisystem_abst_has_reached_final_propagation_time)
  end function multisystem_abst_has_reached_final_propagation_time

  ! ---------------------------------------------------------
  recursive subroutine multisystem_abst_propagation_step_finish(this, iteration)
    class(multisystem_abst_t),      intent(inout) :: this
    integer,                   intent(in)    :: iteration

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_propagation_step_finish)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      if (system%propagation_step_is_done()) then
        call system%propagation_step_finish(iteration)
      end if
    end do

    POP_SUB(multisystem_abst_propagation_step_finish)
  end subroutine multisystem_abst_propagation_step_finish

  ! ---------------------------------------------------------------------------------------
  recursive logical function multisystem_abst_propagation_step_is_done(this)
    class(multisystem_abst_t),      intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_propagation_step_is_done)

    multisystem_abst_propagation_step_is_done = .false.
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      multisystem_abst_propagation_step_is_done = multisystem_abst_propagation_step_is_done .or. system%propagation_step_is_done()
    end do

    POP_SUB(multisystem_abst_propagation_step_is_done)
  end function multisystem_abst_propagation_step_is_done

  ! ---------------------------------------------------------
  recursive subroutine multisystem_abst_init_all_interactions(this)
    class(multisystem_abst_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter_i
    class(interaction_t), pointer :: interaction
    type(system_iterator_t) :: iter_s
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_init_all_interactions)

    ! Initialize interactions directly owned by the multisystem
    call iter_i%start(this%interactions)
    do while (iter_i%has_next())
      interaction => iter_i%get_next()
      select type (interaction)
      type is (ghost_interaction_t)
        ! Skip the ghost interactions
      class default
        call this%init_interaction(interaction)
      end select
    end do

    ! Initialize interactions owned by the subsystems
    call iter_s%start(this%list)
    do while (iter_s%has_next())
      system => iter_s%get_next()
      call system%init_all_interactions()
    end do

    POP_SUB(multisystem_abst_init_all_interactions)
  end subroutine multisystem_abst_init_all_interactions

  ! ---------------------------------------------------------
  subroutine multisystem_abst_init_interaction(this, interaction)
    class(multisystem_abst_t), target, intent(inout) :: this
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(multisystem_abst_init_interaction)

    ! The multitystem class should never know about any specific interaction.
    ! Only classes that extend it can know about specific interactions.
    ! Such classes should override this method to add new supported interactions.
    message(1) = "Trying to initialize an interaction in the multisystem class"
    call messages_fatal(1)

    POP_SUB(multisystem_abst_init_interaction)
  end subroutine multisystem_abst_init_interaction

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_abst_write_interaction_graph(this, iunit)
    class(multisystem_abst_t), intent(in) :: this
    integer,              intent(in) :: iunit

    class(system_t), pointer :: system
    class(interaction_t), pointer :: interaction
    type(system_iterator_t) :: sys_iter
    type(interaction_iterator_t) :: inter_iter

    PUSH_SUB(multisystem_abst_write_interaction_graph)

    ! Loop over all the subsystems
    call sys_iter%start(this%list)
    do while (sys_iter%has_next())
      system => sys_iter%get_next()

      ! Loop over the interactions that this subsystem has
      call inter_iter%start(system%interactions)
      do while (inter_iter%has_next())
        interaction => inter_iter%get_next()

        ! Write interaction to DOT graph if this interaction has a partner
        select type (interaction)
        type is (ghost_interaction_t)
          ! Do not include systems connected by ghost interactions
        class is (interaction_with_partner_t)
          write(iunit, '(2x,a)') '"' + trim(system%namespace%get()) + '" -> "' + trim(interaction%partner%namespace%get()) + &
            '" [label="'+ interaction%label + '"];'
        end select
      end do

      ! If this subsystem is also a multisystem, then we also need to traverse it
      select type (system)
      class is (multisystem_abst_t)
        call system%write_interaction_graph(iunit)
      end select
    end do

    POP_SUB(multisystem_abst_write_interaction_graph)
  end subroutine multisystem_abst_write_interaction_graph

  ! ---------------------------------------------------------
  recursive subroutine multisystem_abst_initial_conditions(this, from_scratch)
    class(multisystem_abst_t), intent(inout) :: this
    logical,              intent(in)    :: from_scratch

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_initial_conditions)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%initial_conditions(from_scratch)
    end do

    POP_SUB(multisystem_abst_initial_conditions)
  end subroutine multisystem_abst_initial_conditions

  ! ---------------------------------------------------------
  recursive subroutine multisystem_abst_do_td_operation(this, operation)
    class(multisystem_abst_t), intent(inout) :: this
    integer,              intent(in)    :: operation

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_do_td_operation)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%do_td_operation(operation)
    end do

    POP_SUB(multisystem_abst_do_td_operation)
  end subroutine multisystem_abst_do_td_operation

  ! ---------------------------------------------------------
  recursive logical function multisystem_abst_is_tolerance_reached(this, tol) result(converged)
    class(multisystem_abst_t), intent(in)    :: this
    FLOAT,                intent(in)    :: tol

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_is_tolerance_reached)

    converged = .true.
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      if (.not. system%is_tolerance_reached(tol)) converged = .false.
    end do

    POP_SUB(multisystem_abst_is_tolerance_reached)
  end function multisystem_abst_is_tolerance_reached

  ! ---------------------------------------------------------
  recursive subroutine multisystem_abst_store_current_status(this)
    class(multisystem_abst_t), intent(inout)    :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_store_current_status)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%store_current_status()
    end do

    POP_SUB(multisystem_abst_store_current_status)
  end subroutine multisystem_abst_store_current_status

   ! ---------------------------------------------------------
  recursive subroutine multisystem_abst_iteration_info(this)
    class(multisystem_abst_t), intent(in) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_iteration_info)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%iteration_info()
    end do

    POP_SUB(multisystem_abst_iteration_info)
  end subroutine multisystem_abst_iteration_info

  ! ---------------------------------------------------------
  subroutine multisystem_abst_output_start(this)
    class(multisystem_abst_t), intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_output_start)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%output_start()
    end do

    POP_SUB(multisystem_abst_output_start)
  end subroutine multisystem_abst_output_start

  ! ---------------------------------------------------------
  recursive subroutine multisystem_abst_output_finish(this)
    class(multisystem_abst_t), intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_output_finish)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%output_finish()
    end do

    POP_SUB(multisystem_abst_output_finish)
  end subroutine multisystem_abst_output_finish

  ! ---------------------------------------------------------
  recursive subroutine multisystem_abst_output_write(this, iter)
    class(multisystem_abst_t), intent(inout) :: this
    integer,              intent(in)    :: iter

    type(system_iterator_t) :: iterator
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_abst_output_write)

    call iterator%start(this%list)
    do while (iterator%has_next())
      system => iterator%get_next()
      call system%output_write(iter)
    end do

    POP_SUB(multisystem_abst_output_write)
  end subroutine multisystem_abst_output_write

  ! ---------------------------------------------------------
  subroutine multisystem_abst_update_quantity(this, iq, requested_time)
    class(multisystem_abst_t), intent(inout) :: this
    integer,              intent(in)    :: iq
    class(clock_t),       intent(in)    :: requested_time

    PUSH_SUB(multisystem_abst_update_quantity)

    ! The multitystem class should never know about any specific quantities.
    ! Only classes that extend it can know about specific quantities.
    ! Such classes should override this method to add new supported quantities.
    message(1) = "Trying to update a quantity in the multisystem class"
    call messages_fatal(1)

    POP_SUB(multisystem_abst_update_quantity)
  end subroutine multisystem_abst_update_quantity

  ! ---------------------------------------------------------
  subroutine multisystem_abst_update_exposed_quantity(partner, iq, requested_time)
    class(multisystem_abst_t), intent(inout) :: partner
    integer,              intent(in)    :: iq
    class(clock_t),       intent(in)    :: requested_time

    PUSH_SUB(multisystem_abst_update_exposed_quantity)

    ! The multitystem class should never know about any specific quantities.
    ! Only classes that extend it can know about specific quantities.
    ! Such classes should override this method to add new supported quantities.
    message(1) = "Trying to update an exposed quantity in the multisystem class"
    call messages_fatal(1)

    POP_SUB(multisystem_abst_update_exposed_quantity)
  end subroutine multisystem_abst_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine multisystem_abst_copy_quantities_to_interaction(partner, interaction)
    class(multisystem_abst_t),         intent(inout) :: partner
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(multisystem_abst_copy_quantities_to_interaction)

    ! The multitystem class should never know about any specific quantities.
    ! Only classes that extend it can know about specific quantities.
    ! Such classes should override this method to add new supported quantities.
    message(1) = "Trying to copy quantities to interaction in the multisystem class"
    call messages_fatal(1)

    POP_SUB(multisystem_abst_copy_quantities_to_interaction)
  end subroutine multisystem_abst_copy_quantities_to_interaction

end module multisystem_abst_oct_m

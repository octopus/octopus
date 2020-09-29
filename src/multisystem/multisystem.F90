!! Copyright (C) 2019-2020 M. Oliveira, Heiko Appel
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

module multisystem_oct_m
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
  use propagator_oct_m
  use system_oct_m
  use system_factory_abst_oct_m
  implicit none

  private
  public ::               &
    multisystem_t,        &
    multisystem_init,     &
    multisystem_end

  type, extends(system_t), abstract :: multisystem_t
    type(system_list_t) :: list
  contains
    procedure :: copy_system => multisystem_copy_system
    procedure :: dt_operation =>  multisystem_dt_operation
    procedure :: init_parallelization => multisystem_init_parallelization
    procedure :: smallest_algo_dt => multisystem_smallest_algo_dt
    procedure :: reset_clocks => multisystem_reset_clocks
    procedure :: init_propagator => multisystem_init_propagator
    procedure :: propagation_start => multisystem_propagation_start
    procedure :: propagation_finish => multisystem_propagation_finish
    procedure :: has_reached_final_propagation_time => multisystem_has_reached_final_propagation_time
    procedure :: init_all_interactions => multisystem_init_all_interactions
    procedure :: init_interaction => multisystem_init_interaction
    procedure :: write_interaction_graph => multisystem_write_interaction_graph
    procedure :: initial_conditions => multisystem_initial_conditions
    procedure :: do_td_operation => multisystem_do_td_operation
    procedure :: iteration_info => multisystem_iteration_info
    procedure :: is_tolerance_reached => multisystem_is_tolerance_reached
    procedure :: store_current_status => multisystem_store_current_status
    procedure :: update_quantity => multisystem_update_quantity
    procedure :: update_exposed_quantity => multisystem_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => multisystem_copy_quantities_to_interaction
    procedure :: process_is_slave => multisystem_process_is_slave
  end type multisystem_t

contains

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_init(this, namespace, factory)
    class(multisystem_t),      intent(inout) :: this
    type(namespace_t),            intent(in) :: namespace
    class(system_factory_abst_t), intent(in) :: factory

    integer :: isys, system_type, ic
    character(len=128) :: system_name
    type(block_t) :: blk

    PUSH_SUB(multisystem_init)

    this%namespace = namespace

    if (parse_block(this%namespace, factory%block_name(), blk) == 0) then

      do isys = 1, parse_block_n(blk)
        ! Parse system name and type
        call parse_block_string(blk, isys - 1, 0, system_name)
        if (len_trim(system_name) == 0) then
          call messages_input_error(this%namespace, factory%block_name(), 'All systems must have a name')
        end if
        do ic = 1, len(parser_varname_excluded_characters)
          if (index(trim(system_name), parser_varname_excluded_characters(ic:ic)) /= 0) then
            call messages_input_error(this%namespace, factory%block_name(), &
              'Illegal character "' // parser_varname_excluded_characters(ic:ic) // '" in system name', row=isys-1, column=0)
          end if
        end do
        call parse_block_integer(blk, isys - 1, 1, system_type)

        call multisystem_create_system(this, system_name, system_type, isys, factory)
      end do
      call parse_block_end(blk)
    else
      message(1) = "Input error while reading block "//trim(this%namespace%get())//"."//trim(factory%block_name())
      call messages_fatal(1, namespace=this%namespace)
    end if

    POP_SUB(multisystem_init)
  end subroutine multisystem_init

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_create_system(this, system_name, system_type, isys, factory)
    class(multisystem_t),      intent(inout) :: this
    character(len=128),           intent(in) :: system_name
    integer,                      intent(in) :: system_type
    integer,                      intent(in) :: isys
    class(system_factory_abst_t), intent(in) :: factory
    class(system_t), pointer       :: dummy

    type(system_iterator_t) :: iter
    class(system_t), pointer :: sys, other

    PUSH_SUB(multisystem_create_system)

    ! Create folder to store system files.
    ! Needs to be done before creating the system as this in turn might create subfolders.
    call io_mkdir(system_name, namespace=this%namespace)

    ! Create system
    sys => factory%create(this%namespace, system_name, system_type)
    if (.not. associated(sys)) then
      call messages_input_error(this%namespace, factory%block_name(), 'Unknown system type.')
    end if

    ! Check that the system is unique
    call iter%start(this%list)
    do while (iter%has_next())
      other => iter%get_next()
      if (sys%namespace == other%namespace) then
        call messages_input_error(this%namespace, factory%block_name(), 'Duplicated system in multisystem', &
          row=isys-1, column=0)
      end if
    end do
    
    if(system_type == 3) then
        dummy => factory%create(this%namespace, "dummy", system_type)

        dummy = sys
    endif
    ! Add system to list of systems
    call this%list%add(sys)

    POP_SUB(multisystem_create_system)
  end subroutine multisystem_create_system

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_init_parallelization(this, grp)
    class(multisystem_t), intent(inout) :: this
    type(mpi_grp_t),      intent(in)    :: grp

    type(system_iterator_t) :: iter
    class(system_t), pointer :: sys
    type(mpi_grp_t) :: sys_grp

    PUSH_SUB(multisystem_init_parallelization)

    call mpi_grp_copy(this%grp, grp)

    ! Now parallelize over systems in this multisystem
    call iter%start(this%list)
    do while (iter%has_next())
      sys => iter%get_next()
      ! for now, duplicate communicator - more complicated parallelization schemes can be implemented here
      call mpi_grp_duplicate(sys_grp, grp)
      call sys%init_parallelization(sys_grp)
    end do

    POP_SUB(multisystem_init_parallelization)
  end subroutine multisystem_init_parallelization

  ! ---------------------------------------------------------------------------------------
  recursive function multisystem_smallest_algo_dt(this) result(smallest_algo_dt)
    class(multisystem_t), intent(inout) :: this
    FLOAT :: smallest_algo_dt

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_smallest_algo_dt)

    smallest_algo_dt = M_HUGE
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      select type (system)
      class is (multisystem_t)
        smallest_algo_dt = min(smallest_algo_dt, system%smallest_algo_dt())
      class default
        smallest_algo_dt = min(smallest_algo_dt, system%prop%dt/system%prop%algo_steps)
      end select
    end do

    POP_SUB(multisystem_smallest_algo_dt)
  end function multisystem_smallest_algo_dt

  ! ---------------------------------------------------------------------------------------
  subroutine multisystem_copy_system(lhs, rhs)
    class(multisystem_t), intent(out)  :: lhs
    class(*), intent(in) :: rhs

    select type (rhs)
    class is (multisystem_t)
      !Do stuff
    end select 
  end subroutine multisystem_copy_system

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_dt_operation(this)
    class(multisystem_t),     intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_dt_operation)

    ! Multisystem
    call system_dt_operation(this)

    ! Subsystems
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%dt_operation()
    end do

    POP_SUB(multisystem_dt_operation)
  end subroutine multisystem_dt_operation

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_reset_clocks(this, accumulated_ticks)
    class(multisystem_t),      intent(inout) :: this
    integer,                   intent(in)    :: accumulated_ticks

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_reset_clocks)

    ! Multisystem clocks
    call system_reset_clocks(this, accumulated_ticks)

    ! Subsystems clocks
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%reset_clocks(accumulated_ticks)
    end do

    POP_SUB(multisystem_reset_clocks)
  end subroutine multisystem_reset_clocks

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_init_propagator(this)
    class(multisystem_t),      intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system
    type(interaction_iterator_t) :: inter_iter
    class(interaction_t), pointer :: interaction

    PUSH_SUB(multisystem_init_propagator)

    ! Now initialized the propagators of the subsystems
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%init_propagator()
    end do

    ! Initialize the propagator of the multisystem
    ! Needs to be done after initializing the subsystems propagators,
    ! as we use the smallest dt of the subsystems.
    this%prop => propagator_t(this%smallest_algo_dt())
    this%interaction_timing = OPTION__INTERACTIONTIMING__TIMING_EXACT
    call this%prop%rewind()

    ! Initialize propagator clock
    this%prop%clock = clock_t(this%namespace%get(), this%prop%dt/this%prop%algo_steps)

    ! Initialize system clock
    this%clock = clock_t(this%namespace%get(), this%prop%dt)

    ! Interaction clocks
    call inter_iter%start(this%interactions)
    do while (inter_iter%has_next())
      interaction => inter_iter%get_next()
      interaction%clock = this%prop%clock - CLOCK_TICK
    end do

    ! Required quantities clocks
    where (this%quantities%required)
      this%quantities%clock = this%prop%clock
    end where

    POP_SUB(multisystem_init_propagator)
  end subroutine multisystem_init_propagator

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_propagation_start(this)
    class(multisystem_t),      intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_propagation_start)

    ! Start the propagation of the multisystem
    call system_propagation_start(this)

    ! Now start the propagation of the subsystems
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%propagation_start()
    end do

    POP_SUB(multisystem_propagation_start)
  end subroutine multisystem_propagation_start

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_propagation_finish(this)
    class(multisystem_t),      intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_propagation_finish)

    ! Finish the propagation of the multisystem
    call system_propagation_finish(this)

    ! Now finish the propagation of the subsystems
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%propagation_finish()
    end do

    POP_SUB(multisystem_propagation_finish)
  end subroutine multisystem_propagation_finish

  ! ---------------------------------------------------------------------------------------
  recursive logical function multisystem_has_reached_final_propagation_time(this, final_time)
    class(multisystem_t),      intent(inout) :: this
    FLOAT,                     intent(in)    :: final_time

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_has_reached_final_propagation_time)

    multisystem_has_reached_final_propagation_time = system_has_reached_final_propagation_time(this, final_time)
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      multisystem_has_reached_final_propagation_time = multisystem_has_reached_final_propagation_time .and. &
        system%has_reached_final_propagation_time(final_time)
    end do

    POP_SUB(multisystem_has_reached_final_propagation_time)
  end function multisystem_has_reached_final_propagation_time

  ! ---------------------------------------------------------
  recursive subroutine multisystem_init_all_interactions(this)
    class(multisystem_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter_i
    class(interaction_t), pointer :: interaction
    type(system_iterator_t) :: iter_s
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_init_all_interactions)

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

    POP_SUB(multisystem_init_all_interactions)
  end subroutine multisystem_init_all_interactions

  ! ---------------------------------------------------------
  subroutine multisystem_init_interaction(this, interaction)
    class(multisystem_t), target, intent(inout) :: this
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(multisystem_init_interaction)

    ! The multitystem class should never know about any specific interaction.
    ! Only classes that extend it can know about specific interactions.
    ! Such classes should override this method to add new supported interactions.
    message(1) = "Trying to initialize an interaction in the multisystem class"
    call messages_fatal(1)

    POP_SUB(multisystem_init_interaction)
  end subroutine multisystem_init_interaction

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_write_interaction_graph(this, iunit)
    class(multisystem_t), intent(in) :: this
    integer,              intent(in) :: iunit

    class(system_t), pointer :: system
    class(interaction_t), pointer :: interaction
    type(system_iterator_t) :: sys_iter
    type(interaction_iterator_t) :: inter_iter

    PUSH_SUB(multisystem_write_interaction_graph)

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
      class is (multisystem_t)
        call system%write_interaction_graph(iunit)
      end select
    end do

    POP_SUB(multisystem_write_interaction_graph)
  end subroutine multisystem_write_interaction_graph

  ! ---------------------------------------------------------
  recursive subroutine multisystem_initial_conditions(this, from_scratch)
    class(multisystem_t), intent(inout) :: this
    logical,              intent(in)    :: from_scratch

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_initial_conditions)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      call system%initial_conditions(from_scratch)
    end do

    POP_SUB(multisystem_initial_conditions)
  end subroutine multisystem_initial_conditions

  ! ---------------------------------------------------------
  subroutine multisystem_do_td_operation(this, operation)
    class(multisystem_t), intent(inout) :: this
    integer,              intent(in)    :: operation

    PUSH_SUB(multisystem_do_td_operation)

    select case (operation)
    case (SKIP)
      ! Nothing to do
    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(multisystem_do_td_operation)
  end subroutine multisystem_do_td_operation

  ! ---------------------------------------------------------
  recursive logical function multisystem_is_tolerance_reached(this, tol) result(converged)
    class(multisystem_t), intent(in)    :: this
    FLOAT,                intent(in)    :: tol

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_is_tolerance_reached)

    converged = .true.
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      if (.not. system%is_tolerance_reached(tol)) converged = .false.
    end do

    POP_SUB(multisystem_is_tolerance_reached)
  end function multisystem_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine multisystem_store_current_status(this)
    class(multisystem_t), intent(inout)    :: this

    PUSH_SUB(multisystem_store_current_status)

    ! The multitystem class does not have any status to be stored.
    ! Classes that extend it might need to override this method in order to
    ! support certain propagators.

    POP_SUB(multisystem_store_current_status)
  end subroutine multisystem_store_current_status

  ! ---------------------------------------------------------
  subroutine multisystem_iteration_info(this)
    class(multisystem_t), intent(in) :: this

    PUSH_SUB(multisystem_iteration_info)

    ! The multitystem class does not output any iteration info.
    ! Classes that extend it and need to output iteration info should override
    ! this method.

    POP_SUB(multisystem_iteration_info)
  end subroutine multisystem_iteration_info

  ! ---------------------------------------------------------
  subroutine multisystem_update_quantity(this, iq)
    class(multisystem_t), intent(inout) :: this
    integer,              intent(in)    :: iq

    PUSH_SUB(multisystem_update_quantity)

    ! The multitystem class should never know about any specific quantities.
    ! Only classes that extend it can know about specific quantities.
    ! Such classes should override this method to add new supported quantities.
    message(1) = "Trying to update a quantity in the multisystem class"
    call messages_fatal(1)

    POP_SUB(multisystem_update_quantity)
  end subroutine multisystem_update_quantity

  ! ---------------------------------------------------------
  subroutine multisystem_update_exposed_quantity(partner, iq)
    class(multisystem_t), intent(inout) :: partner
    integer,              intent(in)    :: iq

    PUSH_SUB(multisystem_update_exposed_quantity)

    ! The multitystem class should never know about any specific quantities.
    ! Only classes that extend it can know about specific quantities.
    ! Such classes should override this method to add new supported quantities.
    message(1) = "Trying to update an exposed quantity in the multisystem class"
    call messages_fatal(1)

    POP_SUB(multisystem_update_exposed_quantity)
  end subroutine multisystem_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine multisystem_copy_quantities_to_interaction(partner, interaction)
    class(multisystem_t),         intent(inout) :: partner
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(multisystem_copy_quantities_to_interaction)

    ! The multitystem class should never know about any specific quantities.
    ! Only classes that extend it can know about specific quantities.
    ! Such classes should override this method to add new supported quantities.
    message(1) = "Trying to copy quantities to interaction in the multisystem class"
    call messages_fatal(1)

    POP_SUB(multisystem_copy_quantities_to_interaction)
  end subroutine multisystem_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  recursive logical function multisystem_process_is_slave(this) result(is_slave)
    class(multisystem_t), intent(in) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_process_is_slave)

    is_slave = .false.
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      if (system%process_is_slave()) is_slave = .true.
    end do

    POP_SUB(multisystem_process_is_slave)
  end function multisystem_process_is_slave

  ! ---------------------------------------------------------
  recursive subroutine multisystem_end(this)
    class(multisystem_t), intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_end)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      SAFE_DEALLOCATE_P(system)
    end do

    call system_end(this)

    POP_SUB(multisystem_end)
  end subroutine multisystem_end

end module multisystem_oct_m

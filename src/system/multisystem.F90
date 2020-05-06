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

module multisystem_oct_m
  use clock_oct_m
  use global_oct_m
  use interaction_abst_oct_m
  use io_oct_m
  use linked_list_oct_m
  use loct_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use system_abst_oct_m
  use system_factory_abst_oct_m
  implicit none

  private
  public ::               &
    multisystem_t

  type, extends(system_abst_t) :: multisystem_t
    type(linked_list_t) :: list
  contains
    procedure :: init_interactions => multisystem_init_interactions
    procedure :: add_interaction_partner => multisystem_add_interaction_partner
    procedure :: has_interaction => multisystem_has_interaction
    procedure :: initial_conditions => multisystem_initial_conditions
    procedure :: do_td_operation => multisystem_do_td_operation
    procedure :: iteration_info => multisystem_iteration_info
    procedure :: output_start => multisystem_output_start
    procedure :: output_write => multisystem_output_write
    procedure :: output_finish => multisystem_output_finish
    procedure :: is_tolerance_reached => multisystem_is_tolerance_reached
    procedure :: store_current_status => multisystem_store_current_status
    procedure :: update_quantity => multisystem_update_quantity
    procedure :: update_exposed_quantity => multisystem_update_exposed_quantity
    procedure :: set_pointers_to_interaction => multisystem_set_pointers_to_interaction
    procedure :: update_interactions_start => multisystem_update_interactions_start
    procedure :: update_interactions_finish => multisystem_update_interactions_finish
    final :: multisystem_finalizer
  end type multisystem_t

  interface multisystem_t
    procedure multisystem_constructor
  end interface multisystem_t

contains

  ! ---------------------------------------------------------------------------------------
  function multisystem_constructor(namespace, factory) result(system)
    type(namespace_t),            intent(in) :: namespace
    class(system_factory_abst_t), intent(in) :: factory
    class(multisystem_t),         pointer    :: system

    integer :: isys, system_type
    character(len=128) :: system_name
    type(block_t) :: blk
    class(*), pointer :: sys
    
    PUSH_SUB(multisystem_constructor)

    SAFE_ALLOCATE(system)

    system%namespace = namespace

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
    !% A maxwell system. (NOT IMPLEMENTED)
    !%Option celestial_body 3
    !% A celestial body. Used for testing purposes only.
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
        sys => factory%create(system%namespace, system_name, system_type)
        if (.not. associated(sys)) then
          call messages_input_error(system%namespace, 'Systems', 'Unknown system type.')
        end if

        ! Add system to list of systems
        call system%list%add(sys)
      end do
      call parse_block_end(blk)
    else
      message(1) = "Input error while reading block Systems."
      call messages_fatal(1, namespace=system%namespace)
    end if

    POP_SUB(multisystem_constructor)
  end function multisystem_constructor

  ! ---------------------------------------------------------------------------------------
  subroutine multisystem_init_interactions(this)
    class(multisystem_t), intent(inout) :: this

    class(system_abst_t), pointer :: sys1, sys2
    type(system_iterator_t) :: iter1, iter2
    type(linked_list_t) :: flat_list
    integer :: iunit_out

    PUSH_SUB(multisystem_init_interactions)

    if (debug%interaction_graph .and. mpi_grp_is_root(mpi_world)) then
      iunit_out = io_open('interaction_graph.dot', this%namespace, action='write')
      write(iunit_out, '(a)') 'digraph {'
    end if

    ! We start by getting a list of all the subsystems that are not multisystems
    call flatten_list(this, flat_list)

    ! Double loop over all the subsystems that can interact
    call iter1%start(flat_list)
    do while (iter1%has_next())
      sys1 => iter1%get_next_system()

      call iter2%start(flat_list)
      do while (iter2%has_next())
        sys2 => iter2%get_next_system()

        !No self interaction
        if(.not.associated(sys1, sys2)) then
          call sys1%add_interaction_partner(sys2)

          !Debug information in form of a DOT graph
          if(debug%interaction_graph .and. mpi_grp_is_root(mpi_world)) then
            write(iunit_out, '(2x,a)') trim(sys1%namespace%get()) + ' -> ' + trim(sys2%namespace%get()) + ';'
          end if
        end if
      end do
    end do

    if (debug%interaction_graph .and. mpi_grp_is_root(mpi_world)) then
      write(iunit_out, '(a)') '}'
      call io_close(iunit_out)
    end if

    POP_SUB(multisystem_init_interactions)
  contains

    recursive subroutine flatten_list(systems, list)
      class(multisystem_t), intent(inout) :: systems
      type(linked_list_t),  intent(inout) :: list

      class(system_abst_t), pointer :: system
      type(system_iterator_t) :: iterator

      call iterator%start(systems%list)
      do while (iterator%has_next())
        system => iterator%get_next_system()

        select type (system)
        class is (multisystem_t)
          call flatten_list(system, list)
        class default
          call list%add(system)
        end select
      end do

    end subroutine flatten_list

  end subroutine multisystem_init_interactions

  ! ---------------------------------------------------------
  subroutine multisystem_add_interaction_partner(this, partner)
    class(multisystem_t), target, intent(inout) :: this
    class(system_abst_t),         intent(inout) :: partner

    type(system_iterator_t) :: iter
    class(system_abst_t), pointer :: system

    PUSH_SUB(multisystem_add_interaction_partner)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next_system()
      call system%add_interaction_partner(partner)
    end do

    POP_SUB(multisystem_add_interaction_partner)
  end subroutine multisystem_add_interaction_partner

  ! ---------------------------------------------------------
  logical function multisystem_has_interaction(this, interaction)
    class(multisystem_t),      intent(in) :: this
    class(interaction_abst_t), intent(in) :: interaction

    type(system_iterator_t) :: iter
    class(system_abst_t), pointer :: system

    PUSH_SUB(multisystem_has_interaction)

    multisystem_has_interaction = .false.
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next_system()
      if (system%has_interaction(interaction)) multisystem_has_interaction = .true.
    end do

    POP_SUB(multisystem_has_interaction)
  end function multisystem_has_interaction

  ! ---------------------------------------------------------
  subroutine multisystem_initial_conditions(this, from_scratch)
    class(multisystem_t), intent(inout) :: this
    logical,              intent(in)    :: from_scratch

    type(system_iterator_t) :: iter
    class(system_abst_t), pointer :: system

    PUSH_SUB(multisystem_initial_conditions)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next_system()
      call system%initial_conditions(from_scratch)
    end do

    POP_SUB(multisystem_initial_conditions)
  end subroutine multisystem_initial_conditions

  ! ---------------------------------------------------------
  subroutine multisystem_do_td_operation(this, operation)
    class(multisystem_t), intent(inout) :: this
    integer,              intent(in)    :: operation

    type(system_iterator_t) :: iter
    class(system_abst_t), pointer :: system

    PUSH_SUB(system_of_system_do_td_operation)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next_system()
      call system%do_td_operation(operation)
    end do

    POP_SUB(system_of_system_do_td_operation)
  end subroutine multisystem_do_td_operation

  ! ---------------------------------------------------------
  logical function multisystem_is_tolerance_reached(this, tol) result(converged)
    class(multisystem_t), intent(in)    :: this
    FLOAT,                intent(in)    :: tol

    type(system_iterator_t) :: iter
    class(system_abst_t), pointer :: system

    PUSH_SUB(multisystem_is_tolerance_reached)

    converged = .true.
    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next_system()
      if (.not. system%is_tolerance_reached(tol)) converged = .false.
    end do

    POP_SUB(multisystem_is_tolerance_reached)
  end function multisystem_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine multisystem_store_current_status(this)
    class(multisystem_t), intent(inout)    :: this

    type(system_iterator_t) :: iter
    class(system_abst_t), pointer :: system

    PUSH_SUB(multisystem_store_current_status)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next_system()
      call system%store_current_status()
    end do

    POP_SUB(multisystem_store_current_status)
  end subroutine multisystem_store_current_status

   ! ---------------------------------------------------------
  subroutine multisystem_iteration_info(this)
    class(multisystem_t), intent(in) :: this

    type(system_iterator_t) :: iter
    class(system_abst_t), pointer :: system

    PUSH_SUB(multisystem_iteration_info)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next_system()
      call system%iteration_info()
    end do

    POP_SUB(multisystem_iteration_info)
  end subroutine multisystem_iteration_info

  ! ---------------------------------------------------------
  subroutine multisystem_output_start(this)
    class(multisystem_t), intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_abst_t), pointer :: system

    PUSH_SUB(multisystem_output_start)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next_system()
      call system%output_start()
    end do

    POP_SUB(multisystem_output_start)
  end subroutine multisystem_output_start

  ! ---------------------------------------------------------
  subroutine multisystem_output_finish(this)
    class(multisystem_t), intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_abst_t), pointer :: system

    PUSH_SUB(multisystem_output_finish)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next_system()
      call system%output_finish()
    end do

    POP_SUB(multisystem_output_finish)
  end subroutine multisystem_output_finish

  ! ---------------------------------------------------------
  subroutine multisystem_output_write(this, iter)
    class(multisystem_t), intent(inout) :: this
    integer,              intent(in)    :: iter

    type(system_iterator_t) :: iterator
    class(system_abst_t), pointer :: system

    PUSH_SUB(multisystem_output_write)

    call iterator%start(this%list)
    do while (iterator%has_next())
      system => iterator%get_next_system()
      call system%output_write(iter)
    end do

    POP_SUB(multisystem_output_write)
  end subroutine multisystem_output_write

  ! ---------------------------------------------------------
  subroutine multisystem_update_quantity(this, iq, clock)
    class(multisystem_t), intent(inout) :: this
    integer,              intent(in)    :: iq
    class(clock_t),       intent(in)    :: clock

    PUSH_SUB(multisystem_update_quantity)

    ! At the moment multitystems cannot expose quantities.
    ! All the quantities are directly exposed by the subsystems
    ASSERT(.false.)

    POP_SUB(multisystem_update_quantity)
  end subroutine multisystem_update_quantity

  ! ---------------------------------------------------------
  logical function multisystem_update_exposed_quantity(this, iq, clock) result(updated)
    class(multisystem_t), intent(inout) :: this
    integer,              intent(in)    :: iq
    class(clock_t),       intent(in)    :: clock

    PUSH_SUB(multisystem_update_exposed_quantity)

    ! At the moment multitystems cannot expose quantities.
    ! All the quantities are directly exposed by the subsystems
    updated = .false.

    POP_SUB(multisystem_update_exposed_quantity)
  end function multisystem_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine multisystem_set_pointers_to_interaction(this, inter)
    class(multisystem_t), target, intent(inout) :: this
    class(interaction_abst_t),    intent(inout) :: inter

    PUSH_SUB(multisystem_set_pointers_to_interaction)

    ! At the moment multitystems cannot have interations.
    ! All the interactions are directly handled by the subsystems
    ASSERT(.false.)

    POP_SUB(multisystem_set_pointers_to_interaction)
  end subroutine multisystem_set_pointers_to_interaction

  ! ---------------------------------------------------------
  subroutine multisystem_update_interactions_start(this)
    class(multisystem_t), intent(inout) :: this

    PUSH_SUB(multisystem_update_interactions_start)

    ! At the moment multitystems cannot have interations
    ! All the interactions are directly handled by the subsystems
    ASSERT(.false.)

    POP_SUB(multisystem_update_interactions_start)
  end subroutine multisystem_update_interactions_start

  ! ---------------------------------------------------------
  subroutine multisystem_update_interactions_finish(this)
    class(multisystem_t), intent(inout) :: this

    PUSH_SUB(multisystem_update_interactions_finish)

    ! At the moment multitystems cannot have interations
    ! All the interactions are directly handled by the subsystems
    ASSERT(.false.)

    POP_SUB(multisystem_update_interactions_finish)
  end subroutine multisystem_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine multisystem_finalizer(this)
    type(multisystem_t), intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_abst_t), pointer :: system

    PUSH_SUB(multisystem_finalizer)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next_system()
      SAFE_DEALLOCATE_P(system)
    end do

    POP_SUB(multisystem_finalizer)
  end subroutine multisystem_finalizer

end module multisystem_oct_m

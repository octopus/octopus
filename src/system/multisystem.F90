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
  use global_oct_m
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
  public ::                        &
    multisystem_init,              &
    multisystem_init_interactions, &
    multisystem_end
  
contains

  ! ---------------------------------------------------------------------------------------
  subroutine multisystem_init(systems, namespace, factory)
    type(linked_list_t),          intent(inout) :: systems
    type(namespace_t),            intent(in)    :: namespace
    class(system_factory_abst_t), intent(in)    :: factory

    integer :: isys, system_type
    character(len=128) :: system_name
    type(block_t) :: blk
    class(*), pointer :: sys
    
    PUSH_SUB(multisystem_init)

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
    if (parse_block(namespace, 'Systems', blk) == 0) then

      do isys = 1, parse_block_n(blk)
        call parse_block_string(blk, isys - 1, 0, system_name)
        if (len_trim(system_name) == 0) then
          call messages_input_error(namespace, 'Systems', 'All systems must have a name.')
        end if
        call parse_block_integer(blk, isys - 1, 1, system_type)

        sys => factory%create(namespace, system_name, system_type)
        if (.not. associated(sys)) then
          call messages_input_error(namespace, 'Systems', 'Unknown system type.')
        end if
        call systems%add(sys)
      end do
      call parse_block_end(blk)
    else
      message(1) = "Input error while reading block Systems."
      call messages_fatal(1, namespace=namespace)
    end if

    POP_SUB(multisystem_init)
  end subroutine multisystem_init


  ! ---------------------------------------------------------------------------------------
  subroutine multisystem_init_interactions(systems, namespace)
    type(linked_list_t), intent(inout) :: systems
    type(namespace_t),   intent(in)    :: namespace

    class(system_abst_t), pointer :: sys1, sys2
    type(system_iterator_t) :: iter1, iter2
    integer :: iunit_out

    PUSH_SUB(multisystem_init_interactions)

    if (debug%interaction_graph .and. mpi_grp_is_root(mpi_world)) then
      iunit_out = io_open('interaction_graph.dot', namespace, action='write')
      write(iunit_out, '(a)') 'digraph {'
    end if

    call iter1%start(systems)
    do while (iter1%has_next())
      sys1 => iter1%get_next_system()
      call iter2%start(systems)
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

  end subroutine multisystem_init_interactions

  ! ---------------------------------------------------------------------------------------
  subroutine multisystem_end(systems)
    type(linked_list_t), intent(inout) :: systems

    type(system_iterator_t) :: iter
    class(system_abst_t), pointer :: sys

    PUSH_SUB(multisystem_end)

    call iter%start(systems)
    do while (iter%has_next())
      sys => iter%get_next_system()
      SAFE_DEALLOCATE_P(sys)
    end do

    POP_SUB(multisystem_end)
  end subroutine multisystem_end

end module multisystem_oct_m

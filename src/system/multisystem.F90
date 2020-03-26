!! Copyright (C) 2019 M. Oliveira
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
  use celestial_body_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use system_oct_m
  use system_abst_oct_m
  implicit none

  private
  public ::                        &
    multisystem_init,              &
    multisystem_init_interactions, &
    multisystem_end

  integer, parameter ::         &
    SYSTEM_ELECTRONIC     = 1,  &
    SYSTEM_MAXWELL        = 2,  &
    SYSTEM_CELESTIAL_BODY = 3
  
contains

  ! ---------------------------------------------------------------------------------------
  subroutine multisystem_init(systems, global_namespace)
    type(linked_list_t), intent(inout) :: systems
    type(namespace_t),   intent(in)  :: global_namespace

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
    !% An electronic system.
    !%Option maxwell 2
    !% A maxwell system.
    !%Option celestial_body 3
    !% A celestial body. Used for testing purposes only.
    !%End
    if(parse_block(global_namespace, 'Systems', blk) == 0) then

      do isys = 1, parse_block_n(blk)
        call parse_block_string(blk, isys - 1, 0, system_name)
        call parse_block_integer(blk, isys - 1, 1, system_type)

        select case (system_type)
        case (SYSTEM_ELECTRONIC)
          sys => system_init(namespace_t(system_name))
          call systems%add(sys)
        case (SYSTEM_CELESTIAL_BODY)
          sys => celestial_body_t(namespace_t(system_name))
          call systems%add(sys)
        case default
          call messages_input_error('Systems')
        end select
      end do
      call parse_block_end(blk)
    else
      sys => system_init(global_namespace)
      call systems%add(sys)
    end if

    POP_SUB(multisystem_init)
  end subroutine multisystem_init


  ! ---------------------------------------------------------------------------------------
  subroutine multisystem_init_interactions(systems)
    type(linked_list_t), intent(inout) :: systems

    class(system_abst_t), pointer :: sys1, sys2
    type(system_iterator_t) :: iter1, iter2

    PUSH_SUB(multisystem_init_interactions)

    call iter1%start(systems)
    do while (iter1%has_next())
      sys1 => iter1%get_next_system()
      call iter2%start(systems)
      do while (iter2%has_next())
        sys2 => iter2%get_next_system()
        !No self interaction
        if(.not.associated(sys1, sys2)) then
          call sys1%add_interaction_partner(sys2)
        end if
      end do
    end do

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

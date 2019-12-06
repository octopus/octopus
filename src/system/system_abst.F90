!! Copyright (C) 2019 N. Tancogne-Dejean
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
  use messages_oct_m
  use linked_list_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m

  implicit none

  private
  public ::               &
    system_abst_t

  integer, public, parameter ::        &
    FINISHED                     = 0,  &
    VERLET_UPDATE_POS            = 1,  &
    VERLET_COMPUTE_ACC           = 2,  &
    VERLET_COMPUTE_VEL           = 3,  &
    VERLET_SYNC_DT               = 4,  &
    UPDATE_INTERACTION           = 5

  integer, public, parameter ::        &
    TOTAL_CURRENT                = 1,  &
    FORCE                        = 2


  type, abstract :: system_abst_t
    private

    type(linked_list_t), public :: interactions
    integer, public :: nb_partners = 0

  contains
    procedure(system_do_td_op), deferred :: do_td_operation
    procedure(system_pull_interaction), deferred :: pull_interaction
    procedure(system_get_needed_quantity), deferred :: get_needed_quantity
    procedure :: add_interaction_partner
    procedure(system_set_propagator), deferred :: set_propagator
    procedure(system_alloc_receiver), deferred :: allocate_receiv_structure
    procedure :: system_dt
  end type system_abst_t

  abstract interface
    subroutine system_do_td_op(this, operation)
      import system_abst_t
      class(system_abst_t), intent(inout) :: this
      integer             , intent(in)    :: operation
    end subroutine system_do_td_op

    subroutine system_pull_interaction(sys, remote, interaction, partner_index)
      import system_abst_t
      class(system_abst_t),     intent(inout) :: sys
      class(system_abst_t),     intent(inout) :: remote
      integer,                  intent(in)    :: interaction
      integer,                  intent(in)    :: partner_index
    end subroutine system_pull_interaction

    integer function system_get_needed_quantity(sys)
      import system_abst_t
      class(system_abst_t),     intent(in) :: sys
    end function system_get_needed_quantity

    subroutine system_set_propagator(sys, prop)
      import system_abst_t
      import propagator_abst_t
      class(system_abst_t),             intent(inout) :: sys
      class(propagator_abst_t), target, intent(in)    :: prop
    end subroutine system_set_propagator

    subroutine system_alloc_receiver(this)
      import system_abst_t
      class(system_abst_t), intent(inout) :: this
    end subroutine system_alloc_receiver
  end interface

contains

  subroutine add_interaction_partner(sys, remote)
    class(system_abst_t),     intent(inout) :: sys
    class(system_abst_t),     intent(in)    :: remote

    PUSH_SUB(add_interaction_partner)

    call sys%interactions%add_node(remote)
    sys%nb_partners = sys%nb_partners+1

    POP_SUB(add_interaction_partner)
  end subroutine add_interaction_partner

  subroutine system_dt(this, prop)
    class(system_abst_t),     intent(inout) :: this
    class(propagator_abst_t), intent(inout) :: prop

    integer :: tdop, inter, counter
    class(*), pointer :: sys

    PUSH_SUB(system_dt)

    tdop = prop%get_td_operation()
    select case(tdop)
    case(FINISHED)
     ! call propagator_finished(prop)
      !DO OUTPUT HERE AND BROADCAST NEEDED QUANTITIES
      !ONLY IF WE ARE NOT YET FINISHED

    case(UPDATE_INTERACTION)

      inter = this%get_needed_quantity()

      ! Loop over systems
      counter = 1
      call this%interactions%rewind()
      do while (this%interactions%has_more_values())
        sys => this%interactions%current()
        select type (sys)
        class is (system_abst_t)
          call this%pull_interaction(sys, inter, counter)
        class default
          message(1) = "Unknow system type."
          call messages_fatal(1)
        end select
        call this%interactions%next()
        counter = counter + 1
      end do

    case default
      call this%do_td_operation(tdop)
    end select

    POP_SUB(system_dt)
  end subroutine system_dt

end module system_abst_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

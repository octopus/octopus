!! Copyright (C)  2019 N. Tancogne-Dejean
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

module propagator_abst_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                            &
    propagator_abst_t

  type, extends(linked_list_t) :: propagator_abst_t
    private

    type(list_iterator_t) :: iter
    integer               :: current_ops

    FLOAT, public   :: internal_time
    FLOAT, public   :: dt
    integer, public :: algo_steps
    logical :: step_done

  contains
    !Below are the list of operations that needs to be implemented
    procedure :: get_td_operation => propagator_get_tdop
    procedure :: step_is_done => propagator_step_is_done
    procedure :: next => propagator_next
    procedure :: rewind => propagator_rewind
    procedure :: finished => propagator_finished
  end type propagator_abst_t

  ! Known propagation operations
  integer, public, parameter ::        &
    FINISHED                     = 0,  &
    VERLET_UPDATE_POS            = 1,  &
    VERLET_COMPUTE_ACC           = 2,  &
    VERLET_COMPUTE_VEL           = 3,  &
    VERLET_SYNC_DT               = 4,  &
    UPDATE_INTERACTIONS          = 5

contains

  subroutine propagator_rewind(this)
    class(propagator_abst_t), intent(inout) :: this

    PUSH_SUB(propagator_rewind)

    call this%iter%start(this)
    call this%next()
    this%step_done = .false.

    POP_SUB(propagator_rewind)
  end subroutine propagator_rewind

  subroutine propagator_finished(this)
    class(propagator_abst_t), intent(inout) :: this

    PUSH_SUB(propagator_finished)

    this%step_done = .true.

    POP_SUB(propagator_finished)
  end subroutine propagator_finished

  subroutine propagator_next(this)
    class(propagator_abst_t), intent(inout) :: this

    PUSH_SUB(propagator_next)

    select type(next_ops => this%iter%get_next())
    type is (integer)
      this%current_ops = next_ops
    class default
      message(1) = "Corrupted list."
      call messages_fatal(1)
    end select

    POP_SUB(propagator_next)
  end subroutine propagator_next

  integer function propagator_get_tdop(this) result(tdop)
    class(propagator_abst_t), intent(in) :: this

    PUSH_SUB(propagator_get_tdop)

    tdop = this%current_ops

    POP_SUB(propagator_get_tdop)
  end function propagator_get_tdop

  logical pure function propagator_step_is_done(this) result(step_is_done)
    class(propagator_abst_t), intent(in) :: this

    step_is_done = this%step_done

  end function propagator_step_is_done

end module propagator_abst_oct_m


!!o, Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

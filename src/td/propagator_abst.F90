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
  use clock_oct_m
  use global_oct_m
  use linked_list_oct_m
  use list_node_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                            &
    propagator_abst_t

  type, extends(linked_list_t) :: propagator_abst_t
    private

    type(list_iterator_t) :: iter
    type(list_iterator_t) :: scf_start
    integer               :: current_ops

    integer, public :: algo_steps
    FLOAT, public   :: dt

    !< Options related to predictor-corrector propagators
    logical, public :: predictor_corrector = .false.
    integer, public :: scf_count
    integer, public :: max_scf_count
    FLOAT, public   :: scf_tol

    logical :: step_done
    integer, public :: last_step_done_tick = -1

    type(clock_t), public :: clock

  contains
    !Below are the list of operations that needs to be implemented
    procedure :: get_td_operation => propagator_get_tdop
    procedure :: step_is_done => propagator_step_is_done
    procedure :: next => propagator_next
    procedure :: rewind => propagator_rewind
    procedure :: finished => propagator_finished
    procedure :: save_scf_start => propagator_save_scf_start
    procedure :: rewind_scf_loop => propagator_rewind_scf_loop
  end type propagator_abst_t

  ! Known propagation operations
  integer, public, parameter ::        &
    FINISHED                     =  0,  &
    VERLET_UPDATE_POS            =  1,  &
    VERLET_COMPUTE_ACC           =  2,  &
    VERLET_COMPUTE_VEL           =  3,  &
    SYNC                         =  4,  &
    UPDATE_INTERACTIONS          =  5,  &
    START_SCF_LOOP               =  6,  &
    END_SCF_LOOP                 =  7,  &
    STORE_CURRENT_STATUS         =  8,  &
    BEEMAN_PREDICT_POS           =  9,  &
    BEEMAN_PREDICT_VEL           = 10,  &
    BEEMAN_CORRECT_POS           = 11,  &
    BEEMAN_CORRECT_VEL           = 12

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
    this%last_step_done_tick = this%clock%get_tick()

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

  subroutine propagator_save_scf_start(this)
    class(propagator_abst_t), intent(inout) :: this
    
    PUSH_SUB(propagator_save_scf_start)

    !Save the current iteration state (START_SCF_LOOP) and move to next step
    this%scf_start = this%iter
    call this%next()
    this%scf_count = 0

    POP_SUB(propagator_save_scf_start) 

  end subroutine propagator_save_scf_start

  subroutine propagator_rewind_scf_loop(this)
    class(propagator_abst_t), intent(inout) :: this

    PUSH_SUB(propagator_rewind_scf_loop)

    !Reset the iteration state to the beginning of the loop (START_SCF_LOOP) and move to next step
    this%iter = this%scf_start
    call this%next()
    this%scf_count = this%scf_count + 1

    POP_SUB(propagator_rewind_scf_loop)

  end subroutine propagator_rewind_scf_loop

end module propagator_abst_oct_m


!!o, Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

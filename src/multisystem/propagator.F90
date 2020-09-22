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

module propagator_oct_m
  use clock_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                       &
    propagator_t,                 &
    propagator_step_debug_message

  type, extends(integer_list_t) :: propagator_t
    private

    type(integer_iterator_t) :: iter
    type(integer_iterator_t) :: scf_start
    integer               :: current_ops

    integer, public       :: start_step
    integer, public       :: final_step

    integer, public :: algo_steps
    FLOAT, public   :: dt

    !< Options related to predictor-corrector propagators
    logical, public :: predictor_corrector = .false.
    integer, public :: scf_count
    integer, public :: max_scf_count
    FLOAT, public   :: scf_tol

    logical :: step_done
    logical, public :: inside_scf = .false.

    type(clock_t), public :: clock

  contains
    ! Below are the list of operations that needs to be implemented
    procedure :: get_td_operation => propagator_get_tdop
    procedure :: step_is_done => propagator_step_is_done
    procedure :: next => propagator_next
    procedure :: rewind => propagator_rewind
    procedure :: finished => propagator_finished
    procedure :: save_scf_start => propagator_save_scf_start
    procedure :: rewind_scf_loop => propagator_rewind_scf_loop
  end type propagator_t

  ! Known propagation operations
  integer, public, parameter ::         &
    SKIP                         = -1,  &
    FINISHED                     =  0,  &
    VERLET_START                 =  1,  &
    VERLET_FINISH                =  2,  &
    VERLET_UPDATE_POS            =  3,  &
    VERLET_COMPUTE_ACC           =  4,  &
    VERLET_COMPUTE_VEL           =  5,  &
    SYNC                         =  6,  &
    UPDATE_INTERACTIONS          =  7,  &
    START_SCF_LOOP               =  8,  &
    END_SCF_LOOP                 =  9,  &
    STORE_CURRENT_STATUS         = 10,  &
    BEEMAN_START                 = 11,  &
    BEEMAN_FINISH                = 12,  &
    BEEMAN_COMPUTE_ACC           = 13,  &
    BEEMAN_PREDICT_POS           = 14,  &
    BEEMAN_PREDICT_VEL           = 15,  &
    BEEMAN_CORRECT_POS           = 16,  &
    BEEMAN_CORRECT_VEL           = 17,  &
    EXPMID_START                 = 18,  &
    EXPMID_FINISH                = 19,  &
    EXPMID_PREDICT_DT_2          = 20,  &
    EXPMID_PREDICT_DT            = 21,  &
    EXPMID_CORRECT_DT_2          = 22,  &
    UPDATE_HAMILTONIAN           = 23

  ! Known multisystem propagators
  integer, public, parameter ::        &
    PROP_VERLET                  = 1,  &
    PROP_BEEMAN                  = 2,  &
    PROP_BEEMAN_SCF              = 3,  &
    PROP_EXPMID                  = 4,  &
    PROP_EXPMID_SCF              = 5

  interface propagator_t
    procedure propagator_constructor
  end interface propagator_t

contains

  ! ---------------------------------------------------------
  function propagator_constructor(dt) result(this)
    FLOAT,              intent(in) :: dt
    type(propagator_t), pointer    :: this

    PUSH_SUB(propagator_constructor)

    SAFE_ALLOCATE(this)

    this%start_step = SKIP
    this%final_step = SKIP

    call this%add(UPDATE_INTERACTIONS)
    call this%add(FINISHED)

    this%algo_steps = 1

    this%dt = dt

    POP_SUB(propagator_constructor)
  end function propagator_constructor

  ! ---------------------------------------------------------
  subroutine propagator_rewind(this)
    class(propagator_t), intent(inout) :: this

    PUSH_SUB(propagator_rewind)

    call this%iter%start(this)
    call this%next()
    this%step_done = .false.

    POP_SUB(propagator_rewind)
  end subroutine propagator_rewind

  ! ---------------------------------------------------------
  subroutine propagator_finished(this)
    class(propagator_t), intent(inout) :: this

    PUSH_SUB(propagator_finished)

    this%step_done = .true.

    POP_SUB(propagator_finished)
  end subroutine propagator_finished

  ! ---------------------------------------------------------
  subroutine propagator_next(this)
    class(propagator_t), intent(inout) :: this

    PUSH_SUB(propagator_next)

    this%current_ops = this%iter%get_next()

    POP_SUB(propagator_next)
  end subroutine propagator_next

  ! ---------------------------------------------------------
  integer function propagator_get_tdop(this) result(tdop)
    class(propagator_t), intent(in) :: this

    PUSH_SUB(propagator_get_tdop)

    tdop = this%current_ops

    POP_SUB(propagator_get_tdop)
  end function propagator_get_tdop

  ! ---------------------------------------------------------
  logical pure function propagator_step_is_done(this) result(step_is_done)
    class(propagator_t), intent(in) :: this

    step_is_done = this%step_done

  end function propagator_step_is_done

  ! ---------------------------------------------------------
  subroutine propagator_save_scf_start(this)
    class(propagator_t), intent(inout) :: this
    
    PUSH_SUB(propagator_save_scf_start)

    ! Save the current iteration state (START_SCF_LOOP) and move to next step
    this%scf_start = this%iter
    call this%next()
    this%scf_count = 0

    POP_SUB(propagator_save_scf_start) 

  end subroutine propagator_save_scf_start

  ! ---------------------------------------------------------
  subroutine propagator_rewind_scf_loop(this)
    class(propagator_t), intent(inout) :: this

    PUSH_SUB(propagator_rewind_scf_loop)

    ! Reset the iteration state to the beginning of the loop (START_SCF_LOOP) and move to next step
    this%iter = this%scf_start
    call this%next()
    this%scf_count = this%scf_count + 1

    POP_SUB(propagator_rewind_scf_loop)

  end subroutine propagator_rewind_scf_loop

  ! ---------------------------------------------------------
  function propagator_step_debug_message(tdop) result(description)
    integer, intent(in) :: tdop
    character(len=100) :: description

    PUSH_SUB(propagator_step_debug_message)

    select case (tdop)
    case (SKIP)
      description = "Skipping propagation step for"
    case (FINISHED)
      description = "Propagation step finished for"
    case (UPDATE_INTERACTIONS)
      description = "Updating interactions for"
    case (START_SCF_LOOP)
      description = "Starting SCF loop for"
    case (END_SCF_LOOP)
      description = "End of SCF iteration for"
    case (STORE_CURRENT_STATUS)
      description = "Storing the current status for"
    case (VERLET_UPDATE_POS)
      description = "Propagation step - Updating positions for"
    case (VERLET_COMPUTE_ACC)
      description = "Propagation step - Computing acceleration for"
    case (VERLET_COMPUTE_VEL)
      description = "Propagation step - Computing velocity for"
    case (BEEMAN_COMPUTE_ACC)
      description = "Propagation step - Computing acceleration for"
    case (BEEMAN_PREDICT_POS)
      description = "Prediction step  - Computing position for"
    case (BEEMAN_PREDICT_VEL)
      description = "Prediction step  - Computing velocity for"
    case (BEEMAN_CORRECT_POS)
      description = "Correction step  - Computing position for"
    case (BEEMAN_CORRECT_VEL)
      description = "Correction step  - Computing velocity for"
    case (EXPMID_PREDICT_DT_2)
      description = "Predict state at dt/2 for"
    case (EXPMID_PREDICT_DT)
      description = "Predict state at dt for"
    case (EXPMID_CORRECT_DT_2)
      description = "Correct state at dt/2 for"
    case (UPDATE_HAMILTONIAN)
      description = "Updating Hamiltonian for"
    case default
      description = "Unknown step for "
    end select

    POP_SUB(propagator_step_debug_message)
  end function propagator_step_debug_message

end module propagator_oct_m


!!o, Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

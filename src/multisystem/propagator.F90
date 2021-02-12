!! Copyright (C)  2019 N. Tancogne-Dejean
!! Copyright (C)  2020 M. Oliveira
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
  use algorithm_oct_m
  use clock_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                       &
    propagator_t

  type, extends(algorithm_t) :: propagator_t
    private

    type(algorithm_iterator_t) :: iter
    type(algorithm_iterator_t) :: scf_start
    type(algorithmic_operation_t) :: current_ops

    type(algorithmic_operation_t), public       :: start_step
    type(algorithmic_operation_t), public       :: final_step

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

  !# doc_start general_propagation_operations
  ! Known propagation operations
  character(len=ALGO_LABEL_LEN), public, parameter ::   &
    SKIP                 = 'SKIP',                      &
    FINISHED             = 'FINISHED',                  &
    UPDATE_INTERACTIONS  = 'UPDATE_INTERACTIONS',       &
    START_SCF_LOOP       = 'START_SCF_LOOP',            &
    END_SCF_LOOP         = 'END_SCF_LOOP',              &
    STORE_CURRENT_STATUS = 'STORE_CURRENT_STATUS'

  type(algorithmic_operation_t), public, parameter :: &
    OP_SKIP                 = algorithmic_operation_t(SKIP, 'Skipping propagation step'), &
    OP_FINISHED             = algorithmic_operation_t(FINISHED,             'Propagation step finished'), &
    OP_UPDATE_INTERACTIONS  = algorithmic_operation_t(UPDATE_INTERACTIONS,  'Updating interactions'),     &
    OP_START_SCF_LOOP       = algorithmic_operation_t(START_SCF_LOOP,       'Starting SCF loop'),         &
    OP_END_SCF_LOOP         = algorithmic_operation_t(END_SCF_LOOP,         'End of SCF iteration'),      &
    OP_STORE_CURRENT_STATUS = algorithmic_operation_t(STORE_CURRENT_STATUS, '')
  !# doc_end

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

    this%start_step = OP_SKIP
    this%final_step = OP_SKIP

    call this%add_operation(OP_UPDATE_INTERACTIONS)
    call this%add_operation(OP_FINISHED)

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
  type(algorithmic_operation_t) function propagator_get_tdop(this) result(tdop)
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

end module propagator_oct_m


!!o, Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

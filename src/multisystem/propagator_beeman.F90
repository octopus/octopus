!! Copyright (C)  2020 N. Tancogne-Dejean
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

module propagator_beeman_oct_m
  use algorithm_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use propagator_oct_m

  implicit none

  private
  public ::                            &
    propagator_beeman_t

  type, extends(propagator_t) :: propagator_beeman_t
    private
  end type propagator_beeman_t

  interface propagator_beeman_t
    procedure propagator_beeman_constructor
  end interface propagator_beeman_t

  ! Specific beeman propagation operations identifiers
  character(len=30), public, parameter ::      &
    BEEMAN_START       = 'BEEMAN_START',        &
    BEEMAN_FINISH      = 'BEEMAN_FINISH',       &
    BEEMAN_COMPUTE_ACC = 'BEEMAN_COMPUTE_ACC',  &
    BEEMAN_PREDICT_POS = 'BEEMAN_PREDICT_POS',  &
    BEEMAN_PREDICT_VEL = 'BEEMAN_PREDICT_VEL',  &
    BEEMAN_CORRECT_POS = 'BEEMAN_CORRECT_POS',  &
    BEEMAN_CORRECT_VEL = 'BEEMAN_CORRECT_VEL'

  ! Specific beeman propagation operations
  type(algorithmic_operation_t), public, parameter :: &
    OP_BEEMAN_START       = algorithmic_operation_t(BEEMAN_START,       'Starting Beeman propagation'),               &
    OP_BEEMAN_FINISH      = algorithmic_operation_t(BEEMAN_FINISH,      'Finishing Beeman propagation'),              &
    OP_BEEMAN_COMPUTE_ACC = algorithmic_operation_t(BEEMAN_COMPUTE_ACC, 'Propagation step - Computing acceleration'), &
    OP_BEEMAN_PREDICT_POS = algorithmic_operation_t(BEEMAN_PREDICT_POS, 'Prediction step  - Computing position'),     &
    OP_BEEMAN_PREDICT_VEL = algorithmic_operation_t(BEEMAN_PREDICT_VEL, 'Prediction step  - Computing velocity'),     &
    OP_BEEMAN_CORRECT_POS = algorithmic_operation_t(BEEMAN_CORRECT_POS, 'Correction step  - Computing position'),     &
    OP_BEEMAN_CORRECT_VEL = algorithmic_operation_t(BEEMAN_CORRECT_VEL, 'Correction step  - Computing velocity')

contains

  ! ---------------------------------------------------------
  function propagator_beeman_constructor(dt, predictor_corrector) result(this)
    FLOAT,               intent(in)    :: dt
    logical,             intent(in)    :: predictor_corrector
    type(propagator_beeman_t), pointer :: this

    PUSH_SUB(propagator_beeman_constructor)

    SAFE_ALLOCATE(this)

    this%predictor_corrector = predictor_corrector

    this%start_step = OP_BEEMAN_START
    this%final_step = OP_BEEMAN_FINISH

    if(predictor_corrector) then

      call this%add_operation(OP_STORE_CURRENT_STATUS)
      call this%add_operation(OP_BEEMAN_PREDICT_POS)
      call this%add_operation(OP_START_SCF_LOOP)
      call this%add_operation(OP_UPDATE_INTERACTIONS)
      call this%add_operation(OP_BEEMAN_COMPUTE_ACC)
      call this%add_operation(OP_BEEMAN_CORRECT_POS)
      call this%add_operation(OP_BEEMAN_CORRECT_VEL)
      call this%add_operation(OP_END_SCF_LOOP)
      call this%add_operation(OP_FINISHED)

      this%max_scf_count = 2 !From Wikipedia
      this%scf_tol = CNST(1e-6) !At the moment arbitrary
 
    else

      call this%add_operation(OP_BEEMAN_PREDICT_POS)
      call this%add_operation(OP_UPDATE_INTERACTIONS)
      call this%add_operation(OP_BEEMAN_COMPUTE_ACC)
      call this%add_operation(OP_BEEMAN_PREDICT_VEL)
      call this%add_operation(OP_FINISHED)

    end if

    ! Beeman has only one algorithmic step
    this%algo_steps = 1

    this%dt = dt

    POP_SUB(propagator_beeman_constructor)
  end function propagator_beeman_constructor

end module propagator_beeman_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

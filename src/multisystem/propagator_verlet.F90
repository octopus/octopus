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

module propagator_verlet_oct_m
  use algorithm_oct_m
  use clock_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use propagator_oct_m

  implicit none

  private
  public ::                            &
    propagator_verlet_t

  type, extends(propagator_t) :: propagator_verlet_t
    private
  end type propagator_verlet_t

  interface propagator_verlet_t
    procedure propagator_verlet_constructor
  end interface propagator_verlet_t

  ! Specific verlet propagation operations identifiers
  character(len=30), public, parameter ::      &
    VERLET_START       = 'VERLET_START',       &
    VERLET_FINISH      = 'VERLET_FINISH',      &
    VERLET_UPDATE_POS  = 'VERLET_UPDATE_POS',  &
    VERLET_COMPUTE_ACC = 'VERLET_COMPUTE_ACC', &
    VERLET_COMPUTE_VEL = 'VERLET_COMPUTE_VEL'

  ! Specific verlet propagation operations
  type(algorithmic_operation_t), public, parameter :: &
    OP_VERLET_START       = algorithmic_operation_t(VERLET_START,       'Starting Verlet propagation'),               &
    OP_VERLET_FINISH      = algorithmic_operation_t(VERLET_FINISH,      'Finishing Verlet propagation'),              &
    OP_VERLET_UPDATE_POS  = algorithmic_operation_t(VERLET_UPDATE_POS,  'Propagation step - Updating positions'),     &
    OP_VERLET_COMPUTE_ACC = algorithmic_operation_t(VERLET_COMPUTE_ACC, 'Propagation step - Computing acceleration'), &
    OP_VERLET_COMPUTE_VEL = algorithmic_operation_t(VERLET_COMPUTE_VEL, 'Propagation step - Computing velocity')

contains

  ! ---------------------------------------------------------
  function propagator_verlet_constructor(dt) result(this)
    FLOAT,                     intent(in) :: dt
    type(propagator_verlet_t), pointer    :: this

    PUSH_SUB(propagator_verlet_constructor)

    SAFE_ALLOCATE(this)

    this%start_step = OP_VERLET_START
    this%final_step = OP_VERLET_FINISH

    call this%add_operation(OP_VERLET_UPDATE_POS)
    call this%add_operation(OP_UPDATE_INTERACTIONS)
    call this%add_operation(OP_VERLET_COMPUTE_ACC)
    call this%add_operation(OP_VERLET_COMPUTE_VEL)
    call this%add_operation(OP_FINISHED)

    ! Verlet has only one algorithmic step
    this%algo_steps = 1

    this%dt = dt

    POP_SUB(propagator_verlet_constructor)
  end function propagator_verlet_constructor

end module propagator_verlet_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

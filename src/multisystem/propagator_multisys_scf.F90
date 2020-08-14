!! Copyright (C) 2020 Heiko Appel
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

module propagator_multisys_scf_oct_m
  use clock_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use propagator_oct_m

  implicit none

  private
  public ::                            &
    propagator_multisys_scf_t

  type, extends(propagator_t) :: propagator_multisys_scf_t
    private
  end type propagator_multisys_scf_t

  interface propagator_multisys_scf_t
    procedure propagator_multisys_scf_constructor
  end interface propagator_multisys_scf_t

contains

  ! ---------------------------------------------------------
  function propagator_multisys_scf_constructor(dt, predictor_corrector) result(this)
    FLOAT,                     intent(in) :: dt
    logical,                   intent(in) :: predictor_corrector
    type(propagator_multisys_scf_t), pointer   :: this

    PUSH_SUB(propagator_multisys_scf_constructor)

    SAFE_ALLOCATE(this)

    this%predictor_corrector = .true.    ! by definition this propagator is always self-consistent
    this%start_step = MULTISYS_SCF_START
    this%final_step = MULTISYS_SCF_FINISH

    call this%add(STORE_CURRENT_STATUS)
    call this%add(MULTISYS_SCF_PROPAGATE_SYSTEMS)
    call this%add(FINISHED)

    this%max_scf_count = 10
    this%scf_tol = CNST(1e-6) ! At the moment arbitrary. This is system specific and should be adapted.

    ! The SCF multi-system propagator has only one algorithmic step
    this%algo_steps = 2

    this%dt = dt

    POP_SUB(propagator_multisys_scf_constructor)
  end function propagator_multisys_scf_constructor

end module propagator_multisys_scf_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

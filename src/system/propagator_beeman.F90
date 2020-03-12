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
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m

  implicit none

  private
  public ::                            &
    propagator_beeman_t

  type, extends(propagator_abst_t) :: propagator_beeman_t
    private
  end type propagator_beeman_t

  interface propagator_beeman_t
    procedure propagator_beeman_init
  end interface propagator_beeman_t

contains

  ! ---------------------------------------------------------
  function propagator_beeman_init(namespace, predictor_corrector) result(this)
    type(namespace_t),   intent(in)    :: namespace
    logical,             intent(in)    :: predictor_corrector
    type(propagator_beeman_t), pointer :: this

    PUSH_SUB(propagator_beeman_init)

    SAFE_ALLOCATE(this)

    this%predictor_corrector = predictor_corrector

    if(predictor_corrector) then

      call this%add(STORE_CURRENT_STATUS)
      call this%add(BEEMAN_PREDICT_POS)
      call this%add(START_SCF_LOOP)
      call this%add(UPDATE_INTERACTIONS)
      call this%add(VERLET_COMPUTE_ACC)
      call this%add(BEEMAN_CORRECT_POS)
      call this%add(BEEMAN_CORRECT_VEL)
      call this%add(END_SCF_LOOP)
      call this%add(FINISHED)

      this%max_scf_count = 2 !From Wikipedia
      this%scf_tol = CNST(1e-6) !At the moment arbitrary
 
    else

      call this%add(BEEMAN_PREDICT_POS)
      call this%add(UPDATE_INTERACTIONS)
      call this%add(VERLET_COMPUTE_ACC)
      call this%add(BEEMAN_PREDICT_VEL)
      call this%add(FINISHED)

    end if

    ! Beeman has only one algorithmic step
    this%algo_steps = 1

    call this%parse_td_variables(namespace)

    POP_SUB(propagator_beeman_init)
  end function propagator_beeman_init

end module propagator_beeman_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

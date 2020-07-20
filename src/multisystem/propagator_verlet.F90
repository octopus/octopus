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
  use clock_oct_m
  use global_oct_m
  use gauge_field_oct_m
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

contains

  ! ---------------------------------------------------------
  function propagator_verlet_constructor(namespace) result(this)
    type(namespace_t),         intent(in) :: namespace
    type(propagator_verlet_t), pointer    :: this

    PUSH_SUB(propagator_verlet_constructor)

    SAFE_ALLOCATE(this)

    this%start_step = VERLET_START
    this%final_step = VERLET_FINISH

    call this%add(VERLET_UPDATE_POS)
    call this%add(UPDATE_INTERACTIONS)
    call this%add(VERLET_COMPUTE_ACC)
    call this%add(VERLET_COMPUTE_VEL)
    call this%add(FINISHED)

    ! Verlet has only one algorithmic step
    this%algo_steps = 1

    call this%parse_td_variables(namespace)

    POP_SUB(propagator_verlet_constructor)
  end function propagator_verlet_constructor

end module propagator_verlet_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

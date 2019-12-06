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
  use global_oct_m
  use gauge_field_oct_m
  use messages_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m
  use system_abst_oct_m

  implicit none

  private
  public ::                            &
    propagator_verlet_t

  type, extends(propagator_abst_t) :: propagator_verlet_t
    private

    contains
    procedure init => propagator_verlet_init
  end type propagator_verlet_t

contains

  ! ---------------------------------------------------------
  subroutine propagator_verlet_init(prop, time, dt)
    class(propagator_verlet_t),  intent(inout) :: prop
    FLOAT,                       intent(in)    :: time
    FLOAT,                       intent(in)    :: dt

    PUSH_SUB(propagator_verlet_init)

    call prop%list%add_node(VERLET_UPDATE_POS)
    call prop%list%add_node(VERLET_SYNC_DT)
    call prop%list%add_node(UPDATE_INTERACTION)
    call prop%list%add_node(VERLET_COMPUTE_ACC)
    call prop%list%add_node(VERLET_COMPUTE_VEL)
    call prop%list%add_node(FINISHED)

    prop%internal_time = time
    prop%dt = dt

    POP_SUB(propagator_verlet_init)
  end subroutine propagator_verlet_init

end module propagator_verlet_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

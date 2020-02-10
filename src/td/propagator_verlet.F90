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

  type, abstract, extends(propagator_abst_t) :: propagator_verlet_t
    private
    class(system_abst_t), pointer, public :: system
  contains
    procedure :: init => propagator_verlet_init
    procedure :: do_td_op => propagator_verlet_do_td
    procedure(system_sync_dt), deferred :: sync_dt
    procedure(system_update_pos), deferred :: update_pos
    procedure(system_compute_acc), deferred :: compute_acc
    procedure(system_compute_vel), deferred :: compute_vel
  end type propagator_verlet_t

  abstract interface
    subroutine system_sync_dt(this)
      import propagator_verlet_t
      class(propagator_verlet_t), intent(inout) :: this
    end subroutine system_sync_dt

    subroutine system_update_pos(this)
      import propagator_verlet_t
      class(propagator_verlet_t), intent(inout) :: this
    end subroutine system_update_pos

    subroutine system_compute_acc(this)
      import propagator_verlet_t
      class(propagator_verlet_t), intent(inout) :: this
    end subroutine system_compute_acc

    subroutine system_compute_vel(this)
      import propagator_verlet_t
      class(propagator_verlet_t), intent(inout) :: this
    end subroutine system_compute_vel
  end interface

contains

  ! ---------------------------------------------------------
  subroutine propagator_verlet_init(this, time, dt, system)
    class(propagator_verlet_t), intent(inout) :: this
    FLOAT, intent(in)    :: time
    FLOAT, intent(in)    :: dt
    class(system_abst_t), target, intent(inout) :: system

    PUSH_SUB(propagator_verlet_init)

    call this%list%add_node(VERLET_UPDATE_POS)
    call this%list%add_node(VERLET_SYNC_DT)
    call this%list%add_node(UPDATE_INTERACTIONS)
    call this%list%add_node(VERLET_COMPUTE_ACC)
    call this%list%add_node(VERLET_COMPUTE_VEL)
    call this%list%add_node(FINISHED)

    this%internal_time = time
    this%dt = dt
    this%system => system

    POP_SUB(propagator_verlet_init)
  end subroutine propagator_verlet_init

  ! ---------------------------------------------------------
  subroutine propagator_verlet_do_td(this, tdop)
    class(propagator_verlet_t),  intent(inout) :: this
    integer,               intent(in)    :: tdop

    integer :: iint

    PUSH_SUB(propagator_verlet_do_td)

    select case(tdop)
    case(VERLET_SYNC_DT)
      if (debug%info) then
        message(1) = "Debug: Propagation step - Synchronizing time for " + trim(this%system%namespace%get())
        call messages_info(1)
      end if

      call this%sync_dt()
      call this%list%next()

    case(VERLET_UPDATE_POS)
      if (debug%info) then
        message(1) = "Debug: Propagation step - Updating positions for " + trim(this%system%namespace%get())
        call messages_info(1)
      end if

      call this%update_pos()
      call this%list%next()

    case(VERLET_COMPUTE_ACC)
      if (debug%info) then
        message(1) = "Debug: Propagation step - Computing acceleration for " + trim(this%system%namespace%get())
        call messages_info(1)
      end if

      call this%compute_acc()
      call this%list%next()

    case(VERLET_COMPUTE_VEL)
      if (debug%info) then
        message(1) = "Debug: Propagation step - Computing velocity for " + trim(this%system%namespace%get())
        call messages_info(1)
      end if

      call this%compute_vel()
      call this%list%next()

    case default
      message(1) = "Unsupported TD tdop."
      call messages_fatal(1, namespace=this%system%namespace)
    end select

    POP_SUB(propagator_verlet_do_td)
  end subroutine propagator_verlet_do_td
end module propagator_verlet_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

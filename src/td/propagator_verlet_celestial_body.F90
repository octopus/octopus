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

module propagator_verlet_celestial_oct_m
  use celestial_body_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m
  use propagator_verlet_oct_m

  implicit none

  private
  public ::                            &
    propagator_verlet_celestial_t

  type, extends(propagator_verlet_t) :: propagator_verlet_celestial_t
    private
    type(celestial_body_t), pointer :: system
  contains
    procedure :: do_td_op => celestial_body_do_td
  end type propagator_verlet_celestial_t

  interface propagator_verlet_celestial_t
    procedure propagator_verlet_celestial_init
  end interface propagator_verlet_celestial_t

contains

  ! ---------------------------------------------------------
  type(propagator_verlet_celestial_t) function propagator_verlet_celestial_init(time, dt, system) result(this)
    FLOAT, intent(in)    :: time
    FLOAT, intent(in)    :: dt
    type(celestial_body_t), target, intent(inout) :: system

    PUSH_SUB(propagator_verlet_celestial_init)

    call this%init_steps()

    this%internal_time = time
    this%dt = dt
    this%system => system

    POP_SUB(propagator_verlet_celestial_init)
  end function propagator_verlet_celestial_init

  ! ---------------------------------------------------------
  subroutine celestial_body_do_td(this, tdop)
    class(propagator_verlet_celestial_t),  intent(inout) :: this
    integer,               intent(in)    :: tdop

    integer :: iint

    PUSH_SUB(celestial_body_do_td)

    associate(system => this%system)

    select case(tdop)
    case(VERLET_SYNC_DT)
      if (debug%info) then
        message(1) = "Debug: Propagation step - Synchronizing time for " + trim(system%namespace%get())
        call messages_info(1)
      end if

      this%internal_time = this%internal_time + this%dt
      call this%list%next()

    case(VERLET_UPDATE_POS)
      if (debug%info) then
        message(1) = "Debug: Propagation step - Updating positions for " + trim(system%namespace%get())
        call messages_info(1)
      end if

      system%acc(1:system%space%dim) = system%tot_force(1:system%space%dim)
      system%pos(1:system%space%dim) = system%pos(1:system%space%dim) + this%dt * system%vel(1:system%space%dim) &
                                   + M_HALF * this%dt**2 * system%tot_force(1:system%space%dim)
      call this%list%next()

    case(VERLET_COMPUTE_ACC)
      if (debug%info) then
        message(1) = "Debug: Propagation step - Computing acceleration for " + trim(system%namespace%get())
        call messages_info(1)
      end if

      call system%compute_total_force()
      call this%list%next()

    case(VERLET_COMPUTE_VEL)
      if (debug%info) then
        message(1) = "Debug: Propagation step - Computing velocity for " + trim(system%namespace%get())
        call messages_info(1)
      end if

      system%vel(1:system%space%dim) = system%vel(1:system%space%dim) + &
         M_HALF * this%dt * (system%acc(1:system%space%dim) + system%tot_force(1:system%space%dim))
      call this%list%next()

    case default
      message(1) = "Unsupported TD tdop."
      call messages_fatal(1, namespace=system%namespace)
    end select
    end associate

    POP_SUB(celestial_body_do_td)
  end subroutine celestial_body_do_td

end module propagator_verlet_celestial_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

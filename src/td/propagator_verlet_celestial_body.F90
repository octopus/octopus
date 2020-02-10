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
  use propagator_verlet_oct_m

  implicit none

  private
  public ::                            &
    propagator_verlet_celestial_t

  type, extends(propagator_verlet_t) :: propagator_verlet_celestial_t
    private
    type(celestial_body_t), pointer :: celestial_body
  contains
    procedure :: sync_dt => celestial_body_sync_dt
    procedure :: update_pos => celestial_body_update_pos
    procedure :: compute_acc => celestial_body_compute_acc
    procedure :: compute_vel => celestial_body_compute_vel
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

    call this%init(time, dt, system)

    ! also save a pointer to the subclass to avoid many type select statements
    this%celestial_body => system

    POP_SUB(propagator_verlet_celestial_init)
  end function propagator_verlet_celestial_init

  ! ---------------------------------------------------------
  subroutine celestial_body_sync_dt(this)
    class(propagator_verlet_celestial_t), intent(inout) :: this

    PUSH_SUB(celestial_body_sync_dt)

    this%internal_time = this%internal_time + this%dt

    POP_SUB(celestial_body_sync_dt)
  end subroutine celestial_body_sync_dt

  ! ---------------------------------------------------------
  subroutine celestial_body_update_pos(this)
    class(propagator_verlet_celestial_t), intent(inout) :: this

    PUSH_SUB(celestial_body_update_pos)

    associate(system => this%celestial_body)
      system%acc(1:system%space%dim) = system%tot_force(1:system%space%dim)
      system%pos(1:system%space%dim) = system%pos(1:system%space%dim) + this%dt * system%vel(1:system%space%dim) &
                                   + M_HALF * this%dt**2 * system%tot_force(1:system%space%dim)
    end associate

    POP_SUB(celestial_body_update_pos)
  end subroutine celestial_body_update_pos

  ! ---------------------------------------------------------
  subroutine celestial_body_compute_acc(this)
    class(propagator_verlet_celestial_t), intent(inout) :: this

    PUSH_SUB(celestial_body_compute_acc)

    call this%celestial_body%compute_total_force()

    POP_SUB(celestial_body_compute_acc)
  end subroutine celestial_body_compute_acc

  ! ---------------------------------------------------------
  subroutine celestial_body_compute_vel(this)
    class(propagator_verlet_celestial_t), intent(inout) :: this

    PUSH_SUB(celestial_body_compute_vel)

    associate(system => this%celestial_body)
      system%vel(1:system%space%dim) = system%vel(1:system%space%dim) + &
         M_HALF * this%dt * (system%acc(1:system%space%dim) + system%tot_force(1:system%space%dim))
    end associate

    POP_SUB(celestial_body_compute_vel)
  end subroutine celestial_body_compute_vel

end module propagator_verlet_celestial_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

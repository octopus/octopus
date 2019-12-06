!! Copyright (C) 2019 N. Tancogne-Dejean
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

module celestial_body_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m
  use system_abst_oct_m

  implicit none

  private
  public ::               &
    celestial_body_t

  type, extends(system_abst_t) :: celestial_body_t
    private

    FLOAT, public :: mass
    FLOAT, public :: pos(1:MAX_DIM)
    FLOAT, public :: vel(1:MAX_DIM)
    FLOAT, public :: acc(1:MAX_DIM)
    FLOAT, public :: tot_force(1:MAX_DIM)

    FLOAT, allocatable :: forces(:,:)
    integer :: ndim

    class(propagator_abst_t), pointer :: prop

  contains
    procedure :: do_td_operation => celestial_body_do_td
    procedure :: pull_interaction => celestial_body_pull
    procedure :: get_needed_quantity => celestial_body_needed_quantity
    procedure :: init => celestial_body_init
    procedure :: set_propagator => celestial_body_set_prop
    procedure :: allocate_receiv_structure => celestial_body_alloc_receiver
  end type celestial_body_t

contains

  subroutine celestial_body_init(sys)
    class(celestial_body_t),  intent(inout) :: sys

    PUSH_SUB(celestial_body_init)

    sys%ndim = 2

    POP_SUB(celestial_body_init)

  end subroutine celestial_body_init

  ! ---------------------------------------------------------
  integer function celestial_body_needed_quantity(sys)
    class(celestial_body_t),  intent(in) :: sys

    PUSH_SUB(celestial_body_needed_quantity)

    celestial_body_needed_quantity = FORCE

    POP_SUB(celestial_body_needed_quantity)

  end function celestial_body_needed_quantity

  ! ---------------------------------------------------------
  subroutine celestial_body_set_prop(sys, prop)
    class(celestial_body_t),          intent(inout) :: sys
    class(propagator_abst_t), target, intent(in)    :: prop

    PUSH_SUB(celestial_body_set_prop)

    sys%prop => prop

    POP_SUB(celestial_body_set_prop)

  end subroutine celestial_body_set_prop

  ! ---------------------------------------------------------
  subroutine celestial_body_alloc_receiver(this)
    class(celestial_body_t),          intent(inout) :: this

    PUSH_SUB(celestial_body_alloc_receiver)

    SAFE_ALLOCATE(this%forces(1:this%ndim, 1:this%nb_partners)) 

    POP_SUB(celestial_body_alloc_receiver)

  end subroutine celestial_body_alloc_receiver


  ! ---------------------------------------------------------
  subroutine celestial_body_do_td(this, operation)
    class(celestial_body_t),  intent(inout) :: this
    integer,               intent(in)    :: operation

    integer :: ip

    PUSH_SUB(celestial_body_do_td)

    select case(operation)
    case(VERLET_SYNC_DT)

      this%prop%internal_time = this%prop%internal_time + this%prop%dt
      call this%prop%list%next()

    case(VERLET_UPDATE_POS)

      this%acc(1:this%ndim) = this%tot_force(1:this%ndim)
      this%pos(1:this%ndim) = this%pos(1:this%ndim) + this%prop%dt * this%vel(1:this%ndim) &
                                   + M_HALF * this%prop%dt**2 * this%tot_force(1:this%ndim)
      call this%prop%list%next()

    case(VERLET_COMPUTE_ACC)

      !We sum the forces from the different partners
      this%tot_force(1:this%ndim) = M_ZERO
      do ip = 1, this%nb_partners
        this%tot_force(1:this%ndim) = this%tot_force(1:this%ndim) + this%forces(1:this%ndim, ip)
      end do
      call this%prop%list%next()

    case(VERLET_COMPUTE_VEL)

      this%vel(1:this%ndim) = this%vel(1:this%ndim) + &
         M_HALF * this%prop%dt * (this%acc(1:this%ndim) + this%tot_force(1:this%ndim))
      call this%prop%list%next()

    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1)
    end select

   POP_SUB(celestial_body_do_td)

  end subroutine celestial_body_do_td

   ! ---------------------------------------------------------
  subroutine celestial_body_pull(sys, remote, interaction, partner_index)
    class(celestial_body_t),  intent(inout) :: sys
    class(system_abst_t),  intent(inout) :: remote
    integer,               intent(in)    :: interaction
    integer,               intent(in)    :: partner_index

    FLOAT :: GG
    FLOAT :: dist

    PUSH_SUB(celestial_body_pull)

    GG = CNST(1.0)

    select case(interaction)
    case(FORCE)

      select type(remote)
      class is(celestial_body_t)

        dist = sqrt(sum((remote%pos(1:sys%ndim)-sys%pos(1:sys%ndim))**2))
        sys%forces(1:sys%ndim, partner_index) = (remote%pos(1:sys%ndim)-sys%pos(1:sys%ndim)) &
                                                  * GG * sys%mass * remote%mass / dist**3

      class default
        message(1) = "Unsupported partner class for force interaction"
        call messages_fatal(1)
      end select

    case default
      message(1) = "Unsupported pull interaction."
      call messages_fatal(1)
    end select

   POP_SUB(celestial_body_pull)

  end subroutine celestial_body_pull


end module celestial_body_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

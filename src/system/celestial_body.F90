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
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m
  use space_oct_m
  use system_abst_oct_m

  implicit none

  private
  public ::               &
    celestial_body_t

  type, extends(system_abst_t) :: celestial_body_t
    private

    FLOAT :: mass
    FLOAT, public :: pos(1:MAX_DIM)
    FLOAT, public :: vel(1:MAX_DIM)
    FLOAT, public :: acc(1:MAX_DIM)
    FLOAT, public :: tot_force(1:MAX_DIM)

    FLOAT, allocatable :: forces(:,:)
    type(space_t) :: space

    class(propagator_abst_t), pointer :: prop

  contains
    procedure :: do_td_operation => celestial_body_do_td
    procedure :: pull_interaction => celestial_body_pull
    procedure :: get_needed_quantity => celestial_body_needed_quantity
    procedure :: init => celestial_body_init
    procedure :: set_propagator => celestial_body_set_prop
    procedure :: allocate_receiv_structure => celestial_body_alloc_receiver
    procedure :: gravitational_interaction => celestial_body_gravitational_interaction
  end type celestial_body_t

contains

  subroutine celestial_body_init(sys, namespace)
    class(celestial_body_t),  intent(inout) :: sys
    type(namespace_t),        intent(in)    :: namespace

    integer :: n_rows, idir
    type(block_t) :: blk

    PUSH_SUB(celestial_body_init)

    sys%namespace = namespace

    call space_init(sys%space, namespace)

    !%Variable CelestialBodyMass
    !%Type float
    !%Section CelestialDynamics
    !%Description
    !% Mass of celestial body.
    !%End
    call parse_variable(namespace, 'CelestialBodyMass', M_ONE, sys%mass)

    !%Variable CelestialBodyInitialPosition
    !%Type block
    !%Section CelestialDynamics
    !%Description
    !% Initial position of celestial body.
    !%End
    sys%pos = M_ZERO
    if (parse_block(namespace, 'CelestialBodyInitialPosition', blk) == 0) then
      n_rows = parse_block_n(blk)
      if (n_rows > 1) call  messages_input_error('CelestialBodyInitialPosition')

      do idir = 1, sys%space%dim
        call parse_block_float(blk, 0, idir - 1, sys%pos(idir))
      end do
      call parse_block_end(blk)
    end if

    !%Variable CelestialBodyInitialVelocity
    !%Type block
    !%Section CelestialDynamics
    !%Description
    !% Initial velocity of celestial body.
    !%End
    sys%vel = M_ZERO
    if (parse_block(namespace, 'CelestialBodyInitialVelocity', blk) == 0) then
      n_rows = parse_block_n(blk)
      if (n_rows > 1) call  messages_input_error('CelestialBodyInitialVelocity')
      do idir = 1, sys%space%dim
        call parse_block_float(blk, 0, idir - 1, sys%vel(idir))
      end do
      call parse_block_end(blk)
    end if

    sys%acc = M_ZERO
    sys%tot_force = M_ZERO

    POP_SUB(celestial_body_init)
  end subroutine celestial_body_init

  ! ---------------------------------------------------------
  integer function celestial_body_needed_quantity(this)
    class(celestial_body_t),  intent(in) :: this

    PUSH_SUB(celestial_body_needed_quantity)

    celestial_body_needed_quantity = FORCE

    POP_SUB(celestial_body_needed_quantity)
  end function celestial_body_needed_quantity

  ! ---------------------------------------------------------
  subroutine celestial_body_set_prop(this, prop)
    class(celestial_body_t),          intent(inout) :: this
    class(propagator_abst_t), target, intent(in)    :: prop

    PUSH_SUB(celestial_body_set_prop)

    this%prop => prop

    POP_SUB(celestial_body_set_prop)
  end subroutine celestial_body_set_prop

  ! ---------------------------------------------------------
  subroutine celestial_body_alloc_receiver(this)
    class(celestial_body_t),          intent(inout) :: this

    PUSH_SUB(celestial_body_alloc_receiver)

    SAFE_ALLOCATE(this%forces(1:this%space%dim, 1:this%nb_partners))
    this%forces = M_ZERO

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

      this%acc(1:this%space%dim) = this%tot_force(1:this%space%dim)
      this%pos(1:this%space%dim) = this%pos(1:this%space%dim) + this%prop%dt * this%vel(1:this%space%dim) &
                                   + M_HALF * this%prop%dt**2 * this%tot_force(1:this%space%dim)
      call this%prop%list%next()

    case(VERLET_COMPUTE_ACC)

      !We sum the forces from the different partners
      this%tot_force(1:this%space%dim) = M_ZERO
      do ip = 1, this%nb_partners
        this%tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim) + this%forces(1:this%space%dim, ip)
      end do
      this%tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim) / this%mass
      call this%prop%list%next()

    case(VERLET_COMPUTE_VEL)

      this%vel(1:this%space%dim) = this%vel(1:this%space%dim) + &
         M_HALF * this%prop%dt * (this%acc(1:this%space%dim) + this%tot_force(1:this%space%dim))
      call this%prop%list%next()

    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

   POP_SUB(celestial_body_do_td)

  end subroutine celestial_body_do_td

   ! ---------------------------------------------------------
  subroutine celestial_body_pull(this, remote, interaction, partner_index)
    class(celestial_body_t), intent(inout) :: this
    class(system_abst_t),    intent(inout) :: remote
    integer,                 intent(in)    :: interaction
    integer,                 intent(in)    :: partner_index

    PUSH_SUB(celestial_body_pull)

    select case(interaction)
    case(FORCE)

      select type(remote)
      class is(celestial_body_t)
        call remote%gravitational_interaction(this%pos, this%mass, this%forces(1:this%space%dim, partner_index))
      class default
        message(1) = "Unsupported partner class for force interaction"
        call messages_fatal(1, namespace=this%namespace)
      end select

    case default
      message(1) = "Unsupported pull interaction."
      call messages_fatal(1, namespace=this%namespace)
    end select

    POP_SUB(celestial_body_pull)
  end subroutine celestial_body_pull

  ! ---------------------------------------------------------
  subroutine celestial_body_gravitational_interaction(this, position, mass, force)
    class(celestial_body_t), intent(in)  :: this
    FLOAT,                   intent(in)  :: position(this%space%dim)
    FLOAT,                   intent(in)  :: mass
    FLOAT,                   intent(out) :: force(this%space%dim)

    FLOAT, parameter :: GG = CNST(6.67430e-11)
    FLOAT :: dist3

    PUSH_SUB(celestial_body_gravitational_interaction)

    dist3 = sum((this%pos(1:this%space%dim) - position(1:this%space%dim))**2)**(M_THREE/M_TWO)

    force(1:this%space%dim) = (this%pos(1:this%space%dim) - position(1:this%space%dim)) / dist3 * (GG * mass * this%mass)

    POP_SUB(celestial_body_force_interaction)
  end subroutine celestial_body_gravitational_interaction

end module celestial_body_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

!! Copyright (C) 2019 N. Tancogne-Dejean
!! Copyright (C) 2020 M. Oliveira
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
  use interaction_abst_oct_m
  use interaction_gravity_oct_m
  use linked_list_oct_m
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
    type(linked_list_t) :: interactions

    type(space_t) :: space

    class(propagator_abst_t), pointer :: prop

  contains
    procedure :: add_interaction_partner => celestial_body_add_interaction_partner
    procedure :: has_interaction => celestial_body_has_interaction
    procedure :: do_td_operation => celestial_body_do_td
    procedure :: update_interactions => celestial_body_update_interactions
    procedure :: update_interaction_as_partner => celestial_body_update_interaction_as_partner
    procedure :: set_propagator => celestial_body_set_prop
    procedure :: write_td_info => celestial_body_write_td_info
    final :: celestial_body_finalize
  end type celestial_body_t

  interface celestial_body_t
    procedure celestial_body_init
  end interface celestial_body_t

contains

  ! ---------------------------------------------------------
  function celestial_body_init(namespace) result(sys)
    class(celestial_body_t), pointer    :: sys
    type(namespace_t),       intent(in) :: namespace

    integer :: n_rows, idir
    type(block_t) :: blk

    PUSH_SUB(celestial_body_init)

    SAFE_ALLOCATE(sys)

    sys%namespace = namespace

    call messages_print_stress(stdout, "Celestial Body", namespace=namespace)

    call space_init(sys%space, namespace)

    !%Variable CelestialBodyMass
    !%Type float
    !%Section CelestialDynamics
    !%Description
    !% Mass of celestial body.
    !%End
    call parse_variable(namespace, 'CelestialBodyMass', M_ONE, sys%mass)
    call messages_print_var_value(stdout, 'CelestialBodyMass', sys%mass)

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
    call messages_print_var_value(stdout, 'CelestialBodyInitialPosition', sys%pos(1:sys%space%dim))

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
    call messages_print_var_value(stdout, 'CelestialBodyInitialVelocity', sys%vel(1:sys%space%dim))

    sys%acc = M_ZERO
    sys%tot_force = M_ZERO

    call messages_print_stress(stdout, namespace=namespace)

    POP_SUB(celestial_body_init)
  end function celestial_body_init

  ! ---------------------------------------------------------
  subroutine celestial_body_add_interaction_partner(this, partner)
    class(celestial_body_t), intent(inout) :: this
    class(system_abst_t),    intent(in)    :: partner

    class(interaction_gravity_t), pointer :: gravity
    type(interaction_gravity_t) :: gravity_t

    PUSH_SUB(celestial_body_add_interaction_partner)

    if (partner%has_interaction(gravity_t)) then
      gravity => interaction_gravity_t(this%space%dim, partner)
      call this%interactions%add(gravity)
    end if

    POP_SUB(celestial_body_add_interaction_partner)
  end subroutine celestial_body_add_interaction_partner

  ! ---------------------------------------------------------
  logical function celestial_body_has_interaction(this, interaction)
    class(celestial_body_t),   intent(in) :: this
    class(interaction_abst_t), intent(in) :: interaction

    PUSH_SUB(celestial_body_has_interaction)

    select type (interaction)
    type is (interaction_gravity_t)
      celestial_body_has_interaction = .true.
    class default
      celestial_body_has_interaction = .false.
    end select

    POP_SUB(celestial_body_has_interaction)
  end function celestial_body_has_interaction

  ! ---------------------------------------------------------
  subroutine celestial_body_set_prop(this, prop)
    class(celestial_body_t),          intent(inout) :: this
    class(propagator_abst_t), target, intent(in)    :: prop

    PUSH_SUB(celestial_body_set_prop)

    this%prop => prop

    POP_SUB(celestial_body_set_prop)
  end subroutine celestial_body_set_prop

  ! ---------------------------------------------------------
  subroutine celestial_body_do_td(this, operation)
    class(celestial_body_t),  intent(inout) :: this
    integer,               intent(in)    :: operation

    class(*), pointer :: interaction

    PUSH_SUB(celestial_body_do_td)

    select case(operation)
    case(VERLET_SYNC_DT)
      if (debug%info) then
        message(1) = "Debug: Propagation step - Synchronizing time for " + trim(this%namespace%get())
        call messages_info(1)
      end if

      this%prop%internal_time = this%prop%internal_time + this%prop%dt
      call this%prop%list%next()

    case(VERLET_UPDATE_POS)
      if (debug%info) then
        message(1) = "Debug: Propagation step - Updating positions for " + trim(this%namespace%get())
        call messages_info(1)
      end if

      this%acc(1:this%space%dim) = this%tot_force(1:this%space%dim)
      this%pos(1:this%space%dim) = this%pos(1:this%space%dim) + this%prop%dt * this%vel(1:this%space%dim) &
                                   + M_HALF * this%prop%dt**2 * this%tot_force(1:this%space%dim)
      call this%prop%list%next()

    case(VERLET_COMPUTE_ACC)
      if (debug%info) then
        message(1) = "Debug: Propagation step - Computing acceleration for " + trim(this%namespace%get())
        call messages_info(1)
      end if

      !We sum the forces from the different partners
      this%tot_force(1:this%space%dim) = M_ZERO
      call this%interactions%rewind()
      do while (this%interactions%has_more_values())
        interaction => this%interactions%current()
        select type (interaction)
        type is (interaction_gravity_t)
          this%tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim) + interaction%force(1:this%space%dim)
        class default
          message(1) = "Unknown interaction by the celestial body " + this%namespace%get()
          call messages_fatal(1)
        end select
        call this%interactions%next()
      end do
      this%tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim) / this%mass
      call this%prop%list%next()

    case(VERLET_COMPUTE_VEL)
      if (debug%info) then
        message(1) = "Debug: Propagation step - Computing velocity for " + trim(this%namespace%get())
        call messages_info(1)
      end if

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
  subroutine celestial_body_update_interactions(this)
    class(celestial_body_t), intent(inout) :: this

    class(*), pointer :: interaction

    PUSH_SUB(celestial_body_update_interactions)

    call this%interactions%rewind()
    do while (this%interactions%has_more_values())
      interaction => this%interactions%current()
      select type (interaction)
      type is (interaction_gravity_t)
        call interaction%update(this%mass, this%pos)
      class default
        message(1) = "Unknown interaction by the celestial body " + this%namespace%get()
        call messages_fatal(1)
      end select
      call this%interactions%next()
    end do

    POP_SUB(celestial_body_update_interactions)
  end subroutine celestial_body_update_interactions

  ! ---------------------------------------------------------
  subroutine celestial_body_update_interaction_as_partner(this, interaction)
    class(celestial_body_t),   intent(in)    :: this
    class(interaction_abst_t), intent(inout) :: interaction

    PUSH_SUB(celestial_body_update_interaction_as_partner)

    select type (interaction)
    type is (interaction_gravity_t)
      interaction%partner_mass = this%mass
      interaction%partner_pos = this%pos

    class default
      message(1) = "Unknown interaction by the celestial body " + this%namespace%get()
      call messages_fatal(1)
    end select

    POP_SUB(celestial_body_update_interaction_as_partner)
  end subroutine celestial_body_update_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine celestial_body_write_td_info(this)
    class(celestial_body_t), intent(in) :: this

    integer :: idir
    character(len=20) :: fmt

    PUSH_SUB(celestial_body_write_td_info)

    write(message(1),'(2X,A,1X,A)') "Celestial body:", trim(this%namespace%get())

    write(fmt,'("(4X,A,1X,",I2,"e14.6)")') this%space%dim
    write(message(2),fmt) "Coordinates: ", (this%pos(idir), idir = 1, this%space%dim)
    write(message(3),fmt) "Velocity:    ", (this%vel(idir), idir = 1, this%space%dim)
    write(message(4),fmt) "Acceleration:", (this%acc(idir), idir = 1, this%space%dim)
    call messages_info(4)

    POP_SUB(celestial_body_write_td_info)
  end subroutine celestial_body_write_td_info

  ! ---------------------------------------------------------
  subroutine celestial_body_finalize(this)
    type(celestial_body_t), intent(inout) :: this

    class(*), pointer :: interaction

    PUSH_SUB(celestial_body_finalize)

    call this%interactions%rewind()
    do while (this%interactions%has_more_values())
      interaction => this%interactions%current()
      SAFE_DEALLOCATE_P(interaction)
      call this%interactions%next()
    end do

    POP_SUB(celestial_body_finalize)
  end subroutine celestial_body_finalize

end module celestial_body_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

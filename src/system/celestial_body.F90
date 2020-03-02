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
  use io_oct_m
  use iso_c_binding  
  use linked_list_oct_m
  use list_node_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m
  use space_oct_m
  use system_abst_oct_m
  use write_iter_oct_m

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

    type(c_ptr) :: output_handle
  contains
    procedure :: add_interaction_partner => celestial_body_add_interaction_partner
    procedure :: has_interaction => celestial_body_has_interaction
    procedure :: do_td_operation => celestial_body_do_td
    procedure :: update_interactions => celestial_body_update_interactions
    procedure :: update_interaction_as_partner => celestial_body_update_interaction_as_partner
    procedure :: write_td_info => celestial_body_write_td_info
    procedure :: td_write_init => celestial_body_td_write_init
    procedure :: td_write_iter => celestial_body_td_write_iter
    procedure :: td_write_end => celestial_body_td_write_end
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
    !% Mass of celestial body in Kg.
    !%End
    call parse_variable(namespace, 'CelestialBodyMass', M_ONE, sys%mass)
    call messages_print_var_value(stdout, 'CelestialBodyMass', sys%mass)

    !%Variable CelestialBodyInitialPosition
    !%Type block
    !%Section CelestialDynamics
    !%Description
    !% Initial position of celestial body, in Km.
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
    !% Initial velocity of celestial body in Km/s.
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
  subroutine celestial_body_do_td(this, operation)
    class(celestial_body_t),  intent(inout) :: this
    integer,               intent(in)    :: operation

    type(list_node_t), pointer :: iint
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
      nullify(iint)
      do while (this%interactions%iterate(iint, interaction))
        select type (interaction)
        type is (interaction_gravity_t)
          this%tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim) + interaction%force(1:this%space%dim)
        class default
          message(1) = "Unknown interaction by the celestial body " + this%namespace%get()
          call messages_fatal(1)
        end select
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
    type(list_node_t), pointer :: iint

    PUSH_SUB(celestial_body_update_interactions)

    nullify(iint)
    do while (this%interactions%iterate(iint, interaction))
      select type (interaction)
      type is (interaction_gravity_t)
        call interaction%update(this%mass, this%pos)
      class default
        message(1) = "Unknown interaction by the celestial body " + this%namespace%get()
        call messages_fatal(1)
      end select
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
    write(message(5),'(4x,A,I8.7)') 'Clock tick: ', this%clock%get_tick()
    call messages_info(5)

    POP_SUB(celestial_body_write_td_info)
  end subroutine celestial_body_write_td_info

  ! ---------------------------------------------------------
  subroutine celestial_body_td_write_init(this, dt)
    class(celestial_body_t), intent(inout) :: this
    FLOAT,                   intent(in)    :: dt

    PUSH_SUB(celestial_body_td_write_init)

    call io_mkdir('td.general', this%namespace)
    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_init(this%output_handle, 0, dt, trim(io_workpath("td.general/coordinates", this%namespace)))
    end if

    POP_SUB(celestial_body_td_write_init)
  end subroutine celestial_body_td_write_init

  ! ---------------------------------------------------------
  subroutine celestial_body_td_write_iter(this, iter)
    class(celestial_body_t), intent(inout) :: this
    integer,                 intent(in)    :: iter

    integer :: idir
    character(len=50) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(celestial_body_td_write_iter)

    if(iter == 0) then
      ! header
      call write_iter_clear(this%output_handle)
      call write_iter_string(this%output_handle,'################################################################################')
      call write_iter_nl(this%output_handle)
      call write_iter_string(this%output_handle,'# HEADER')
      call write_iter_nl(this%output_handle)

      ! first line: column names
      call write_iter_header_start(this%output_handle)

      do idir = 1, this%space%dim
        write(aux, '(a2,i3,a1)') 'x(', idir, ')'
        call write_iter_header(this%output_handle, aux)
      end do
      do idir = 1, this%space%dim
        write(aux, '(a2,i3,a1)') 'v(', idir, ')'
        call write_iter_header(this%output_handle, aux)
      end do
      call write_iter_nl(this%output_handle)

      ! second line: units
      call write_iter_string(this%output_handle, '#[Iter n.]')
      call write_iter_header(this%output_handle, '[seconds]')
      call write_iter_string(this%output_handle, 'Positions in Km, Velocities in Km/s')
      call write_iter_nl(this%output_handle)

      call write_iter_string(this%output_handle,'################################################################################')
      call write_iter_nl(this%output_handle)
    end if

    call write_iter_start(this%output_handle)

    call write_iter_double(this%output_handle, this%pos, this%space%dim)
    call write_iter_double(this%output_handle, this%vel, this%space%dim)
    call write_iter_nl(this%output_handle)
    
    POP_SUB(celestial_body_td_write_iter)
  end subroutine celestial_body_td_write_iter

  ! ---------------------------------------------------------
  subroutine celestial_body_td_write_end(this)
    class(celestial_body_t), intent(inout) :: this

    PUSH_SUB(celestial_body_td_write_end)

    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_end(this%output_handle)
    end if

    POP_SUB(celestial_body_td_write_end)
  end subroutine celestial_body_td_write_end

  ! ---------------------------------------------------------
  subroutine celestial_body_finalize(this)
    type(celestial_body_t), intent(inout) :: this

    type(list_node_t), pointer :: iint
    class(*), pointer :: interaction

    PUSH_SUB(celestial_body_finalize)

    nullify(iint)
    do while (this%interactions%iterate(iint, interaction))
      SAFE_DEALLOCATE_P(interaction)
    end do

    POP_SUB(celestial_body_finalize)
  end subroutine celestial_body_finalize

end module celestial_body_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

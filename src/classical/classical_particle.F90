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

module classical_particle_oct_m
  use algorithm_oct_m
  use classical_particles_oct_m
  use clock_oct_m
  use global_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use gravity_oct_m
  use io_oct_m
  use iso_c_binding
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use space_oct_m
  use system_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::                  &
    classical_particle_t,    &
    classical_particle_init

   type, extends(classical_particles_t) :: classical_particle_t
    type(c_ptr) :: output_handle
  contains
    procedure :: init_interaction => classical_particle_init_interaction
    procedure :: initial_conditions => classical_particle_initial_conditions
    procedure :: iteration_info => classical_particle_iteration_info
    procedure :: output_start => classical_particle_output_start
    procedure :: output_write => classical_particle_output_write
    procedure :: output_finish => classical_particle_output_finish
    procedure :: update_quantity => classical_particle_update_quantity
    procedure :: update_exposed_quantity => classical_particle_update_exposed_quantity
    procedure :: init_interaction_as_partner => classical_particle_init_interaction_as_partner
    procedure :: copy_quantities_to_interaction => classical_particle_copy_quantities_to_interaction
    final :: classical_particle_finalize
  end type classical_particle_t

  interface classical_particle_t
    procedure classical_particle_constructor
  end interface classical_particle_t

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function classical_particle_constructor(namespace) result(sys)
    class(classical_particle_t), pointer    :: sys
    type(namespace_t),           intent(in) :: namespace

    PUSH_SUB(classical_particle_constructor)

    SAFE_ALLOCATE(sys)

    call classical_particle_init(sys, namespace)

    POP_SUB(classical_particle_constructor)
  end function classical_particle_constructor

  ! ---------------------------------------------------------
  !> The init routine is a module level procedure
  !! This has the advantage that different classes can have different
  !! signatures for the initialization routines because they are not
  !! type-bound and thus also not inherited.
  ! ---------------------------------------------------------
  subroutine classical_particle_init(this, namespace)
    class(classical_particle_t), intent(inout) :: this
    type(namespace_t),           intent(in)    :: namespace

    PUSH_SUB(classical_particle_init)

    this%namespace = namespace

    call space_init(this%space, namespace)
    if (this%space%periodic_dim > 0) then
      call messages_not_implemented('Classical particle for periodic systems')
    end if

    call messages_print_stress(stdout, "Classical Particle", namespace=namespace)

    call classical_particles_init(this, 1)

    !%Variable ParticleMass
    !%Type float
    !%Section ClassicalParticles
    !%Description
    !% Mass of classical particle in Kg.
    !%End
    call parse_variable(namespace, 'ParticleMass', M_ONE, this%mass(1))
    call messages_print_var_value(stdout, 'ParticleMass', this%mass(1))

    call this%supported_interactions%add(GRAVITY)
    call this%supported_interactions_as_partner%add(GRAVITY)

    call messages_print_stress(stdout, namespace=namespace)

    POP_SUB(classical_particle_init)
  end subroutine classical_particle_init

  ! ---------------------------------------------------------
  subroutine classical_particle_init_interaction(this, interaction)
    class(classical_particle_t), target, intent(inout) :: this
    class(interaction_t),                intent(inout) :: interaction

    PUSH_SUB(classical_particle_init_interaction)

    select type (interaction)
    type is (gravity_t)
      call interaction%init(this%space%dim, 1, this%quantities, this%mass, this%pos)
    class default
      call classical_particles_init_interaction(this, interaction)
    end select

    POP_SUB(classical_particle_init_interaction)
  end subroutine classical_particle_init_interaction

  ! ---------------------------------------------------------
  subroutine classical_particle_initial_conditions(this, from_scratch)
    class(classical_particle_t), intent(inout) :: this
    logical,                 intent(in)    :: from_scratch

    integer :: n_rows, idir
    type(block_t) :: blk

    PUSH_SUB(classical_particle_initial_conditions)

    if (.not. from_scratch) then
      message(1) = "Classical particle propagation is currently only allowed from scratch"
      call messages_fatal(1, namespace=this%namespace)
    end if

    !%Variable ParticleInitialPosition
    !%Type block
    !%Section ClassicalParticles
    !%Description
    !% Initial position of classical particle, in Km.
    !%End
    this%pos = M_ZERO
    if (parse_block(this%namespace, 'ParticleInitialPosition', blk) == 0) then
      n_rows = parse_block_n(blk)
      if (n_rows > 1) call  messages_input_error(this%namespace, 'ParticleInitialPosition')

      do idir = 1, this%space%dim
        call parse_block_float(blk, 0, idir - 1, this%pos(idir, 1))
      end do
      call parse_block_end(blk)
    end if
    call messages_print_var_value(stdout, 'ParticleInitialPosition', this%pos(1:this%space%dim, 1))

    !%Variable ParticleInitialVelocity
    !%Type block
    !%Section ClassicalParticles
    !%Description
    !% Initial velocity of classical particle in Km/s.
    !%End
    this%vel = M_ZERO
    if (parse_block(this%namespace, 'ParticleInitialVelocity', blk) == 0) then
      n_rows = parse_block_n(blk)
      if (n_rows > 1) call  messages_input_error(this%namespace, 'ParticleInitialVelocity')
      do idir = 1, this%space%dim
        call parse_block_float(blk, 0, idir - 1, this%vel(idir, 1))
      end do
      call parse_block_end(blk)
    end if
    call messages_print_var_value(stdout, 'ParticleInitialVelocity', this%vel(1:this%space%dim, 1))

    POP_SUB(classical_particle_initial_conditions)
  end subroutine classical_particle_initial_conditions

  ! ---------------------------------------------------------
  subroutine classical_particle_iteration_info(this)
    class(classical_particle_t), intent(in) :: this

    integer :: idir
    character(len=20) :: fmt

    PUSH_SUB(classical_particle_iteration_info)

    write(message(1),'(2X,A,1X,A)') "Classical particle:", trim(this%namespace%get())

    write(fmt,'("(4X,A,1X,",I2,"e14.6)")') this%space%dim
    write(message(2),fmt) "Coordinates: ", (this%pos(idir, 1), idir = 1, this%space%dim)
    write(message(3),fmt) "Velocity:    ", (this%vel(idir, 1), idir = 1, this%space%dim)
    write(message(4),fmt) "Force:       ", (this%tot_force(idir, 1), idir = 1, this%space%dim)
    write(message(5),'(4x,A,I8.7)')  'Clock tick:      ', this%clock%get_tick()
    write(message(6),'(4x,A,e14.6)') 'Simulation time: ', this%clock%time()
    call messages_info(6)

    POP_SUB(classical_particle_iteration_info)
  end subroutine classical_particle_iteration_info

  ! ---------------------------------------------------------
  subroutine classical_particle_output_start(this)
    class(classical_particle_t), intent(inout) :: this

    PUSH_SUB(classical_particle_output_start)

    ! Create output handle
    call io_mkdir('td.general', this%namespace)
    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_init(this%output_handle, 0, this%prop%dt, trim(io_workpath("td.general/coordinates", this%namespace)))
    end if

    ! Output info for first iteration
    call this%output_write()

    POP_SUB(classical_particle_output_start)
  end subroutine classical_particle_output_start

  ! ---------------------------------------------------------
  subroutine classical_particle_output_finish(this)
    class(classical_particle_t), intent(inout) :: this

    PUSH_SUB(classical_particle_output_finish)

    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_end(this%output_handle)
    end if

    POP_SUB(classical_particle_output_finish)
  end subroutine classical_particle_output_finish

  ! ---------------------------------------------------------
  subroutine classical_particle_output_write(this)
    class(classical_particle_t), intent(inout) :: this

    integer :: idir
    character(len=50) :: aux
    FLOAT :: tmp(MAX_DIM)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(classical_particle_output_write)

    if (this%clock%get_tick() == 0) then
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
      do idir = 1, this%space%dim
        write(aux, '(a2,i3,a1)') 'f(', idir, ')'
        call write_iter_header(this%output_handle, aux)
      end do
      call write_iter_nl(this%output_handle)

      ! second line: units
      call write_iter_string(this%output_handle, '#[Iter n.]')
      call write_iter_header(this%output_handle, '[' // trim(units_abbrev(units_out%time)) // ']')
      call write_iter_string(this%output_handle, &
        'Position in '   // trim(units_abbrev(units_out%length))   //   &
        ', Velocity in '// trim(units_abbrev(units_out%velocity)) //   &
        ', Force in '    // trim(units_abbrev(units_out%force)))
      call write_iter_nl(this%output_handle)

      call write_iter_string(this%output_handle,'################################################################################')
      call write_iter_nl(this%output_handle)
    end if

    call write_iter_start(this%output_handle)

    ! Position
    tmp(1:this%space%dim) = units_from_atomic(units_out%length, this%pos(1:this%space%dim, 1))
    call write_iter_double(this%output_handle, tmp, this%space%dim)
    ! Velocity
    tmp(1:this%space%dim) = units_from_atomic(units_out%velocity, this%vel(1:this%space%dim, 1))
    call write_iter_double(this%output_handle, tmp, this%space%dim)
    ! Force
    tmp(1:this%space%dim) = units_from_atomic(units_out%force, this%tot_force(1:this%space%dim, 1))
    call write_iter_double(this%output_handle, tmp, this%space%dim)

    call write_iter_nl(this%output_handle)

    POP_SUB(classical_particle_output_write)
  end subroutine classical_particle_output_write

  ! ---------------------------------------------------------
  subroutine classical_particle_update_quantity(this, iq)
    class(classical_particle_t), intent(inout) :: this
    integer,                     intent(in)    :: iq

    PUSH_SUB(classical_particle_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case default
      ! Other quantities should be handled by the parent class
      call classical_particles_update_quantity(this, iq)
    end select

    POP_SUB(classical_particle_update_quantity)
  end subroutine classical_particle_update_quantity

  ! ---------------------------------------------------------
  subroutine classical_particle_update_exposed_quantity(partner, iq)
    class(classical_particle_t), intent(inout) :: partner
    integer,                     intent(in)    :: iq

    PUSH_SUB(classical_particle_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case default
      ! Other quantities should be handled by the parent class
      call classical_particles_update_exposed_quantity(partner, iq)
    end select

    POP_SUB(classical_particle_update_exposed_quantity)
  end subroutine classical_particle_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine classical_particle_init_interaction_as_partner(partner, interaction)
    class(classical_particle_t), intent(in)    :: partner
    class(interaction_t),        intent(inout) :: interaction

    PUSH_SUB(classical_particle_init_interaction_as_partner)

    select type (interaction)
    type is (gravity_t)
      interaction%partner_np = 1
      SAFE_ALLOCATE(interaction%partner_mass(1))
      SAFE_ALLOCATE(interaction%partner_pos(partner%space%dim, 1))

    class default
      ! Other interactions should be handled by the parent class
      call classical_particles_init_interaction_as_partner(partner, interaction)
    end select

    POP_SUB(classical_particle_init_interaction_as_partner)
  end subroutine classical_particle_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine classical_particle_copy_quantities_to_interaction(partner, interaction)
    class(classical_particle_t),          intent(inout) :: partner
    class(interaction_t),                 intent(inout) :: interaction

    PUSH_SUB(classical_particle_copy_quantities_to_interaction)

    select type (interaction)
    type is (gravity_t)
      interaction%partner_mass(1) = partner%mass(1)
      interaction%partner_pos(:,1) = partner%pos(:,1)

    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(classical_particle_copy_quantities_to_interaction)
  end subroutine classical_particle_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine classical_particle_finalize(this)
    type(classical_particle_t), intent(inout) :: this

    PUSH_SUB(classical_particle_finalize)

    call classical_particles_end(this)

    POP_SUB(classical_particle_finalize)
  end subroutine classical_particle_finalize

end module classical_particle_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

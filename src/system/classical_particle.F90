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
  use clock_oct_m
  use global_oct_m
  use interaction_abst_oct_m
  use interaction_gravity_oct_m
  use io_oct_m
  use iso_c_binding
  use linked_list_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m
  use quantity_oct_m
  use space_oct_m
  use system_abst_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::                 &
    classical_particle_t,   &
    classical_particle_init

  type, extends(system_abst_t) :: classical_particle_t
    private

    FLOAT :: mass
    FLOAT, public :: pos(1:MAX_DIM)
    FLOAT, public :: vel(1:MAX_DIM)
    FLOAT, public :: acc(1:MAX_DIM)
    FLOAT, public :: prev_acc(1:MAX_DIM,1) !< A storage of the prior times. At the moment we only can store one time 
    FLOAT, public :: save_pos(1:MAX_DIM) !< A storage for the SCF loops
    FLOAT, public :: save_vel(1:MAX_DIM) !< A storage for the SCF loops
    FLOAT, public :: tot_force(1:MAX_DIM)
    FLOAT, public :: prev_tot_force(1:MAX_DIM) !< Used for the SCF convergence criterium

    type(c_ptr) :: output_handle
  contains
    procedure :: add_interaction_partner => classical_particle_add_interaction_partner
    procedure :: has_interaction => classical_particle_has_interaction
    procedure :: do_td_operation => classical_particle_do_td
    procedure :: write_td_info => classical_particle_write_td_info
    procedure :: td_write_init => classical_particle_td_write_init
    procedure :: td_write_iter => classical_particle_td_write_iter
    procedure :: td_write_end => classical_particle_td_write_end
    procedure :: is_tolerance_reached => classical_particle_is_tolerance_reached
    procedure :: store_current_status => classical_particle_store_current_status
    procedure :: update_quantity => classical_particle_update_quantity
    procedure :: update_exposed_quantity => classical_particle_update_exposed_quantity
    procedure :: set_pointers_to_interaction => classical_set_pointers_to_interaction
    final :: classical_particle_finalize
  end type classical_particle_t

  interface classical_particle_t
    procedure classical_particle_constructor
  end interface classical_particle_t

contains
  ! ---------------------------------------------------------
  function classical_particle_constructor(namespace) result(sys)
    class(classical_particle_t), pointer    :: sys
    type(namespace_t),       intent(in) :: namespace

    PUSH_SUB(classical_particle_constructor)

    SAFE_ALLOCATE(sys)

    call classical_particle_init(sys, namespace)

    POP_SUB(classical_particle_constructor)
  end function classical_particle_constructor

  ! ---------------------------------------------------------
  subroutine classical_particle_init(this, namespace)
    class(classical_particle_t), intent(inout) :: this
    type(namespace_t),           intent(in)    :: namespace

    integer :: n_rows, idir
    type(block_t) :: blk

    PUSH_SUB(classical_particle_init)

    this%namespace = namespace

    call messages_print_stress(stdout, "Classicle Particle", namespace=namespace)

    call space_init(this%space, namespace)

    !%Variable ClassicleParticleMass
    !%Type float
    !%Section CelestialDynamics
    !%Description
    !% Mass of classical body in Kg.
    !%End
    call parse_variable(namespace, 'ClassicleParticleMass', M_ONE, this%mass)
    call messages_print_var_value(stdout, 'ClassicleParticleMass', this%mass)

    !%Variable ClassicleParticleInitialPosition
    !%Type block
    !%Section CelestialDynamics
    !%Description
    !% Initial position of classical body, in Km.
    !%End
    this%pos = M_ZERO
    if (parse_block(namespace, 'ClassicleParticleInitialPosition', blk) == 0) then
      n_rows = parse_block_n(blk)
      if (n_rows > 1) call  messages_input_error('ClassicleParticleInitialPosition')

      do idir = 1, this%space%dim
        call parse_block_float(blk, 0, idir - 1, this%pos(idir))
      end do
      call parse_block_end(blk)
    end if
    call messages_print_var_value(stdout, 'ClassicleParticleInitialPosition', this%pos(1:this%space%dim))

    !%Variable ClassicleParticleInitialVelocity
    !%Type block
    !%Section CelestialDynamics
    !%Description
    !% Initial velocity of classical body in Km/s.
    !%End
    this%vel = M_ZERO
    if (parse_block(namespace, 'ClassicleParticleInitialVelocity', blk) == 0) then
      n_rows = parse_block_n(blk)
      if (n_rows > 1) call  messages_input_error('ClassicleParticleInitialVelocity')
      do idir = 1, this%space%dim
        call parse_block_float(blk, 0, idir - 1, this%vel(idir))
      end do
      call parse_block_end(blk)
    end if
    call messages_print_var_value(stdout, 'ClassicleParticleInitialVelocity', this%vel(1:this%space%dim))

    this%acc = M_ZERO
    this%prev_acc = M_ZERO
    this%tot_force = M_ZERO

    this%quantities(POSITION)%internal = .true.
    this%quantities(VELOCITY)%internal = .true.

    call messages_print_stress(stdout, namespace=namespace)

    POP_SUB(classical_particle_init)
  end subroutine classical_particle_init

  ! ---------------------------------------------------------
  subroutine classical_particle_add_interaction_partner(this, partner)
    class(classical_particle_t), target, intent(inout) :: this
    class(system_abst_t),                intent(in)    :: partner

    class(interaction_gravity_t), pointer :: gravity
    type(interaction_gravity_t) :: gravity_t

    PUSH_SUB(classical_particle_add_interaction_partner)

    if (partner%has_interaction(gravity_t)) then
      gravity => interaction_gravity_t(this%space%dim, partner)
      gravity%system_mass => this%mass
      gravity%system_pos  => this%pos
      call this%interactions%add(gravity)
    end if

    POP_SUB(classical_particle_add_interaction_partner)
  end subroutine classical_particle_add_interaction_partner

  ! ---------------------------------------------------------
  logical function classical_particle_has_interaction(this, interaction)
    class(classical_particle_t),   intent(in) :: this
    class(interaction_abst_t),     intent(in) :: interaction

    PUSH_SUB(classical_particle_has_interaction)

    select type (interaction)
    type is (interaction_gravity_t)
      classical_particle_has_interaction = .true.
    class default
      classical_particle_has_interaction = .false.
    end select

    POP_SUB(classical_particle_has_interaction)
  end function classical_particle_has_interaction

  ! ---------------------------------------------------------
  subroutine classical_particle_do_td(this, operation)
    class(classical_particle_t), intent(inout) :: this
    integer,                     intent(in)    :: operation

    type(interaction_iterator_t) :: iter

    PUSH_SUB(classical_particle_do_td)

    select case(operation)
    case (VERLET_UPDATE_POS)
      this%acc(1:this%space%dim) = this%tot_force(1:this%space%dim)
      this%pos(1:this%space%dim) = this%pos(1:this%space%dim) + this%prop%dt * this%vel(1:this%space%dim) &
                                 + M_HALF * this%prop%dt**2 * this%tot_force(1:this%space%dim)

      call this%quantities(POSITION)%clock%increment()

    case (VERLET_COMPUTE_ACC)
      !Use as SCF criterium
      this%prev_tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim)

      !We sum the forces from the different partners
      this%tot_force(1:this%space%dim) = M_ZERO
      call iter%start(this%interactions)
      do while (iter%has_next())
        select type (interaction => iter%get_next_interaction())
        type is (interaction_gravity_t)
          this%tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim) + interaction%force(1:this%space%dim)
        class default
          message(1) = "Unknown interaction by the classical body " + this%namespace%get()
          call messages_fatal(1)
        end select
      end do
      this%tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim) / this%mass

    case (VERLET_COMPUTE_VEL)
      this%vel(1:this%space%dim) = this%vel(1:this%space%dim) &
                                 + M_HALF * this%prop%dt * (this%acc(1:this%space%dim) + this%tot_force(1:this%space%dim))

      call this%quantities(VELOCITY)%clock%increment()

    case (BEEMAN_PREDICT_POS)
      this%pos(1:this%space%dim) = this%pos(1:this%space%dim) + this%prop%dt * this%vel(1:this%space%dim) &
                                 + M_ONE/CNST(6.0) * this%prop%dt**2  &
                                 * (M_FOUR*this%acc(1:this%space%dim) - this%prev_acc(1:this%space%dim, 1))
      this%prev_acc(1:this%space%dim, 1) = this%acc(1:this%space%dim)
      this%acc(1:this%space%dim) = this%tot_force(1:this%space%dim)

      if (.not. this%prop%predictor_corrector) then
        call this%quantities(POSITION)%clock%increment()
      end if

    case (BEEMAN_PREDICT_VEL)
      this%vel(1:this%space%dim) = this%vel(1:this%space%dim)  &
                                 + M_ONE/CNST(6.0) * this%prop%dt * (this%acc(1:this%space%dim) &
                                 + M_TWO * this%tot_force(1:this%space%dim) - this%prev_acc(1:this%space%dim, 1))

      call this%quantities(VELOCITY)%clock%increment()

    case( BEEMAN_CORRECT_POS)
      this%pos(1:this%space%dim) = this%save_pos(1:this%space%dim) + this%prop%dt * this%save_vel(1:this%space%dim) &
                                 + M_ONE/CNST(6.0) * this%prop%dt**2  &
                                 * (M_TWO * this%acc(1:this%space%dim) + this%tot_force(1:this%space%dim))

      !We set it to the propagation time to avoid double increment
      call this%quantities(POSITION)%clock%set_time(this%prop%clock)

    case (BEEMAN_CORRECT_VEL)
      this%vel(1:this%space%dim) = this%save_vel(1:this%space%dim) &
                                 + M_HALF * this%prop%dt * (this%acc(1:this%space%dim) + this%tot_force(1:this%space%dim))

      !We set it to the propagation time to avoid double increment
      call this%quantities(VELOCITY)%clock%set_time(this%prop%clock)

    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

   POP_SUB(classical_particle_do_td)
  end subroutine classical_particle_do_td

  ! ---------------------------------------------------------
  logical function classical_particle_is_tolerance_reached(this, tol) result(converged)
    class(classical_particle_t),   intent(in)    :: this
    FLOAT,                         intent(in)    :: tol

    PUSH_SUB(classical_particle_is_tolerance_reached)

    !Here we put the criterium that acceleration change is below the tolerance
    converged = .false.
    if(sum((this%prev_tot_force(1:this%space%dim) -this%tot_force(1:this%space%dim))**2) < tol**2) then
      converged = .true.
    end if 

    if (debug%info) then
      write(message(1), '(a, e12.6, a, e12.6)') "Debug: -- Change in acceleration  ", &
          sqrt(sum((this%prev_tot_force(1:this%space%dim) - this%tot_force(1:this%space%dim))**2)), " and tolerance ", tol
      call messages_info(1)
    end if

    POP_SUB(classical_particle_is_tolerance_reached)
   end function classical_particle_is_tolerance_reached

   ! ---------------------------------------------------------
   subroutine classical_particle_store_current_status(this)
     class(classical_particle_t),   intent(inout)    :: this

     PUSH_SUB(classical_particle_store_current_status) 

     this%save_pos(1:this%space%dim) = this%pos(1:this%space%dim)
     this%save_vel(1:this%space%dim) = this%vel(1:this%space%dim)

     POP_SUB(classical_particle_store_current_status)
   end subroutine classical_particle_store_current_status

  ! ---------------------------------------------------------
  subroutine classical_particle_write_td_info(this)
    class(classical_particle_t), intent(in) :: this

    integer :: idir
    character(len=20) :: fmt

    PUSH_SUB(classical_particle_write_td_info)

    write(message(1),'(2X,A,1X,A)') "Celestial body:", trim(this%namespace%get())

    write(fmt,'("(4X,A,1X,",I2,"e14.6)")') this%space%dim
    write(message(2),fmt) "Coordinates: ", (this%pos(idir), idir = 1, this%space%dim)
    write(message(3),fmt) "Velocity:    ", (this%vel(idir), idir = 1, this%space%dim)
    write(message(4),fmt) "Acceleration:", (this%acc(idir), idir = 1, this%space%dim)
    write(message(5),'(4x,A,I8.7)') 'Clock tick: ', this%clock%get_tick()
    call messages_info(5)

    POP_SUB(classical_particle_write_td_info)
  end subroutine classical_particle_write_td_info

  ! ---------------------------------------------------------
  subroutine classical_particle_td_write_init(this, dt)
    class(classical_particle_t), intent(inout) :: this
    FLOAT,                       intent(in)    :: dt

    PUSH_SUB(classical_particle_td_write_init)

    call io_mkdir('td.general', this%namespace)
    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_init(this%output_handle, 0, dt, trim(io_workpath("td.general/coordinates", this%namespace)))
    end if

    POP_SUB(classical_particle_td_write_init)
  end subroutine classical_particle_td_write_init

  ! ---------------------------------------------------------
  subroutine classical_particle_td_write_iter(this, iter)
    class(classical_particle_t), intent(inout) :: this
    integer,                 intent(in)    :: iter

    integer :: idir
    character(len=50) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(classical_particle_td_write_iter)

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
    
    POP_SUB(classical_particle_td_write_iter)
  end subroutine classical_particle_td_write_iter

  ! ---------------------------------------------------------
  subroutine classical_particle_td_write_end(this)
    class(classical_particle_t), intent(inout) :: this

    PUSH_SUB(classical_particle_td_write_end)

    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_end(this%output_handle)
    end if

    POP_SUB(classical_particle_td_write_end)
  end subroutine classical_particle_td_write_end

  ! ---------------------------------------------------------
  logical function classical_particle_update_quantity(this, iq, clock) result(updated)
    class(classical_particle_t),   intent(inout) :: this
    integer,                       intent(in)    :: iq
    class(clock_t),                intent(in)    :: clock

    PUSH_SUB(classical_particle_update_quantity)

    if(this%quantities(iq)%clock > clock) then
      call this%quantities(iq)%clock%print()
      call clock%print()
      message(1) = "The system quantity is in advance compared to the requested clock."
      call messages_fatal(1)
    end if

    if (this%quantities(iq)%internal) then
      !Don`t do anything, this is a protected quantity. The propagator update it
      !If I have (system) a SCF propagator, this is not a problem here, as I handle the self-concistency.
      updated = .true.
    else
      select case (iq)
      case (MASS)
        !The classical body has a mass, but it does not require any update, as it does not change with time.
        updated = .true.
      case default
        message(1) = "Incompatible quantity."
        call messages_fatal(1)
      end select
    end if

    POP_SUB(classical_particle_update_quantity)
  end function classical_particle_update_quantity

 ! ---------------------------------------------------------
 logical function classical_particle_update_exposed_quantity(this, iq, clock) result(updated)
    class(classical_particle_t),   intent(inout) :: this
    integer,                       intent(in)    :: iq
    class(clock_t),                intent(in)    :: clock

    PUSH_SUB(classical_particle_update_exposed_quantity)

    if (this%quantities(iq)%clock > clock) then
      message(1) = "The partner quantity is in advance compared to the requested clock."
      call messages_fatal(1)
    end if

    if (this%quantities(iq)%internal) then
      !Don`t do anything, this is a protected quantity. The propagator update it.
      !However, it can only be used if the predictor-corrector step is done.
      if (this%prop%predictor_corrector) then
        updated = this%prop%last_step_done_tick == clock%get_tick()
      else
        updated = .true.
      end if
    else
      select case (iq)
      case (MASS)
        !The classical body has a mass, but it does not require any update, as it does not change with time.
        updated = .true.
      case default
        message(1) = "Incompatible quantity."
        call messages_fatal(1)
      end select
    end if

    POP_SUB(classical_particle_update_exposed_quantity)
  end function classical_particle_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine classical_set_pointers_to_interaction(this, inter)
    class(classical_particle_t), target,  intent(in)   :: this
    class(interaction_abst_t),           intent(inout) :: inter

    PUSH_SUB(classical_set_pointers_to_interaction)

    select type(inter)
    type is(interaction_gravity_t)
      inter%partner_mass => this%mass
      inter%partner_pos => this%pos
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(classical_set_pointers_to_interaction)
  end subroutine classical_set_pointers_to_interaction

  ! ---------------------------------------------------------
  subroutine classical_particle_finalize(this)
    type(classical_particle_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter
    class(interaction_abst_t), pointer :: interaction

    PUSH_SUB(classical_particle_finalize)

    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next_interaction()
      SAFE_DEALLOCATE_P(interaction)
    end do

    POP_SUB(classical_particle_finalize)
  end subroutine classical_particle_finalize

end module classical_particle_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

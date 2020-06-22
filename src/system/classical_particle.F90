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
  use interaction_lorentz_force_oct_m
  use interactions_factory_oct_m
  use ghost_interaction_oct_m
  use io_oct_m
  use iso_c_binding
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
  public ::                  &
    classical_particle_t,    &
    classical_particle_init

   type, extends(system_abst_t) :: classical_particle_t
    FLOAT :: mass
    FLOAT :: pos(1:MAX_DIM)
    FLOAT :: vel(1:MAX_DIM)
    FLOAT :: acc(1:MAX_DIM)
    FLOAT, allocatable :: prev_acc(:,:) !< A storage of the prior times.
    FLOAT :: save_pos(1:MAX_DIM)   !< A storage for the SCF loops
    FLOAT :: save_vel(1:MAX_DIM)   !< A storage for the SCF loops
    FLOAT :: tot_force(1:MAX_DIM)
    FLOAT :: prev_tot_force(1:MAX_DIM) !< Used for the SCF convergence criterium
    FLOAT, allocatable :: prev_pos(:, :) !< Used for extrapolation
    FLOAT, allocatable :: prev_vel(:, :) !< Used for extrapolation
    FLOAT :: hamiltonian_elements(1:MAX_DIM)

    type(c_ptr) :: output_handle
  contains
    procedure :: init_interaction => classical_particle_init_interaction
    procedure :: initial_conditions => classical_particle_initial_conditions
    procedure :: do_td_operation => classical_particle_do_td
    procedure :: iteration_info => classical_particle_iteration_info
    procedure :: output_start => classical_particle_output_start
    procedure :: output_write => classical_particle_output_write
    procedure :: output_finish => classical_particle_output_finish
    procedure :: is_tolerance_reached => classical_particle_is_tolerance_reached
    procedure :: store_current_status => classical_particle_store_current_status
    procedure :: update_quantity => classical_particle_update_quantity
    procedure :: update_exposed_quantity => classical_particle_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => classical_particle_copy_quantities_to_interaction
    procedure :: update_interactions_start => classical_particle_update_interactions_start
    procedure :: update_interactions_finish => classical_particle_update_interactions_finish
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

    call messages_print_stress(stdout, "Classical Particle", namespace=namespace)

    call space_init(this%space, namespace)

    !%Variable ClassicalParticleMass
    !%Type float
    !%Section ClassicalParticles
    !%Description
    !% Mass of classical particle in Kg.
    !%End
    call parse_variable(namespace, 'ClassicalParticleMass', M_ONE, this%mass)
    call messages_print_var_value(stdout, 'ClassicalParticleMass', this%mass)

    this%quantities(POSITION)%required = .true.
    this%quantities(VELOCITY)%required = .true.
    this%quantities(POSITION)%protected = .true.
    this%quantities(VELOCITY)%protected = .true.
    call this%supported_interactions%add(GRAVITY)
    call this%supported_interactions_as_partner%add(GRAVITY)

    call messages_print_stress(stdout, namespace=namespace)

    POP_SUB(classical_particle_init)
  end subroutine classical_particle_init

  ! ---------------------------------------------------------
  subroutine classical_particle_init_interaction(this, interaction)
    class(classical_particle_t), target, intent(inout) :: this
    class(interaction_abst_t),           intent(inout) :: interaction

    PUSH_SUB(classical_particle_init_interaction)

    select type (interaction)
    type is (interaction_gravity_t)
      call interaction%init(this%space%dim, this%quantities, this%mass, this%pos)
    class default
      message(1) = "Trying to initialize an unsupported interaction by classical particles."
      call messages_fatal(1)
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

    !%Variable ClassicalParticleInitialPosition
    !%Type block
    !%Section ClassicalParticles
    !%Description
    !% Initial position of classical particle, in Km.
    !%End
    this%pos = M_ZERO
    if (parse_block(this%namespace, 'ClassicalParticleInitialPosition', blk) == 0) then
      n_rows = parse_block_n(blk)
      if (n_rows > 1) call  messages_input_error(this%namespace, 'ClassicalParticleInitialPosition')

      do idir = 1, this%space%dim
        call parse_block_float(blk, 0, idir - 1, this%pos(idir))
      end do
      call parse_block_end(blk)
    end if
    call messages_print_var_value(stdout, 'ClassicalParticleInitialPosition', this%pos(1:this%space%dim))

    !%Variable ClassicalParticleInitialVelocity
    !%Type block
    !%Section ClassicalParticles
    !%Description
    !% Initial velocity of classical particle in Km/s.
    !%End
    this%vel = M_ZERO
    if (parse_block(this%namespace, 'ClassicalParticleInitialVelocity', blk) == 0) then
      n_rows = parse_block_n(blk)
      if (n_rows > 1) call  messages_input_error(this%namespace, 'ClassicalParticleInitialVelocity')
      do idir = 1, this%space%dim
        call parse_block_float(blk, 0, idir - 1, this%vel(idir))
      end do
      call parse_block_end(blk)
    end if
    call messages_print_var_value(stdout, 'ClassicalParticleInitialVelocity', this%vel(1:this%space%dim))

    POP_SUB(classical_particle_initial_conditions)
  end subroutine classical_particle_initial_conditions

  ! ---------------------------------------------------------
  subroutine classical_particle_do_td(this, operation)
    class(classical_particle_t), intent(inout) :: this
    integer,                     intent(in)    :: operation

    integer :: ii, sdim
    FLOAT, allocatable :: tmp_pos(:, :), tmp_vel(:, :)
    FLOAT :: factor

    PUSH_SUB(classical_particle_do_td)

    sdim = this%space%dim

    select case(operation)
    case (VERLET_START)
      SAFE_ALLOCATE(this%prev_acc(1:this%space%dim, 1))
      this%acc(1:this%space%dim) = this%tot_force(1:this%space%dim) / this%mass

    case (VERLET_FINISH, BEEMAN_FINISH)
      SAFE_DEALLOCATE_A(this%prev_acc)

    case (VERLET_UPDATE_POS)
      this%pos(1:this%space%dim) = this%pos(1:this%space%dim) + this%prop%dt * this%vel(1:this%space%dim) &
                                 + M_HALF * this%prop%dt**2 * this%acc(1:this%space%dim)

      call this%quantities(POSITION)%clock%increment()

    case (VERLET_COMPUTE_ACC)
      do ii = size(this%prev_acc, dim=2) - 1, 1, -1
        this%prev_acc(1:this%space%dim, ii + 1) = this%prev_acc(1:this%space%dim, ii)
      end do
      this%prev_acc(1:this%space%dim, 1) = this%acc(1:this%space%dim)
      this%acc(1:this%space%dim) = this%tot_force(1:this%space%dim) / this%mass

    case (VERLET_COMPUTE_VEL)
      this%vel(1:this%space%dim) = this%vel(1:this%space%dim) &
        + M_HALF * this%prop%dt * (this%prev_acc(1:this%space%dim, 1) + this%acc(1:this%space%dim))

      call this%quantities(VELOCITY)%clock%increment()


    case (BEEMAN_START)
      SAFE_ALLOCATE(this%prev_acc(1:this%space%dim, 2))
      this%acc(1:this%space%dim) = this%tot_force(1:this%space%dim) / this%mass
      this%prev_acc(1:this%space%dim, 1) = this%acc(1:this%space%dim)

    case (BEEMAN_PREDICT_POS)
      this%pos(1:this%space%dim) = this%pos(1:this%space%dim) + this%prop%dt * this%vel(1:this%space%dim) &
                                 + M_ONE/CNST(6.0) * this%prop%dt**2  &
                                 * (M_FOUR*this%acc(1:this%space%dim) - this%prev_acc(1:this%space%dim, 1))

      if (.not. this%prop%predictor_corrector) then
        call this%quantities(POSITION)%clock%increment()
      end if

    case (BEEMAN_PREDICT_VEL)
      this%vel(1:this%space%dim) = this%vel(1:this%space%dim)  &
        + M_ONE/CNST(6.0) * this%prop%dt * ( M_TWO * this%acc(1:this%space%dim) + &
        CNST(5.0) * this%prev_acc(1:this%space%dim, 1) - this%prev_acc(1:this%space%dim, 2))

      call this%quantities(VELOCITY)%clock%increment()

    case (BEEMAN_CORRECT_POS)
      this%pos(1:this%space%dim) = this%save_pos(1:this%space%dim) + this%prop%dt * this%save_vel(1:this%space%dim) &
                                 + M_ONE/CNST(6.0) * this%prop%dt**2  &
                                 * (this%acc(1:this%space%dim) + M_TWO * this%prev_acc(1:this%space%dim, 1))

      ! We set it to the propagation time to avoid double increment
      call this%quantities(POSITION)%clock%set_time(this%prop%clock)

    case (BEEMAN_CORRECT_VEL)
      this%vel(1:this%space%dim) = this%save_vel(1:this%space%dim) &
                                 + M_HALF * this%prop%dt * (this%acc(1:this%space%dim) + this%prev_acc(1:this%space%dim, 1))

      ! We set it to the propagation time to avoid double increment
      call this%quantities(VELOCITY)%clock%set_time(this%prop%clock)

    case (EXPMID_START)
      SAFE_ALLOCATE(this%prev_pos(1:sdim, 1))
      SAFE_ALLOCATE(this%prev_vel(1:sdim, 1))
      this%prev_pos(1:sdim, 1) = this%pos(1:sdim)
      this%prev_vel(1:sdim, 1) = this%vel(1:sdim)

    case (EXPMID_FINISH)
      SAFE_DEALLOCATE_A(this%prev_pos)
      SAFE_DEALLOCATE_A(this%prev_vel)

    case (EXPMID_PREDICT_DT_2)
      this%pos(1:sdim) = CNST(1.5)*this%save_pos(1:sdim) - &
                         CNST(0.5)*this%prev_pos(1:sdim, 1)
      this%vel(1:sdim) = CNST(1.5)*this%save_vel(1:sdim) - &
                         CNST(0.5)*this%prev_vel(1:sdim, 1)
      this%prev_pos(1:sdim, 1) = this%save_pos(1:sdim)
      this%prev_vel(1:sdim, 1) = this%save_vel(1:sdim)
      call this%quantities(POSITION)%clock%increment()
      call this%quantities(VELOCITY)%clock%increment()

    case (UPDATE_HAMILTONIAN)
      this%hamiltonian_elements(1:sdim) = this%tot_force(1:sdim) / (this%mass * this%pos(1:sdim))

    case (EXPMID_PREDICT_DT)
      SAFE_ALLOCATE(tmp_pos(1:sdim, 2))
      SAFE_ALLOCATE(tmp_vel(1:sdim, 2))
      ! apply exponential - at some point this could use the machinery of
      !   exponential_apply (but this would require a lot of boilerplate code
      !   like a Hamiltonian class etc)
      ! save_pos/vel contain the state at t - this is the state we want to
      !   apply the Hamiltonian to
      tmp_pos(1:sdim, 1) = this%save_pos(1:sdim)
      tmp_vel(1:sdim, 1) = this%save_vel(1:sdim)
      this%pos(1:sdim) = this%save_pos(1:sdim)
      this%vel(1:sdim) = this%save_vel(1:sdim)
      ! compute exponential with Taylor expansion
      factor = M_ONE
      do ii = 1, 4
        factor = factor * this%prop%dt / ii
        ! apply hamiltonian
        tmp_pos(1:sdim, 2) = tmp_vel(1:sdim, 1)
        tmp_vel(1:sdim, 2) = this%hamiltonian_elements(1:sdim) * tmp_pos(1:sdim, 1)
        ! swap temporary variables
        tmp_pos(1:sdim, 1) = tmp_pos(1:sdim, 2)
        tmp_vel(1:sdim, 1) = tmp_vel(1:sdim, 2)
        ! accumulate components of Taylor expansion
        this%pos(1:sdim) = this%pos(1:sdim) + factor * tmp_pos(1:sdim, 1)
        this%vel(1:sdim) = this%vel(1:sdim) + factor * tmp_vel(1:sdim, 1)
      end do
      SAFE_DEALLOCATE_A(tmp_pos)
      SAFE_DEALLOCATE_A(tmp_vel)
      call this%quantities(POSITION)%clock%increment()
      call this%quantities(VELOCITY)%clock%increment()

    case (EXPMID_CORRECT_DT_2)
      ! only correct for dt/2 if not converged yet
      if(.not. this%is_tolerance_reached(this%prop%scf_tol)) then
        this%pos(1:sdim) = CNST(0.5)*(this%pos(1:sdim) + &
                                      this%save_pos(1:sdim))
        this%vel(1:sdim) = CNST(0.5)*(this%vel(1:sdim) + &
                                      this%save_vel(1:sdim))
        call this%quantities(POSITION)%clock%increment()
        call this%quantities(VELOCITY)%clock%increment()
      end if

    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

   POP_SUB(classical_particle_do_td)
  end subroutine classical_particle_do_td

  ! ---------------------------------------------------------
  logical function classical_particle_is_tolerance_reached(this, tol) result(converged)
    class(classical_particle_t),   intent(in)    :: this
    FLOAT,                     intent(in)    :: tol

    PUSH_SUB(classical_particle_is_tolerance_reached)

    ! Here we put the criterion that acceleration change is below the tolerance
    converged = .false.
    if ( (sum((this%prev_tot_force(1:this%space%dim) - this%tot_force(1:this%space%dim))**2)/ this%mass) < tol**2) then
      converged = .true.
    end if

    if (debug%info) then
      write(message(1), '(a, e12.6, a, e12.6)') "Debug: -- Change in acceleration  ", &
        sqrt(sum((this%prev_tot_force(1:this%space%dim) - this%tot_force(1:this%space%dim))**2))/this%mass, " and tolerance ", tol
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
  subroutine classical_particle_iteration_info(this)
    class(classical_particle_t), intent(in) :: this

    integer :: idir
    character(len=20) :: fmt

    PUSH_SUB(classical_particle_iteration_info)

    write(message(1),'(2X,A,1X,A)') "Classical particle:", trim(this%namespace%get())

    write(fmt,'("(4X,A,1X,",I2,"e14.6)")') this%space%dim
    write(message(2),fmt) "Coordinates: ", (this%pos(idir), idir = 1, this%space%dim)
    write(message(3),fmt) "Velocity:    ", (this%vel(idir), idir = 1, this%space%dim)
    write(message(4),fmt) "Acceleration:", (this%acc(idir), idir = 1, this%space%dim)
    write(message(5),fmt) "Force:       ", (this%tot_force(idir), idir = 1, this%space%dim)
    write(message(6),'(4x,A,I8.7)')  'Clock tick:      ', this%clock%get_tick()
    write(message(7),'(4x,A,e14.6)') 'Simulation time: ', this%clock%get_sim_time()
    call messages_info(7)

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
    call this%output_write(0)

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
  subroutine classical_particle_output_write(this, iter)
    class(classical_particle_t), intent(inout) :: this
    integer,                 intent(in)    :: iter

    integer :: idir
    character(len=50) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(classical_particle_output_write)

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
      do idir = 1, this%space%dim
        write(aux, '(a2,i3,a1)') 'f(', idir, ')'
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
    call write_iter_double(this%output_handle, this%tot_force, this%space%dim)
    call write_iter_nl(this%output_handle)

    POP_SUB(classical_particle_output_write)
  end subroutine classical_particle_output_write

  ! ---------------------------------------------------------
  subroutine classical_particle_update_quantity(this, iq, requested_time)
    class(classical_particle_t),   intent(inout) :: this
    integer,                   intent(in)    :: iq
    class(clock_t),            intent(in)    :: requested_time

    PUSH_SUB(classical_particle_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case (MASS)
      ! The classical particle has a mass, but it is not necessary to update it, as it does not change with time.
      call this%quantities(iq)%clock%set_time(requested_time)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(classical_particle_update_quantity)
  end subroutine classical_particle_update_quantity

  ! ---------------------------------------------------------
  subroutine classical_particle_update_exposed_quantity(partner, iq, requested_time)
    class(classical_particle_t),   intent(inout) :: partner
    integer,                   intent(in)    :: iq
    class(clock_t),            intent(in)    :: requested_time

    PUSH_SUB(classical_particle_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case (MASS)
      ! The classical particle has a mass, but it does not require any update, as it does not change with time.
      call partner%quantities(iq)%clock%set_time(requested_time)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(classical_particle_update_exposed_quantity)
  end subroutine classical_particle_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine classical_particle_copy_quantities_to_interaction(partner, interaction)
    class(classical_particle_t),          intent(inout) :: partner
    class(interaction_abst_t),        intent(inout) :: interaction

    PUSH_SUB(classical_particle_copy_quantities_to_interaction)

    select type (interaction)
    type is (interaction_gravity_t)
      interaction%partner_mass = partner%mass
      interaction%partner_pos = partner%pos
    type is (interaction_lorentz_force_t)
      ! Nothing to copy
    type is (ghost_interaction_t)
      ! Nothing to copy
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(classical_particle_copy_quantities_to_interaction)
  end subroutine classical_particle_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine classical_particle_update_interactions_start(this)
    class(classical_particle_t), intent(inout) :: this

    PUSH_SUB(classical_particle_update_interactions_start)

    ! Store previous force, as it is used as SCF criterium
    this%prev_tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim)

    POP_SUB(classical_particle_update_interactions_start)
  end subroutine classical_particle_update_interactions_start

  ! ---------------------------------------------------------
  subroutine classical_particle_update_interactions_finish(this)
    class(classical_particle_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter

    PUSH_SUB(classical_particle_update_interactions_finish)

    ! Compute the total force acting on the classical particle
    this%tot_force(1:this%space%dim) = M_ZERO
    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      type is (interaction_lorentz_force_t)
        this%tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim) + interaction%force(1:this%space%dim)
      type is (interaction_gravity_t)
        this%tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim) + interaction%force(1:this%space%dim)
      type is (ghost_interaction_t)
        ! Nothing to do
      class default
        message(1) = "Unknown interaction by the classical particle " + this%namespace%get()
        call messages_fatal(1)
      end select
    end do

    POP_SUB(classical_particle_update_interactions_finish)
  end subroutine classical_particle_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine classical_particle_finalize(this)
    type(classical_particle_t), intent(inout) :: this

    PUSH_SUB(classical_particle_finalize)

    call system_abst_end(this)

    POP_SUB(classical_particle_finalize)
  end subroutine classical_particle_finalize

end module classical_particle_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

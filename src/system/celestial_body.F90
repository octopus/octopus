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
  use clock_oct_m
  use global_oct_m
  use interaction_abst_oct_m
  use interaction_gravity_oct_m
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
  public ::               &
    celestial_body_t

  type, extends(system_abst_t) :: celestial_body_t
    private

    FLOAT :: mass
    FLOAT :: pos(1:MAX_DIM)
    FLOAT :: vel(1:MAX_DIM)
    FLOAT :: acc(1:MAX_DIM)
    FLOAT, allocatable :: prev_acc(:,:) !< A storage of the prior times.
    FLOAT :: save_pos(1:MAX_DIM)   !< A storage for the SCF loops
    FLOAT :: save_vel(1:MAX_DIM)   !< A storage for the SCF loops
    FLOAT :: tot_force(1:MAX_DIM)
    FLOAT :: prev_tot_force(1:MAX_DIM) !< Used for the SCF convergence criterium

    type(c_ptr) :: output_handle
  contains
    procedure :: add_interaction_partner => celestial_body_add_interaction_partner
    procedure :: has_interaction => celestial_body_has_interaction
    procedure :: initial_conditions => celestial_body_initial_conditions
    procedure :: do_td_operation => celestial_body_do_td
    procedure :: iteration_info => celestial_body_iteration_info
    procedure :: output_start => celestial_body_output_start
    procedure :: output_write => celestial_body_output_write
    procedure :: output_finish => celestial_body_output_finish
    procedure :: is_tolerance_reached => celestial_body_is_tolerance_reached
    procedure :: store_current_status => celestial_body_store_current_status
    procedure :: update_quantity => celestial_body_update_quantity
    procedure :: update_exposed_quantity => celestial_body_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => celestial_body_copy_quantities_to_interaction
    procedure :: update_interactions_start => celestial_body_update_interactions_start
    procedure :: update_interactions_finish => celestial_body_update_interactions_finish
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

    sys%quantities(POSITION)%required = .true.
    sys%quantities(VELOCITY)%required = .true.
    sys%quantities(POSITION)%protected = .true.
    sys%quantities(VELOCITY)%protected = .true.

    call messages_print_stress(stdout, namespace=namespace)

    POP_SUB(celestial_body_init)
  end function celestial_body_init

  ! ---------------------------------------------------------
  subroutine celestial_body_add_interaction_partner(this, partner)
    class(celestial_body_t), target, intent(inout) :: this
    class(system_abst_t),            intent(inout) :: partner

    class(interaction_gravity_t), pointer :: gravity
    type(interaction_gravity_t) :: gravity_t

    PUSH_SUB(celestial_body_add_interaction_partner)

    if (partner%has_interaction(gravity_t)) then
      gravity => interaction_gravity_t(this%space%dim, partner)
      this%quantities(POSITION)%required = .true.
      this%quantities(MASS)%required = .true.
      gravity%system_mass => this%mass
      gravity%system_pos  => this%pos
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
  subroutine celestial_body_initial_conditions(this, from_scratch)
    class(celestial_body_t), intent(inout) :: this
    logical,                 intent(in)    :: from_scratch

    integer :: n_rows, idir
    type(block_t) :: blk

    PUSH_SUB(celestial_body_initial_conditions)

    if (.not. from_scratch) then
      message(1) = "Celestial mechanics is currently only allowed from scratch"
      call messages_fatal(1, namespace=this%namespace)
    end if

    !%Variable CelestialBodyInitialPosition
    !%Type block
    !%Section CelestialDynamics
    !%Description
    !% Initial position of celestial body, in Km.
    !%End
    this%pos = M_ZERO
    if (parse_block(this%namespace, 'CelestialBodyInitialPosition', blk) == 0) then
      n_rows = parse_block_n(blk)
      if (n_rows > 1) call  messages_input_error(this%namespace, 'CelestialBodyInitialPosition')

      do idir = 1, this%space%dim
        call parse_block_float(blk, 0, idir - 1, this%pos(idir))
      end do
      call parse_block_end(blk)
    end if
    call messages_print_var_value(stdout, 'CelestialBodyInitialPosition', this%pos(1:this%space%dim))

    !%Variable CelestialBodyInitialVelocity
    !%Type block
    !%Section CelestialDynamics
    !%Description
    !% Initial velocity of celestial body in Km/s.
    !%End
    this%vel = M_ZERO
    if (parse_block(this%namespace, 'CelestialBodyInitialVelocity', blk) == 0) then
      n_rows = parse_block_n(blk)
      if (n_rows > 1) call  messages_input_error(this%namespace, 'CelestialBodyInitialVelocity')
      do idir = 1, this%space%dim
        call parse_block_float(blk, 0, idir - 1, this%vel(idir))
      end do
      call parse_block_end(blk)
    end if
    call messages_print_var_value(stdout, 'CelestialBodyInitialVelocity', this%vel(1:this%space%dim))

    POP_SUB(celestial_body_initial_conditions)
  end subroutine celestial_body_initial_conditions

  ! ---------------------------------------------------------
  subroutine celestial_body_do_td(this, operation)
    class(celestial_body_t), intent(inout) :: this
    integer,                 intent(in)    :: operation

    integer :: ii

    PUSH_SUB(celestial_body_do_td)

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

    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

   POP_SUB(celestial_body_do_td)
  end subroutine celestial_body_do_td

  ! ---------------------------------------------------------
  logical function celestial_body_is_tolerance_reached(this, tol) result(converged)
    class(celestial_body_t),   intent(in)    :: this
    FLOAT,                     intent(in)    :: tol

    PUSH_SUB(celestial_body_is_tolerance_reached)

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

    POP_SUB(celestial_body_is_tolerance_reached)
   end function celestial_body_is_tolerance_reached

   ! ---------------------------------------------------------
   subroutine celestial_body_store_current_status(this)
     class(celestial_body_t),   intent(inout)    :: this

     PUSH_SUB(celestial_body_store_current_status) 

     this%save_pos(1:this%space%dim) = this%pos(1:this%space%dim)
     this%save_vel(1:this%space%dim) = this%vel(1:this%space%dim)

     POP_SUB(celestial_body_store_current_status)
   end subroutine celestial_body_store_current_status

  ! ---------------------------------------------------------
  subroutine celestial_body_iteration_info(this)
    class(celestial_body_t), intent(in) :: this

    integer :: idir
    character(len=20) :: fmt

    PUSH_SUB(celestial_body_iteration_info)

    write(message(1),'(2X,A,1X,A)') "Celestial body:", trim(this%namespace%get())

    write(fmt,'("(4X,A,1X,",I2,"e14.6)")') this%space%dim
    write(message(2),fmt) "Coordinates: ", (this%pos(idir), idir = 1, this%space%dim)
    write(message(3),fmt) "Velocity:    ", (this%vel(idir), idir = 1, this%space%dim)
    write(message(4),fmt) "Acceleration:", (this%acc(idir), idir = 1, this%space%dim)
    write(message(5),fmt) "Force:       ", (this%tot_force(idir), idir = 1, this%space%dim)
    write(message(6),'(4x,A,I8.7)')  'Clock tick:      ', this%clock%get_tick()
    write(message(7),'(4x,A,e14.6)') 'Simulation time: ', this%clock%get_sim_time()
    call messages_info(7)

    POP_SUB(celestial_body_iteration_info)
  end subroutine celestial_body_iteration_info

  ! ---------------------------------------------------------
  subroutine celestial_body_output_start(this)
    class(celestial_body_t), intent(inout) :: this

    PUSH_SUB(celestial_body_output_start)

    ! Create output handle
    call io_mkdir('td.general', this%namespace)
    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_init(this%output_handle, 0, this%prop%dt, trim(io_workpath("td.general/coordinates", this%namespace)))
    end if

    ! Output info for first iteration
    call this%output_write(0)

    POP_SUB(celestial_body_output_start)
  end subroutine celestial_body_output_start

  ! ---------------------------------------------------------
  subroutine celestial_body_output_finish(this)
    class(celestial_body_t), intent(inout) :: this

    PUSH_SUB(celestial_body_output_finish)

    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_end(this%output_handle)
    end if

    POP_SUB(celestial_body_output_finish)
  end subroutine celestial_body_output_finish

  ! ---------------------------------------------------------
  subroutine celestial_body_output_write(this, iter)
    class(celestial_body_t), intent(inout) :: this
    integer,                 intent(in)    :: iter

    integer :: idir
    character(len=50) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(celestial_body_output_write)

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
    
    POP_SUB(celestial_body_output_write)
  end subroutine celestial_body_output_write

  ! ---------------------------------------------------------
  subroutine celestial_body_update_quantity(this, iq, requested_time)
    class(celestial_body_t),   intent(inout) :: this
    integer,                   intent(in)    :: iq
    class(clock_t),            intent(in)    :: requested_time

    PUSH_SUB(celestial_body_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case (MASS)
      ! The celestial body has a mass, but it is not necessary to update it, as it does not change with time.
      call this%quantities(iq)%clock%set_time(requested_time)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(celestial_body_update_quantity)
  end subroutine celestial_body_update_quantity

 ! ---------------------------------------------------------
 subroutine celestial_body_update_exposed_quantity(this, iq, requested_time)
    class(celestial_body_t),   intent(inout) :: this
    integer,                   intent(in)    :: iq
    class(clock_t),            intent(in)    :: requested_time

    PUSH_SUB(celestial_body_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case (MASS)
      ! The celestial body has a mass, but it does not require any update, as it does not change with time.
      call this%quantities(iq)%clock%set_time(this%clock)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(celestial_body_update_exposed_quantity)
  end subroutine celestial_body_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine celestial_body_copy_quantities_to_interaction(this, interaction)
    class(celestial_body_t),          intent(inout) :: this
    class(interaction_abst_t),        intent(inout) :: interaction

    PUSH_SUB(celestial_body_copy_quantities_to_interaction)

    select type (interaction)
    type is (interaction_gravity_t)
      interaction%partner_mass = this%mass
      interaction%partner_pos = this%pos
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(celestial_body_copy_quantities_to_interaction)
  end subroutine celestial_body_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine celestial_body_update_interactions_start(this)
    class(celestial_body_t), intent(inout) :: this

    PUSH_SUB(celestial_body_update_interactions_start)

    ! Store previous force, as it is used as SCF criterium
    this%prev_tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim)

    POP_SUB(celestial_body_update_interactions_start)
  end subroutine celestial_body_update_interactions_start

  ! ---------------------------------------------------------
  subroutine celestial_body_update_interactions_finish(this)
    class(celestial_body_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter

    PUSH_SUB(celestial_body_update_interactions_finish)

    ! Compute the total force acting on the celestial body
    this%tot_force(1:this%space%dim) = M_ZERO
    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next_interaction())
      type is (interaction_gravity_t)
        this%tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim) + interaction%force(1:this%space%dim)
      class default
        message(1) = "Unknown interaction by the celestial body " + this%namespace%get()
        call messages_fatal(1)
      end select
    end do

    POP_SUB(celestial_body_update_interactions_finish)
  end subroutine celestial_body_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine celestial_body_finalize(this)
    type(celestial_body_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter
    class(interaction_abst_t), pointer :: interaction

    PUSH_SUB(celestial_body_finalize)

    deallocate(this%prop)

    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next_interaction()
      SAFE_DEALLOCATE_P(interaction)
    end do

    POP_SUB(celestial_body_finalize)
  end subroutine celestial_body_finalize

end module celestial_body_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

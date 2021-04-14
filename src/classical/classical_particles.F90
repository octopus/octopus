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

module classical_particles_oct_m
  use algorithm_oct_m
  use clock_oct_m
  use force_interaction_oct_m
  use global_oct_m
  use interaction_oct_m
  use gravity_oct_m
  use lorentz_force_oct_m
  use interactions_factory_oct_m
  use io_oct_m
  use iso_c_binding
  use messages_oct_m
  use mpi_oct_m
  use multisystem_debug_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagator_beeman_oct_m
  use propagator_exp_mid_oct_m
  use propagator_oct_m
  use propagator_verlet_oct_m
  use quantity_oct_m
  use space_oct_m
  use system_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::                                      &
    classical_particles_t,                       &
    classical_particles_init,                    &
    classical_particles_end,                     &
    classical_particles_init_interaction,        &
    classical_particles_update_quantity,         &
    classical_particles_update_exposed_quantity

  type, extends(system_t), abstract :: classical_particles_t
    integer :: np
    FLOAT, allocatable :: mass(:)
    FLOAT, allocatable :: pos(:,:)
    FLOAT, allocatable :: vel(:,:)
    FLOAT, allocatable :: acc(:,:)
    FLOAT, allocatable :: prev_acc(:,:,:) !< A storage of the prior times.
    FLOAT, allocatable :: save_pos(:,:)   !< A storage for the SCF loops
    FLOAT, allocatable :: save_vel(:,:)   !< A storage for the SCF loops
    FLOAT, allocatable :: tot_force(:,:)  !< Total force acting on each particle
    FLOAT, allocatable :: prev_tot_force(:,:) !< Used for the SCF convergence criterium
    FLOAT, allocatable :: prev_pos(:,:,:) !< Used for extrapolation
    FLOAT, allocatable :: prev_vel(:,:,:) !< Used for extrapolation
    FLOAT, allocatable :: hamiltonian_elements(:,:)
  contains
    procedure :: do_td_operation => classical_particles_do_td
    procedure :: is_tolerance_reached => classical_particles_is_tolerance_reached
    procedure :: copy_quantities_to_interaction => classical_particles_copy_quantities_to_interaction
    procedure :: update_interactions_start => classical_particles_update_interactions_start
    procedure :: update_interactions_finish => classical_particles_update_interactions_finish
  end type classical_particles_t

contains

  ! ---------------------------------------------------------
  !> The init routine is a module level procedure
  !! This has the advantage that different classes can have different
  !! signatures for the initialization routines because they are not
  !! type-bound and thus also not inherited.
  ! ---------------------------------------------------------
  subroutine classical_particles_init(this, namespace, np)
    class(classical_particles_t), intent(inout) :: this
    type(namespace_t),            intent(in)    :: namespace
    integer,                      intent(in)    :: np !< Number of particles

    PUSH_SUB(classical_particles_init)

    this%namespace = namespace

    call space_init(this%space, namespace)

    this%np = np
    SAFE_ALLOCATE(this%mass(1:np))
    SAFE_ALLOCATE(this%pos(1:this%space%dim, 1:np))
    SAFE_ALLOCATE(this%vel(1:this%space%dim, 1:np))
    SAFE_ALLOCATE(this%acc(1:this%space%dim, 1:np))
    SAFE_ALLOCATE(this%prev_tot_force(1:this%space%dim, 1:np))
    SAFE_ALLOCATE(this%tot_force(1:this%space%dim, 1:np))

    this%quantities(POSITION)%required = .true.
    this%quantities(VELOCITY)%required = .true.
    this%quantities(POSITION)%protected = .true.
    this%quantities(VELOCITY)%protected = .true.

    POP_SUB(classical_particles_init)
  end subroutine classical_particles_init

  ! ---------------------------------------------------------
  subroutine classical_particles_init_interaction(this, interaction)
    class(classical_particles_t), target, intent(inout) :: this
    class(interaction_t),                 intent(inout) :: interaction

    PUSH_SUB(classical_particles_init_interaction)

    select type (interaction)
    class default
      message(1) = "Trying to initialize an unsupported interaction by classical particles."
      call messages_fatal(1)
    end select

    POP_SUB(classical_particles_init_interaction)
  end subroutine classical_particles_init_interaction

  ! ---------------------------------------------------------
  subroutine classical_particles_do_td(this, operation)
    class(classical_particles_t),    intent(inout) :: this
    class(algorithmic_operation_t),  intent(in)    :: operation

    integer :: ii, sdim, ip
    FLOAT, allocatable :: tmp_pos(:,:,:), tmp_vel(:,:,:)
    FLOAT :: factor

    PUSH_SUB(classical_particles_do_td)

    sdim = this%space%dim

    select case (operation%id)
    case (SKIP)
      ! Do nothing
    case (STORE_CURRENT_STATUS)
      this%save_pos(1:sdim, 1:this%np) = this%pos(1:sdim, 1:this%np)
      this%save_vel(1:sdim, 1:this%np) = this%vel(1:sdim, 1:this%np)

    case (VERLET_START)
      SAFE_ALLOCATE(this%prev_acc(1:sdim, 1:this%np, 1))
      do ip = 1, this%np
        this%acc(1:sdim, ip) = this%tot_force(1:sdim, ip) / this%mass(ip)
      end do

    case (VERLET_FINISH)
      SAFE_DEALLOCATE_A(this%prev_acc)

    case (BEEMAN_FINISH)
      SAFE_DEALLOCATE_A(this%prev_acc)
      SAFE_DEALLOCATE_A(this%save_pos)
      SAFE_DEALLOCATE_A(this%save_vel)

    case (VERLET_UPDATE_POS)
      this%pos(1:sdim, 1:this%np) = this%pos(1:sdim, 1:this%np) + &
        this%prop%dt * this%vel(1:sdim, 1:this%np) &
        + M_HALF * this%prop%dt**2 * this%acc(1:sdim, 1:this%np)

      this%quantities(POSITION)%clock = this%quantities(POSITION)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t(clock_name="quantity", & 
                                          clock_detail=QUANTITY_LABEL(POSITION), & 
                                          clock = this%quantities(POSITION)%clock, action="tick") )


    case (VERLET_COMPUTE_ACC, BEEMAN_COMPUTE_ACC)
      do ii = size(this%prev_acc, dim=3) - 1, 1, -1
        this%prev_acc(1:sdim, 1:this%np, ii + 1) = this%prev_acc(1:sdim, 1:this%np, ii)
      end do
      do ip = 1, this%np
        this%prev_acc(1:sdim, ip, 1) = this%acc(1:sdim, ip)
        this%acc(1:sdim, ip) = this%tot_force(1:sdim, ip) / this%mass(ip)
      end do

    case (VERLET_COMPUTE_VEL)
      this%vel(1:sdim, 1:this%np) = this%vel(1:sdim, 1:this%np) &
        + M_HALF * this%prop%dt * (this%prev_acc(1:sdim, 1:this%np, 1) + this%acc(1:sdim, 1:this%np))

      this%quantities(VELOCITY)%clock = this%quantities(VELOCITY)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t(clock_name="quantity", &
                                          clock_detail = QUANTITY_LABEL(VELOCITY), & 
                                          clock = this%quantities(VELOCITY)%clock, action="tick") )

    case (BEEMAN_START)
      if (this%prop%predictor_corrector) then
        SAFE_ALLOCATE(this%save_pos(1:this%space%dim, 1:this%np))
        SAFE_ALLOCATE(this%save_vel(1:this%space%dim, 1:this%np))
      end if
      SAFE_ALLOCATE(this%prev_acc(1:sdim, 1:this%np, 1:2))
      do ip = 1, this%np
        this%acc(1:sdim, ip) = this%tot_force(1:sdim, ip) / this%mass(ip)
        this%prev_acc(1:sdim, ip, 1) = this%acc(1:sdim, ip)
      end do

    case (BEEMAN_PREDICT_POS)
      this%pos(1:sdim, 1:this%np) = this%pos(1:sdim, 1:this%np) + this%prop%dt * this%vel(1:sdim, 1:this%np) + &
        M_ONE/CNST(6.0) * this%prop%dt**2 * (M_FOUR*this%acc(1:sdim, 1:this%np) - this%prev_acc(1:sdim, 1:this%np, 1))

      if (.not. this%prop%predictor_corrector) then
        this%quantities(POSITION)%clock = this%quantities(POSITION)%clock + CLOCK_TICK
        call multisystem_debug_write_marker(this%namespace, event_clock_update_t(clock_name="quantity", &
                                                              clock_detail = QUANTITY_LABEL(POSITION), & 
                                                              clock = this%quantities(POSITION)%clock, action="tick") )
end if

    case (BEEMAN_PREDICT_VEL)
      this%vel(1:sdim, 1:this%np) = this%vel(1:sdim, 1:this%np)  &
        + M_ONE/CNST(6.0) * this%prop%dt * ( M_TWO * this%acc(1:sdim, 1:this%np) + &
        CNST(5.0) * this%prev_acc(1:sdim, 1:this%np, 1) - this%prev_acc(1:sdim, 1:this%np, 2))

      this%quantities(VELOCITY)%clock = this%quantities(VELOCITY)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t(clock_name="quantity", &
                                          clock_detail = QUANTITY_LABEL(VELOCITY), & 
                                          clock = this%quantities(VELOCITY)%clock, action="tick") )

    case (BEEMAN_CORRECT_POS)
      this%pos(1:sdim, 1:this%np) = this%save_pos(1:sdim, 1:this%np) + this%prop%dt * this%save_vel(1:sdim, 1:this%np) &
        + M_ONE/CNST(6.0) * this%prop%dt**2 * (this%acc(1:sdim, 1:this%np) + M_TWO * this%prev_acc(1:sdim, 1:this%np, 1))

      ! We set it to the propagation time to avoid double increment
      call this%quantities(POSITION)%clock%set_time(this%prop%clock)

    case (BEEMAN_CORRECT_VEL)
      this%vel(1:sdim, 1:this%np) = this%save_vel(1:sdim, 1:this%np) &
        + M_HALF * this%prop%dt * (this%acc(1:sdim, 1:this%np) + this%prev_acc(1:sdim, 1:this%np, 1))

      ! We set it to the propagation time to avoid double increment
      call this%quantities(VELOCITY)%clock%set_time(this%prop%clock)
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t(clock_name="quantity", &
                                          clock_detail = QUANTITY_LABEL(VELOCITY), & 
                                          clock = this%quantities(VELOCITY)%clock, action="set") )

    case (EXPMID_START)
      SAFE_ALLOCATE(this%save_pos(1:this%space%dim, 1:this%np))
      SAFE_ALLOCATE(this%save_vel(1:this%space%dim, 1:this%np))
      SAFE_ALLOCATE(this%hamiltonian_elements(1:sdim, 1:this%np))
      SAFE_ALLOCATE(this%prev_pos(1:sdim, 1:this%np, 1))
      SAFE_ALLOCATE(this%prev_vel(1:sdim, 1:this%np, 1))
      this%prev_pos(1:sdim, 1:this%np, 1) = this%pos(1:sdim, 1:this%np)
      this%prev_vel(1:sdim, 1:this%np, 1) = this%vel(1:sdim, 1:this%np)

    case (EXPMID_FINISH)
      SAFE_DEALLOCATE_A(this%save_pos)
      SAFE_DEALLOCATE_A(this%save_vel)
      SAFE_DEALLOCATE_A(this%hamiltonian_elements)
      SAFE_DEALLOCATE_A(this%prev_pos)
      SAFE_DEALLOCATE_A(this%prev_vel)

    case (EXPMID_PREDICT_DT_2)
      this%pos(1:sdim, 1:this%np) = CNST(1.5)*this%save_pos(1:sdim, 1:this%np) - CNST(0.5)*this%prev_pos(1:sdim, 1:this%np, 1)
      this%vel(1:sdim, 1:this%np) = CNST(1.5)*this%save_vel(1:sdim, 1:this%np) - CNST(0.5)*this%prev_vel(1:sdim, 1:this%np, 1)
      this%prev_pos(1:sdim, 1:this%np, 1) = this%save_pos(1:sdim, 1:this%np)
      this%prev_vel(1:sdim, 1:this%np, 1) = this%save_vel(1:sdim, 1:this%np)
      this%quantities(POSITION)%clock = this%quantities(POSITION)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t(clock_name="quantity", &
                                          clock_detail = QUANTITY_LABEL(POSITION), & 
                                          clock = this%quantities(POSITION)%clock, action="tick") )
      this%quantities(VELOCITY)%clock = this%quantities(VELOCITY)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t(clock_name="quantity", &
                                          clock_detail = QUANTITY_LABEL(VELOCITY), & 
                                          clock = this%quantities(VELOCITY)%clock, action="tick") )

    case (UPDATE_HAMILTONIAN)
      do ip = 1, this%np
        this%hamiltonian_elements(1:sdim, ip) = this%tot_force(1:sdim, ip) / (this%mass(ip) * this%pos(1:sdim, ip))
      end do

    case (EXPMID_PREDICT_DT)
      SAFE_ALLOCATE(tmp_pos(1:sdim, 1:this%np, 2))
      SAFE_ALLOCATE(tmp_vel(1:sdim, 1:this%np, 2))
      ! apply exponential - at some point this could use the machinery of
      !   exponential_apply (but this would require a lot of boilerplate code
      !   like a Hamiltonian class etc)
      ! save_pos/vel contain the state at t - this is the state we want to
      !   apply the Hamiltonian to
      tmp_pos(1:sdim, 1:this%np, 1) = this%save_pos(1:sdim, 1:this%np)
      tmp_vel(1:sdim, 1:this%np, 1) = this%save_vel(1:sdim, 1:this%np)
      this%pos(1:sdim, 1:this%np) = this%save_pos(1:sdim, 1:this%np)
      this%vel(1:sdim, 1:this%np) = this%save_vel(1:sdim, 1:this%np)
      ! compute exponential with Taylor expansion
      factor = M_ONE
      do ii = 1, 4
        factor = factor * this%prop%dt / ii
        do ip = 1, this%np          
          ! apply hamiltonian
          tmp_pos(1:sdim, ip, 2) = tmp_vel(1:sdim, ip, 1)
          tmp_vel(1:sdim, ip, 2) = this%hamiltonian_elements(1:sdim, ip) * tmp_pos(1:sdim, ip, 1)
          ! swap temporary variables
          tmp_pos(1:sdim, ip, 1) = tmp_pos(1:sdim, ip, 2)
          tmp_vel(1:sdim, ip, 1) = tmp_vel(1:sdim, ip, 2)
          ! accumulate components of Taylor expansion
          this%pos(1:sdim, ip) = this%pos(1:sdim, ip) + factor * tmp_pos(1:sdim, ip, 1)
          this%vel(1:sdim, ip) = this%vel(1:sdim, ip) + factor * tmp_vel(1:sdim, ip, 1)
        end do
      end do
      SAFE_DEALLOCATE_A(tmp_pos)
      SAFE_DEALLOCATE_A(tmp_vel)
      this%quantities(POSITION)%clock = this%quantities(POSITION)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t(clock_name="quantity", &
                                          clock_detail = QUANTITY_LABEL(POSITION), & 
                                          clock = this%quantities(POSITION)%clock, action="tick") )
      this%quantities(VELOCITY)%clock = this%quantities(VELOCITY)%clock + CLOCK_TICK
      call multisystem_debug_write_marker(this%namespace, event_clock_update_t(clock_name="quantity", &
                                          clock_detail = QUANTITY_LABEL(VELOCITY), & 
                                          clock = this%quantities(VELOCITY)%clock, action="tick") )

    case (EXPMID_CORRECT_DT_2)
      ! only correct for dt/2 if not converged yet
      if (.not. this%is_tolerance_reached(this%prop%scf_tol)) then
        this%pos(1:sdim, 1:this%np) = CNST(0.5)*(this%pos(1:sdim, 1:this%np) + this%save_pos(1:sdim, 1:this%np))
        this%vel(1:sdim, 1:this%np) = CNST(0.5)*(this%vel(1:sdim, 1:this%np) + this%save_vel(1:sdim, 1:this%np))
        this%quantities(POSITION)%clock = this%quantities(POSITION)%clock + CLOCK_TICK
        call multisystem_debug_write_marker(this%namespace, event_clock_update_t(clock_name="quantity", &
                                            clock_detail = QUANTITY_LABEL(POSITION), & 
                                            clock = this%quantities(POSITION)%clock, action="tick") )
        this%quantities(VELOCITY)%clock = this%quantities(VELOCITY)%clock + CLOCK_TICK
        call multisystem_debug_write_marker(this%namespace, event_clock_update_t(clock_name="quantity", &
                                            clock_detail = QUANTITY_LABEL(VELOCITY), & 
                                            clock = this%quantities(VELOCITY)%clock, action="tick") )
end if

    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

   POP_SUB(classical_particles_do_td)
  end subroutine classical_particles_do_td

  ! ---------------------------------------------------------
  logical function classical_particles_is_tolerance_reached(this, tol) result(converged)
    class(classical_particles_t), intent(in)    :: this
    FLOAT,                        intent(in)    :: tol

    integer :: ip
    FLOAT :: change, max_change

    PUSH_SUB(classical_particles_is_tolerance_reached)

    ! Here we put the criterion that maximum acceleration change is below the tolerance
    max_change = M_ZERO
    do ip = 1, this%np
      change = sum((this%prev_tot_force(1:this%space%dim, ip) - this%tot_force(1:this%space%dim, ip))**2)/this%mass(ip)
      if (change > max_change) then
        max_change = change
      end if
    end do
    converged = max_change < tol**2

    if (debug%info) then
      write(message(1), '(a, e12.6, a, e12.6)') "Debug: -- Maximum change in acceleration  ", &
        sqrt(max_change), " and tolerance ", tol
      call messages_info(1)
    end if

    POP_SUB(classical_particles_is_tolerance_reached)
  end function classical_particles_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine classical_particles_update_quantity(this, iq)
    class(classical_particles_t), intent(inout) :: this
    integer,                     intent(in)    :: iq

    PUSH_SUB(classical_particles_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case (MASS)
      ! The classical particles have a mass, but it is not necessary to update them, as they do not change with time.
      ! We still need to set their clock, so we set it to be in sync with the particles position.
      call this%quantities(iq)%clock%set_time(this%quantities(POSITION)%clock)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(classical_particles_update_quantity)
  end subroutine classical_particles_update_quantity

  ! ---------------------------------------------------------
  subroutine classical_particles_update_exposed_quantity(partner, iq)
    class(classical_particles_t), intent(inout) :: partner
    integer,                     intent(in)    :: iq

    PUSH_SUB(classical_particles_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case (MASS)
      ! The classical particles have a mass, but thye do not require any update, as they do not change with time.
      ! We still need to set their clock, so we set it to be in sync with the particles position.
      call partner%quantities(iq)%clock%set_time(partner%quantities(POSITION)%clock)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(classical_particles_update_exposed_quantity)
  end subroutine classical_particles_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine classical_particles_copy_quantities_to_interaction(partner, interaction)
    class(classical_particles_t),         intent(inout) :: partner
    class(interaction_t),                 intent(inout) :: interaction

    PUSH_SUB(classical_particles_copy_quantities_to_interaction)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(classical_particles_copy_quantities_to_interaction)
  end subroutine classical_particles_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine classical_particles_update_interactions_start(this)
    class(classical_particles_t), intent(inout) :: this

    PUSH_SUB(classical_particles_update_interactions_start)

    ! Store previous force, as it is used as SCF criterium
    this%prev_tot_force(1:this%space%dim, 1:this%np) = this%tot_force(1:this%space%dim, 1:this%np)

    POP_SUB(classical_particles_update_interactions_start)
  end subroutine classical_particles_update_interactions_start

  ! ---------------------------------------------------------
  subroutine classical_particles_update_interactions_finish(this)
    class(classical_particles_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter

    PUSH_SUB(classical_particles_update_interactions_finish)

    ! Compute the total force acting on the classical particles
    this%tot_force(1:this%space%dim, 1:this%np) = M_ZERO
    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class is (force_interaction_t)
        this%tot_force(1:this%space%dim, 1:this%np) = this%tot_force(1:this%space%dim, 1:this%np) + &
          interaction%force(1:this%space%dim, 1:this%np)
      end select
    end do

    POP_SUB(classical_particles_update_interactions_finish)
  end subroutine classical_particles_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine classical_particles_end(this)
    class(classical_particles_t), intent(inout) :: this

    PUSH_SUB(classical_particles_end)

    SAFE_DEALLOCATE_A(this%mass)
    SAFE_DEALLOCATE_A(this%pos)
    SAFE_DEALLOCATE_A(this%vel)
    SAFE_DEALLOCATE_A(this%acc)
    SAFE_DEALLOCATE_A(this%prev_tot_force)
    SAFE_DEALLOCATE_A(this%tot_force)
    SAFE_DEALLOCATE_A(this%save_pos)
    SAFE_DEALLOCATE_A(this%save_vel)

    call system_end(this)

    POP_SUB(classical_particles_end)
  end subroutine classical_particles_end

end module classical_particles_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

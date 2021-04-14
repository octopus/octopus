!! Copyright (C) 2020 H. Appel, S. Ohlmann, M. Oliveira, N. Tancogne-Dejean
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

module charged_particle_oct_m
  use algorithm_oct_m
  use classical_particle_oct_m
  use clock_oct_m
  use global_oct_m
  use interaction_oct_m
  use coulomb_force_oct_m
  use lorentz_force_oct_m
  use interactions_factory_oct_m
  use io_oct_m
  use iso_c_binding
  use messages_oct_m
  use multisystem_debug_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use space_oct_m
  use system_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::               &
    charged_particle_t,   &
    charged_particle_init

  type, extends(classical_particle_t) :: charged_particle_t
    private

    FLOAT :: charge(1)

  contains
    procedure :: init_interaction => charged_particle_init_interaction
    procedure :: initial_conditions => charged_particle_initial_conditions
    procedure :: do_td_operation => charged_particle_do_td
    procedure :: iteration_info => charged_particle_iteration_info
    procedure :: is_tolerance_reached => charged_particle_is_tolerance_reached
    procedure :: update_quantity => charged_particle_update_quantity
    procedure :: update_exposed_quantity => charged_particle_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => charged_particle_copy_quantities_to_interaction
  end type charged_particle_t

  interface charged_particle_t
    procedure charged_particle_constructor
  end interface charged_particle_t

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function charged_particle_constructor(namespace) result(sys)
    class(charged_particle_t), pointer  :: sys
    type(namespace_t),       intent(in) :: namespace

    PUSH_SUB(charged_particle_constructor)

    SAFE_ALLOCATE(sys)

    call charged_particle_init(sys, namespace)

    POP_SUB(charged_particle_constructor)
  end function charged_particle_constructor

  ! ---------------------------------------------------------
  !> The init routine is a module level procedure
  !! This has the advantage that different classes can have different
  !! signatures for the initialization routines because they are not
  !! type-bound and thus also not inherited.
  subroutine charged_particle_init(this, namespace)
    class(charged_particle_t), intent(inout) :: this
    type(namespace_t),         intent(in)    :: namespace

    PUSH_SUB(charged_particle_init)

    call classical_particle_init(this%classical_particle_t, namespace)

    !%Variable ParticleCharge
    !%Type float
    !%Section ClassicalParticles
    !%Description
    !% Charge of classical particle
    !%End
    call parse_variable(namespace, 'ParticleCharge', M_ONE, this%charge(1))
    call messages_print_var_value(stdout, 'ParticleCharge', this%charge(1))

    call this%supported_interactions%add(LORENTZ_FORCE)
    call this%supported_interactions%add(COULOMB_FORCE)
    call this%supported_interactions_as_partner%add(COULOMB_FORCE)

    POP_SUB(charged_particle_init)
  end subroutine charged_particle_init

  ! ---------------------------------------------------------
  subroutine charged_particle_init_interaction(this, interaction)
    class(charged_particle_t), target, intent(inout) :: this
    class(interaction_t),              intent(inout) :: interaction

    PUSH_SUB(charged_particle_init_interaction)

    select type (interaction)
    type is (coulomb_force_t)
      call interaction%init(this%space%dim, 1, this%quantities, this%charge, this%pos)
    type is (lorentz_force_t)
      call interaction%init(this%space%dim, 1, this%quantities, this%charge, this%pos, this%vel)
    class default
      call this%classical_particle_t%init_interaction(interaction)
    end select

    POP_SUB(charged_particle_init_interaction)
  end subroutine charged_particle_init_interaction

  ! ---------------------------------------------------------
  subroutine charged_particle_initial_conditions(this, from_scratch)
    class(charged_particle_t), intent(inout) :: this
    logical,                   intent(in)    :: from_scratch

    PUSH_SUB(charged_particle_initial_conditions)

    call this%classical_particle_t%initial_conditions(from_scratch)

    POP_SUB(charged_particle_initial_conditions)
  end subroutine charged_particle_initial_conditions

  ! ---------------------------------------------------------
  subroutine charged_particle_do_td(this, operation)
    class(charged_particle_t),      intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    PUSH_SUB(charged_particle_do_td)

    call this%classical_particle_t%do_td_operation(operation)

    POP_SUB(charged_particle_do_td)
  end subroutine charged_particle_do_td

  ! ---------------------------------------------------------
  logical function charged_particle_is_tolerance_reached(this, tol) result(converged)
    class(charged_particle_t), intent(in) :: this
    FLOAT,                     intent(in) :: tol

    PUSH_SUB(charged_particle_is_tolerance_reached)

    converged = this%classical_particle_t%is_tolerance_reached(tol)

    POP_SUB(charged_particle_is_tolerance_reached)
   end function charged_particle_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine charged_particle_iteration_info(this)
    class(charged_particle_t), intent(in) :: this

    PUSH_SUB(charged_particle_iteration_info)

    call this%classical_particle_t%iteration_info()

    POP_SUB(charged_particle_iteration_info)
  end subroutine charged_particle_iteration_info

  ! ---------------------------------------------------------
  subroutine charged_particle_update_quantity(this, iq)
    class(charged_particle_t), intent(inout) :: this
    integer,                   intent(in)    :: iq

    PUSH_SUB(charged_particle_update_quantity)

    select case (iq)
    case (CHARGE)
      ! The charged particle has a charge, but it is not necessary to update it, as it does not change with time.
      ! We still need to set its clock, so we set it to be in sync with the particle position.
      call this%quantities(iq)%clock%set_time(this%quantities(POSITION)%clock)
    case default
      ! Other quantities should be handled by the parent class
      call this%classical_particle_t%update_quantity(iq)
    end select

    POP_SUB(charged_particle_update_quantity)
  end subroutine charged_particle_update_quantity

 ! ---------------------------------------------------------
 subroutine charged_particle_update_exposed_quantity(partner, iq)
    class(charged_particle_t), intent(inout) :: partner
    integer,                   intent(in)    :: iq

    PUSH_SUB(charged_particle_update_exposed_quantity)

    select case (iq)
    case (CHARGE)
      ! The charged particle has a charge, but it is not necessary to update it, as it does not change with time.
      ! We still need to set its clock, so we set it to be in sync with the particle position.
      call partner%quantities(iq)%clock%set_time(partner%quantities(POSITION)%clock)
      call multisystem_debug_write_marker(partner%namespace, &
        event_clock_update_t(clock_name="quantity", clock_detail=QUANTITY_LABEL(iq), &
                             clock = partner%quantities(iq)%clock, action="set") )

    case default
      ! Other quantities should be handled by the parent class
      call partner%classical_particle_t%update_exposed_quantity(iq)
    end select

    POP_SUB(charged_particle_update_exposed_quantity)
  end subroutine charged_particle_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine charged_particle_copy_quantities_to_interaction(partner, interaction)
    class(charged_particle_t), intent(inout) :: partner
    class(interaction_t),      intent(inout) :: interaction

    PUSH_SUB(charged_particle_copy_quantities_to_interaction)

    select type (interaction)
    type is (coulomb_force_t)
      interaction%partner_np = 1

      if (.not. allocated(interaction%partner_charge)) then
        SAFE_ALLOCATE(interaction%partner_charge(1))
      end if
      interaction%partner_charge(1) = partner%charge(1)

      if (.not. allocated(interaction%partner_pos)) then
        SAFE_ALLOCATE(interaction%partner_pos(partner%space%dim, 1))
      end if
      interaction%partner_pos(:,1) = partner%pos(:, 1)

    class default
      call partner%classical_particle_t%copy_quantities_to_interaction(interaction)
    end select

    POP_SUB(charged_particle_copy_quantities_to_interaction)
  end subroutine charged_particle_copy_quantities_to_interaction

end module charged_particle_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

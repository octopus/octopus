!! Copyright (C) 2020 Heiko Appel
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
  use classical_particle_oct_m
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
  public ::               &
    charged_particle_t

  type, extends(classical_particle_t) :: charged_particle_t
    private

    FLOAT :: charge

  contains
    procedure :: add_interaction_partner => charged_particle_add_interaction_partner
    procedure :: has_interaction => charged_particle_has_interaction
    procedure :: do_td_operation => charged_particle_do_td
    procedure :: write_td_info => charged_particle_write_td_info
    procedure :: td_write_init => charged_particle_td_write_init
    procedure :: td_write_iter => charged_particle_td_write_iter
    procedure :: td_write_end => charged_particle_td_write_end
    procedure :: is_tolerance_reached => charged_particle_is_tolerance_reached
    procedure :: store_current_status => charged_particle_store_current_status
    procedure :: update_quantity => charged_particle_update_quantity
    procedure :: update_exposed_quantity => charged_particle_update_exposed_quantity
    procedure :: set_pointers_to_interaction => classical_set_pointers_to_interaction
    final :: charged_particle_finalize
  end type charged_particle_t

  interface charged_particle_t
    procedure charged_particle_init
  end interface charged_particle_t

contains

  ! ---------------------------------------------------------
  function charged_particle_init(this, namespace) result(sys)
    class(charged_particle_t), target, intent(in) :: this
    class(charged_particle_t), pointer  :: sys
    type(namespace_t),       intent(in) :: namespace

    PUSH_SUB(charged_particle_init)

    ! calling the parent routine here triggers a compiler error
    ! since classical_particle_init is not a member of the 
    ! ‘system_abst_t’ structure
    this%classical_particle_t%classical_particle_init(namespace)

    !%Variable ClassicalParticleCharge
    !%Type float
    !%Section ClassicalDynamics
    !%Description
    !% Charge of classical particle
    !%End
    call parse_variable(namespace, 'ClassicalParticleCharge', M_ONE, sys%charge)
    call messages_print_var_value(stdout, 'ClassicalParticleCharge', sys%charge)

    POP_SUB(charged_particle_init)
  end function charged_particle_init

  ! ---------------------------------------------------------
  subroutine charged_particle_add_interaction_partner(this, partner)
    class(charged_particle_t), target, intent(inout) :: this
    class(system_abst_t),                intent(in)    :: partner

    class(interaction_gravity_t), pointer :: gravity
    type(interaction_gravity_t) :: gravity_t

    PUSH_SUB(charged_particle_add_interaction_partner)


    POP_SUB(charged_particle_add_interaction_partner)
  end subroutine charged_particle_add_interaction_partner

  ! ---------------------------------------------------------
  logical function charged_particle_has_interaction(this, interaction)
    class(charged_particle_t),   intent(in) :: this
    class(interaction_abst_t),     intent(in) :: interaction

    PUSH_SUB(charged_particle_has_interaction)


    POP_SUB(charged_particle_has_interaction)
  end function charged_particle_has_interaction

  ! ---------------------------------------------------------
  subroutine charged_particle_do_td(this, operation)
    class(charged_particle_t), intent(inout) :: this
    integer,                     intent(in)    :: operation

    type(interaction_iterator_t) :: iter

    PUSH_SUB(charged_particle_do_td)


    POP_SUB(charged_particle_do_td)
  end subroutine charged_particle_do_td

  ! ---------------------------------------------------------
  logical function charged_particle_is_tolerance_reached(this, tol) result(converged)
    class(charged_particle_t),   intent(in)    :: this
    FLOAT,                         intent(in)    :: tol

    PUSH_SUB(charged_particle_is_tolerance_reached)

    POP_SUB(charged_particle_is_tolerance_reached)
   end function charged_particle_is_tolerance_reached

   ! ---------------------------------------------------------
   subroutine charged_particle_store_current_status(this)
     class(charged_particle_t),   intent(inout)    :: this

     PUSH_SUB(charged_particle_store_current_status) 

     POP_SUB(charged_particle_store_current_status)
   end subroutine charged_particle_store_current_status

  ! ---------------------------------------------------------
  subroutine charged_particle_write_td_info(this)
    class(charged_particle_t), intent(in) :: this

    integer :: idir
    character(len=20) :: fmt

    PUSH_SUB(charged_particle_write_td_info)

    POP_SUB(charged_particle_write_td_info)
  end subroutine charged_particle_write_td_info

  ! ---------------------------------------------------------
  subroutine charged_particle_td_write_init(this, dt)
    class(charged_particle_t), intent(inout) :: this
    FLOAT,                       intent(in)    :: dt

    PUSH_SUB(charged_particle_td_write_init)

    POP_SUB(charged_particle_td_write_init)
  end subroutine charged_particle_td_write_init

  ! ---------------------------------------------------------
  subroutine charged_particle_td_write_iter(this, iter)
    class(charged_particle_t), intent(inout) :: this
    integer,                 intent(in)    :: iter

    integer :: idir
    character(len=50) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(charged_particle_td_write_iter)

    POP_SUB(charged_particle_td_write_iter)
  end subroutine charged_particle_td_write_iter

  ! ---------------------------------------------------------
  subroutine charged_particle_td_write_end(this)
    class(charged_particle_t), intent(inout) :: this

    PUSH_SUB(charged_particle_td_write_end)

    POP_SUB(charged_particle_td_write_end)
  end subroutine charged_particle_td_write_end

  ! ---------------------------------------------------------
  logical function charged_particle_update_quantity(this, iq, clock) result(updated)
    class(charged_particle_t),   intent(inout) :: this
    integer,                       intent(in)    :: iq
    class(clock_t),                intent(in)    :: clock

    PUSH_SUB(charged_particle_update_quantity)

    POP_SUB(charged_particle_update_quantity)
  end function charged_particle_update_quantity

 ! ---------------------------------------------------------
 logical function charged_particle_update_exposed_quantity(this, iq, clock) result(updated)
    class(charged_particle_t),   intent(inout) :: this
    integer,                       intent(in)    :: iq
    class(clock_t),                intent(in)    :: clock

    PUSH_SUB(charged_particle_update_exposed_quantity)

    POP_SUB(charged_particle_update_exposed_quantity)
  end function charged_particle_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine classical_set_pointers_to_interaction(this, inter)
    class(charged_particle_t), target,  intent(in)   :: this
    class(interaction_abst_t),           intent(inout) :: inter

    PUSH_SUB(classical_set_pointers_to_interaction)

    POP_SUB(classical_set_pointers_to_interaction)
  end subroutine classical_set_pointers_to_interaction

  ! ---------------------------------------------------------
  subroutine charged_particle_finalize(this)
    type(charged_particle_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter
    class(interaction_abst_t), pointer :: interaction

    PUSH_SUB(charged_particle_finalize)

    POP_SUB(charged_particle_finalize)
  end subroutine charged_particle_finalize

end module charged_particle_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

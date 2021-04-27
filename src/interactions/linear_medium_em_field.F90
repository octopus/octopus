!! Copyright (C) 2021 F. Bonafé
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

module linear_medium_em_field_oct_m
  use clock_oct_m
  use ghost_interaction_oct_m
  use global_oct_m
  use iihash_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use interaction_partner_oct_m
  use io_function_oct_m
  use linked_list_oct_m
  use lorentz_force_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  implicit none

  private
  public ::               &
    linear_medium_em_field_t

  type, extends(interaction_partner_t) :: linear_medium_em_field_t
    private

  contains
    procedure :: calculate => linear_medium_em_field_calculate
    procedure :: allocate_memory => linear_medium_em_field_allocate
    procedure :: deallocate_memory => linear_medium_em_field_deallocate
    procedure :: update_exposed_quantities => linear_medium_em_field_update_exposed_quantities
    procedure :: update_exposed_quantity => linear_medium_em_field_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => linear_medium_em_field_copy_quantities_to_interaction
    final :: linear_medium_em_field_finalize
  end type linear_medium_em_field_t


  interface linear_medium_em_field_t
    module procedure linear_medium_em_field_init
  end interface linear_medium_em_field_t

contains

  function linear_medium_em_field_init(namespace) result(this)
    class(linear_medium_em_field_t), pointer :: this
    type(namespace_t), intent(in) :: namespace

    PUSH_SUB(linear_medium_em_field_init)

    SAFE_ALLOCATE(this)

    this%namespace = namespace_t("LinearMediumEMFieldInteraction", parent=namespace)

    POP_SUB(linear_medium_em_field_init)
  end function linear_medium_em_field_init

  ! ---------------------------------------------------------
  subroutine linear_medium_em_field_finalize(this)
    type(linear_medium_em_field_t), intent(inout) :: this

    PUSH_SUB(linear_medium_em_field_finalize)

    call this%deallocate_memory()

    POP_SUB(linear_medium_em_field_finalize)

  end subroutine linear_medium_em_field_finalize

  ! ---------------------------------------------------------
  subroutine linear_medium_em_field_allocate(this, mesh)
    class(linear_medium_em_field_t), intent(inout) :: this
    type(mesh_t),                intent(in)    :: mesh

    PUSH_SUB(linear_medium_em_field_allocate)

    POP_SUB(linear_medium_em_field_allocate)

  end subroutine linear_medium_em_field_allocate

  ! ---------------------------------------------------------
  subroutine linear_medium_em_field_deallocate(this)
    class(linear_medium_em_field_t), intent(inout) :: this

    PUSH_SUB(linear_medium_em_field_deallocate)

    POP_SUB(linear_medium_em_field_deallocate)
  end subroutine linear_medium_em_field_deallocate

  ! ---------------------------------------------------------
  logical function linear_medium_em_field_update_exposed_quantities(partner, requested_time, interaction) &
    result(allowed_to_update)
    class(linear_medium_em_field_t), intent(inout) :: partner
    type(clock_t),               intent(in)    :: requested_time
    class(interaction_t),        intent(inout) :: interaction

    PUSH_SUB(linear_medium_em_field_update_exposed_quantities)

    allowed_to_update = .true.

    call partner%clock%set_time(requested_time)

    select type (interaction)
    class default
      call partner%copy_quantities_to_interaction(interaction)
    end select

    POP_SUB(linear_medium_em_field_update_exposed_quantities)

  end function linear_medium_em_field_update_exposed_quantities

  ! ---------------------------------------------------------
  subroutine linear_medium_em_field_update_exposed_quantity(partner, iq)
    class(linear_medium_em_field_t),      intent(inout) :: partner
    integer,                          intent(in)    :: iq

    PUSH_SUB(linear_medium_em_field_update_exposed_quantities)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(linear_medium_em_field_update_exposed_quantities)

  end subroutine linear_medium_em_field_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine linear_medium_em_field_copy_quantities_to_interaction(partner, interaction)
    class(linear_medium_em_field_t),     intent(inout) :: partner
    class(interaction_t),            intent(inout) :: interaction

    integer :: ip

    PUSH_SUB(linear_medium_em_field_copy_quantities_to_interaction)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(linear_medium_em_field_copy_quantities_to_interaction)
  end subroutine linear_medium_em_field_copy_quantities_to_interaction


  ! ---------------------------------------------------------
  subroutine linear_medium_em_field_calculate(this, namespace)
    class(linear_medium_em_field_t), intent(inout) :: this
    type(namespace_t),           intent(in)    :: namespace

    PUSH_SUB(linear_medium_em_field_calculate)

    POP_SUB(linear_medium_em_field_calculate)
  end subroutine linear_medium_em_field_calculate


end module linear_medium_em_field_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

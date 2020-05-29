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

module interaction_partner_oct_m
  use clock_oct_m
  use global_oct_m
  use interaction_abst_oct_m
  use namespace_oct_m
  use quantity_oct_m
  implicit none

  private
  public ::               &
    interaction_partner_t

  !> Some interactions require a partner. This is usually a system, but it could
  !! also be some external entity, like an external field.
  !! An interaction partner must expose some quantities that the interaction can use.
  type, abstract :: interaction_partner_t
    private
    type(namespace_t), public :: namespace
    type(clock_t),     public :: clock

    type(quantity_t),  public :: quantities(MAX_QUANTITIES) !< Array of all possible quantities.
                                                            !< The elements of the array are accessed using the
                                                            !< quantity`s identifiers.
  contains
    procedure(interaction_partner_update_exposed_quantities),      deferred :: update_exposed_quantities
    procedure(interaction_partner_update_exposed_quantity),        deferred :: update_exposed_quantity
    procedure(interaction_partner_copy_quantities_to_interaction), deferred :: copy_quantities_to_interaction
  end type interaction_partner_t

  abstract interface

    ! ---------------------------------------------------------
    logical function interaction_partner_update_exposed_quantities(this, requested_time, interaction)
      import interaction_partner_t
      import clock_t
      import interaction_abst_t
      class(interaction_partner_t), intent(inout) :: this
      type(clock_t),                intent(in)    :: requested_time
      class(interaction_abst_t),    intent(inout) :: interaction
    end function interaction_partner_update_exposed_quantities

    ! ---------------------------------------------------------
    subroutine interaction_partner_update_exposed_quantity(this, iq, requested_time)
      import interaction_partner_t
      import clock_t
      class(interaction_partner_t),      intent(inout) :: this
      integer,                           intent(in)    :: iq
      class(clock_t),                    intent(in)    :: requested_time
    end subroutine interaction_partner_update_exposed_quantity

    ! ---------------------------------------------------------
    subroutine interaction_partner_copy_quantities_to_interaction(this, interaction)
      import interaction_partner_t
      import interaction_abst_t
      class(interaction_partner_t),     intent(inout) :: this
      class(interaction_abst_t),        intent(inout) :: interaction
    end subroutine interaction_partner_copy_quantities_to_interaction

  end interface

end module interaction_partner_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

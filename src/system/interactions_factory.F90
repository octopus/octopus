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

module interactions_factory_oct_m
  use global_oct_m
  use interaction_abst_oct_m
  use interaction_coulomb_force_oct_m
  use interaction_gravity_oct_m
  use interaction_lorentz_force_oct_m
  use interaction_partner_oct_m
  use interactions_factory_abst_oct_m
  use messages_oct_m
  implicit none

  private
  public :: interactions_factory_t

  integer, parameter, public :: &
    GRAVITY          = 1,       &
    LORENTZ_FORCE    = 2,       &
    COULOMB_FORCE    = 3

  type, extends(interactions_factory_abst_t) :: interactions_factory_t
  contains
    procedure :: create => interactions_factory_create
  end type interactions_factory_t

contains

  ! ---------------------------------------------------------------------------------------
  function interactions_factory_create(this, type, partner) result(interaction)
    class(interactions_factory_t),         intent(in)    :: this
    integer,                               intent(in)    :: type
    class(interaction_partner_t),  target, intent(inout) :: partner
    class(interaction_abst_t),             pointer       :: interaction

    PUSH_SUB(interactions_factory_create)

    select case (type)
    case (GRAVITY)
      interaction => interaction_gravity_t(partner)
    case (COULOMB_FORCE)
      interaction => interaction_coulomb_force_t(partner)
    case (LORENTZ_FORCE)
      interaction => interaction_lorentz_force_t(partner)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(interactions_factory_create)
  end function interactions_factory_create

end module interactions_factory_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

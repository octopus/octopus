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
  use interaction_oct_m
  use coulomb_force_oct_m
  use gravity_oct_m
  use linear_medium_em_field_oct_m
  use lorentz_force_oct_m
  use interaction_partner_oct_m
  use interactions_factory_abst_oct_m
  use messages_oct_m
  implicit none

  private
  public :: interactions_factory_t

  integer, parameter, public :: &
    GRAVITY          = 1,       &
    LORENTZ_FORCE    = 2,       &
    COULOMB_FORCE    = 3,       &
    LINEAR_MEDIUM_EM_FIELD = 4

  type, extends(interactions_factory_abst_t) :: interactions_factory_t
  contains
    procedure :: create => interactions_factory_create
    procedure :: default_mode => interactions_factory_default_mode
    procedure :: block_name => interactions_factory_block_name
  end type interactions_factory_t

contains

  ! ---------------------------------------------------------------------------------------
  function interactions_factory_create(this, type, partner) result(interaction)
    class(interactions_factory_t),         intent(in)    :: this
    integer,                               intent(in)    :: type
    class(interaction_partner_t),  target, intent(inout) :: partner
    class(interaction_t),                  pointer       :: interaction

    PUSH_SUB(interactions_factory_create)

    !%Variable Interactions
    !%Type block
    !%Section System
    !%Description
    !% This input option controls the interactions between systems. It basically
    !% allows to select which systems will interact with another system through
    !% a given interaction type. The format of the block is the following:
    !%
    !%  <br>%</tt>Namespace.Interactions
    !%   <br>&nbsp;&nbsp;interaction_type | interaction_mode | ...
    !%  <br>%</tt>
    !%
    !% Here is an example to better understand how this works:
    !%
    !%  <br>%</tt>SystemA.Interactions
    !%   <br>&nbsp;&nbsp;gravity | all_except | "SystemB"
    !%  <br>%</tt>
    !%
    !% This means that SystemA and all the systems that belong to the same
    !% namespace (i.e., all its subsystems) will interact through gravity with
    !% all interaction partners that are also able to interact through gravity,
    !% except with SystemB. Note that the opposite is not true so, although
    !% clearly unphysical, this will not prevent SystemB from feeling the
    !% gravity from SystemA (in <tt>Octopus</tt> the interactions are always
    !% one-sided).
    !%
    !% NB: Each interaction type should only appear once in the block. Any
    !% further instances beyond the first will be ignored.
    !%
    !% Available modes and interaction types:
    !%Option no_partners -1
    !%  (interaction mode)
    !% Do not interact with any partner.
    !%Option all_partners -2
    !%  (interaction mode)
    !% Interact with all available partners.
    !%Option only_partners -3
    !%  (interaction mode)
    !% Interact only with some specified partners. A list of partner names must
    !% be given.
    !%Option all_except -4
    !%  (interaction mode)
    !% Interact with all available partners except with some specified
    !% partners. A list of partner names to exclude must be given.
    !%Option gravity 1
    !%  (interaction type)
    !% Gravity interaction between two masses.
    !%Option lorenz_force 2
    !%  (interaction type)
    !% Lorenz force resulting from an EM field acting on a moving charge.
    !%Option coulomb_force 3
    !%  (interaction type)
    !% Coulomb force between two charged particles.
    !%Option linear_medium_em_field 4
    !%  (interaction type)
    !% Polarization of linear medium in the presence of a field.
    !%End
    select case (type)
    case (GRAVITY)
      interaction => gravity_t(partner)
    case (COULOMB_FORCE)
      interaction => coulomb_force_t(partner)
    case (LORENTZ_FORCE)
      interaction => lorentz_force_t(partner)
    case (LINEAR_MEDIUM_EM_FIELD)
      interaction => linear_medium_em_field_t(partner)
    case default
      ! This should never happen, as this is handled in
      ! interactions_factory_abst_create_interactions
    end select

    POP_SUB(interactions_factory_create)
  end function interactions_factory_create

  ! ---------------------------------------------------------------------------------------
  integer function interactions_factory_default_mode(this, type) result(mode)
    class(interactions_factory_t), intent(in)    :: this
    integer,                       intent(in)    :: type

    PUSH_SUB(interactions_factory_default_mode)

    select case (type)
    case (GRAVITY)
      mode = NO_PARTNERS
    case (COULOMB_FORCE)
      mode = ALL_PARTNERS
    case (LORENTZ_FORCE)
      mode = ALL_PARTNERS
    case (LINEAR_MEDIUM_EM_FIELD)
      mode = ALL_PARTNERS
    case default
      message(1) = "Unknown interaction type"
      call messages_fatal(1)
    end select

    POP_SUB(interactions_factory_default_mode)
  end function interactions_factory_default_mode

  ! ---------------------------------------------------------------------------------------
  character(len=80) function interactions_factory_block_name(this) result(name)
    class(interactions_factory_t), intent(in)    :: this

    PUSH_SUB(interactions_factory_block_name)

    name = "Interactions"

    POP_SUB(interactions_factory_block_name)
  end function interactions_factory_block_name

end module interactions_factory_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

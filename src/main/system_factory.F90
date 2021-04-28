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

module system_factory_oct_m
  use charged_particle_oct_m
  use classical_particle_oct_m
  use global_oct_m
  use messages_oct_m
  use multisystem_basic_oct_m
  use namespace_oct_m
  use electrons_oct_m
  use system_oct_m
  use system_dftb_oct_m
  use system_factory_abst_oct_m
  use system_linear_medium_oct_m
  use system_mxll_oct_m
  implicit none

  private
  public ::                         &
    system_factory_t

  integer, parameter ::             &
    SYSTEM_ELECTRONIC         = 1,  &
    SYSTEM_MAXWELL            = 2,  &
    SYSTEM_CLASSICAL_PARTICLE = 3,  &
    SYSTEM_CHARGED_PARTICLE   = 4,  &
    SYSTEM_DFTBPLUS           = 5,  &
    SYSTEM_LINEAR_MEDIUM      = 6,  &
    SYSTEM_MULTISYSTEM        = 7


  type, extends(system_factory_abst_t) :: system_factory_t
  contains
    procedure :: create => system_factory_create
    procedure :: block_name => system_factory_block_name
  end type system_factory_t

contains

  ! ---------------------------------------------------------------------------------------
  recursive function system_factory_create(this, namespace, name, type) result(system)
    class(system_factory_t), intent(in) :: this
    type(namespace_t),       intent(in) :: namespace
    character(len=*),        intent(in) :: name
    integer,                 intent(in) :: type
    class(system_t),         pointer    :: system

    PUSH_SUB(system_factory_create)

    !%Variable Systems
    !%Type block
    !%Section System
    !%Description
    !% List of systems that will be treated in the calculation.
    !% The first column should be a string containing the system name.
    !% The second column should be the system type. See below for a list of
    !% available system types.
    !%Option electronic 1
    !% An electronic system. (not fully implemented yet)
    !%Option maxwell 2
    !% A maxwell system.
    !%Option classical_particle 3
    !% A classical particle. Used for testing purposes only.
    !%Option charged_particle 4
    !% A charged classical particle.
    !%Option dftbplus 5
    !% A DFTB+ system
    !%Option linear_medium 6
    !% A linear medium for classical electrodynamics.
    !%Option multisystem 7
    !% A system containing other systems.
    !%End
    select case (type)
    case (SYSTEM_MULTISYSTEM)
      system => multisystem_basic_t(namespace_t(name, parent=namespace), this)
    case (SYSTEM_ELECTRONIC)
      system => electrons_t(namespace_t(name, parent=namespace))
    case (SYSTEM_MAXWELL)
      system => system_mxll_t(namespace_t(name, parent=namespace))
    case (SYSTEM_CLASSICAL_PARTICLE)
      system => classical_particle_t(namespace_t(name, parent=namespace))
    case (SYSTEM_CHARGED_PARTICLE)
      system => charged_particle_t(namespace_t(name, parent=namespace))
    case (SYSTEM_DFTBPLUS)
      system => system_dftb_t(namespace_t(name, parent=namespace))
    case (SYSTEM_LINEAR_MEDIUM)
      system => system_linear_medium_t(namespace_t(name, parent=namespace))
    case default
      system => null()
    end select

    POP_SUB(system_factory_create)
  end function system_factory_create

  ! ---------------------------------------------------------------------------------------
  character(len=80) function system_factory_block_name(this) result(name)
    class(system_factory_t), intent(in)    :: this

    PUSH_SUB(system_factory_block_name)

    name = "Systems"

    POP_SUB(system_factory_block_name)
  end function system_factory_block_name

end module system_factory_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

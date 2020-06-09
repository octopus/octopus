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
  use classical_particle_oct_m
  use global_oct_m
  use messages_oct_m
  use multisystem_oct_m
  use namespace_oct_m
  use system_oct_m
  use system_abst_oct_m
  use system_factory_abst_oct_m
  use system_mxll_oct_m
  implicit none

  private
  public ::                         &
    system_factory_t

  integer, parameter ::             &
    SYSTEM_ELECTRONIC         = 1,  &
    SYSTEM_MAXWELL            = 2,  &
    SYSTEM_CLASSICAL_PARTICLE = 3,  &
    SYSTEM_MULTISYSTEM        = 4

  type, extends(system_factory_abst_t) :: system_factory_t
  contains
    procedure :: create => system_factory_create
  end type system_factory_t

contains

  ! ---------------------------------------------------------------------------------------
  recursive function system_factory_create(this, namespace, name, type) result(system)
    class(system_factory_t), intent(in) :: this
    type(namespace_t),       intent(in) :: namespace
    character(len=*),        intent(in) :: name
    integer,                 intent(in) :: type
    class(system_abst_t),    pointer    :: system

    PUSH_SUB(system_factory_create)

    select case (type)
    case (SYSTEM_MULTISYSTEM)
      system => multisystem_t(namespace_t(name, parent=namespace), this)
    case (SYSTEM_MAXWELL)
      system => system_mxll_t(namespace_t(name, parent=namespace))
    case (SYSTEM_CLASSICAL_PARTICLE)
      system => classical_particle_t(namespace_t(name, parent=namespace))
    case default
      system => null()
    end select

    POP_SUB(system_factory_create)
  end function system_factory_create

end module system_factory_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

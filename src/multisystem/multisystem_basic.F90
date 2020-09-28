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

module multisystem_basic_oct_m
  use global_oct_m
  use messages_oct_m
  use multisystem_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use system_factory_abst_oct_m
  implicit none

  private
  public ::               &
    multisystem_basic_t

  type, extends(multisystem_t) :: multisystem_basic_t

  contains
    final :: multisystem_basic_finalizer
  end type multisystem_basic_t

  interface multisystem_basic_t
    procedure multisystem_basic_constructor
  end interface multisystem_basic_t

contains

  ! ---------------------------------------------------------------------------------------
  recursive function multisystem_basic_constructor(namespace, factory) result(system)
    type(namespace_t),            intent(in) :: namespace
    class(system_factory_abst_t), intent(in) :: factory
    class(multisystem_basic_t),   pointer    :: system

    PUSH_SUB(multisystem_basic_constructor)

    SAFE_ALLOCATE(system)

    call multisystem_basic_init(system, namespace, factory)

    POP_SUB(multisystem_basic_constructor)
  end function multisystem_basic_constructor

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_basic_init(this, namespace, factory)
    class(multisystem_t),      intent(inout) :: this
    type(namespace_t),            intent(in) :: namespace
    class(system_factory_abst_t), intent(in) :: factory

    PUSH_SUB(multisystem_basic_init)

    call multisystem_init(this, namespace, factory)

    POP_SUB(multisystem_basic_init)
  end subroutine multisystem_basic_init

  ! ---------------------------------------------------------
  recursive subroutine multisystem_basic_finalizer(this)
    type(multisystem_basic_t), intent(inout) :: this

    PUSH_SUB(multisystem_basic_finalizer)

    call multisystem_end(this)

    POP_SUB(multisystem_basic_finalizer)
  end subroutine multisystem_basic_finalizer

end module multisystem_basic_oct_m

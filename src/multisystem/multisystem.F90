!! Copyright (C) 2019-2020 M. Oliveira
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

module multisystem_oct_m
  use clock_oct_m
  use global_oct_m
  use ghost_interaction_oct_m
  use interaction_oct_m
  use interaction_with_partner_oct_m
  use io_oct_m
  use loct_oct_m
  use messages_oct_m
  use multisystem_abst_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use system_oct_m
  use system_replica_oct_m
  use system_factory_abst_oct_m
  implicit none

  private
  public ::               &
    multisystem_t

  type, extends(multisystem_abst_t) :: multisystem_t
  contains
    final :: multisystem_finalizer
  end type multisystem_t

  interface multisystem_t
    procedure multisystem_constructor
  end interface multisystem_t

contains

  ! ---------------------------------------------------------------------------------------
  recursive function multisystem_constructor(namespace, factory, system_replica) result(system)
    class(multisystem_t),         pointer    :: system
    type(namespace_t),            intent(in) :: namespace
    class(system_factory_abst_t), intent(in) :: factory
    type(system_replica_t),       intent(inout) :: system_replica

    PUSH_SUB(multisystem_constructor)

    SAFE_ALLOCATE(system)

    call multisystem_init(system, namespace, factory, system_replica)

    POP_SUB(multisystem_constructor)
  end function multisystem_constructor

  ! ---------------------------------------------------------------------------------------
  recursive subroutine multisystem_init(this, namespace, factory, system_replica)
    class(multisystem_t),      intent(inout) :: this
    type(namespace_t),            intent(in) :: namespace
    class(system_factory_abst_t), intent(in) :: factory
    type(system_replica_t),       intent(inout) :: system_replica

    PUSH_SUB(multisystem_init)

    call multisystem_abst_init(this, namespace, factory, system_replica)

    POP_SUB(multisystem_init)
  end subroutine multisystem_init


  ! ---------------------------------------------------------
  recursive subroutine multisystem_finalizer(this)
    type(multisystem_t), intent(inout) :: this

    type(system_iterator_t) :: iter
    class(system_t), pointer :: system

    PUSH_SUB(multisystem_finalizer)

    call iter%start(this%list)
    do while (iter%has_next())
      system => iter%get_next()
      SAFE_DEALLOCATE_P(system)
    end do

    call system_end(this)

    POP_SUB(multisystem_finalizer)
  end subroutine multisystem_finalizer

end module multisystem_oct_m

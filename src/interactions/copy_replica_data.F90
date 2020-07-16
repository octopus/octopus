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

module copy_replica_data_oct_m
  use global_oct_m
  use interaction_with_partner_oct_m
  use interaction_partner_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use quantity_oct_m

  implicit none

  private
  public ::                &
    copy_replica_data_t

  type, extends(interaction_with_partner_t) :: copy_replica_data_t
    private
    integer :: dim

    integer :: num_replicas

    FLOAT, pointer, public :: system_pos(:)
    FLOAT, pointer, public :: system_vel(:)
    FLOAT, pointer, public :: system_pos2(:)
    FLOAT, pointer, public :: system_vel2(:)
    FLOAT, pointer, public :: system_varpos(:)
    FLOAT, pointer, public :: system_varvel(:)


    FLOAT, allocatable, public :: partner_pos(:)
    FLOAT, allocatable, public :: partner_vel(:)

  contains
    procedure :: init => copy_replica_data_init
    procedure :: calculate => copy_replica_data_calculate
    final :: copy_replica_data_finalize
  end type copy_replica_data_t

  interface copy_replica_data_t
    module procedure copy_replica_data_constructor
  end interface copy_replica_data_t

contains

  ! ---------------------------------------------------------
  function copy_replica_data_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(copy_replica_data_t),           pointer       :: this

    integer :: num_replicas_default = 0

    PUSH_SUB(copy_replica_data_constructor)

    SAFE_ALLOCATE(this)

    this%label = "copy_replica_data"
    this%partner => partner

    num_replicas_default = 0
    call parse_variable(partner%namespace, 'SystemReplicas', num_replicas_default, this%num_replicas)

    ! The copy_replica_data interaction needs the position, velocity of the system
    this%n_system_quantities = 2
    SAFE_ALLOCATE(this%system_quantities(this%n_system_quantities))
    this%system_quantities(1) = POSITION
    this%system_quantities(2) = VELOCITY
    nullify(this%system_pos)
    nullify(this%system_vel)
    nullify(this%system_pos2)
    nullify(this%system_vel2)
    nullify(this%system_varpos)
    nullify(this%system_varvel)

    ! The copy_replica_data interaction partner needs the particle position and velocity
    this%n_partner_quantities = 2
    SAFE_ALLOCATE(this%partner_quantities(this%n_partner_quantities))
    this%partner_quantities(1) = POSITION
    this%partner_quantities(2) = VELOCITY

    POP_SUB(copy_replica_data_constructor)
  end function copy_replica_data_constructor

  ! ---------------------------------------------------------
  subroutine copy_replica_data_init(this, dim, system_quantities, system_pos, system_vel, &
          system_pos2, system_vel2, system_varpos, system_varvel)
    class(copy_replica_data_t),           intent(inout) :: this
    integer,                              intent(in)    :: dim
    type(quantity_t),                     intent(inout) :: system_quantities(:)
    FLOAT,                        target, intent(in)    :: system_pos(:)
    FLOAT,                        target, intent(in)    :: system_vel(:)
    FLOAT,                        target, intent(in)    :: system_pos2(:)
    FLOAT,                        target, intent(in)    :: system_vel2(:)
    FLOAT,                        target, intent(in)    :: system_varpos(:)
    FLOAT,                        target, intent(in)    :: system_varvel(:)

    PUSH_SUB(copy_replica_data_init)

    this%dim = dim
    SAFE_ALLOCATE(this%partner_pos(dim))
    SAFE_ALLOCATE(this%partner_vel(dim))

    system_quantities(POSITION)%required = .true.
    system_quantities(VELOCITY)%required = .true.
    this%system_pos => system_pos
    this%system_vel => system_vel
    this%system_pos2 => system_pos2
    this%system_vel2 => system_vel2
    this%system_varpos => system_varpos
    this%system_varvel => system_varvel

    POP_SUB(copy_replica_data_init)
  end subroutine copy_replica_data_init

  ! ---------------------------------------------------------
  subroutine copy_replica_data_calculate(this, namespace)
    class(copy_replica_data_t),             intent(inout) :: this
    type(namespace_t),                  intent(in)    :: namespace

    PUSH_SUB(copy_replica_data_calculate)

    this%system_pos = this%system_pos + this%partner_pos/this%num_replicas
    this%system_vel = this%system_vel + this%partner_vel/this%num_replicas
    this%system_pos2 = this%system_pos2 + this%partner_pos**2/this%num_replicas
    this%system_vel2 = this%system_vel2 + this%partner_vel**2/this%num_replicas

    POP_SUB(copy_replica_data_calculate)
  end subroutine copy_replica_data_calculate

  ! ---------------------------------------------------------
  subroutine copy_replica_data_finalize(this)
    type(copy_replica_data_t), intent(inout) :: this

    PUSH_SUB(copy_replica_data_finalize)

    nullify(this%system_pos)
    nullify(this%system_vel)
    nullify(this%system_pos2)
    nullify(this%system_vel2)
    nullify(this%system_varpos)
    nullify(this%system_varvel)

    call interaction_with_partner_end(this)

    POP_SUB(copy_replica_data_finalize)
  end subroutine copy_replica_data_finalize

end module copy_replica_data_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

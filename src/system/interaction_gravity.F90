!! Copyright (C) 2020 M. Oliveira, Heiko Appel
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

module interaction_gravity_oct_m
  use global_oct_m
  use interaction_with_partner_oct_m
  use interaction_partner_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use quantity_oct_m

  implicit none

  private
  public ::                &
    interaction_gravity_t

  type, extends(interaction_with_partner_t) :: interaction_gravity_t
    integer :: dim

    FLOAT :: force(MAX_DIM)

    FLOAT, pointer :: system_mass
    FLOAT, pointer :: system_pos(:)

    FLOAT :: partner_mass
    FLOAT, allocatable :: partner_pos(:)

  contains
    procedure :: init => interaction_gravity_init
    procedure :: calculate => interaction_gravity_calculate
    final :: interaction_gravity_finalize
  end type interaction_gravity_t

  interface interaction_gravity_t
    module procedure interaction_gravity_constructor
  end interface interaction_gravity_t

contains
  ! ---------------------------------------------------------
  function interaction_gravity_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(interaction_gravity_t),         pointer       :: this

    PUSH_SUB(interaction_gravity_constructor)

    SAFE_ALLOCATE(this)

    this%label = "gravity"

    this%partner => partner

    ! Gravity interaction needs two quantities from each system: the position and the mass
    ! From the sytem:
    this%n_system_quantities = 2
    SAFE_ALLOCATE(this%system_quantities(this%n_system_quantities))
    this%system_quantities(1) = POSITION
    this%system_quantities(2) = MASS
    nullify(this%system_mass)
    nullify(this%system_pos)

    ! From the partner:
    this%n_partner_quantities = 2
    SAFE_ALLOCATE(this%partner_quantities(this%n_partner_quantities))
    this%partner_quantities(1) = POSITION
    this%partner_quantities(2) = MASS
    this%partner%quantities(POSITION)%required = .true.
    this%partner%quantities(MASS)%required = .true.

    POP_SUB(interaction_gravity_constructor)
  end function interaction_gravity_constructor

  ! ---------------------------------------------------------
  subroutine interaction_gravity_init(this, dim, system_quantities, system_mass, system_pos)
    class(interaction_gravity_t),         intent(inout) :: this
    integer,                              intent(in)    :: dim
    type(quantity_t),                     intent(inout) :: system_quantities(:)
    FLOAT,                        target, intent(in)    :: system_mass
    FLOAT,                        target, intent(in)    :: system_pos(:)

    PUSH_SUB(interaction_gravity_init)

    this%dim = dim
    SAFE_ALLOCATE(this%partner_pos(dim))

    system_quantities(POSITION)%required = .true.
    system_quantities(MASS)%required = .true.
    this%system_mass => system_mass
    this%system_pos => system_pos

    POP_SUB(interaction_gravity_init)
  end subroutine interaction_gravity_init

  ! ---------------------------------------------------------
  subroutine interaction_gravity_calculate(this, namespace)
    class(interaction_gravity_t), intent(inout) :: this
    type(namespace_t),            intent(in)    :: namespace

    FLOAT, parameter :: GG = CNST(6.67430e-11)
    FLOAT :: dist3

    PUSH_SUB(interaction_gravity_calculate)

    dist3 = sum((this%partner_pos(1:this%dim) - this%system_pos(1:this%dim))**2)**(M_THREE/M_TWO)

    this%force(1:this%dim) = (this%partner_pos(1:this%dim) - this%system_pos(1:this%dim)) &
      / dist3 * (GG * this%system_mass * this%partner_mass)

    POP_SUB(interaction_gravity_calculate)
  end subroutine interaction_gravity_calculate

  ! ---------------------------------------------------------
  subroutine interaction_gravity_finalize(this)
    type(interaction_gravity_t), intent(inout) :: this

    PUSH_SUB(interaction_gravity_finalize)

    this%force = M_ZERO
    nullify(this%system_mass)
    nullify(this%system_pos)
    SAFE_DEALLOCATE_A(this%partner_pos)

    call interaction_with_partner_end(this)

    POP_SUB(interaction_gravity_finalize)
  end subroutine interaction_gravity_finalize

end module interaction_gravity_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

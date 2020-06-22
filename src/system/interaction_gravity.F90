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
    procedure :: calculate => interaction_gravity_calculate
    procedure :: end => interaction_gravity_end
    final :: interaction_gravity_finalize
  end type interaction_gravity_t

  interface interaction_gravity_t
    module procedure interaction_gravity_init
  end interface interaction_gravity_t

contains

  ! ---------------------------------------------------------

  function interaction_gravity_init(dim, partner) result(this)
    integer,                              intent(in)    :: dim
    class(interaction_partner_t), target, intent(inout) :: partner
    class(interaction_gravity_t),         pointer       :: this

    PUSH_SUB(interaction_gravity_init)

    SAFE_ALLOCATE(this)

    this%dim = dim
    this%partner => partner

    ! Gravity interaction needs two quantities from each system: the position and the mass
    ! From the sytem:
    this%n_system_quantities = 2
    SAFE_ALLOCATE(this%system_quantities(this%n_system_quantities))
    this%system_quantities(1) = POSITION
    this%system_quantities(2) = MASS
    ! From the partner:
    this%n_partner_quantities = 2
    SAFE_ALLOCATE(this%partner_quantities(this%n_partner_quantities))
    this%partner_quantities(1) = POSITION
    this%partner_quantities(2) = MASS
    this%partner%quantities(POSITION)%required = .true.
    this%partner%quantities(MASS)%required = .true.
    SAFE_ALLOCATE(this%partner_pos(dim))

    POP_SUB(interaction_gravity_init)
  end function interaction_gravity_init

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
  subroutine interaction_gravity_end(this)
    class(interaction_gravity_t), intent(inout) :: this

    PUSH_SUB(interaction_gravity_end)

    this%force = M_ZERO
    nullify(this%system_mass)
    nullify(this%system_pos)
    SAFE_DEALLOCATE_A(this%partner_pos)

    call interaction_with_partner_end(this)

    POP_SUB(interaction_gravity_end)
  end subroutine interaction_gravity_end

  ! ---------------------------------------------------------
  subroutine interaction_gravity_finalize(this)
    type(interaction_gravity_t), intent(inout) :: this

    PUSH_SUB(interaction_gravity_finalize)

    call this%end()

    POP_SUB(interaction_gravity_finalize)
  end subroutine interaction_gravity_finalize

end module interaction_gravity_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

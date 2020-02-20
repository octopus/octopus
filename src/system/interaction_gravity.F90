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

module interaction_gravity_oct_m
  use global_oct_m
  use interaction_abst_oct_m
  use messages_oct_m
  use observable_oct_m
  use profiling_oct_m
  use system_abst_oct_m

  implicit none

  private
  public ::                &
    interaction_gravity_t

  type, extends(interaction_abst_t) :: interaction_gravity_t
    integer :: dim

    FLOAT :: force(MAX_DIM)

    class(system_abst_t), pointer :: partner

    FLOAT, pointer :: system_mass
    FLOAT, pointer :: system_pos(:)

    FLOAT, pointer :: partner_mass
    FLOAT, pointer :: partner_pos(:)

  contains
    procedure :: update => interaction_gravity_update
    procedure :: end => interaction_gravity_end
    final :: interaction_gravity_finalize
  end type interaction_gravity_t

  interface interaction_gravity_t
    module procedure interaction_gravity_init
  end interface interaction_gravity_t

contains

  ! ---------------------------------------------------------

  function interaction_gravity_init(dim, partner) result(this)
    integer,                      intent(in) :: dim
    class(system_abst_t), target, intent(in) :: partner
    class(interaction_gravity_t), pointer    :: this

    PUSH_SUB(interaction_gravity_init)

    SAFE_ALLOCATE(this)

    this%dim = dim
    this%partner => partner

    !Gravity interaction needs only one observable from each system, which is the position
    this%n_system_observables = 1
    this%n_partner_observables = 1
    SAFE_ALLOCATE(this%system_observables(this%n_system_observables))
    SAFE_ALLOCATE(this%partner_observables(this%n_partner_observables))
    this%system_observables(1) = POSITION
    this%partner_observables(1) = POSITION

    call partner%set_pointers_to_interaction(this)

    POP_SUB(interaction_gravity_init)
  end function interaction_gravity_init

  ! ---------------------------------------------------------
  subroutine interaction_gravity_update(this)
    class(interaction_gravity_t), intent(inout) :: this

    FLOAT, parameter :: GG = CNST(6.67430e-11)
    FLOAT :: dist3

    PUSH_SUB(interaction_gravity_update)

    ASSERT(associated(this%partner_pos))
    ASSERT(associated(this%system_pos))
    ASSERT(associated(this%partner_mass))
    ASSERT(associated(this%system_mass))

    ! Now calculate the gravitational force
    dist3 = sum((this%partner_pos(1:this%dim) - this%system_pos(1:this%dim))**2)**(M_THREE/M_TWO)

    this%force(1:this%dim) = (this%partner_pos(1:this%dim) - this%system_pos(1:this%dim)) &
                                      / dist3 * (GG * this%system_mass * this%partner_mass)

    POP_SUB(interaction_gravity_update)
  end subroutine interaction_gravity_update

  ! ---------------------------------------------------------
  subroutine interaction_gravity_end(this)
    class(interaction_gravity_t), intent(inout) :: this

    PUSH_SUB(interaction_gravity_end)

    nullify(this%partner)
    this%force = M_ZERO
    nullify(this%system_mass)
    nullify(this%system_pos)
    nullify(this%partner_mass)
    nullify(this%partner_pos)

    call interaction_abst_end(this)

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

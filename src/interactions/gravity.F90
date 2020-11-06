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

module gravity_oct_m
  use force_interaction_oct_m
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
    gravity_t

  !> Gravity interaction between two particles. This should be used
  !! for testing purposes only. Note that this interaction assumes all
  !! quantities are in S.I. units instead of atomic units.
  type, extends(force_interaction_t) :: gravity_t
    private
    integer :: dim

    FLOAT, pointer :: system_mass
    FLOAT, pointer :: system_pos(:)

    FLOAT, public :: partner_mass
    FLOAT, allocatable, public :: partner_pos(:)

  contains
    procedure :: init => gravity_init
    procedure :: calculate => gravity_calculate
    final :: gravity_finalize
  end type gravity_t

  interface gravity_t
    module procedure gravity_constructor
  end interface gravity_t

contains
  ! ---------------------------------------------------------
  function gravity_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(gravity_t),                     pointer       :: this

    PUSH_SUB(gravity_constructor)

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

    POP_SUB(gravity_constructor)
  end function gravity_constructor

  ! ---------------------------------------------------------
  subroutine gravity_init(this, dim, system_quantities, system_mass, system_pos)
    class(gravity_t),                     intent(inout) :: this
    integer,                              intent(in)    :: dim
    type(quantity_t),                     intent(inout) :: system_quantities(:)
    FLOAT,                        target, intent(in)    :: system_mass
    FLOAT,                        target, intent(in)    :: system_pos(:)

    PUSH_SUB(gravity_init)

    this%dim = dim
    SAFE_ALLOCATE(this%partner_pos(dim))

    system_quantities(POSITION)%required = .true.
    system_quantities(MASS)%required = .true.
    this%system_mass => system_mass
    this%system_pos => system_pos

    POP_SUB(gravity_init)
  end subroutine gravity_init

  ! ---------------------------------------------------------
  subroutine gravity_calculate(this, namespace)
    class(gravity_t),             intent(inout) :: this
    type(namespace_t),            intent(in)    :: namespace

    FLOAT, parameter :: GG = CNST(6.67430e-11) ! In S.I. units!
    FLOAT :: dist3

    PUSH_SUB(gravity_calculate)

    dist3 = sum((this%partner_pos(1:this%dim) - this%system_pos(1:this%dim))**2)**(M_THREE/M_TWO)

    this%force(1:this%dim) = (this%partner_pos(1:this%dim) - this%system_pos(1:this%dim)) &
      / dist3 * (GG * this%system_mass * this%partner_mass)

    POP_SUB(gravity_calculate)
  end subroutine gravity_calculate

  ! ---------------------------------------------------------
  subroutine gravity_finalize(this)
    type(gravity_t), intent(inout) :: this

    PUSH_SUB(gravity_finalize)

    this%force = M_ZERO
    nullify(this%system_mass)
    nullify(this%system_pos)
    SAFE_DEALLOCATE_A(this%partner_pos)

    call interaction_with_partner_end(this)

    POP_SUB(gravity_finalize)
  end subroutine gravity_finalize

end module gravity_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

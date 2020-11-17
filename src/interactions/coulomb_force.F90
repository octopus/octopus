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

module coulomb_force_oct_m
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
    coulomb_force_t

  type, extends(force_interaction_t) :: coulomb_force_t
    private
    integer :: dim

    FLOAT, pointer :: system_charge
    FLOAT, pointer :: system_pos(:)

    FLOAT, public :: partner_charge
    FLOAT, allocatable, public :: partner_pos(:)

  contains
    procedure :: init => coulomb_force_init
    procedure :: calculate => coulomb_force_calculate
    final :: coulomb_force_finalize
  end type coulomb_force_t

  interface coulomb_force_t
    module procedure coulomb_force_constructor
  end interface coulomb_force_t

contains

  ! ---------------------------------------------------------
  function coulomb_force_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(coulomb_force_t),               pointer       :: this

    PUSH_SUB(coulomb_force_constructor)

    SAFE_ALLOCATE(this)

    this%label = "coulomb_force"

    this%partner => partner

    ! The Coulomb force needs the position and charge of the system
    this%n_system_quantities = 2
    SAFE_ALLOCATE(this%system_quantities(this%n_system_quantities))
    this%system_quantities(1) = POSITION
    this%system_quantities(2) = CHARGE
    nullify(this%system_pos)
    nullify(this%system_charge)

    ! The Coulomb force needs the position and charge of the partner
    this%n_partner_quantities = 2
    SAFE_ALLOCATE(this%partner_quantities(this%n_partner_quantities))
    this%partner_quantities(1) = POSITION
    this%partner_quantities(2) = CHARGE
    this%partner%quantities(POSITION)%required = .true.
    this%partner%quantities(CHARGE)%required = .true.

    POP_SUB(coulomb_force_constructor)
  end function coulomb_force_constructor

  ! ---------------------------------------------------------
  subroutine coulomb_force_init(this, dim, system_quantities, system_charge, system_pos)
    class(coulomb_force_t),             intent(inout) :: this
    integer,                            intent(in)    :: dim
    type(quantity_t),                   intent(inout) :: system_quantities(:)
    FLOAT,                      target, intent(in)    :: system_charge
    FLOAT,                      target, intent(in)    :: system_pos(:)

    PUSH_SUB(coulomb_force_init)

    this%dim = dim
    SAFE_ALLOCATE(this%partner_pos(dim))

    system_quantities(POSITION)%required = .true.
    system_quantities(CHARGE)%required = .true.
    this%system_charge => system_charge
    this%system_pos => system_pos

    POP_SUB(coulomb_force_init)
  end subroutine coulomb_force_init

  ! ---------------------------------------------------------
  subroutine coulomb_force_calculate(this, namespace)
    class(coulomb_force_t),             intent(inout) :: this
    type(namespace_t),                  intent(in)    :: namespace

    FLOAT, parameter :: COULCONST = M_ONE ! Coulomb constant in atomic units
    FLOAT :: dist3

    PUSH_SUB(coulomb_force_calculate)

    dist3 = sum((this%partner_pos(1:this%dim) - this%system_pos(1:this%dim))**2)**(M_THREE/M_TWO)

    this%force(1:this%dim) = -(this%partner_pos(1:this%dim) - this%system_pos(1:this%dim)) &
      / (dist3 + M_EPSILON) * (COULCONST * this%system_charge * this%partner_charge)


    POP_SUB(coulomb_force_calculate)
  end subroutine coulomb_force_calculate

  ! ---------------------------------------------------------
  subroutine coulomb_force_finalize(this)
    type(coulomb_force_t), intent(inout) :: this

    PUSH_SUB(coulomb_force_finalize)

    this%force = M_ZERO
    nullify(this%system_charge)
    nullify(this%system_pos)
    SAFE_DEALLOCATE_A(this%partner_pos)

    call interaction_with_partner_end(this)

    POP_SUB(coulomb_force_finalize)
  end subroutine coulomb_force_finalize

end module coulomb_force_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

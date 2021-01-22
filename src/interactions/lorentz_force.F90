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

module lorentz_force_oct_m
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
    lorentz_force_t

  type, extends(force_interaction_t) :: lorentz_force_t
    private
    integer :: dim

    FLOAT, pointer :: system_charge
    FLOAT, pointer, public :: system_pos(:)
    FLOAT, pointer :: system_vel(:)

    FLOAT, allocatable, public :: partner_E_field(:)
    FLOAT, allocatable, public :: partner_B_field(:)

  contains
    procedure :: init => lorentz_force_init
    procedure :: calculate => lorentz_force_calculate
    final :: lorentz_force_finalize
  end type lorentz_force_t

  interface lorentz_force_t
    module procedure lorentz_force_constructor
  end interface lorentz_force_t

contains

  ! ---------------------------------------------------------
  function lorentz_force_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(lorentz_force_t),               pointer       :: this

    PUSH_SUB(lorentz_force_constructor)

    SAFE_ALLOCATE(this)

    this%label = "lorenz_force"

    this%partner => partner

    ! The Lorentz force needs the position, velocity and charge of the system
    this%n_system_quantities = 3
    SAFE_ALLOCATE(this%system_quantities(this%n_system_quantities))
    this%system_quantities(1) = POSITION
    this%system_quantities(2) = VELOCITY
    this%system_quantities(3) = CHARGE
    nullify(this%system_pos)
    nullify(this%system_vel)

    ! The Lorenz force needs the E and B field of the interaction partner at the particle position
    this%n_partner_quantities = 2
    SAFE_ALLOCATE(this%partner_quantities(this%n_partner_quantities))
    this%partner_quantities(1) = E_FIELD
    this%partner_quantities(2) = B_FIELD

    POP_SUB(lorentz_force_constructor)
  end function lorentz_force_constructor

  ! ---------------------------------------------------------
  subroutine lorentz_force_init(this, dim, system_quantities, system_charge, system_pos, system_vel)
    class(lorentz_force_t),               intent(inout) :: this
    integer,                              intent(in)    :: dim
    type(quantity_t),                     intent(inout) :: system_quantities(:)
    FLOAT,                        target, intent(in)    :: system_charge
    FLOAT,                        target, intent(in)    :: system_pos(:)
    FLOAT,                        target, intent(in)    :: system_vel(:)

    PUSH_SUB(lorentz_force_init)

    this%dim = dim
    SAFE_ALLOCATE(this%partner_E_field(dim))
    SAFE_ALLOCATE(this%partner_B_field(dim))
  
    system_quantities(POSITION)%required = .true.
    system_quantities(VELOCITY)%required = .true.
    system_quantities(CHARGE)%required = .true.
    this%system_charge => system_charge
    this%system_pos => system_pos
    this%system_vel => system_vel

    POP_SUB(lorentz_force_init)
  end subroutine lorentz_force_init

  ! ---------------------------------------------------------
  subroutine lorentz_force_calculate(this)
    class(lorentz_force_t),             intent(inout) :: this

    PUSH_SUB(lorentz_force_calculate)

    ! Curl is defined only in 3D
    ASSERT(this%dim == 3)

    this%force(1) = this%system_charge * (this%partner_E_field(1) + &
      this%system_vel(2) * this%partner_B_field(3) - this%system_vel(3) * this%partner_B_field(2))
    this%force(2) = this%system_charge * (this%partner_E_field(2) + &
      this%system_vel(3) * this%partner_B_field(1) - this%system_vel(1) * this%partner_B_field(3))
    this%force(3) = this%system_charge * (this%partner_E_field(3) + &
      this%system_vel(1) * this%partner_B_field(2) - this%system_vel(2) * this%partner_B_field(1))

    POP_SUB(lorentz_force_calculate)
  end subroutine lorentz_force_calculate

  ! ---------------------------------------------------------
  subroutine lorentz_force_finalize(this)
    type(lorentz_force_t), intent(inout) :: this

    PUSH_SUB(lorentz_force_finalize)

    this%force = M_ZERO
    nullify(this%system_charge)
    nullify(this%system_pos)
    nullify(this%system_vel)
    SAFE_DEALLOCATE_A(this%partner_E_field)
    SAFE_DEALLOCATE_A(this%partner_B_field)

    call interaction_with_partner_end(this)

    POP_SUB(lorentz_force_finalize)
  end subroutine lorentz_force_finalize

end module lorentz_force_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

!! Copyright (C) 2021 F. Bonafé
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

module linear_medium_em_field_oct_m
  use clock_oct_m
  use global_oct_m
  use grid_oct_m
  use interaction_with_partner_oct_m
  use interaction_partner_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use quantity_oct_m

  implicit none

  private
  public ::               &
    linear_medium_em_field_t

  type, extends(interaction_with_partner_t) :: linear_medium_em_field_t
    private

    type(grid_t), pointer, public    :: system_gr !< pointer to grid of the Maxwell system

    FLOAT, allocatable, public       :: partner_ep(:) !< permitivity of the linear medium
    FLOAT, allocatable, public       :: partner_mu(:) !< permeability of the linear medium
    FLOAT, allocatable, public       :: partner_c(:) !< speed of light in the linear medium
    FLOAT, allocatable, public       :: partner_sigma_e(:) !< permeability of the linear medium
    FLOAT, allocatable, public       :: partner_sigma_m(:) !< permeability of the linear medium
    integer, public                  :: partner_points_number !< number of points of linear medium
    integer, allocatable, public     :: partner_points_map(:) !< points map of linear medium
    FLOAT, allocatable, public       :: partner_aux_ep(:,:) !< auxiliary array for storing the epsilon derivative profile
    FLOAT, allocatable, public       :: partner_aux_mu(:,:) !< auxiliary array for storing the softened mu profile

    logical, public :: allocated_partner_arrays

  contains
    procedure :: init => linear_medium_em_field_init
    procedure :: calculate => linear_medium_em_field_calculate
    final :: linear_medium_em_field_finalize
  end type linear_medium_em_field_t


  interface linear_medium_em_field_t
    module procedure linear_medium_em_field_constructor
  end interface linear_medium_em_field_t

contains

  ! ---------------------------------------------------------
  function linear_medium_em_field_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(linear_medium_em_field_t),               pointer       :: this

    PUSH_SUB(linear_medium_em_field_constructor)

    SAFE_ALLOCATE(this)

    this%label = "linear_medium_em_field"
    this%partner => partner

    this%n_system_quantities = 0
    nullify(this%system_gr)

    this%n_partner_quantities = 4
    SAFE_ALLOCATE(this%partner_quantities(this%n_partner_quantities))
    this%partner_quantities(1) = PERMITTIVITY
    this%partner_quantities(2) = PERMEABILITY
    this%partner_quantities(3) = E_CONDUCTIVITY
    this%partner_quantities(4) = M_CONDUCTIVITY

    this%allocated_partner_arrays = .false.

    POP_SUB(linear_medium_em_field_constructor)
  end function linear_medium_em_field_constructor


  subroutine linear_medium_em_field_init(this, gr)
    class(linear_medium_em_field_t), intent(inout) :: this
    type(grid_t), target, intent(in)         :: gr

    PUSH_SUB(linear_medium_em_field_init)

    this%system_gr => gr

    POP_SUB(linear_medium_em_field_init)
  end subroutine linear_medium_em_field_init

  ! ---------------------------------------------------------
  subroutine linear_medium_em_field_finalize(this)
    type(linear_medium_em_field_t), intent(inout) :: this

    PUSH_SUB(linear_medium_em_field_finalize)

    POP_SUB(linear_medium_em_field_finalize)

  end subroutine linear_medium_em_field_finalize


  ! ---------------------------------------------------------
  subroutine linear_medium_em_field_calculate(this)
    class(linear_medium_em_field_t), intent(inout) :: this

    PUSH_SUB(linear_medium_em_field_calculate)

    POP_SUB(linear_medium_em_field_calculate)
  end subroutine linear_medium_em_field_calculate


end module linear_medium_em_field_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

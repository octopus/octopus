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
  public ::                    &
    linear_medium_em_field_t,  &
    single_medium_box_t,       &
    single_medium_box_allocate,&
    single_medium_box_end

  type single_medium_box_t
    FLOAT, allocatable            :: ep(:) !< permitivity of the linear media
    FLOAT, allocatable            :: mu(:) !< permeability of the linear media
    FLOAT, allocatable            :: c(:) !< speed of light in the linear media
    FLOAT, allocatable            :: sigma_e(:) !< electric conductivy of (lossy) medium
    FLOAT, allocatable            :: sigma_m(:) !< magnetic conductivy of (lossy) medium
    integer                       :: points_number
    integer                       :: global_points_number
    integer, allocatable          :: points_map(:)
    integer                       :: bdry_number
    integer, allocatable          :: bdry_map(:)
    FLOAT, allocatable            :: aux_ep(:,:) !< auxiliary array for the epsilon derivative profile
    FLOAT, allocatable            :: aux_mu(:,:) !< auxiliary array for the softened mu profile
  end type single_medium_box_t

  type, extends(interaction_with_partner_t) :: linear_medium_em_field_t
    private

    type(grid_t), pointer, public    :: system_gr !< pointer to grid of the Maxwell system

    type(single_medium_box_t), public :: medium_box

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


  ! ---------------------------------------------------------
  !> Allocation of medium_box components
  subroutine single_medium_box_allocate(medium_box, n_points)
    type(single_medium_box_t),   intent(inout)    :: medium_box
    integer,                     intent(in)       :: n_points

    type(profile_t), save :: prof

    PUSH_SUB(medium_box_allocate)

    call profiling_in(prof, 'MEDIUM_BOX_ALLOC')
    SAFE_ALLOCATE(medium_box%aux_ep(n_points,1:3))
    SAFE_ALLOCATE(medium_box%aux_mu(n_points,1:3))
    SAFE_ALLOCATE(medium_box%c(n_points))
    SAFE_ALLOCATE(medium_box%ep(n_points))
    SAFE_ALLOCATE(medium_box%mu(n_points))
    SAFE_ALLOCATE(medium_box%sigma_e(n_points))
    SAFE_ALLOCATE(medium_box%sigma_m(n_points))
    SAFE_ALLOCATE(medium_box%points_map(n_points))
    medium_box%points_map = 0
    call profiling_out(prof)

    POP_SUB(medium_box_allocate)

  end subroutine single_medium_box_allocate

  ! ---------------------------------------------------------
  !> Deallocation of medium_box components
  subroutine single_medium_box_end(medium_box)
    type(single_medium_box_t),   intent(inout)    :: medium_box

    type(profile_t), save :: prof

    PUSH_SUB(medium_box_end)

    call profiling_in(prof, 'MEDIUM_BOX_END')

    SAFE_DEALLOCATE_A(medium_box%points_map)
    SAFE_DEALLOCATE_A(medium_box%bdry_map)
    SAFE_DEALLOCATE_A(medium_box%aux_ep)
    SAFE_DEALLOCATE_A(medium_box%aux_mu)
    SAFE_DEALLOCATE_A(medium_box%c)
    SAFE_DEALLOCATE_A(medium_box%ep)
    SAFE_DEALLOCATE_A(medium_box%mu)
    SAFE_DEALLOCATE_A(medium_box%sigma_e)
    SAFE_DEALLOCATE_A(medium_box%sigma_m)

    call profiling_out(prof)

    POP_SUB(medium_box_end)

  end subroutine single_medium_box_end

end module linear_medium_em_field_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

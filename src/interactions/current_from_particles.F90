!! Copyright (C) 2020 Heiko Appel, Franco Bonafé
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

module current_from_particles_oct_m
  use global_oct_m
  use grid_oct_m
  use interaction_with_partner_oct_m
  use interaction_partner_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use states_mxll_oct_m

  implicit none

  private
  public ::                &
    current_from_particles_t

  type, extends(interaction_with_partner_t) :: current_from_particles_t
    private
    integer :: dim

    type(grid_t), pointer, public :: system_gr    !< the mesh of the Maxwell system

    FLOAT, public :: partner_charge
    FLOAT, public :: partner_charge_smearing
    FLOAT, allocatable, public :: partner_pos(:)
    FLOAT, allocatable, public :: partner_vel(:)

    FLOAT, allocatable :: current_density(:,:)
    CMPLX, allocatable, public :: current_rs_state(:,:)

  contains
    procedure :: init => current_from_particles_init
    procedure :: calculate => current_from_particles_calculate
    final :: current_from_particles_finalize
  end type current_from_particles_t

  interface current_from_particles_t
    module procedure current_from_particles_constructor
  end interface current_from_particles_t

contains

  ! ---------------------------------------------------------
  function current_from_particles_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(current_from_particles_t),      pointer       :: this

    PUSH_SUB(current_from_particles_constructor)

    SAFE_ALLOCATE(this)

    this%label = "classical_current_from_charged_particles"

    this%partner => partner

    ! The current interaction needs no quantities of the Maxwell system
    this%n_system_quantities = 0

    ! The current interaction need the charge, charge smearing, velocity, and position of the charged particle
    this%n_partner_quantities = 4
    SAFE_ALLOCATE(this%partner_quantities(this%n_partner_quantities))
    this%partner_quantities(1) = POSITION
    this%partner_quantities(2) = VELOCITY
    this%partner_quantities(3) = CHARGE
    this%partner_quantities(4) = CHARGE_SMEARING

    POP_SUB(current_from_particles_constructor)
  end function current_from_particles_constructor

  ! ---------------------------------------------------------
  subroutine current_from_particles_init(this, dim, system_gr)
    class(current_from_particles_t), intent(inout) :: this
    integer,                         intent(in)    :: dim
    type(grid_t),            target, intent(in)    :: system_gr

    PUSH_SUB(current_from_particles_init)

    this%dim = dim
    SAFE_ALLOCATE(this%partner_pos(dim))
    SAFE_ALLOCATE(this%partner_vel(dim))

    this%system_gr => system_gr

    SAFE_ALLOCATE(this%current_density(1:this%system_gr%mesh%np, 1:this%system_gr%mesh%sb%dim))
    SAFE_ALLOCATE(this%current_rs_state(1:this%system_gr%mesh%np, 1:this%system_gr%mesh%sb%dim))

    POP_SUB(current_from_particles_init)
  end subroutine current_from_particles_init

  ! ---------------------------------------------------------
  subroutine current_from_particles_calculate(this, namespace)
    class(current_from_particles_t),    intent(inout) :: this
    type(namespace_t),                  intent(in)    :: namespace

    FLOAT :: rr(MAX_DIM)
    integer :: ii, ip

    PUSH_SUB(current_from_particles_calculate)

    ! Calculate classical current: j(r,t) = \sum_j z_j*v_j(t)*\delta(r - r_j(t))
    ! Note, that we use a Gaussian to represent the delta function in the expression above
    do ip = 1, this%system_gr%mesh%np
      rr(1:this%system_gr%mesh%sb%dim) = this%system_gr%mesh%x(ip, 1:this%system_gr%mesh%sb%dim)
      do ii = 1, this%system_gr%mesh%sb%dim
        this%current_density(ip, ii) = this%partner_charge * this%partner_vel(ii)  &
                * exp((rr(ii) - this%partner_pos(ii))**2/this%partner_charge_smearing)
      end do
    end do
    call build_rs_current_state(this%current_density, this%system_gr%mesh, this%current_rs_state)

    POP_SUB(current_from_particles_calculate)
  end subroutine current_from_particles_calculate

  ! ---------------------------------------------------------
  subroutine current_from_particles_finalize(this)
    type(current_from_particles_t), intent(inout) :: this

    PUSH_SUB(current_from_particles_finalize)

    SAFE_DEALLOCATE_A(this%current_density)
    SAFE_DEALLOCATE_A(this%current_rs_state)
    SAFE_DEALLOCATE_A(this%partner_pos)
    SAFE_DEALLOCATE_A(this%partner_vel)

    call interaction_with_partner_end(this)

    POP_SUB(current_from_particles_finalize)
  end subroutine current_from_particles_finalize

end module current_from_particles_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

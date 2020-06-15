!! Copyright (C) 2020 F. Bonafe, R. Jestaedt, H. Appel, N. Tancogne-Dejean
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

module energy_mxll_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::               &
    energy_mxll_t,        &
    energy_mxll_nullify,  &
    energy_mxll_allocate, &
    energy_mxll_end


  type energy_mxll_t
    ! Components are public by default
    ! Energies
    FLOAT                          :: energy
    FLOAT                          :: boundaries
    FLOAT                          :: e_energy
    FLOAT                          :: b_energy
    FLOAT                          :: energy_plane_waves
    FLOAT                          :: e_energy_plane_waves
    FLOAT                          :: b_energy_plane_waves

    FLOAT, allocatable             :: energy_density(:)
    FLOAT, allocatable             :: energy_density_plane_waves(:)
    FLOAT, allocatable             :: e_energy_density(:)
    FLOAT, allocatable             :: b_energy_density(:)
    FLOAT                          :: energy_trans
    FLOAT                          :: energy_long
    FLOAT                          :: e_energy_trans
    FLOAT                          :: b_energy_trans
    FLOAT                          :: energy_incident_waves
  end type energy_mxll_t

contains

  subroutine energy_mxll_nullify(this)
    type(energy_mxll_t), intent(inout) :: this

    PUSH_SUB(energy_mxll_nullify)

    this%energy = M_ZERO
    this%boundaries = M_ZERO
    this%e_energy = M_ZERO
    this%b_energy = M_ZERO
    this%energy_plane_waves = M_ZERO
    this%e_energy_plane_waves = M_ZERO
    this%b_energy_plane_waves = M_ZERO
    this%energy_trans = M_ZERO
    this%energy_long = M_ZERO
    this%e_energy_trans = M_ZERO
    this%b_energy_trans = M_ZERO
    this%energy_incident_waves = M_ZERO

    POP_SUB(energy_mxll_nullify)
  end subroutine energy_mxll_nullify


  subroutine energy_mxll_allocate(this, np)
    type(energy_mxll_t), intent(inout) :: this
    integer,             intent(in)    :: np

    PUSH_SUB(energy_mxll_allocate)

    SAFE_ALLOCATE(this%energy_density(1:np))
    SAFE_ALLOCATE(this%energy_density_plane_waves(1:np))
    SAFE_ALLOCATE(this%e_energy_density(1:np))
    SAFE_ALLOCATE(this%b_energy_density(1:np))

    this%energy_density(:) = M_ZERO
    this%e_energy_density(:) = M_ZERO
    this%b_energy_density(:) = M_ZERO
    this%energy_density_plane_waves(:) = M_ZERO

    POP_SUB(energy_mxll_allocate)
  end subroutine energy_mxll_allocate


  subroutine energy_mxll_end(this)
    type(energy_mxll_t), intent(inout) :: this

    PUSH_SUB(energy_mxll_end)

    SAFE_DEALLOCATE_A(this%energy_density)
    SAFE_DEALLOCATE_A(this%energy_density_plane_waves)
    SAFE_DEALLOCATE_A(this%e_energy_density)
    SAFE_DEALLOCATE_A(this%b_energy_density)

    POP_SUB(energy_mxll_end)
  end subroutine energy_mxll_end

end module energy_mxll_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

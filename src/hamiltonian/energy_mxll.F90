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

  implicit none

  private

  public ::         &
    energy_mxll_t,  &
    energy_mxll_nullify  

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

    FLOAT, pointer                 :: energy_density(:)
    FLOAT, pointer                 :: energy_density_plane_waves(:)
    FLOAT, pointer                 :: e_energy_density(:)
    FLOAT, pointer                 :: b_energy_density(:)
    FLOAT                          :: energy_trans
    FLOAT                          :: energy_long
    FLOAT                          :: e_energy_trans
    FLOAT                          :: b_energy_trans
    FLOAT                          :: energy_incident_waves
  end type energy_mxll_t

contains

  subroutine energy_mxll_nullify(this)
    type(energy_mxll_t), intent(out) :: this

    PUSH_SUB(energy_nullify)

     this%energy = M_ZERO
     this%boundaries = M_ZERO
     this%e_energy = M_ZERO
     this%b_energy = M_ZERO
     this%energy_plane_waves = M_ZERO
     this%e_energy_plane_waves = M_ZERO
     this%b_energy_plane_waves = M_ZERO

     nullify(this%energy_density)
     nullify(this%energy_density_plane_waves)
     nullify(this%e_energy_density)
     nullify(this%b_energy_density)
     this%energy_trans = M_ZERO
     this%energy_long = M_ZERO
     this%e_energy_trans = M_ZERO
     this%b_energy_trans = M_ZERO
     this%energy_incident_waves = M_ZERO

    POP_SUB(energy_nullify)
  end subroutine energy_mxll_nullify

end module energy_mxll_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

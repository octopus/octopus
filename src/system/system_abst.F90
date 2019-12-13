!! Copyright (C) 2019 F. Bonafe
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
!! along with st program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
#include "global.h"

module system_abst_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use namespace_oct_m
  use multicomm_oct_m
  use output_oct_m
  use simul_box_oct_m
  use space_oct_m

  implicit none

  private

  public ::                           &
    system_abst_t

  type :: system_abst_t
    type(space_t)                :: space
    type(geometry_t)             :: geo
    type(grid_t),        pointer :: gr    !< the mesh
    type(output_t)               :: outp  !< the output
    type(multicomm_t)            :: mc    !< index and domain communicators
    type(namespace_t)            :: namespace
  end type system_abst_t

end module system_abst_oct_m

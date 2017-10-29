!! Copyright (C) 2016 N. Tancogne-Dejean 
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
!! $Id$

#include "global.h"

module atomic_orbital_oct_m
  use geometry_oct_m
  use global_oct_m
  use mesh_oct_m
  use messages_oct_m
  use orbital_set_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use species_oct_m
  use species_pot_oct_m
  use submesh_oct_m
 
  implicit none

  private

  public ::                             &
           dget_atomic_orbital,         &
           zget_atomic_orbital

contains

#include "undef.F90"
#include "real.F90"
#include "atomic_orbital_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "atomic_orbital_inc.F90"

end module atomic_orbital_oct_m

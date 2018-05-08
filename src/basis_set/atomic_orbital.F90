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
  use loct_math_oct_m
  use math_oct_m
  use mesh_oct_m
  use messages_oct_m
  use orbitalset_oct_m
  use profiling_oct_m
  use ps_oct_m
  use simul_box_oct_m
  use species_oct_m
  use splines_oct_m
  use submesh_oct_m
 
  implicit none

  private

  public ::                             &
           atomic_orbital_get_radius,   &
           datomic_orbital_get_submesh, &
           zatomic_orbital_get_submesh, &
           dget_atomic_orbital,         &
           zget_atomic_orbital,         &
           l_notation

  character(len=1), parameter :: &
    l_notation(0:3) = (/ 's', 'p', 'd', 'f' /)

  integer, public, parameter :: MAX_L = 4

contains

  ! ---------------------------------------------------------
  FLOAT function atomic_orbital_get_radius(geo, mesh, ia, iorb, ispin, threshold) result(radius)
    type(geometry_t), target, intent(in)   :: geo
    type(mesh_t),             intent(in)   :: mesh
    integer,                  intent(in)   :: ia, iorb, ispin
    FLOAT,                    intent(in)   :: threshold

    type(species_t), pointer :: spec
    integer :: ii, ll, mm

    PUSH_SUB(atomic_orbital_get_radius)

    spec => geo%atom(ia)%species
    call species_iwf_ilm(spec, iorb, ispin, ii, ll, mm)

    radius = species_get_iwf_radius(spec, ii, ispin, threshold)

    POP_SUB(atomic_orbital_get_radius)
  end function atomic_orbital_get_radius


#include "undef.F90"
#include "real.F90"
#include "atomic_orbital_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "atomic_orbital_inc.F90"

end module atomic_orbital_oct_m

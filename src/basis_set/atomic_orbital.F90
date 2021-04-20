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

#include "global.h"

module atomic_orbital_oct_m
  use box_cylinder_oct_m
  use box_minimum_oct_m
  use box_sphere_oct_m
  use global_oct_m
  use ions_oct_m
  use lalg_basic_oct_m
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

  public ::                                  &
           atomic_orbital_get_radius,        &
           datomic_orbital_get_submesh,      &
           zatomic_orbital_get_submesh,      &
           datomic_orbital_get_submesh_safe, &
           zatomic_orbital_get_submesh_safe, &
           dget_atomic_orbital,              &
           zget_atomic_orbital,              &
           l_notation

  character(len=1), parameter :: &
    l_notation(0:3) = (/ 's', 'p', 'd', 'f' /)

  integer, public, parameter :: MAX_L = 4

contains

  ! ---------------------------------------------------------
  FLOAT function atomic_orbital_get_radius(ions, mesh, ia, iorb, ispin, truncation, threshold) result(radius)
    type(ions_t),     target, intent(in)   :: ions
    type(mesh_t),             intent(in)   :: mesh
    integer,                  intent(in)   :: ia, iorb, ispin
    integer(8),               intent(in)   :: truncation
    FLOAT,                    intent(in)   :: threshold

    type(species_t), pointer :: spec
    integer :: ii, ll, mm

    PUSH_SUB(atomic_orbital_get_radius)

    spec => ions%atom(ia)%species
    call species_iwf_ilm(spec, iorb, ispin, ii, ll, mm)

    if(truncation == OPTION__AOTRUNCATION__AO_FULL) then
      radius = species_get_iwf_radius(spec, ii, ispin, threshold)
    else
      radius = species_get_iwf_radius(spec, ii, ispin)

      if(truncation == OPTION__AOTRUNCATION__AO_BOX) then
        ! if the orbital is larger than the size of the box, we restrict it to this size, 
        ! otherwise the orbital will overlap more than one time with the simulation box.
        ! This would induces phase problem if the complete mesh is used instead of the sphere
        radius = min(radius, minval(mesh%sb%lsize(1:mesh%sb%dim)-mesh%spacing(1:mesh%sb%dim)*CNST(1.01)))
      else
        !If asked, we truncate the orbital to the radius on the projector spheres 
        !of the NL part of the pseudopotential.
        !This is a way to garanty no overlap between orbitals of different atoms.
        if(species_is_ps(spec)) &
          radius = min(radius,species_get_ps_radius(spec))
        end if

      end if
      ! make sure that if the spacing is too large, the orbitals fit in a few points at least
      radius = max(radius, CNST(2.0)*maxval(mesh%spacing(1:mesh%sb%dim)))

    POP_SUB(atomic_orbital_get_radius)
  end function atomic_orbital_get_radius


#include "undef.F90"
#include "real.F90"
#include "atomic_orbital_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "atomic_orbital_inc.F90"

end module atomic_orbital_oct_m

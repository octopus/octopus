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

module orbitalset_utils_oct_m
  use atomic_orbital_oct_m
  use geometry_oct_m
  use global_oct_m
  use mesh_oct_m
  use messages_oct_m
  use orbitalset_oct_m
  use species_oct_m
  use submesh_oct_m
 
  implicit none

  private

  public ::                            &
       orbitalset_utils_count,         &
       dorbitalset_utils_getorbitals,  &
       zorbitalset_utils_getorbitals

contains

  integer function orbitalset_utils_count(geo, ia, iselect) result(norb)
    type(geometry_t),     intent(in) :: geo
    integer,              intent(in) :: ia
    integer, optional,    intent(in) :: iselect

    integer :: iorb, ii, ll, mm

    !We count the number of orbital sets we have for a given atom
    !If iselect is present, this routine return instead the number of orbital for a given
    !value of i
    norb = 0
    do iorb = 1, species_niwfs(geo%atom(ia)%species)
      call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm)
      if(present(iselect)) then
        if(ii == iselect) norb = norb + 1
      else
        norb = max(norb, ii)
      end if
    end do
  end function orbitalset_utils_count

#include "undef.F90"
#include "real.F90"
#include "orbitalset_utils_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "orbitalset_utils_inc.F90"

end module orbitalset_utils_oct_m

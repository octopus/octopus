!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.

#include "global.h"

program centergeom
  use global
  use units
  use geometry
  use xyz_adjust

  implicit none

  FLOAT :: dummy
  type(geometry_type) :: geo

  call global_init()                       ! initialize

  call units_init()
  call geometry_init_xyz(geo)              ! we need the geometry
  call geometry_init_species(geo, dummy)   ! we also need the masses

  call xyz_adjust_it(geo)
  call atom_write_xyz(".", "adjusted", geo)

  call geometry_end(geo)                   ! clean up
  call global_end()
end program centergeom

!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!!
!! $Id$

#include "global.h"

program centergeom
  use global_m
  use messages_m
  use datasets_m
  use parser_m
  use io_m
  use units_m
  use geometry_m
  use xyz_adjust_m

  implicit none

  type(geometry_t) :: geo

  call global_init()                       ! initialize
  call parser_init()
  call parse_integer('DebugLevel', 0, conf%debug_level)
  if(conf%debug_level>0) then
    in_debug_mode = .true.
  else
    in_debug_mode = .false.
  end if

  call datasets_init(1)
  call io_init()
  call units_init()

  call geometry_init_xyz(geo)              ! we need the geometry
  call geometry_init_species(geo)          ! we also need the masses

  call xyz_adjust_it(geo)
  call atom_write_xyz(".", "adjusted", geo, 3)

  call geometry_end(geo)                   ! clean up

  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()
end program centergeom

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

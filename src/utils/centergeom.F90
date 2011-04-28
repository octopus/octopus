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
  use command_line_m
  use datasets_m
  use geometry_m
  use global_m
  use io_m
  use messages_m
  use parser_m
  use simul_box_m
  use space_m
  use unit_m
  use unit_system_m
  use xyz_adjust_m

  implicit none

  integer :: ierr
  type(simul_box_t) :: sb
  type(geometry_t)  :: geo
  type(space_t)     :: space

  call global_init()                       ! initialize

  call getopt_init(ierr)
  if(ierr.eq.0) call getopt_center_geom
  call getopt_end()

  call parser_init()
  call messages_init()

  call datasets_init(1)
  call io_init()
  call unit_system_init()

  call space_init(space)
  call geometry_init(geo, space)
  call simul_box_init(sb, geo, space)

  call xyz_adjust_it(geo)
  call atom_write_xyz(".", "adjusted", geo, 3)

  call simul_box_end(sb)
  call geometry_end(geo)
  call space_end(space)

  call io_end()
  call datasets_end()
  call messages_end()
  call parser_end()
  call global_end()
end program centergeom

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

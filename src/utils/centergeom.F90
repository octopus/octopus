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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

program centergeom
  use command_line_oct_m
  use geometry_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use simul_box_oct_m
  use space_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use xyz_adjust_oct_m

  implicit none

  integer :: ierr
  type(simul_box_t) :: sb
  type(geometry_t)  :: geo
  type(space_t)     :: space
  type(namespace_t) :: default_namespace

  call global_init(is_serial = .true.)

  call getopt_init(ierr)
  if(ierr == 0) call getopt_center_geom()
  call getopt_end()

  call parser_init()
  default_namespace = namespace_t("")

  call messages_init(default_namespace)

  call io_init(default_namespace)
  call unit_system_init(default_namespace)

  call space_init(space, default_namespace)
  call geometry_init(geo, default_namespace, space)
  call simul_box_init(sb, default_namespace, geo, space)

  call xyz_adjust_it(geo, default_namespace)
  call geometry_write_xyz(geo, './adjusted', default_namespace)

  call simul_box_end(sb)
  call geometry_end(geo)
  call space_end(space)

  call io_end()
  call messages_end()

  call parser_end()

  call global_end()
end program centergeom

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

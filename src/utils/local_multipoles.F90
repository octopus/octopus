!! Copyright (C) 2014 M. Oliveira
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
!! $Id: dielectric_function.F90 11591 2013-12-17 09:10:32Z joseba $

#include "global.h"

program local_multipoles
  use batch_m
  use datasets_m
  use command_line_m
  use geometry_m
  use global_m
  use io_m
  use lalg_adv_m
  use loct_m
  use messages_m
  use parser_m
  use profiling_m
  use space_m
  use simul_box_m
  use unit_m
  use unit_system_m
  use utils_m

  implicit none

  integer :: ierr
  character*256 :: config_str
  type(space_t)     :: space
  type(geometry_t)  :: geo
  type(simul_box_t) :: sb

  ! Initialize stuff
  call global_init(is_serial = .true.)

  config_str = trim(get_config_opts()) // trim(get_optional_libraries())
  if(ierr  ==  0) call getopt_octopus(config_str)
  call getopt_end()

  call parser_init()
  call messages_init()

  call parse_logical('ExperimentalFeatures', .false., conf%devel_version)

  call datasets_init(1)
  call io_init()

  call unit_system_init()

  call space_init(space)
  call geometry_init(geo, space)
  call simul_box_init(sb, geo, space)

  


  call simul_box_end(sb)
  call geometry_end(geo)
  call space_end(space)
  call io_end()
  call datasets_end()
  call messages_end()
  call parser_end()
  call global_end()

end program local_multipoles

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

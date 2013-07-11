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
!! $Id$

#include "global.h"

program xyzanim
  use command_line_m
  use global_m
  use messages_m
  use datasets_m
  use loct_m
  use parser_m
  use io_m
  use unit_m
  use unit_system_m
  use simul_box_m
  use space_m
  use geometry_m

  implicit none

  character(len=256) :: coords_file, comment
  integer :: ierr, sampling, i, coords_unit, iter, j, record_length
  FLOAT :: time
  type(simul_box_t) :: sb
  type(geometry_t)  :: geo
  type(space_t)     :: space

  ! Initialize stuff
  call global_init()

  call getopt_init(ierr)
  if(ierr == 0) call getopt_xyz_anim
  call getopt_end()

  call parser_init()
  call messages_init()
  call datasets_init(1)
  call io_init()
  call unit_system_init()

  ! Sets the filenames
  coords_file = 'td.general/coordinates'

  !%Variable AnimationSampling
  !%Type integer
  !%Default 100
  !%Section Utilities::oct-xyz-anim
  !%Description
  !% Sampling rate of the animation. The animation will be constructed using
  !% the iteration numbers that are multiples of <tt>AnimationSampling<tt>.
  !%End
  call parse_integer(datasets_check('AnimationSampling'), 100, sampling)
  if(sampling < 1) then
    message(1) = 'Sampling rate (AnimationSampling) should be bigger than 0'
    call messages_fatal(1)
  end if

  call space_init(space)
  call geometry_init(geo, space)
  call simul_box_init(sb, geo, space)

  record_length = 100 + geo%space%dim*geo%natoms*3*20

  ! Opens the coordinates file
  coords_unit = io_open(coords_file, action='read', recl = record_length)

  call io_skip_header(coords_unit)
  ierr = 0
  do while(ierr == 0)
    read(unit = coords_unit, iostat = ierr, fmt = *) iter, time, &
      ((geo%atom(i)%x(j), j = 1, geo%space%dim), i = 1, geo%natoms)
    if(mod(iter, sampling) == 0) then
      write(comment, '(i10,f20.6)') iter, time
      call io_mkdir('td.general')
      call geometry_write_xyz(geo, 'td.general/movie', append = .true., comment = trim(comment))
    end if
  end do

  call simul_box_end(sb)
  call geometry_end(geo)
  call space_end(space)

  call io_close(coords_unit)

  call io_end()
  call datasets_end()
  call messages_end()
  call parser_end()
  call global_end()

end program xyzanim

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

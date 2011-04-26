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

  character(len=256) :: coords_file, xyzfile
  integer :: ierr, sampling, i, coords_unit, xyz_unit, iter, j, record_length
  FLOAT :: time
  type(simul_box_t) :: sb
  type(geometry_t)  :: geo
  type(space_t)     :: space

  ! Initialize stuff
  call global_init()
  call getopt_init(ierr)
  if(ierr.eq.0) call getopt_xyz_anim
  call parser_init()
  call messages_init()
  call datasets_init(1)
  call io_init()
  call unit_system_init()

  ! Sets the filenames
  coords_file = 'td.general/coordinates'
  xyzfile = 'td.general/movie.xyz'

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

  record_length = 100 + 3*geo%natoms*3*20

  ! Opens the coordinates file
  coords_unit = io_open(coords_file, action='read', recl = record_length)

  ! Opens the xyz file
  xyz_unit = io_open(xyzfile, action='write')

  call io_skip_header(coords_unit)
  ierr = 0
  do while(ierr == 0)
    read(unit = coords_unit, iostat = ierr, fmt = *) iter, time, &
      ((geo%atom(i)%x(j), j = 1, 3), i = 1, geo%natoms)
    if(mod(iter, sampling) == 0) then
      call write_xyz()
    end if
  end do

  call simul_box_end(sb)
  call geometry_end(geo)
  call space_end(space)

  call io_close(coords_unit); call io_close(xyz_unit)

  call io_end()
  call datasets_end()
  call messages_end()
  call parser_end()
  call global_end()

contains

  ! ---------------------------------------------------------
  subroutine write_xyz
    integer :: i
    ! xyz format
    write(xyz_unit, '(i4)') geo%natoms
    write(xyz_unit, '(i10,f20.6)') iter, time
    do i = 1, geo%natoms
      write(xyz_unit, '(6x,a,2x,3f12.6)') geo%atom(i)%label, geo%atom(i)%x(:)
    end do
  end subroutine write_xyz

end program xyzanim

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

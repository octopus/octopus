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
!!
!! $Id$

#include "global.h"

program xyzanim
  use global
  use messages
  use syslabels
  use lib_oct
  use lib_oct_parser
  use io
  use units
  use geometry

  implicit none

  character(len=256) :: nbofile, xyzfile
  integer :: ierr, sampling, i, nbo_unit, xyz_unit, iter, j
  FLOAT :: dump

  type(geometry_type) :: geo

  ! Initialize stuff
  call global_init()
  call parser_init()
  call io_init()
  call syslabels_init(1)
  if(in_debug_mode) then
     call io_mkdir('debug')
  endif
  call units_init()

  ! Sets the filenames
  nbofile = 'td.general/coordinates'
  xyzfile = 'td.general/movie.xyz'

  ! how often do we sample?
  call loct_parse_int(check_inp('AnimationSampling'), 100, sampling)
  if(sampling < 1) then
    message(1) = 'Sampling rate (AnimationSampling) should be bigger than 0'
    call write_fatal(1)
  end if

  ! Initializes the atom system
  call geometry_init_xyz(geo)

  ! Opens the nbo file
  nbo_unit = io_open(nbofile, action='read')

  ! Opens the xyz file
  xyz_unit = io_open(xyzfile, action='write')

  ! Skips the header
  rewind(unit = nbo_unit)
  read(unit = nbo_unit, fmt = *); read(unit = nbo_unit, fmt = *)

  ierr = 0
  do while(ierr == 0)
     read(unit = nbo_unit, iostat = ierr, fmt = *) iter, dump, dump, dump, dump, &
         ((geo%atom(i)%x(j), j = 1, 3), i = 1, geo%natoms)
     if(mod(iter, sampling) == 0) then
       call write_xyz()
     endif
  enddo

  call io_close(nbo_unit); call io_close(xyz_unit)

  call syslabels_end()
  call io_end()
  call parser_end()
  call global_end()

contains

  subroutine write_xyz
    integer :: i
    ! xyz format
    write(xyz_unit, '(i4)') geo%natoms
    write(xyz_unit, '(i10)') iter
    do i = 1, geo%natoms
      write(xyz_unit, '(6x,a,2x,3f12.6)') geo%atom(i)%label, geo%atom(i)%x(:)
    end do
  end subroutine write_xyz

end program xyzanim

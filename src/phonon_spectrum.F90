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

program phonon_spectrum
  use global
  use messages
  use syslabels
  use io
  use lib_oct_parser
  use units

  implicit none

  integer :: iunit

  ! Initialize stuff
  call global_init()
  call parser_init()
  call io_init()
  call syslabels_init(1)
  if(in_debug_mode) then
     call io_mkdir('debug')
  endif
  call units_init()

  ! Opens the coordinates files.
  iunit = io_open('td.general/coordinates', action='read', status='old')

  ! Closes file and exits
  call io_close(iunit)

  call syslabels_end()
  call io_end()
  call parser_end()
  call global_end()

end program

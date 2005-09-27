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

program rotational_strength
  use global
  use messages
  use syslabels
  use lib_oct_parser
  use io
  use units
  use spectrum

  implicit none

  character(len=100) :: txt
  type(spec_type) :: s
  integer :: i, in_file, out_file
  integer(POINTER_SIZE) :: blk

  ! Initialize stuff
  call global_init()
  call parser_init()
  call io_init()
  call syslabels_init(1)
  if(in_debug_mode) then
     call io_mkdir('debug')
  endif
  call units_init()

  call spectrum_init(s)

  in_file = io_open('angular', action='read', status='old', die=.false.)
  if(in_file < 0) in_file = io_open('td.general/angular', action='read', status='old', die=.false.)
  if(in_file < 0) then
     write(message(1),'(a)') 'No "angular" or "td.general/angular" file found. At least one of those'
     write(message(2),'(a)') 'should be visible.'
     call write_fatal(2)
  endif
  out_file = io_open('rotatory_strength', action='write')

  call spectrum_rotatory_strength(in_file, out_file, s)
  call io_close(in_file); call io_close(out_file)

  call syslabels_end()
  call io_end()
  call parser_end()
  call global_end()
end program rotational_strength

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

program rotational_strength
  use command_line_m
  use datasets_m
  use global_m
  use io_m
  use messages_m
  use parser_m
  use spectrum_m
  use unit_m
  use unit_system_m

  implicit none

  type(spec_t) :: spectrum
  integer :: in_file, out_file, ierr

  ! Initialize stuff
  call global_init()

  call getopt_init(ierr)
  if(ierr.eq.0) call getopt_rotatory_strength
  call getopt_end()

  call parser_init()
  call messages_init()
  call datasets_init(1)
  call io_init()
  call unit_system_init()

  call spectrum_init(spectrum)

  in_file = io_open('angular', action='read', status='old', die=.false.)
  if(in_file < 0) in_file = io_open('td.general/angular', action='read', status='old', die=.false.)
  if(in_file < 0) then
    write(message(1),'(a)') 'No "angular" or "'//trim(io_workpath('td.general/angular'))//'" file found. At least one of those'
    write(message(2),'(a)') 'should be visible.'
    call messages_fatal(2)
  end if
  out_file = io_open('rotatory_strength', action='write')

  call spectrum_rotatory_strength(in_file, out_file, spectrum)
  call io_close(in_file); call io_close(out_file)

  call io_end()
  call datasets_end()
  call messages_end()
  call parser_end()
  call global_end()
end program rotational_strength

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

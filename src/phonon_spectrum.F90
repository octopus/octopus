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

#include "global.h"

program phonon_spectrum
  use global
  use io
  use units

  implicit none

  integer :: ierr, iunit, time_steps

  ! Initialize stuff
  call global_init()
  call units_init()

  ! Opens the coordinates files.
  call io_assign(iunit)
  open(iunit, file='td.general/coordinates', status='old', action='read', iostat=ierr)
  if(ierr .ne. 0) then
    message(1) = 'Error opening file td.general/coordinates'
    call write_fatal(1)
  end if

  ! Counts the number of lines in the file ( 2 is the number of lines in the header)
  time_steps = number_of_lines('td.general/coordinates') - 2

  ! Closes file and exits
  call io_close(iunit)

end program

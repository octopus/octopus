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

#include "config_F90.h"

program octopus
  use global
  use liboct
  use run_prog

  implicit none

  integer :: ierr, val(8)

#ifdef HAVE_MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpiv%node, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiv%numprocs, ierr)
  write(stdout,'(a,i4,a,i4,a)') 'Process ', mpiv%node, ' of ', mpiv%numprocs, ' is alive'  
#else
  mpiv%node = 0
  mpiv%numprocs = 1
#endif

  ! init some of the stuff
  ierr = oct_parse_init(C_string('inp'), C_string('out.oct'))
  if(ierr .ne. 0) then
    ierr = oct_parse_init(C_string("-"), C_string('out.oct'))
    if(ierr .ne. 0) then
      message(1) = "Error initializing liboct"
      call write_fatal(1)
    end if
  end if
  
  call oct_parse_int(C_string('verbose'), 30, conf%verbose)
  if(conf%verbose > 999 .and. mpiv%node == 0) then
    message(1) = 'Entering DEBUG mode'
    call write_warning(1)
  end if
  
  ! Let us print our logo
  if(conf%verbose > 20 .and. mpiv%node == 0) &
       ierr = print_file(C_string(SHARE_OCTOPUS//'/logo'))

  ! print date
  call date_and_time(values=val)
  write(message(1),'(a,i4,a1,i2.2,a1,i2.2,a,i2.2,a1,i2.2,a1,i2.2)') &
       "Info: Calculation started on ", val(1), "/", val(2), "/", val(3), &
       " at ", val(5), ":", val(6), ":", val(7)
  call write_info(1)

  ! Sets the dimensionaliy of the problem.
  call oct_parse_int(C_string('Dimensions'), 3, conf%dim)
  if(conf%dim<1 .or. conf%dim>3) then
    message(1) = 'Dimensions must be either 1, 2, or 3'
    call write_fatal(1)
  end if
  write(message(1), '(a,i1,a)') 'Octupus will run in ', conf%dim, ' dimension(s)'

  ! create temporary dir (is always necessary)
  call oct_mkdir(C_string("tmp"))

  ! now we really start
  call run()

  ! print date
  call date_and_time(values=val)
  write(message(1),'(a,i4,a1,i2.2,a1,i2.2,a,i2.2,a1,i2.2,a1,i2.2)') &
       "Info: Calculation ended on ", val(1), "/", val(2), "/", val(3), &
       " at ", val(5), ":", val(6), ":", val(7)
  call write_info(1)

#ifdef HAVE_MPI
  call MPI_FINALIZE(ierr)
#endif
  
  call oct_parse_end()
  stop
end program octopus

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

program octopus
  use global
  use lib_oct
  use run_prog

  implicit none

  integer :: ierr, val(8)
  
  call global_init()

  ! Let us print our logo
  if(mpiv%node == 0) then
    ! Let us print our logo
    if(conf%verbose > 20) ierr = loct_print_file(SHARE_OCTOPUS+'/logo')
#ifdef DEBUG
    if(conf%verbose > 999) write(stderr, '(5a)') "# ", " A ", "Time", "Mem", "Call"
#endif
  end if

  ! Let us print the version
  message(1) = ""
  message(2) = str_center(trim("Running octopus, version " + OCTOPUS_VERSION), 70)
  message(3) = str_center(trim("(build time - " + BUILD_TIME + ")"), 70)
  message(4) = ""
  call write_info(4)

  ! Let us print where we are running
  if(conf%verbose > 20) then
    call loct_sysname(message(1))
    write(stdout, '(a)') str_center("The octopus is swimming in "+trim(message(1)), 70)
  end if
#if defined(HAVE_MPI)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

  ! print date
  call date_and_time(values=val)
  write(message(1),'(a,i4,a1,i2.2,a1,i2.2,a,i2.2,a1,i2.2,a1,i2.2)') &
       "Info: Calculation started on ", val(1), "/", val(2), "/", val(3), &
       " at ", val(5), ":", val(6), ":", val(7)
  call write_info(1)

  ! create temporary dir (we will need it)
  call loct_mkdir("tmp")

  ! now we really start
  call run()

  ! print date
  call date_and_time(values=val)
  write(message(1),'(a,i4,a1,i2.2,a1,i2.2,a,i2.2,a1,i2.2,a1,i2.2)') &
       "Info: Calculation ended on ", val(1), "/", val(2), "/", val(3), &
       " at ", val(5), ":", val(6), ":", val(7)
  call write_info(1)

#if defined(HAVE_MPI)
  ! wait for all processors to finish
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

  call global_end()

  stop
end program octopus

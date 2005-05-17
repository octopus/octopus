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

program octopus
  use string
  use global
  use lib_oct
  use run_prog

  implicit none

  integer :: ierr, ns, val(8)

  call global_init()

  subsystems: do ns = 1,no_syslabels

     ! set system label
     current_label = trim(subsys_label(subsys_run_order(ns)))
     current_subsystem = subsys_run_order(ns)
     tmpdir = trim(current_label)//'tmp/'
     calc_mode = subsys_runmode(subsys_run_order(ns))
     
     ! Let us print our logo
     if(mpiv%node == 0) then
        ! Let us print our logo
        if(conf%verbose >= VERBOSE_NORMAL) call io_dump_file(stdout, trim(trim(conf%share) // '/logo'))
#ifdef DEBUG
        if(conf%verbose >= VERBOSE_DEBUG) write(stderr, '(5a)') "# ", " A ", "Time", "Mem", "Call"
#endif
     end if
     
     ! Let us print the version
     message(1) = ""
     message(2) = str_center("Running octopus, version " // trim(conf%version), 70)
     message(3) = str_center("(build time - " // trim(conf%build_time) // ")", 70)
     message(4) = ""
     call write_info(4)
     
     ! Let us print where we are running
     call loct_sysname(message(1))
     write(message(1), '(a)') str_center("The octopus is swimming in " // trim(message(1)), 70)
     message(2) = ""
     call write_info(2)

#if defined(HAVE_MPI)
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
     
     ! print date
     call date_and_time(values=val)
     write(message(2),'(a,i4,a1,i2.2,a1,i2.2,a,i2.2,a1,i2.2,a1,i2.2)') &
          "Calculation started on ", val(1), "/", val(2), "/", val(3), &
          " at ", val(5), ":", val(6), ":", val(7)
     message(1) = str_center(trim(message(2)), 70)
     message(2) = ""
     call write_info(2)

     if(no_syslabels > 1) then
        message(1) = 'Info: Multi-Subsystem Mode'
        message(2) = 'Info: Running '//current_label
        call write_info(2, stress = .true.)
     endif

     ! create temporary dir (we will need it)
     call io_mkdir(tmpdir)
     
     ! create debug directory if in debugging mode
     if(conf%verbose>=VERBOSE_DEBUG) call io_mkdir(trim(current_label)//'debug')
     
     ! now we really start
     call run_init()
     call run()
     call run_end()

     ! print date
     call date_and_time(values=val)
     message(1) = ""
     write(message(3),'(a,i4,a1,i2.2,a1,i2.2,a,i2.2,a1,i2.2,a1,i2.2)') &
          "Calculation ended on ", val(1), "/", val(2), "/", val(3), &
          " at ", val(5), ":", val(6), ":", val(7)
     message(2) = str_center(trim(message(3)), 70)
     message(3) = ""
     call write_info(3)
     
#if defined(HAVE_MPI)
     ! wait for all processors to finish
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
     
  enddo subsystems

  call global_end()

  stop
end program octopus

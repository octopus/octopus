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
  
  call global_init()

  ! Let us print our logo
  if(conf%verbose > 20 .and. mpiv%node == 0) &
       ierr = print_file(SHARE_OCTOPUS//'/logo')

  ! Let us print the version
  message(1) = ""
  message(2) = str_center(trim("Running octopus, version " + OCTOPUS_VERSION), 70)
  message(3) = str_center(trim("(build time - " + BUILD_TIME + ")"), 70)
  message(4) = ""
  call write_info(4)

  ! print date
  call date_and_time(values=val)
  write(message(1),'(a,i4,a1,i2.2,a1,i2.2,a,i2.2,a1,i2.2,a1,i2.2)') &
       "Info: Calculation started on ", val(1), "/", val(2), "/", val(3), &
       " at ", val(5), ":", val(6), ":", val(7)
  call write_info(1)

  ! create temporary dir (we will need it)
  call oct_mkdir("tmp")

  ! now we really start
  call run()

  ! print date
  call date_and_time(values=val)
  write(message(1),'(a,i4,a1,i2.2,a1,i2.2,a,i2.2,a1,i2.2,a1,i2.2)') &
       "Info: Calculation ended on ", val(1), "/", val(2), "/", val(3), &
       " at ", val(5), ":", val(6), ":", val(7)
  call write_info(1)

  call global_end()

  stop
end program octopus

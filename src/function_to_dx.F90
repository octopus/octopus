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

program strength_function
  use global
  use liboct
  use states
  use system
  use mix
  use dx

  implicit none

  character(len=20) :: sysname
  character(len=30) :: filename
  integer :: ierr, i, iunit_data, iunit_header
  type(system_type) :: sys

  ! init liboct
  ierr = oct_parse_init(C_string('inp'), C_string('out.oct'))
  if(ierr .ne. 0) then
    message(1) = "Error initializing liboct"
    call write_fatal(1)
  end if
  call oct_parse_int(C_string("verbose"), 30, conf%verbose)

  conf%verbose = 0
  conf%dim     = 3

  call units_init()

  call system_init(sys)

  allocate(sys%st%dpsi(0:sys%m%np, sys%st%dim, sys%st%nst, sys%st%nik))

  call oct_parse_str('SystemName', 'system', sysname)
  filename = trim(sysname)//'.restart'

  if(dstates_load_restart(trim(sys%sysname)//".restart", sys%m, sys%st)) then
     call R_FUNC(calcdens)(sys%st, sys%m%np, sys%st%rho)
  else
     message(1) = "Error opening restart file"
     call write_fatal(1)
  endif

  call io_assign(iunit_data)
  call io_assign(iunit_header)
  open(unit=iunit_data, file=trim(sys%sysname)//".dxdata", form='formatted')
  open(unit=iunit_header, file=trim(sys%sysname)//".general", form='formatted')
  call create_dx_file(iunit_data, trim(sys%sysname)//".dxdata", &
                      iunit_header, sys%m, sys%st%rho(1:sys%m%np, 1) )
  call io_close(iunit_data)
  call io_close(iunit_header)

  stop  
end program strength_function

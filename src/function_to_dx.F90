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

program function_to_dx
#ifdef THREE_D
  use global
  use liboct
  use states
  use system
  use mix
  use dx

  implicit none

  type orbital_type
    integer :: n, k, spin
  end type

  character(len=20) :: sysname
  character(len=30) :: filename, str
  character(len=5)  :: nstr1, nstr2, nstr3
  integer :: ierr, i, j, iunit_data, iunit_header, norbitals, unocc_states
  real(r8), pointer :: f(:)
  type(system_type) :: sys
  type(orbital_type), allocatable :: orbital(:)

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

  call oct_parse_int(C_string("UnoccNumberStates"), 0, unocc_states)
  sys%st%nst = sys%st%nst + unocc_states
  sys%st%st_start = 1
  sys%st%st_end   = sys%st%nst

  allocate(sys%st%dpsi(0:sys%m%np, sys%st%dim, sys%st%nst, sys%st%nik))

  call oct_parse_str('SystemName', 'system', sysname)
  filename = trim(sysname)+'.occ_restart'

  if(dstates_load_restart(trim(sys%sysname)+".occ_restart", sys%m, sys%st)) then
     call R_FUNC(calcdens)(sys%st, sys%m%np, sys%st%rho)
  else
     message(1) = "Error opening restart file"
     call write_fatal(1)
  endif

  str = C_string('DXOrbitalToPlot')
  norbitals = oct_parse_block_n(trim(str))
  allocate(orbital(norbitals))
  if(norbitals > 0) then
    do j = 1, norbitals
       call oct_parse_block_int(str, j-1, 0, orbital(j)%n)
       call oct_parse_block_int(str, j-1, 1, orbital(j)%k)
       call oct_parse_block_int(str, j-1, 2, orbital(j)%spin)
       call io_assign(iunit_data); call io_assign(iunit_header)
       write(nstr1,*) orbital(j)%n;    nstr1 = adjustl(nstr1)
       write(nstr2,*) orbital(j)%k;    nstr2 = adjustl(nstr2)
       write(nstr3,*) orbital(j)%spin; nstr3 = adjustl(nstr3)
       filename = trim(sys%sysname)+"."+trim(nstr1)+"."+trim(nstr2)+"."+trim(nstr3)+".general"
       open(unit=iunit_header, file=filename, form='formatted')
       filename = trim(sys%sysname)+"."+trim(nstr1)+"."+trim(nstr2)+"."+trim(nstr3)+".dxdata"
       open(unit=iunit_data, file=filename, form='formatted')
       call create_dx_file(iunit_data, filename, iunit_header, sys%m, &
                           (sys%st%dpsi(1:sys%m%np, orbital(j)%spin, orbital(j)%n, orbital(j)%k))**2 )
       call io_close(iunit_data) ;call io_close(iunit_header)
       enddo
  else
    call io_assign(iunit_data); call io_assign(iunit_header)
    open(unit=iunit_data, file=trim(sys%sysname)+".dxdata", form='formatted')
    open(unit=iunit_header, file=trim(sys%sysname)+".general", form='formatted')
    call create_dx_file(iunit_data, trim(sys%sysname)+".dxdata", &
                      iunit_header, sys%m, sys%st%rho(1:sys%m%np, 1) )
    call io_close(iunit_data) ;call io_close(iunit_header)
  endif

  stop
#endif
end program function_to_dx

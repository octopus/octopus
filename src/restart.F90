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

module restart
use lib_oct_parser
use global
use io
use states
use mesh
use output

implicit none

integer, parameter :: RESTART_PLAIN  = 1, &
                      RESTART_NETCDF = 2
integer :: restart_format

contains

subroutine restart_init
  integer :: i
  
  ! read restart format information
  call loct_parse_int('RestartFileFormat', RESTART_PLAIN, i)
  if (i<RESTART_PLAIN .or. i>RESTART_NETCDF) then
    write(message(1),'(a,i4,a)') "Input: '", i,"' is not a valid RestartFileFormat"
    message(2) = '(RestartFileFormat = plain | netcdf)'
    call write_fatal(2)
  end if

  ! Fix the restart format...
  restart_format = output_fill_how("Plain")
#if defined(HAVE_NETCDF)
  if(i == RESTART_NETCDF) then
    restart_format = output_fill_how("NETCDF")
  end if
#endif

end subroutine restart_init


#include "undef.F90"
#include "real.F90"
#include "restart_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "restart_inc.F90"

end module restart

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

module wave_matching
  use global
  use messages
  use syslabels
  use system
  use hamiltonian
  use io

  implicit none

  private
  public :: &
    wave_matching_run


contains

  ! ---------------------------------------------------------
  integer function wave_matching_run() result(ierr)
!!$    type(system_type)       :: sys
!!$    type(hamiltonian_type)  :: h

    call push_sub('wave_matching.wave_matching_run')

    ierr = 0

    call check_params()

    call pop_sub()


  contains

    ! ---------------------------------------------------------
    subroutine check_params()

      if( current_subsystem-1 .lt. 1 ) then
        message(1) = 'Error: Missing left neighbor'
        message(2) = 'Please correct your input file'
        call write_fatal(2)
      end if

      if( current_subsystem+1 .gt. no_syslabels ) then
        message(1) = 'Error: Missing right neighbor'
        message(2) = 'Please correct your input file'
        call write_fatal(2)
      end if

      if( subsys_runmode(current_subsystem-1) .eq. 10 ) then
        message(1) = 'Error: Left neighbor cannot be in runmode wave_matching'
        message(2) = 'Please correct your input file'
        call write_fatal(2)
      end if

      if( subsys_runmode(current_subsystem+1) .eq. 10 ) then
        message(1) = 'Error: Right neighbor cannot be in runmode wave_matching'
        message(2) = 'Please correct your input file'
        call write_fatal(2)
      end if

      message(1) = 'Info: Starting Wave-Matching'
      message(2) = 'Info: We are                    : '//subsys_label(current_subsystem)
      message(3) = 'Info: Our left neighbor is      : '//subsys_label(current_subsystem-1)
      message(4) = 'Info: And our right neighbor is : '//subsys_label(current_subsystem+1)
      call write_info(4, stress=.true.)

    end subroutine check_params

  end function wave_matching_run


end module wave_matching

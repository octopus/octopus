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

module wave_matching_m
  use global_m
  use messages_m
  use datasets_m
  use system_m
  use hamiltonian_m
  use io_m

  implicit none

  private
  public :: &
    wave_matching_run


contains

  ! ---------------------------------------------------------
  subroutine wave_matching_run()
!!$    type(system_t)       :: sys
!!$    type(hamiltonian_t)  :: h

    call push_sub('wave_matching.wave_matching_run')

    call check_params()

    call pop_sub()


  contains

    ! ---------------------------------------------------------
    subroutine check_params()

      if( current_dataset-1 .lt. 1 ) then
        message(1) = 'Error: Missing left neighbor'
        message(2) = 'Please correct your input file'
        call write_fatal(2)
      end if

      if( current_dataset+1 .gt. no_datasets ) then
        message(1) = 'Error: Missing right neighbor'
        message(2) = 'Please correct your input file'
        call write_fatal(2)
      end if

      if( dataset_runmode(current_dataset-1) .eq. 10 ) then
        message(1) = 'Error: Left neighbor cannot be in runmode wave_matching'
        message(2) = 'Please correct your input file'
        call write_fatal(2)
      end if

      if( dataset_runmode(current_dataset+1) .eq. 10 ) then
        message(1) = 'Error: Right neighbor cannot be in runmode wave_matching'
        message(2) = 'Please correct your input file'
        call write_fatal(2)
      end if

      message(1) = 'Info: Starting Wave-Matching'
      message(2) = 'Info: We are                    : '//dataset_label(current_dataset)
      message(3) = 'Info: Our left neighbor is      : '//dataset_label(current_dataset-1)
      message(4) = 'Info: And our right neighbor is : '//dataset_label(current_dataset+1)
      call write_info(4, stress=.true.)

    end subroutine check_params

  end subroutine wave_matching_run


end module wave_matching_m

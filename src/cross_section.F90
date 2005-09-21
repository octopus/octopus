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

program cross_section
  use global
  use messages
  use syslabels
  use lib_oct_parser
  use io
  use units
  use spectrum

  implicit none

  integer :: in_file(3), out_file(3), eq_axis, i, j, nspin, lmax, time_steps
  FLOAT :: dt
  integer(POINTER_SIZE) :: blk
  logical :: calculate_tensor
  type(spec_type) :: s
  type(kick_type) :: kick

  ! Initialize stuff
  call global_init()
  call parser_init()
  call io_init()
  call syslabels_init(1)
  current_label = trim(subsys_label(subsys_run_order(1)))
  call units_init()

  call spectrum_init(s)
  call read_files()

  if(.not.calculate_tensor) then

     out_file(1) = io_open('cross_section_vector', action='write')
     call spectrum_cross_section(in_file(1), out_file(1), s)
     call io_close(in_file(1)); call io_close(out_file(1))

  else

      select case(eq_axis)

        case(0, 1)

          out_file(1) = io_open('cross_section_vector.1', action='write')
          out_file(2) = io_open('cross_section_vector.2', action='write')
          out_file(3) = io_open('cross_section_vector.2', action='write')

          call spectrum_cross_section(in_file(1), out_file(1), s)
          call io_close(in_file(1)); call io_close(out_file(1))
          call spectrum_cross_section(in_file(2), out_file(2), s)
          call io_close(in_file(2)); call io_close(out_file(2))
          call spectrum_cross_section(in_file(3), out_file(3), s)
          call io_close(in_file(3)); call io_close(out_file(3))

          ! And now we build the cross section tensor
          in_file(1)  = io_open('cross_section_vector.1', action='read', status='old')
          in_file(2)  = io_open('cross_section_vector.2', action='read', status='old')
          in_file(3)  = io_open('cross_section_vector.2', action='read', status='old')
          out_file(1) = io_open('cross_section_tensor', action='write')

          call spectrum_cross_section_tensor(s, out_file(1), in_file(1), in_file(2), in_file(3))
          call io_close(in_file(1)); call io_close(in_file(2)); call io_close(in_file(3)); call io_close(out_file(1))    

        case(2)

          out_file(1) = io_open('cross_section_vector.1', action='write')
          out_file(2) = io_open('cross_section_vector.2', action='write')

          call spectrum_cross_section(in_file(1), out_file(1), s)
          call io_close(in_file(1)); call io_close(out_file(1))
          call spectrum_cross_section(in_file(2), out_file(2), s)
          call io_close(in_file(2)); call io_close(out_file(2))

          ! And now we build the cross section tensor
          in_file(1)  = io_open('cross_section_vector.1', action='read', status='old')
          in_file(2)  = io_open('cross_section_vector.2', action='read', status='old')
          out_file(1) = io_open('cross_section_tensor', action='write')

          call spectrum_cross_section_tensor(s, out_file(1), in_file(1), in_file(2))
          call io_close(in_file(1)); call io_close(in_file(2)); call io_close(out_file(1))

        case(3)

          out_file(1) = io_open('cross_section_vector', action='write')

          call spectrum_cross_section(in_file(1), out_file(1), s)
          call io_close(in_file(1)); call io_close(out_file(1))

          ! And now we build the cross section tensor
          in_file(1)  = io_open('cross_section_vector', action='read', status='old')
          out_file(1) = io_open('cross_section_tensor', action='write')

          ! The following routine should now build the tensor...
          call spectrum_cross_section_tensor(s, out_file(1), in_file(1))

          call io_close(in_file(1)); call io_close(out_file(1))

      end select

  endif

  call io_end()
  call syslabels_end()
  call parser_end()
  call global_end()

  contains

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here we start the search for the file(s) that we must process. In future 
! version of octopus, I would like this to be controlled by a command line 
! option. The routine sets the eq_axis and calculate_tensor variables, as well
! as the in_file unit numbers.
!
! The options are:
! (i)  A file called "multipoles" is found. In this case, no other file is 
!      considered.
!      (i.1) If this file signals three equivalent axis, the full tensor is
!            calculated, and placed in "cross_section_tensor".
!      (i.2) If this file signals less than three equivalent axis, the full
!            tensor cannot be calculated, and instead a "cross_section_vector"
!            will be generated.
! (ii) A file called "multipoles.1" is found. In this case, the program will
!      always try to generate the full tensor (calculate_tensor = .true.).
!      The file "cross_esection_tensor" should be generated at the end.
!      Other files are searched for, depending on the equivalent axis that
!      are written in the "multipoles.1" file.
!      (ii.1) Three equivalent axis. No other file is searched for.
!      (ii.2) Two equivalent axis. File "multipoles.2" is searched for; the
!             program ends if it is not found.
!      (ii.3) No equivalent axis. Files "multipoles.2" and "multipoles.3" are
!             searched for; the program ends if they are not found.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/!
subroutine read_files()

  in_file(1) = io_open('multipoles', action='read', status='old', die=.false.)
  if(in_file(1) < 0) in_file(1) = io_open('td.general/multipoles', action='read', status='old', die=.false.)
  if(in_file(1) >= 0) then
    write(message(1),'(a)') 'File "multipoles" found. This will be the only file to be processed.'
    write(message(2),'(a)') '(If more than one file is to be used, the files should be called'
    write(message(3),'(a)') '"multipoles.1", "multipoles.2", etc.'
    write(message(4),'(a)')
    call write_info(4)

    ! OK, so we only have one file. Now we have to see what it has inside.
    call spectrum_mult_info(in_file(1), nspin, lmax, kick, time_steps, dt)
    eq_axis = kick%pol_equiv_axis
    if(eq_axis == 3) then
       calculate_tensor = .true.
       write(message(1),'(a)') 'The file "multipoles" tells me that the system has three equivalent axis.'
       write(message(2),'(a)') 'I will calculate the full tensor, written in file "cross_section_tensor".'
       call write_info(2)
    elseif(eq_axis == 2) then
       write(message(1),'(a)') 'The file "multipoles" tells me that the system has two equivalent axis.'
       write(message(2),'(a)') 'However, I am only using this file; cannot calculate the full tensor.'
       write(message(3),'(a)') 'A file "cross_section_vector" will be generated instead.'
       call write_warning(3)
       calculate_tensor = .false.
    else
       write(message(1),'(a)') 'The file "multipoles" tells me that the system has no usable symmetry. '
       write(message(2),'(a)') 'However, I am only using this file; cannot calculate the full tensor.'
       write(message(3),'(a)') 'A file "cross_section_vector" will be generated instead.'
       call write_warning(3)
       calculate_tensor = .false.
    endif

  else  ! We will try to load more multipoles.1 files...

       ! In this case, we will always want the full tensor
       calculate_tensor = .false.

       in_file(1) = io_open('multipoles.1', action='read', status='old', die=.false.)
       if(in_file(1) < 0) in_file(1) = io_open('td.general/multipoles.1', action='read', status='old', die=.false.)
       if(in_file(1) < 0) then ! Could not find proper files. Die and complain.
           write(message(1),'(a)') 'No "multipoles" or "multipoles.1" file found. At least one of those'
           write(message(2),'(a)') 'should be visible.'
           call write_fatal(2)
       endif

       call spectrum_mult_info(in_file(1), nspin, lmax, kick, time_steps, dt)
       eq_axis = kick%pol_equiv_axis

       if(eq_axis == 3) then
          write(message(1),'(a)') 'The file "multipoles.1" tells me that the system has three equivalent axis.'
          write(message(2),'(a)') 'I will calculate the full tensor, written in file "cross_section_tensor".'
          call write_info(2)
       elseif(eq_axis == 2) then
          in_file(2) = io_open('multipoles.2', action='read', status='old', die=.false.)
          if(in_file(2) < 0) in_file(2) = io_open('td.general/multipoles.2', action='read', status='old', die=.false.)
          if(in_file(2) < 0) then
             write(message(1),'(a)') 'The file "multipoles.1" tells me that the system has two equivalent axis,'
             write(message(2),'(a)') 'but I cannot find a "multipoles.2".'
             call write_fatal(2)
          endif
       else ! No equivalent axis
          in_file(2) = io_open('multipoles.2', action='read', status='old', die=.false.)
          if(in_file(2) < 0) in_file(2) = io_open('td.general/multipoles.2', action='read', status='old', die=.false.)
          if(in_file(2) < 0) then
             write(message(1),'(a)') 'The file "multipoles.1" tells me that the system has three equivalent axis,'
             write(message(2),'(a)') 'but I cannot find a "multipoles.2".'
             call write_fatal(2)
          endif
          in_file(3) = io_open('multipoles.3', action='read', status='old', die=.false.)
          if(in_file(2) < 0) in_file(2) = io_open('td.general/multipoles.3', action='read', status='old', die=.false.)
          if(in_file(2) < 0) then
             write(message(1),'(a)') 'The file "multipoles.1" tells me that the system has two equivalent axis,'
             write(message(2),'(a)') 'but I cannot find a "multipoles.3".'
             call write_fatal(2)
          endif
       endif

  endif

  end subroutine read_files

end program cross_section

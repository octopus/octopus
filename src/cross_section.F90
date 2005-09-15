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

  integer :: in_file_1, in_file_2, out_file_1, out_file_2, eq_axis, i, j
  integer(POINTER_SIZE) :: blk
  logical :: calculate_tensor
  type(spec_type) :: s

  ! Initialize stuff
  call global_init()
  call parser_init()
  call io_init()
  call syslabels_init(1)
  current_label = trim(subsys_label(subsys_run_order(1)))
  call units_init()

  call spectrum_init(s)

  ! In the future, it would be nice if these two were entered through the command line,
  ! so that no input file is needed for this utility. Also, the file names (multipoles.x)
  ! should be entered in the command line.
  call loct_parse_logical(check_inp('SpectrumCalculateTensor'), .false., calculate_tensor)
  call loct_parse_int(check_inp('TDPolarizationEquivAxis'), 0, eq_axis)
  select case(eq_axis)
    case(0, 1)
      write(message(1),'(a)') 'No symmetry information will be used.'
      write(message(2),'(a)') 'Sorry, this case is not implemnted yet.'
      call write_fatal(2)
    case(2)
      write(message(1),'(a)') 'The system contains two equivalent axis.'
    case(3)
      write(message(1),'(a)') 'The system contains three equivalent axis.'
    case default
      write(message(1),'(a)') 'The number of equivalent axis must be a number between zero and three.'
      call write_fatal(1)
  end select
  call write_info(1)

  if(.not.calculate_tensor) then

     call io_assign(in_file_1)
     in_file_1 = io_open('multipoles', action='read', status='old', die=.false.)
     if(in_file_1 < 0) then
       in_file_1 = io_open('td.general/multipoles', action='read', status='old')
     end if

     call io_assign(out_file_1)
     out_file_1 = io_open('cross_section_vector', action='write')

     call spectrum_cross_section(in_file_1, out_file_1, s)

     call io_close(in_file_1)
     call io_close(out_file_1)

  else



      select case(eq_axis)

        case(2)

          call io_assign(in_file_1)
          in_file_1 = io_open('multipoles.1', action='read', status='old', die=.false.)
          if(in_file_1 < 0) then
             in_file_1 = io_open('td.general/multipoles.1', action='read', status='old')
          end if
          call io_assign(in_file_2)
          in_file_2 = io_open('multipoles.2', action='read', status='old', die=.false.)
          if(in_file_2 < 0) then
             in_file_2 = io_open('td.general/multipoles.2', action='read', status='old')
          end if
          call io_assign(out_file_1)
          out_file_1 = io_open('cross_section_vector.1', action='write')
          call io_assign(out_file_2)
          out_file_2 = io_open('cross_section_vector.2', action='write')

          call spectrum_cross_section(in_file_1, out_file_1, s)
          call io_close(in_file_1)
          call io_close(out_file_1)
          call spectrum_cross_section(in_file_2, out_file_2, s)
          call io_close(in_file_2)
          call io_close(out_file_2)

          ! And now we build the cross section tensor
          in_file_1  = io_open('cross_section_vector.1', action='read', status='old')
          in_file_2  = io_open('cross_section_vector.2', action='read', status='old')
          out_file_1 = io_open('cross_section_tensor', action='write')

          !call spectrum_cross_section(s, out_file, in_file)
          call io_close(in_file_1)
          call io_close(out_file_1)

        case(3)

          call io_assign(in_file_1)
          in_file_1 = io_open('multipoles', action='read', status='old', die=.false.)
          if(in_file_1 < 0) then
             in_file_1 = io_open('td.general/multipoles', action='read', status='old')
          end if
          call io_assign(out_file_1)
          out_file_1 = io_open('cross_section_vector', action='write')

          call spectrum_cross_section(in_file_1, out_file_1, s)
          call io_close(in_file_1)
          call io_close(out_file_1)

          ! And now we build the cross section tensor
          in_file_1  = io_open('cross_section_vector', action='read', status='old')
          out_file_1 = io_open('cross_section_tensor', action='write')

          ! The following routine should now build the tensor...
          call spectrum_cross_section_tensor(s, out_file_1, in_file_1)

          call io_close(in_file_1)
          call io_close(out_file_1)

        case default

      end select

  endif

  call syslabels_end()
  call io_end()
  call parser_end()
  call global_end()
end program cross_section

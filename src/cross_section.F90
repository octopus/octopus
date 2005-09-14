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

  integer :: in_file, out_file, eq_axis, i, j
  integer(POINTER_SIZE) :: blk
  logical :: calculate_tensor
  character(len=100) :: txt
  type(spec_type) :: s

  ! Initialize stuff
  call global_init()
  call parser_init()
  call io_init()
  call syslabels_init(1)
  current_label = trim(subsys_label(subsys_run_order(1)))
  call units_init()

  call spectrum_init(s)

  call loct_parse_logical(check_inp('SpectrumCalculateTensor'), .false., calculate_tensor)

  ! Hard coded for the moment...
  eq_axis = 3

  if(.not.calculate_tensor) then

     call io_assign(in_file)
     in_file = io_open('multipoles', action='read', status='old', die=.false.)
     if(in_file < 0) then
       in_file = io_open('td.general/multipoles', action='read', status='old')
     end if

     call io_assign(out_file)
     out_file = io_open('cross_section_vector', action='write')

     call spectrum_cross_section(in_file, out_file, s)

     call io_close(in_file)
     call io_close(out_file)

  else



      select case(eq_axis)

        case(2)

        case(3)

          call io_assign(in_file)
          in_file = io_open('multipoles', action='read', status='old', die=.false.)
          if(in_file < 0) then
             in_file = io_open('td.general/multipoles', action='read', status='old')
          end if
          call io_assign(out_file)
          out_file = io_open('cross_section_vector', action='write')

          call spectrum_cross_section(in_file, out_file, s)
          call io_close(in_file)
          call io_close(out_file)

          ! And now we build the cross section tensor
          in_file  = io_open('cross_section_vector', action='read', status='old')
          out_file = io_open('cross_section_tensor', action='write')

          ! The following routine should now build the tensor...
          call spectrum_cross_section_tensor(in_file, out_file, s)

          call io_close(in_file)
          call io_close(out_file)

        case default

      end select

  endif

  call syslabels_end()
  call io_end()
  call parser_end()
  call global_end()
end program cross_section

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
  FLOAT, allocatable :: basis(:, :)
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


  allocate(basis(3, 3))
  ! WARNING : this is broken for 2D or 3D
  if(loct_parse_block(check_inp('SpectrumBasis'), blk) == 0) then
    i = loct_parse_block_n(blk)
    if(i.ne.3) then
       write(message(1),'(a,i1,a)') 'The number of vectors in the SpectrumBasis block is not ',i,', as it should.'
       call write_fatal(1)
    endif
    do j = 1, 3
      call loct_parse_block_float(blk, j-1, 0, basis(1, j))
      call loct_parse_block_float(blk, j-1, 1, basis(2, j))
      call loct_parse_block_float(blk, j-1, 2, basis(3, j))
    enddo
  else
    basis(1:3, 1) = (/ CNST(1.0), CNST(0.0), CNST(0.0) /)
    basis(1:3, 2) = (/ CNST(0.0), CNST(1.0), CNST(0.0) /)
    basis(1:3, 3) = (/ CNST(0.0), CNST(0.0), CNST(1.0) /)
  endif

  ! Normalize the input vectors.
  do j = 1, 3
     basis(1:3, j) = basis(1:3, j)/sqrt(dot_product(basis(1:3, j),basis(1:3, j)))
  enddo

  call loct_parse_logical(check_inp('SpectrumCalculateTensor'), .false., calculate_tensor)
  call loct_parse_int(check_inp('SpectrumEquivalentAxis'), 1, eq_axis)

  if(.not.calculate_tensor) then

     call io_assign(in_file)
     in_file = io_open('multipoles', action='read', status='old', die=.false.)
     if(in_file < 0) then
       in_file = io_open('td.general/multipoles', action='read', status='old')
     end if

     call io_assign(out_file)
     out_file = io_open('cross_section_vector', action='write')

     call spectrum_cross_section(in_file, out_file, s, basis)

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

          ! Here one should read the polarization vector of in_file, and check that it
          ! corresponds with basis(1:3, 1)
          call spectrum_cross_section(in_file, out_file, s, basis)
          call io_close(in_file)
          call io_close(out_file)

          ! And now we build the cross section tensor
          in_file  = io_open('cross_section_vector', action='read', status='old')
          out_file = io_open('cross_section_tensor', action='write')

          ! The following routine should now build the tensor...
          !call spectrum_cross_section_tensor(in_file, out_file, basis)


        case default

      end select

  endif

  call syslabels_end()
  call io_end()
  call parser_end()
  call global_end()
end program cross_section

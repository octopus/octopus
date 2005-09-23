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

program rotational_strength
  use global
  use messages
  use syslabels
  use lib_oct_parser
  use io
  use units
  use spectrum

  implicit none

  character(len=100) :: txt
  type(spec_type) :: s
  type(spec_rsf) :: rsf

  integer :: i
  integer(POINTER_SIZE) :: blk

  ! Initialize stuff
  call global_init()
  call parser_init()
  call io_init()
  call syslabels_init(1)
  if(in_debug_mode) then
     call io_mkdir('debug')
  endif
  call units_init()

  call spectrum_init(s)

  call loct_parse_float(check_inp('TDDeltaStrength'), CNST(0.05), rsf%delta_strength)
  rsf%delta_strength = rsf%delta_strength / units_inp%length%factor

  !!! read in the default direction for the polarization
  rsf%pol(:) = M_ZERO
  if(loct_parse_block(check_inp('TDPolarization'), blk)==0) then
    do i = 1, 3
      call loct_parse_block_float(blk, 0, i-1, rsf%pol(i))
    end do
    call loct_parse_block_end(blk)
  else  !default along the x-direction
    rsf%pol(1) = M_ONE
  endif

  call spectrum_rotatory_strength('rotatory_strength', s, rsf, .true.)

  call syslabels_end()
  call io_end()
  call parser_end()
  call global_end()

  stop
end program rotational_strength

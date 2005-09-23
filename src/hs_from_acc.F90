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

program hs_from_acc
  use global
  use messages
  use syslabels
  use io
  use units
  use lib_oct_parser
  use spectrum

  character(len=100) :: txt
  type(spec_type) :: s
  type(spec_sh) :: sh

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

  call loct_parse_string(check_inp('HSPolarization'), 'z', txt)
  sh%pol = txt(1:1)
  if(sh%pol.ne.'x' .and. sh%pol.ne.'y' .and. sh%pol.ne.'z' .and. &
       sh%pol.ne.'+' .and. sh%pol.ne.'-') then
    message(1) = "HSPolarization has an invalid value"
    message(2) = "Valid values are ('x' | 'y' | 'z' | '+' | '-')"
    call write_fatal(2)
  end if

  call spectrum_hs_from_acc('hs-acc', s, sh)

  deallocate(sh%sp)

  call io_end()
  call syslabels_end()
  call parser_end()
  call global_end()

  stop
end program hs_from_acc

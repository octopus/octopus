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

program hs_from_mult
  use global
  use units
  use lib_oct_parser
  use spectrum

  character(len=100) :: txt
  type(spec_type) :: s
  type(spec_sh) :: sh

  ! Initialize stuff
  call global_init()
  call units_init()

  call loct_parse_float(check_inp('SpecStartTime'),  M_ZERO,     s%start_time)
  call loct_parse_float(check_inp('SpecEndTime'),   -M_ONE,      s%end_time)
  call loct_parse_float(check_inp('SpecEnergyStep'), CNST(0.05), s%energy_step)
  call loct_parse_float(check_inp('SpecMinEnergy'),  M_ZERO,     s%min_energy)
  call loct_parse_float(check_inp('SpecMaxEnergy'),  CNST(20.0), s%max_energy)

  ! adjust units
  s%start_time  = s%start_time  * units_inp%time%factor
  s%end_time    = s%end_time    * units_inp%time%factor
  s%energy_step = s%energy_step * units_inp%energy%factor
  s%min_energy  = s%min_energy  * units_inp%energy%factor
  s%max_energy  = s%max_energy  * units_inp%energy%factor

  call loct_parse_string(check_inp('HSPolarization'), 'z', txt)
  sh%pol = txt(1:1)
  if(sh%pol.ne.'x' .and. sh%pol.ne.'y' .and. sh%pol.ne.'z' .and. &
       sh%pol.ne.'+' .and. sh%pol.ne.'-') then
    message(1) = "HSPolarization has an invalid value"
    message(2) = "Valid values are ('x' | 'y' | 'z' | '+' | '-')"
    call write_fatal(2)
  end if

  call spectrum_hs_from_mult('hs-mult', s, sh, .true.)

  deallocate(sh%sp)
  stop  
end program hs_from_mult

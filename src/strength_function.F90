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

program strength_function
  use global
  use messages
  use syslabels
  use lib_oct_parser
  use units
  use spectrum

  character(len=100) :: txt
  type(spec_type) :: s
  type(spec_sf) :: sf

  ! Initialize stuff
  call global_init()
  call parser_init()
  call io_init()
  call syslabels_init(1)
  call units_init()

  call loct_parse_string(check_inp('SpecTransformMode'), "sin", txt)
  select case(txt(1:3))
  case('sin')
    sf%transform = 1
  case('cos')
    sf%transform = 2
  case default
    write(message(1), '(2a)') trim(txt), ' is not a valid transform mode'
    message(2) = "SpecTransformMode = sin | cos"
    call write_fatal(2)
  end select
  
  call loct_parse_string(check_inp('SpecDampMode'), "exp", txt)
  select case(txt(1:3))
  case('exp')
    sf%damp = SPECTRUM_DAMP_LORENTZIAN
  case('pol')
    sf%damp = SPECTRUM_DAMP_POLYNOMIAL
  case('gau')
    sf%damp = SPECTRUM_DAMP_GAUSSIAN
  case default
    sf%damp = SPECTRUM_DAMP_NONE
  end select

  call loct_parse_float(check_inp('SpecDampFactor'), CNST(0.15), sf%damp_factor)
  call loct_parse_float(check_inp('SpecStartTime'),  M_ZERO,      s%start_time)
  call loct_parse_float(check_inp('SpecEndTime'),   -M_ONE,       s%end_time)
  call loct_parse_float(check_inp('SpecEnergyStep'), CNST(0.05),  s%energy_step)
  call loct_parse_float(check_inp('SpecMaxEnergy'),  CNST(20.0),  s%max_energy)
  call loct_parse_float(check_inp('SpecMinEnergy'),  M_ZERO,      s%min_energy)
  call loct_parse_float(check_inp('TDDeltaStrength'),CNST(0.05), sf%delta_strength)

  ! adjust units
  sf%damp_factor    = sf%damp_factor    / units_inp%time%factor
  s%start_time      = s%start_time      * units_inp%time%factor
  s%end_time        = s%end_time        * units_inp%time%factor
  s%energy_step     = s%energy_step     * units_inp%energy%factor
  s%max_energy      = s%max_energy      * units_inp%energy%factor
  s%min_energy      = s%min_energy      * units_inp%energy%factor
  sf%delta_strength = sf%delta_strength / units_inp%length%factor

  call spectrum_strength_function('spectrum', s, sf, .true.)

  deallocate(sf%sp)

  call syslabels_end()
  call io_end()
  call parser_end()
  call global_end()

  stop  
end program strength_function

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

program strength_function
  use spectrum

  integer :: ierr
  character(len=100) :: txt
  type(spec_type) :: s
  type(spec_sf) :: sf

  ! Initialize stuff
  call global_init()
  call units_init()

  call oct_parse_string("SpecTransformMode", "sin", txt)
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
  
  call oct_parse_string("SpecDampMode", "exp", txt)
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

  call oct_parse_double("SpecDampFactor", CNST(0.15), sf%damp_factor)
  call oct_parse_double("SpecStartTime",  M_ZERO,      s%start_time)
  call oct_parse_double("SpecEndTime",   -M_ONE,       s%end_time)
  call oct_parse_double("SpecEnergyStep", CNST(0.05),  s%energy_step)
  call oct_parse_double("SpecMaxEnergy",  CNST(20.0),  s%max_energy)
  call oct_parse_double("SpecMinEnergy",  M_ZERO,      s%min_energy)
  call oct_parse_double("TDDeltaStrength",CNST(0.05), sf%delta_strength)

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
  stop  
end program strength_function

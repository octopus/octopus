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

#include "config_F90.h"

program strength_function
  use global
  use spectrum
  use liboct

  integer :: ierr
  character(len=100) :: txt
  type(spec_type) :: s
  type(spec_sf) :: sf

  ! init liboct
  ierr = oct_parse_init(C_string('inp'), C_string('out.oct'))
  if(ierr .ne. 0) then
    message(1) = "Error initializing liboct"
    call write_fatal(1)
  end if
  call oct_parse_int(C_string("verbose"), 30, conf%verbose)
  call units_init()

  call oct_parse_str("SpecTransformMode", "sin", txt)
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
  
  call oct_parse_str("SpecDampMode", "exp", txt)
  select case(txt(1:3))
  case('exp')
    sf%damp = 1
  case('pol')
    sf%damp = 2
  case default
    sf%damp = 0
  end select

  call oct_parse_double(C_string("SpecDampFactor"), 0.15_r8, sf%damp_factor)
  call oct_parse_double(C_string("SpecStartTime"), 0._r8, s%start_time)
  call oct_parse_double(C_string("SpecEndTime"), -1._r8, s%end_time)
  call oct_parse_double(C_string("SpecEnergyStep"), 0.05_r8, s%energy_step)
  call oct_parse_double(C_string("SpecMaxEnergy"), 20._r8, s%max_energy)
  call oct_parse_double(C_string("TDDeltaStrength"), 0.05_r8, sf%delta_strength)

  ! adjust units
  sf%damp_factor    = sf%damp_factor    / units_inp%time%factor
  s%start_time      = s%start_time      * units_inp%time%factor
  s%end_time        = s%end_time        * units_inp%time%factor
  s%energy_step     = s%energy_step     * units_inp%energy%factor
  s%max_energy      = s%max_energy      * units_inp%energy%factor
  sf%delta_strength = sf%delta_strength / units_inp%length%factor

  call oct_parse_str('SystemName', 'system', txt)
  call spectrum_strength_function(trim(txt), trim(txt)//'.spectrum', s, sf, .true.)

  deallocate(sf%sp)
  stop  
end program strength_function

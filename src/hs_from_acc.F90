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

program hs_from_acc
  use global
  use spectrum
  use liboct

  integer :: ierr
  character(len=100) :: txt
  type(spec_type) :: s
  type(spec_sh) :: sh

  ! init liboct
  ierr = oct_parse_init(C_string('inp'), C_string('out.oct'))
  if(ierr .ne. 0) then
    ierr = oct_parse_init(C_string("-"), C_string('out.oct'))
    if(ierr .ne. 0) then
      message(1) = "Error initializing liboct"
      call write_fatal(1)
    end if
  end if

  call oct_parse_int(C_string("verbose"), 30, conf%verbose)
  call units_init()

  call oct_parse_double(C_string("SpecStartTime"), 0._r8, s%start_time)
  call oct_parse_double(C_string("SpecEndTime"), -1._r8, s%end_time)
  call oct_parse_double(C_string("SpecEnergyStep"), 0.05_r8, s%energy_step)
  call oct_parse_double(C_string("SpecMinEnergy"),  0._r8, s%min_energy)
  call oct_parse_double(C_string("SpecMaxEnergy"), 20._r8, s%max_energy)

  ! adjust units
  s%start_time  = s%start_time  * units_inp%time%factor
  s%end_time    = s%end_time    * units_inp%time%factor
  s%energy_step = s%energy_step * units_inp%energy%factor
  s%min_energy  = s%min_energy  * units_inp%energy%factor
  s%max_energy  = s%max_energy  * units_inp%energy%factor

  call oct_parse_str('HSPolarization', 'z', txt)
  sh%pol = txt(1:1)
  if(sh%pol.ne.'x' .and. sh%pol.ne.'y' .and. sh%pol.ne.'z' .and. &
       sh%pol.ne.'+' .and. sh%pol.ne.'-') then
    message(1) = "HSPolarization has an invalid value"
    message(2) = "Valid values are ('x' | 'y' | 'z' | '+' | '-')"
    call write_fatal(2)
  end if

  call spectrum_hs_from_acc('hs-acc', s, sh, .true.)

  deallocate(sh%sp)
  stop  
end program hs_from_acc

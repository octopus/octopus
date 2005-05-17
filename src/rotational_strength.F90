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
  use lib_oct_parser
  use units
  use spectrum

  implicit none


!!$  integer :: ierr
  character(len=100) :: txt
  type(spec_type) :: s
  type(spec_rsf) :: rsf

  integer :: i
  integer(POINTER_SIZE) :: blk

  ! Initialize stuff
  call global_init()
  call units_init()

  call loct_parse_string(check_inp('SpecDampMode'), "exp", txt)
  select case(txt(1:3))
  case('exp')
    rsf%damp = SPECTRUM_DAMP_LORENTZIAN
  case('pol')
    rsf%damp = SPECTRUM_DAMP_POLYNOMIAL
  case('gau')
    rsf%damp = SPECTRUM_DAMP_GAUSSIAN
  case default
    rsf%damp = SPECTRUM_DAMP_NONE
  end select

  call loct_parse_float(check_inp('SpecDampFactor'),  CNST(0.15), rsf%damp_factor)
  call loct_parse_float(check_inp('SpecStartTime'),   M_ZERO,       s%start_time)
  call loct_parse_float(check_inp('SpecEndTime'),    -M_ONE,        s%end_time)
  call loct_parse_float(check_inp('SpecEnergyStep'),  CNST(0.05),   s%energy_step)
  call loct_parse_float(check_inp('SpecMaxEnergy'),   CNST(20.0),   s%max_energy)
  call loct_parse_float(check_inp('SpecMinEnergy'),   M_ZERO,       s%min_energy)
  call loct_parse_float(check_inp('TDDeltaStrength'), CNST(0.05), rsf%delta_strength)
  !!! read in the default direction for the polarization
  rsf%pol(:) = M_ZERO
  if(loct_parse_block(check_inp('TDPolarization'), blk)==0) then
    do i = 1, conf%dim
      call loct_parse_block_float(blk, 0, i-1, rsf%pol(i))
    end do
    call loct_parse_block_end(blk)
  else  !default along the x-direction
    rsf%pol(1) = M_ONE
  endif

  ! adjust units
  rsf%damp_factor    = rsf%damp_factor    / units_inp%time%factor
  s%start_time      = s%start_time      * units_inp%time%factor
  s%end_time        = s%end_time        * units_inp%time%factor
  s%energy_step     = s%energy_step     * units_inp%energy%factor
  s%max_energy      = s%max_energy      * units_inp%energy%factor
  s%min_energy      = s%min_energy      * units_inp%energy%factor
  rsf%delta_strength = rsf%delta_strength / units_inp%length%factor

  call spectrum_rotatory_strength('rotatory_strength', s, rsf, .true.)

  stop  
end program rotational_strength

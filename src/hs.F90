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

program harmonic_spectrum
  use global_m
  use messages_m
  use datasets_m
  use io_m
  use units_m
  use lib_oct_parser_m
  use spectrum_m
  use varinfo_m

  implicit none

  character(len=100) :: txt
  integer :: mode
  type(spec_t) :: s
  type(spec_sh)   :: sh

  integer, parameter :: &
    HS_FROM_MULT = 1,   &
    HS_FROM_ACC  = 2

  ! Initialize stuff
  call global_init()
  call parser_init()
  call datasets_init(1)
  call io_init()
  if(in_debug_mode) then
    call io_mkdir('debug')
  end if
  call units_init()

  call spectrum_init(s)

  !%Variable HarmonicSpectrumPolarization
  !%Type string
  !%Default "z"
  !%Section Utilities::oct-harmonic-spectrum
  !%Description
  !% The oct-harmonic-spectrum utility program needs to know the direction along
  !% which the emission raidiation is considered to be polarized. It may be
  !% linearly polarized or circularly polarized.
  !%Option "x"
  !% Linearly polarized field in the x direction.
  !%Option "y"
  !% Linearly polarized field in the y direction.
  !%Option "z"
  !% Linearly polarized field in the z direction.
  !%Option "+"
  !% Circularly polarized field, counter clock-wise.
  !%Option "-"
  !% Circularly polarized field, clock-wise.
  !%End
  call loct_parse_string(check_inp('HarmonicSpectrumPolarization'), 'z', txt)
  sh%pol = txt(1:1)
  if(sh%pol.ne.'x' .and. sh%pol.ne.'y' .and. sh%pol.ne.'z' .and. &
    sh%pol.ne.'+' .and. sh%pol.ne.'-') then
    call input_error('HarmonicSpectrumPolarization')
  end if

  !%Variable HarmonicSpectrumMode
  !%Type integer
  !%Default hs_from_dipole
  !%Section Utilities::oct-harmonic-spectrum
  !%Description
  !% The oct-harmonic-spectrum may calculate the spectrum in two alternative ways,
  !% mathematically equivalent but numerically diferent: by reading the dipole
  !% moment (from the multipoles file) and calculating the accelaratio numerically
  !% from it, or by reading directly the acceleration from the acceleration file,
  !% which may also be generated during a time-dependent run of octopus.
  !%Option hs_from_dipole 1
  !% Calculate the harmonic spectrum by numerically differentiating the multipoles file.
  !%Option hs_from_acceleration 2
  !% Calculate the harmonic spectrum by reading the acceleration file.
  !%End
  call loct_parse_int(check_inp('HarmonicSpectrumMode'), HS_FROM_MULT, mode)
  if(.not.varinfo_valid_option('HarmonicSpectrumMode', mode)) call input_error('HarmonicSpectrumMode')

  select case(mode)
  case(HS_FROM_MULT)
    call spectrum_hs_from_mult('hs-mult', s, sh)
  case(HS_FROM_ACC)
    call spectrum_hs_from_acc('hs-acc', s, sh)
  end select

  deallocate(sh%sp)

  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()
end program harmonic_spectrum

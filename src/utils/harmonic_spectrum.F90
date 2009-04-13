!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
  use command_line_m
  use global_m
  use messages_m
  use datasets_m
  use io_m
  use units_m
  use loct_parser_m
  use spectrum_m
  use varinfo_m

  implicit none

  integer :: mode, ierr
  FLOAT :: w0
  type(spec_t) :: s
  character :: pol
  logical :: get_maxima

  integer, parameter :: &
    HS_FROM_MULT = 1,   &
    HS_FROM_ACC  = 2

  call getopt_init(ierr)
  if(ierr.ne.0) then
    write(stderr, '(a)') "Your fortran compiler doesn't support command line arguments;"
    write(stderr, '(a)') "the oct-harmonic-spectrum command is not available."
    stop
  end if

  ! These are the default values.
  get_maxima = .true.
  w0 = M_ZERO
  pol = 'x'
  mode = 1
  call getopt_harmonic_spectrum(w0, mode, pol)
  if(w0 <= M_ZERO) get_maxima = .false.

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

  call obsolete_variable('HarmonicSpectrumPolarization')
  call obsolete_variable('HarmonicSpectrumMode')

  if( (pol.ne.'x') .and. &
      (pol.ne.'y') .and. &
      (pol.ne.'z') .and. &
      (pol.ne.'+') .and. &
      (pol.ne.'-') ) then
    message(1) = 'The polarization direction given in the command line is not valid.'
    call write_fatal(1)
  end if
  if( (mode.ne.HS_FROM_MULT) .and. &
      (mode .ne.HS_FROM_ACC) ) then
    message(1) = 'The harmonic spectrum mode given in the command line is not valid.'
    call write_fatal(1)
  end if

  select case(mode)
  case(HS_FROM_MULT)
    if(get_maxima) then
      call spectrum_hs_from_mult('hs-mult-maxima', s, pol, w0)
    else
      call spectrum_hs_from_mult('hs-mult', s, pol)
    end if
  case(HS_FROM_ACC)
    if(get_maxima) then
      call spectrum_hs_from_acc('hs-acc-maxima', s, pol, w0)
    else
      call spectrum_hs_from_acc('hs-acc', s, pol)
    end if
  end select


  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()
end program harmonic_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

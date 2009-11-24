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

program broad
  use datasets_m
  use global_m
  use io_m
  use messages_m
  use parser_m
  use profiling_m
  use unit_m
  use unit_system_m

  implicit none

  type broad_t
    FLOAT :: br, energy_step, min_energy, max_energy
  end type broad_t

  type(broad_t) :: br

  ! Initialize stuff
  call global_init()
  call parser_init()
  call datasets_init(1)
  call io_init()
  if(in_debug_mode) then
     call io_mkdir('debug')
  end if
  call unit_system_init()

  !%Variable LinBroadening
  !%Type float
  !%Default 0.005
  !%Section Utilities::oct-broad
  !%Description
  !% Width of the Lorentzian used to broaden the excitations.
  !%End
  call parse_float(datasets_check('LinBroadening'), CNST(0.005), br%br, units_inp%energy)

  !%Variable LinEnergyStep
  !%Type float
  !%Default 0.001
  !%Section Utilities::oct-broad
  !%Description
  !% Sampling rate for the spectrum. 
  !%End
  call parse_float(datasets_check('LinEnergyStep'), CNST(0.001), br%energy_step, units_inp%energy)

  !%Variable LinMinEnergy
  !%Type float
  !%Default 0.0
  !%Section Utilities::oct-broad
  !%Description
  !% The broadening is done for energies greater than LinMinEnergy.
  !%End
  call parse_float(datasets_check('LinMinEnergy'), M_ZERO, br%min_energy, units_inp%energy)

  !%Variable LinMaxEnergy
  !%Type float
  !%Default 1.0
  !%Section Utilities::oct-broad
  !%Description
  !% The broadening is done for energies smaller than LinMaxEnergy.
  !%End
  call parse_float(datasets_check('LinMaxEnergy'), M_ONE, br%max_energy, units_inp%energy)

  call calc_broad(br, CASIDA_DIR, 'eps-diff', .true.)
  call calc_broad(br, CASIDA_DIR, 'petersilka', .true.)
  call calc_broad(br, CASIDA_DIR, 'casida', .false.)

  call io_end()
  call datasets_end()
  call parser_end()
  call global_end()

contains

  ! ---------------------------------------------------------
  subroutine calc_broad(br, dir, fname, extracols)
    type(broad_t),    intent(in) :: br
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: fname
    logical,          intent(in) :: extracols

    FLOAT, allocatable :: spectrum(:,:)
    FLOAT :: omega, energy, ff(4)
    integer :: nsteps, iunit, j1, j2, ii

    nsteps = (br%max_energy - br%min_energy) / br%energy_step
    SAFE_ALLOCATE(spectrum(1:4, 1:nsteps))
    spectrum = M_ZERO

    iunit = io_open(trim(dir)//"/"// fname, action='read', status='old', die = .false.)
    if(iunit < 0) return

    read(iunit, *) ! skip header
    do
      if(extracols) then
        read(iunit, *, end=100) j1, j2, energy, ff
      else
        read(iunit, *, end=100) energy, ff
      end if

      energy = units_to_atomic(units_out%energy, energy)

      do j1 = 1, nsteps
        omega = br%min_energy + real(j1-1, REAL_PRECISION)*br%energy_step
        spectrum(1:4, j1) = spectrum(1:4, j1) + ff(1:4)*br%br/((omega-energy)**2 + br%br**2)/M_PI ! Lorentzian
      end do
    end do
100   continue
    call io_close(iunit)

    ! print spectra
    iunit = io_open(trim(dir)//"/spectrum."//fname, action='write')
    do j1 = 1, nsteps
      write(iunit, '(5es14.6)') units_from_atomic(units_out%energy, br%min_energy + real(j1 - 1, REAL_PRECISION) &
        *br%energy_step), (units_from_atomic(units_out%energy, spectrum(ii, j1)), ii = 1, 4)
    end do

    call io_close(iunit)

    SAFE_DEALLOCATE_A(spectrum)
  end subroutine calc_broad

end program broad

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

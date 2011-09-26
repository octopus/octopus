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

program casida_spectrum
  use command_line_m
  use datasets_m
  use global_m
  use io_m
  use messages_m
  use parser_m
  use profiling_m
  use unit_m
  use unit_system_m

  implicit none

  type casida_spectrum_t
    FLOAT :: br, energy_step, min_energy, max_energy
  end type casida_spectrum_t

  integer :: ierr
  type(casida_spectrum_t) :: cs

  ! Initialize stuff
  call global_init()

  call getopt_init(ierr)
  if(ierr.eq.0) call getopt_casida_spectrum
  call getopt_end()

  call parser_init()
  call messages_init()
  call datasets_init(1)
  call io_init()
  call unit_system_init()

  !%Variable CasidaSpectrumBroadening
  !%Type float
  !%Default 0.005
  !%Section Utilities::oct-casida_spectrum
  !%Description
  !% Width of the Lorentzian used to broaden the excitations.
  !%End
  call parse_float(datasets_check('CasidaSpectrumBroadening'), CNST(0.005), cs%br, units_inp%energy)

  call messages_print_var_value(stdout, datasets_check('CasidaSpectrumBroadening'), cs%br, unit = units_out%energy)

  !%Variable CasidaSpectrumEnergyStep
  !%Type float
  !%Default 0.001
  !%Section Utilities::oct-casida_spectrum
  !%Description
  !% Sampling rate for the spectrum. 
  !%End
  call parse_float(datasets_check('CasidaSpectrumEnergyStep'), CNST(0.001), cs%energy_step, units_inp%energy)

  call messages_print_var_value(stdout, datasets_check('CasidaSpectrumEnergyStep'), cs%energy_step, unit = units_out%energy)

  !%Variable CasidaSpectrumMinEnergy
  !%Type float
  !%Default 0.0
  !%Section Utilities::oct-casida_spectrum
  !%Description
  !% The broadening is done for energies greater than <tt>CasidaSpectrumMinEnergy</tt>.
  !%End
  call parse_float(datasets_check('CasidaSpectrumMinEnergy'), M_ZERO, cs%min_energy, units_inp%energy)

  call messages_print_var_value(stdout, datasets_check('CasidaSpectrumMinEnergy'), cs%min_energy, unit = units_out%energy)

  !%Variable CasidaSpectrumMaxEnergy
  !%Type float
  !%Default 1.0
  !%Section Utilities::oct-casida_spectrum
  !%Description
  !% The broadening is done for energies smaller than <tt>CasidaSpectrumMaxEnergy</tt>.
  !%End
  call parse_float(datasets_check('CasidaSpectrumMaxEnergy'), M_ONE, cs%max_energy, units_inp%energy)

  call messages_print_var_value(stdout, datasets_check('CasidaSpectrumMaxEnergy'), cs%max_energy, unit = units_out%energy)

  call calc_broad(cs, CASIDA_DIR, 'eps_diff', .true.)
  call calc_broad(cs, CASIDA_DIR, 'petersilka', .true.)
  call calc_broad(cs, CASIDA_DIR, 'casida', .false.)

  call io_end()
  call datasets_end()
  call messages_end()
  call parser_end()
  call global_end()

contains

  ! ---------------------------------------------------------
  subroutine calc_broad(cs, dir, fname, extracols)
    type(casida_spectrum_t), intent(in) :: cs
    character(len=*),        intent(in) :: dir
    character(len=*),        intent(in) :: fname
    logical,                 intent(in) :: extracols

    FLOAT, allocatable :: spectrum(:,:)
    FLOAT :: omega, energy, ff(4)
    integer :: nsteps, iunit, j1, j2, ii

    nsteps = (cs%max_energy - cs%min_energy) / cs%energy_step
    SAFE_ALLOCATE(spectrum(1:4, 1:nsteps))
    spectrum = M_ZERO

    iunit = io_open(trim(dir)//"/"// fname, action='read', status='old', die = .false.)

    if(iunit < 0) then
      message(1) = 'Cannot open file "'//trim(dir)//'/'//trim(fname)//'".'
      message(2) = 'The '//trim(fname)//' spectrum was not generated.'
      call messages_warning(2)
      return
    end if

    read(iunit, *) ! skip header
    do
      if(extracols) then
        read(iunit, *, end=100) j1, j2, energy, ff
      else
        read(iunit, *, end=100) energy, ff
      end if

      energy = units_to_atomic(units_out%energy, energy)

      do j1 = 1, nsteps
        omega = cs%min_energy + real(j1-1, REAL_PRECISION)*cs%energy_step
        spectrum(1:4, j1) = spectrum(1:4, j1) + ff(1:4)*cs%br/((omega-energy)**2 + cs%br**2)/M_PI ! Lorentzian
      end do
    end do
100   continue
    call io_close(iunit)

    ! print spectra
    iunit = io_open(trim(dir)//"/spectrum."//fname, action='write')
    do j1 = 1, nsteps
      write(iunit, '(5es14.6)') units_from_atomic(units_out%energy, cs%min_energy + real(j1 - 1, REAL_PRECISION) &
        *cs%energy_step), (units_from_atomic(units_out%energy, spectrum(ii, j1)), ii = 1, 4)
    end do

    call io_close(iunit)

    SAFE_DEALLOCATE_A(spectrum)
  end subroutine calc_broad

end program casida_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

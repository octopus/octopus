!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2012-2013 D. Strubbe
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
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
  use space_m
  use unit_m
  use unit_system_m
  use utils_m
  use varinfo_m

  implicit none

  type casida_spectrum_t
    FLOAT :: br, energy_step, min_energy, max_energy
    type(space_t) :: space
    integer :: ispin
  end type casida_spectrum_t

  integer :: ierr
  type(casida_spectrum_t) :: cs

  ! Initialize stuff
  call global_init()

  call getopt_init(ierr)
  if(ierr == 0) call getopt_casida_spectrum
  call getopt_end()

  call parser_init()
  call messages_init()
  call datasets_init(1)
  call io_init()
  call unit_system_init()
  call space_init(cs%space)

  ! Reads the spin components. This is read here, as well as in states_init.
  call parse_integer(datasets_check('SpinComponents'), 1, cs%ispin)
  if(.not.varinfo_valid_option('SpinComponents', cs%ispin)) call input_error('SpinComponents')
  cs%ispin = min(2, cs%ispin)

  !%Variable CasidaSpectrumBroadening
  !%Type float
  !%Default 0.005 Ha
  !%Section Utilities::oct-casida_spectrum
  !%Description
  !% Width of the Lorentzian used to broaden the excitations.
  !%End
  call parse_float(datasets_check('CasidaSpectrumBroadening'), CNST(0.005), cs%br, units_inp%energy)

  call messages_print_var_value(stdout, datasets_check('CasidaSpectrumBroadening'), cs%br, unit = units_out%energy)

  !%Variable CasidaSpectrumEnergyStep
  !%Type float
  !%Default 0.001 Ha
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
  !%Default 1.0 Ha
  !%Section Utilities::oct-casida_spectrum
  !%Description
  !% The broadening is done for energies smaller than <tt>CasidaSpectrumMaxEnergy</tt>.
  !%End
  call parse_float(datasets_check('CasidaSpectrumMaxEnergy'), M_ONE, cs%max_energy, units_inp%energy)

  call messages_print_var_value(stdout, datasets_check('CasidaSpectrumMaxEnergy'), cs%max_energy, unit = units_out%energy)

  call calc_broad(cs, CASIDA_DIR, 'eps_diff', .true.)
  call calc_broad(cs, CASIDA_DIR, 'petersilka', .false.)
  call calc_broad(cs, CASIDA_DIR, 'tamm_dancoff', .false.)
  call calc_broad(cs, CASIDA_DIR, 'variational', .false.)
  call calc_broad(cs, CASIDA_DIR, 'casida', .false.)

  call space_end(cs%space)
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
    FLOAT :: omega, energy, tm(MAX_DIM), ff(MAX_DIM+1)
    integer :: istep, nsteps, iunit, trash(3), ii, ncols

    nsteps = (cs%max_energy - cs%min_energy) / cs%energy_step
    SAFE_ALLOCATE(spectrum(1:cs%space%dim+1, 1:nsteps))
    spectrum = M_ZERO

    iunit = io_open(trim(dir)// fname, action='read', status='old', die = .false.)

    if(iunit < 0) then
      message(1) = 'Cannot open file "'//trim(dir)//trim(fname)//'".'
      message(2) = 'The '//trim(fname)//' spectrum was not generated.'
      call messages_warning(2)
      return
    end if

    ! For Casida, CV(2), Tamm-Dancoff, Petersilka: first column is the index of the excitation
    ! For eps_diff: first two columns are occ and unocc states, then spin if spin-polarized
    ncols = 1
    if(extracols) then
      ncols = ncols + cs%ispin
    endif

    read(iunit, *) ! skip header
    do
      read(iunit, *, end=100) trash(1:ncols), energy, tm(1:cs%space%dim), ff(cs%space%dim+1)

      energy = units_to_atomic(units_out%energy, energy)

      do istep = 1, nsteps
        omega = cs%min_energy + real(istep-1, REAL_PRECISION)*cs%energy_step

        ! transition matrix elements by themselves are dependent on gauge in degenerate subspaces
        ! make into oscillator strengths, as in casida_inc.F90 X(oscillator_strengths), and like the last column
        ff(1:cs%space%dim) = (M_TWO / cs%space%dim) * energy * (tm(1:cs%space%dim))**2
        spectrum(1:cs%space%dim+1, istep) = spectrum(1:cs%space%dim+1, istep) + &
          ff(1:cs%space%dim+1)*cs%br/((omega-energy)**2 + cs%br**2)/M_PI ! Lorentzian
      end do
    end do
100 continue
    call io_close(iunit)

    ! print spectra
    iunit = io_open(trim(dir)//"/spectrum."//fname, action='write')

    write(iunit, '(a2,a12)', advance = 'no') '# ', 'E [' // trim(units_abbrev(units_out%energy)) // ']'
    do ii = 1, cs%space%dim
      write(iunit, '(a14)', advance = 'no') '<' // index2axis(ii) // '>^2'
    enddo
    write(iunit, '(a14)') '<f>'

    do istep = 1, nsteps
      write(iunit, '(99es14.6)') units_from_atomic(units_out%energy, cs%min_energy + real(istep - 1, REAL_PRECISION) &
        *cs%energy_step), (units_from_atomic(unit_one/units_out%energy, spectrum(ii, istep)), ii = 1, cs%space%dim+1)
    end do

    call io_close(iunit)

    SAFE_DEALLOCATE_A(spectrum)
  end subroutine calc_broad

end program casida_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

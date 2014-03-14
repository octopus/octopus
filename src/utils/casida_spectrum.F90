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
  use geometry_m
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

  integer :: ierr, idir, jdir, iatom
  type(casida_spectrum_t) :: cs
  FLOAT :: rotation(MAX_DIM, MAX_DIM), rot2(MAX_DIM, MAX_DIM), identity(MAX_DIM, MAX_DIM), coord(MAX_DIM)
  type(block_t) :: blk
  type(geometry_t) :: geo

  ! Initialize stuff
  call global_init(is_serial = .true.)

  call getopt_init(ierr)
  if(ierr == 0) call getopt_casida_spectrum()
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

  identity = M_ZERO
  do idir = 1, cs%space%dim
    identity(idir, idir) = M_ONE
  enddo

  !%Variable CasidaSpectrumRotationMatrix
  !%Type block
  !%Default identity
  !%Section Utilities::oct-casida_spectrum
  !%Description
  !% Supply a rotation matrix to apply to the transition dipoles in generating the spectrum. The rotated atomic structure
  !% will also be output. Size of matrix must be <tt>Dimensions</tt>.
  !%End

  if (parse_block(datasets_check('CasidaSpectrumRotationMatrix'), blk) == 0) then 
    rotation(:,:) = M_ZERO
    do idir = 1, cs%space%dim
      do jdir = 1, cs%space%dim
        call parse_block_float(blk, idir - 1,  jdir - 1, rotation(idir, jdir))
      end do
    end do
    call parse_block_end(blk)

    message(1) = "Info: Applying rotation matrix"
    call messages_info(1)
    call output_tensor(stdout, rotation, cs%space%dim, unit_one, write_average = .false.)

    ! allowing inversions is fine
    rot2(:,:) = abs(matmul(transpose(rotation), rotation))
    if(any(abs(rot2(:,:) - identity(:,:)) > CNST(1e-6))) then
      write(message(1),'(a,es12.6)') "Rotation matrix is not orthogonal. max discrepancy in product = ", &
        maxval(abs(rot2(:,:) - identity(:,:)))
      call messages_warning(1)
    endif

    ! apply rotation to geometry
    call geometry_init(geo, cs%space)
    do iatom = 1, geo%natoms
      coord(1:cs%space%dim) = geo%atom(iatom)%x(1:cs%space%dim)
      geo%atom(iatom)%x(1:cs%space%dim) = matmul(rotation(1:cs%space%dim, 1:cs%space%dim), coord(1:cs%space%dim))
    enddo
    call geometry_write_xyz(geo, trim(CASIDA_DIR)//'rotated')
    call geometry_end(geo)
  else
    rotation(:,:) = identity(:,:)
  end if

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
    FLOAT :: omega, energy, re_tm(MAX_DIM), im_tm(MAX_DIM), ff(MAX_DIM+1), tm_sq(MAX_DIM)
    integer :: istep, nsteps, iunit, trash(3), idir, ncols, ios
    character(len=256) :: string
    logical :: is_complex

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
    
    read(iunit, '(a)') string
    ! complex has this many columns; real lacks the imaginary columns (im_tm)
    read(string, *, iostat = ios) trash(1:ncols), energy, (re_tm(idir), im_tm(idir), idir = 1, cs%space%dim), ff(cs%space%dim+1)
    is_complex = (ios == 0)
    ! put it back to read the data in the loop
    backspace(iunit)

    do
      if(is_complex) then
        read(iunit, *, iostat = ios) trash(1:ncols), energy, (re_tm(idir), im_tm(idir), idir = 1, cs%space%dim), ff(cs%space%dim+1)
        im_tm(:) = units_to_atomic(units_out%length, im_tm(:))
      else
        read(iunit, *, iostat = ios) trash(1:ncols), energy, (re_tm(idir),              idir = 1, cs%space%dim), ff(cs%space%dim+1)
      endif
      re_tm(:) = units_to_atomic(units_out%length, re_tm(:))
      ! ff, the last column, is a dimensionless number

      if(ios < 0) then
        exit ! end of file
      else if(ios > 0) then
        message(1) = "Error parsing file " // trim(fname)
        call messages_fatal(1)
      endif

      energy = units_to_atomic(units_out%energy, energy)
      ! transition matrix elements by themselves are dependent on gauge in degenerate subspaces
      ! make into oscillator strengths, as in casida_inc.F90 X(oscillator_strengths), and like the last column
      tm_sq(1:cs%space%dim) = (matmul(rotation(1:cs%space%dim, 1:cs%space%dim), re_tm(1:cs%space%dim)))**2
      if(is_complex) tm_sq(1:cs%space%dim) = tm_sq(1:cs%space%dim) + &
        (matmul(rotation(1:cs%space%dim, 1:cs%space%dim), im_tm(1:cs%space%dim)))**2
      ff(1:cs%space%dim) = M_TWO * energy * tm_sq(1:cs%space%dim)

      do istep = 1, nsteps
        omega = cs%min_energy + real(istep-1, REAL_PRECISION)*cs%energy_step
        spectrum(1:cs%space%dim+1, istep) = spectrum(1:cs%space%dim+1, istep) + &
          ff(1:cs%space%dim+1)*cs%br/((omega-energy)**2 + cs%br**2)/M_PI ! Lorentzian
      end do
    end do
    call io_close(iunit)

    ! print spectra
    iunit = io_open(trim(dir)//"/spectrum."//fname, action='write')

    write(iunit, '(a2,a12)', advance = 'no') '# ', 'E [' // trim(units_abbrev(units_out%energy)) // ']'
    do idir = 1, cs%space%dim
      write(iunit, '(a14)', advance = 'no') '<' // index2axis(idir) // index2axis(idir) // '>'
    enddo
    write(iunit, '(a14)') '<f>'

    do istep = 1, nsteps
      write(iunit, '(99es14.6)') units_from_atomic(units_out%energy, cs%min_energy + real(istep - 1, REAL_PRECISION) &
        *cs%energy_step), (units_from_atomic(unit_one/units_out%energy, spectrum(idir, istep)), idir = 1, cs%space%dim+1)
    end do

    call io_close(iunit)

    SAFE_DEALLOCATE_A(spectrum)
  end subroutine calc_broad

end program casida_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

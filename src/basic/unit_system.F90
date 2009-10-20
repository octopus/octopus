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

!
! Atomic weights should be read in "atomic mass units" (u) (not to
! be confused with mass in "atomic units"), that is, it should be given
! the relative atomic weight). 1 u is roughly the mass of the proton,
! and exactly one twelfth of mass of the ^{12}C isotope. The relation of the
! atomic mass unit and the atomic unit of mass, au_[mass], is:
!
! 1 au_[mass] = 5.485799110e-4 u
!
! The atomic unit of mass is the mass of the electron. Unfortunately, the
! code uses units of mass of (eV/A^2)(h/(2pieV))^2, which are related to
! atomic units through 1 cu_[mass] = 7.619963358 au_[mass] . So:
!
! 1 u = (1/5.485799110e-4) au_[mass] = (1/5.485799110e-4) *
!      (1/7.619963358) cu_[mass] = 239.225360 cu_[mass].

#include "global.h"

module unit_system_m
  use datasets_m
  use global_m
  use parser_m
  use messages_m
  use io_m
  use unit_m
  use varinfo_m

  implicit none

  private
  public ::                  &
    unit_system_t,           &
    unit_system_init,        &
    unit_system_get,         &
    unit_system_from_file

  type unit_system_t
    type(unit_t) :: length
    type(unit_t) :: energy
    type(unit_t) :: time
    type(unit_t) :: velocity
    type(unit_t) :: mass
    type(unit_t) :: force
    type(unit_t) :: acceleration
    type(unit_t) :: polarizability
    type(unit_t) :: hyperpolarizability
  end type unit_system_t

  type(unit_t),        public :: unit_one, unit_debye, unit_invcm, unit_ppm, unit_susc_ppm_cgs
  type(unit_system_t), public :: units_inp, units_out

  integer, parameter, public :: UNITS_ATOMIC = 1, UNITS_EVA = 2

contains


  ! ---------------------------------------------------------
  subroutine unit_system_init()
    integer :: c, cinp, cout

    call push_sub('unit_system.unit_system_init')

    !%Variable Units
    !%Type integer
    !%Default atomic
    !%Section Execution::Units
    !%Description
    !% Atomic units seem to be the preferred system in the atomic and
    !% molecular physics community. Internally, the code works in
    !% atomic units. However, for input or output, some people like
    !% to use a system based in eV for energies and <math>\AA</math>
    !% for length. The default is atomic units.
    !%
    !% Warning 1: All files read on input will also be treated using
    !% these units, including XYZ geometry files.
    !%
    !% Warning 2: Some values are treated in their most common units,
    !% for example atomic masses (a.m.u.), vibrational frequencies
    !% (cm<sup>-1</sup>) or temperatures (Kelvin). The unit of charge is always
    !% the electronic charge <i>e</i>.
    !%
    !%Option atomic        1
    !% Atomic units
    !%Option ev_angstrom   2
    !% Electronvolts for energy, Angstrom for length, the rest of the
    !% units are derived from these and <math>hbar=1</math>.
    !%End

    !%Variable UnitsInput
    !%Type integer
    !%Default atomic
    !%Section Execution::Units
    !%Description
    !% Same as "Units", but only refers to input values.
    !%Option atomic        1
    !% Atomic units
    !%Option ev_angstrom   2
    !% Electronvolts for energy, Angstroms for length, and the rest of the
    !% units are derived from these and <math>hbar=1</math>.
    !%End

    !%Variable UnitsOutput
    !%Type integer
    !%Default atomic
    !%Section Execution::Units
    !%Description
    !% Same as "Units", but only refers to output values.
    !%Option atomic        1
    !% Atomic units
    !%Option ev_angstrom   2
    !% Electronvolts for energy, Angstrom for length, the rest of the
    !% units are derived from these and <math>hbar=1</math>.
    !%End

    if(parse_isdef(datasets_check('Units')).ne.0) then
      call parse_integer(datasets_check('Units'), UNITS_ATOMIC, c)
      if(.not.varinfo_valid_option('Units', c)) call input_error('Units')
      cinp = c
      cout = c
    else
      call parse_integer(datasets_check('UnitsInput'), UNITS_ATOMIC, c)
      if(.not.varinfo_valid_option('UnitsInput', c)) call input_error('UnitsInput')
      cinp = c
      call parse_integer(datasets_check('UnitsOutput'), UNITS_ATOMIC, c)
      if(.not.varinfo_valid_option('UnitsOutput', c)) call input_error('UnitsOutput')
      cout = c
    end if

    unit_one%factor = M_ONE
    unit_one%abbrev = '1'
    unit_one%name   = 'one'

    unit_ppm%factor = CNST(1e-6)
    unit_ppm%abbrev = 'ppm a.u.'
    unit_ppm%name   = 'parts per million atomic units'

    unit_susc_ppm_cgs%factor = CNST(8.9238878e-2)*CNST(1e-6)
    unit_susc_ppm_cgs%abbrev = 'ppm cgs/mol'
    unit_susc_ppm_cgs%name   = 'magnetic susceptibility parts per million atomic units'

    unit_debye%factor = M_ONE/CNST(2.5417462)
    unit_debye%abbrev = 'Debye'
    unit_debye%name   = 'Debye'

    unit_invcm%factor = M_ONE/CNST(219474.63)
    unit_invcm%abbrev = 'cm^-1'
    unit_invcm%name   = 'h times c over centimeters'

    call unit_system_get(units_inp, cinp)
    call unit_system_get(units_out, cout)

    call pop_sub()

  end subroutine unit_system_init

  ! ---------------------------------------------------------
  subroutine unit_system_get(u, c)
    type(unit_system_t), intent(out) :: u
    integer,             intent(in)  :: c

    select case(c)
    case (UNITS_ATOMIC)
      call unit_system_init_atomic(u)
    case (UNITS_EVA)
      call unit_system_init_eV_Ang(u)
    case default
      call input_error('Units')
    end select
  end subroutine unit_system_get


  ! these routines output the unit conversions factors, defined by
  ! [a.u.] = <input>*u.unit
  ! <output> = [a.u.]/u.unit

  ! ---------------------------------------------------------
  subroutine unit_system_init_atomic(u)
    type(unit_system_t), intent(out) :: u

    u%length%abbrev = "b"
    u%length%name   = "Bohr"
    u%length%factor = M_ONE

    u%energy%abbrev = "H"
    u%energy%name   = "Hartree"
    u%energy%factor = M_ONE

    u%time%abbrev = "hbar/H"
    u%time%name   = "hbar/Hartree"
    u%time%factor = M_ONE/u%energy%factor

    u%velocity%abbrev = "bH(2pi/h)"
    u%velocity%name   = "Bohr times Hartree over hbar"
    u%velocity%factor = M_ONE

    u%mass%abbrev   = "u"
    u%mass%name     = "1/12 of the mass of C^12"
    u%mass%factor   = M_ONE/CNST(5.485799110e-4)

    u%force%abbrev  = "H/b"
    u%force%name    = "Hartree/Bohr"
    u%force%factor  = M_ONE

    u%acceleration%abbrev = "bH(2pi/h)^2"
    u%acceleration%name   = "Bohr times (Hartree over h bar) squared"
    u%acceleration%factor = M_ONE

    u%polarizability%abbrev  = "b^3"
    u%polarizability%name    = "Bohr^3"
    u%polarizability%factor  = M_ONE
    ! By convention, this unit appears more commonly than the
    ! equivalent b^2/H. It does not depend on the dimensionality
    ! of the system, despite analogies to a volume.

    u%hyperpolarizability%abbrev  = "b^5"
    u%hyperpolarizability%name    = "Bohr^5"
    u%hyperpolarizability%factor  = M_ONE
  end subroutine unit_system_init_atomic


  ! ---------------------------------------------------------
  subroutine unit_system_init_eV_Ang(u)
    type(unit_system_t), intent(out) :: u

    u%length%abbrev = "A"
    u%length%name   = "Angstrom"
    u%length%factor = P_Ang  ! 1 a.u. = 0.529 A

    u%energy%abbrev = "eV"
    u%energy%name   = "electronvolt"
    u%energy%factor = M_ONE/(M_TWO*P_Ry)   ! 1 a.u. = 27.2 eV

    u%time%abbrev = "hbar/eV"
    u%time%name   = "hbar/electronvolt"
    u%time%factor = M_ONE/u%energy%factor

    u%velocity%abbrev = "AeV(2pi/h)"
    u%velocity%name   = "Angstrom times electronvolts over hbar"
    u%velocity%factor = u%length%factor*u%energy%factor

    u%mass%abbrev   = "u"
    u%mass%name     = "1/12 of the mass of C^12"
    u%mass%factor   = M_ONE/CNST(5.485799110e-4)

    u%force%abbrev  = "eV/A"
    u%force%name    = "electronvolt/Angstrom"
    u%force%factor  = u%energy%factor/u%length%factor

    u%acceleration%abbrev = "AeV(2pi/h)^2"
    u%acceleration%name   = "Angstrom times (electronvolt over hbar) squared"
    u%acceleration%factor = u%length%factor/u%time%factor**2

    u%polarizability      = u%length**2 / u%energy
    u%hyperpolarizability = u%length**3 / u%energy**2
  end subroutine unit_system_init_eV_Ang


  ! ---------------------------------------------------------
  ! This is a very primitive procedure that attempts to find out
  ! which units were used to write one octopus file, be it a
  ! "multipoles", a "cross_section_tensor", etc. 
  ! TODO: although it seems to work in most cases, it is obviously
  ! a very weak code.
  ! ---------------------------------------------------------
  subroutine unit_system_from_file(u, fname, ierr)
    type(unit_system_t), intent(inout) :: u
    character(len=*),    intent(in)    :: fname
    integer,             intent(inout) :: ierr

    integer            :: iunit, ios
    character(len=256) :: line

    call push_sub('unit_system.unit_system_from_file')

    iunit = io_open(file = trim(fname), action = 'read', status = 'old', die = .false.)
    if(iunit < 0) then
      ierr = -2
      call pop_sub(); return
    end if

    ierr = 0
    do
      read(iunit, '(a)', iostat = ios) line
      if(ios.ne.0) exit
      if(index(line,'[A]').ne.0  .or.  index(line,'eV').ne.0) then
        call unit_system_get(u, UNITS_EVA)
        call pop_sub(); return
      elseif(index(line,'[b]').ne.0) then
        call unit_system_get(u, UNITS_ATOMIC)
        call pop_sub(); return
      end if
    end do

    ierr = -1

    call pop_sub()
  end subroutine unit_system_from_file

end module unit_system_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

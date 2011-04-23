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
!
!> Atomic weights should be read in "atomic mass units" (u) (not to
!! be confused with mass in "atomic units"), that is, it should be given
!! the relative atomic weight). 1 u is roughly the mass of the proton,
!! and exactly one twelfth of mass of the ^{12}C isotope. The relation of the
!! atomic mass unit and the atomic unit of mass, au_[mass], is:
!!
!! 1 au_[mass] = 5.485799110e-4 u
!!
module unit_system_m
  use datasets_m
  use global_m
  use io_m
  use messages_m
  use parser_m
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

  ! the units systems for reading and writing
  type(unit_system_t), public :: units_inp, units_out

  ! some special units required for particular quantities
  type(unit_t),        public :: unit_one           !< For unitless quantities and arithmetics with units.
  type(unit_t),        public :: unit_ppm           !< Parts per million.
  type(unit_t),        public :: unit_debye         !< For dipoles.
  type(unit_t),        public :: unit_invcm         !< For vibrational frequencies.
  type(unit_t),        public :: unit_susc_ppm_cgs  !< Some magnetic stuff.
  type(unit_t),        public :: unit_kelvin        !< For converting energies into temperatures.
  type(unit_t),        public :: unit_femtosecond   !< Time in femtoseconds.

  integer, parameter, public :: UNITS_ATOMIC = 0, UNITS_EVA = 1, UNITS_FS = 2

contains


  ! ---------------------------------------------------------
  subroutine unit_system_init()
    integer :: cc, cinp, cout

    PUSH_SUB(unit_system_init)

    !%Variable Units
    !%Type integer
    !%Default atomic
    !%Section Execution::Units
    !%Description
    !% This variable selects the units that Octopus use for
    !% input and output.
    !%
    !% Atomic units seem to be the preferred system in the atomic and
    !% molecular physics community. Internally, the code works in
    !% atomic units. However, for input or output, some people like
    !% to use a system based in eV for energies and <math>\AA</math>
    !% for length. The default is atomic units.
    !%
    !% Normally time units are derived from energy and length units,
    !% so it is measured in <math>\hbar/Hartree</math> or
    !% <math>\hbar/electronvolt</math>. Alternatively you can tell
    !% Octopus to use femtoseconds as the time unit by adding the
    !% value <tt>femtoseconds</tt> (Note that no other unit will be 
    !% based on femtoseconds). So for example you can use:
    !%
    !% <tt>Units = femtoseconds</tt>
    !%
    !% or
    !%
    !% <tt>Units = ev_angstrom + femtoseconds</tt>
    !%
    !% You can use different unit systems for input and output by
    !% setting the <tt>UnitsInput</tt> and <tt>UnitsOutput</tt>.
    !%
    !% Warning 1: All files read on input will also be treated using
    !% these units, including XYZ geometry files.
    !%
    !% Warning 2: Some values are treated in their most common units,
    !% for example atomic masses (a.m.u.), electron effective masses
    !% (electron mass), vibrational frequencies
    !% (cm<sup>-1</sup>) or temperatures (Kelvin). The unit of charge is always
    !% the electronic charge <i>e</i>.
    !%
    !%Option atomic        0
    !% Atomic units.
    !%Option ev_angstrom   1
    !% Electronvolts for energy, Angstroms for length, the rest of the
    !% units are derived from these and <math>hbar=1</math>.
    !%Option femtoseconds  2
    !% (Experimental) If you add this value to the other options,
    !% Octopus will treat time in femtoseconds units.
    !%End

    !%Variable UnitsInput
    !%Type integer
    !%Default atomic
    !%Section Execution::Units
    !%Description
    !% Same as <tt>Units</tt>, but only refers to input values.
    !%End

    !%Variable UnitsOutput
    !%Type integer
    !%Default atomic
    !%Section Execution::Units
    !%Description
    !% Same as <tt>Units</tt>, but only refers to output values.
    !%End

    if(parse_isdef(datasets_check('Units')).ne.0) then
      call parse_integer(datasets_check('Units'), UNITS_ATOMIC, cc)
      if(.not.varinfo_valid_option('Units', cc, is_flag = .true.)) call input_error('Units')
      cinp = cc
      cout = cc
    else
      ! note that we check the value is valid for the 'Units' variable
      call parse_integer(datasets_check('UnitsInput'), UNITS_ATOMIC, cc)
      if(.not.varinfo_valid_option('Units', cc, is_flag = .true.)) call input_error('UnitsInput')
      cinp = cc
      call parse_integer(datasets_check('UnitsOutput'), UNITS_ATOMIC, cc)
      if(.not.varinfo_valid_option('Units', cc, is_flag = .true.)) call input_error('UnitsOutput')
      cout = cc
    end if

    unit_one%factor = M_ONE
    unit_one%abbrev = '1'
    unit_one%name   = 'one'

    unit_ppm%factor = CNST(1e-6)
    unit_ppm%abbrev = 'ppm a.u.'
    unit_ppm%name   = 'parts per million'

    unit_susc_ppm_cgs%factor = CNST(1e-6)/CNST(8.9238878e-2)
    unit_susc_ppm_cgs%abbrev = 'ppm cgs/mol'
    unit_susc_ppm_cgs%name   = 'magnetic susceptibility parts per million cgs'

    unit_debye%factor = M_ONE/CNST(2.5417462)
    unit_debye%abbrev = 'Debye'
    unit_debye%name   = 'Debye'

    unit_invcm%factor = M_ONE/CNST(219474.63)
    unit_invcm%abbrev = 'cm^-1'
    unit_invcm%name   = 'h times c over centimeters'

    unit_kelvin%factor = P_KB
    unit_kelvin%abbrev = 'K'
    unit_kelvin%name   = 'degrees Kelvin'

    unit_femtosecond%factor = CNST(1.0)/CNST(0.024188843)
    unit_femtosecond%abbrev = 'fs'
    unit_femtosecond%name   = 'femtoseconds'

    call unit_system_get(units_inp, mod(cinp, 2))
    call unit_system_get(units_out, mod(cout, 2))

    if(cinp/2 == 1 .or. cout/2 == 1) then
      call messages_experimental('Femtosecond units')
    end if

    if(cinp/2 == 1) units_inp%time = unit_femtosecond
    if(cout/2 == 1) units_out%time = unit_femtosecond

    POP_SUB(unit_system_init)

  end subroutine unit_system_init

  ! ---------------------------------------------------------
  subroutine unit_system_get(uu, cc)
    type(unit_system_t), intent(out) :: uu
    integer,             intent(in)  :: cc

    PUSH_SUB(unit_system_get)

    select case(cc)
    case (UNITS_ATOMIC)
      call unit_system_init_atomic(uu)
    case (UNITS_EVA)
      call unit_system_init_eV_Ang(uu)
    case default
      call input_error('Units')
    end select

    POP_SUB(unit_system_get)
  end subroutine unit_system_get


  ! ---------------------------------------------------------
  !> These routines output the unit-conversion factors, defined by
  !! [a.u.] = <input>*u.unit
  !! <output> = [a.u.]/u.unit
  ! ---------------------------------------------------------
  subroutine unit_system_init_atomic(uu)
    type(unit_system_t), intent(out) :: uu

    PUSH_SUB(unit_system_init_atomic)

    uu%length%abbrev = "b"
    uu%length%name   = "Bohr"
    uu%length%factor = M_ONE

    uu%energy%abbrev = "H"
    uu%energy%name   = "Hartree"
    uu%energy%factor = M_ONE

    uu%time%abbrev = "hbar/H"
    uu%time%name   = "hbar/Hartree"
    uu%time%factor = M_ONE/uu%energy%factor

    uu%velocity%abbrev = "bH(2pi/h)"
    uu%velocity%name   = "Bohr times Hartree over hbar"
    uu%velocity%factor = M_ONE

    uu%mass%abbrev   = "u"
    uu%mass%name     = "1/12 of the mass of C^12"
    uu%mass%factor   = M_ONE/CNST(5.485799110e-4)

    uu%force%abbrev  = "H/b"
    uu%force%name    = "Hartree/Bohr"
    uu%force%factor  = M_ONE

    uu%acceleration%abbrev = "bH(2pi/h)^2"
    uu%acceleration%name   = "Bohr times (Hartree over h bar) squared"
    uu%acceleration%factor = M_ONE

    uu%polarizability%abbrev  = "b^3"
    uu%polarizability%name    = "Bohr^3"
    uu%polarizability%factor  = M_ONE
    ! By convention, this unit appears more commonly than the
    ! equivalent b^2/H. It does not depend on the dimensionality
    ! of the system, despite analogies to a volume.

    uu%hyperpolarizability%abbrev  = "b^5"
    uu%hyperpolarizability%name    = "Bohr^5"
    uu%hyperpolarizability%factor  = M_ONE

    POP_SUB(unit_system_init_atomic)
  end subroutine unit_system_init_atomic


  ! ---------------------------------------------------------
  subroutine unit_system_init_eV_Ang(uu)
    type(unit_system_t), intent(out) :: uu

    PUSH_SUB(unit_system_init_eV_Ang)

    uu%length%abbrev = "A"
    uu%length%name   = "Angstrom"
    uu%length%factor = P_Ang  ! 1 a.u. = 0.529 A

    uu%energy%abbrev = "eV"
    uu%energy%name   = "electronvolt"
    uu%energy%factor = M_ONE/(M_TWO*P_Ry)   ! 1 a.u. = 27.2 eV

    uu%time%abbrev = "hbar/eV"
    uu%time%name   = "hbar/electronvolt"
    uu%time%factor = M_ONE/uu%energy%factor

    uu%velocity%abbrev = "AeV(2pi/h)"
    uu%velocity%name   = "Angstrom times electronvolts over hbar"
    uu%velocity%factor = uu%length%factor*uu%energy%factor

    uu%mass%abbrev   = "u"
    uu%mass%name     = "1/12 of the mass of C^12"
    uu%mass%factor   = M_ONE/CNST(5.485799110e-4)

    uu%force%abbrev  = "eV/A"
    uu%force%name    = "electronvolt/Angstrom"
    uu%force%factor  = uu%energy%factor/uu%length%factor

    uu%acceleration%abbrev = "AeV(2pi/h)^2"
    uu%acceleration%name   = "Angstrom times (electronvolt over hbar) squared"
    uu%acceleration%factor = uu%length%factor/uu%time%factor**2

    uu%polarizability      = uu%length**3
    uu%hyperpolarizability = uu%length**5

    POP_SUB(unit_system_init_eV_Ang)
  end subroutine unit_system_init_eV_Ang


  ! ---------------------------------------------------------
  !> This is a very primitive procedure that attempts to find out
  !! which units were used to write an octopus file, whether
  !! "multipoles", "cross_section_tensor", etc.
  !! \todo  Although it seems to work in most cases, it is obviously
  !! a very weak code.
  ! ---------------------------------------------------------
  subroutine unit_system_from_file(uu, fname, ierr)
    type(unit_system_t), intent(inout) :: uu
    character(len=*),    intent(in)    :: fname
    integer,             intent(inout) :: ierr

    integer            :: iunit, ios
    character(len=256) :: line

    PUSH_SUB(unit_system_from_file)

    iunit = io_open(file = trim(fname), action = 'read', status = 'old', die = .false.)
    if(iunit < 0) then
      ierr = -2
      POP_SUB(unit_system_from_file)
      return
    end if

    ierr = 0
    do
      read(iunit, '(a)', iostat = ios) line
      if(ios.ne.0) exit
      if(index(line,'[A]').ne.0  .or.  index(line,'eV').ne.0) then
        call unit_system_get(uu, UNITS_EVA)
        POP_SUB(unit_system_from_file)
        return
      elseif(index(line,'[b]').ne.0) then
        call unit_system_get(uu, UNITS_ATOMIC)
        POP_SUB(unit_system_from_file)
        return
      end if
    end do

    ierr = -1

    POP_SUB(unit_system_from_file)
  end subroutine unit_system_from_file

end module unit_system_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

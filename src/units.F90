!***********************************************************************
!                                Comments
!***********************************************************************
! Atomic weights should be read in [atomic mass units] (u) (not to mistake
! with [atomic units]), that is, it should be given the relative
! atomic weight). 1u is roughly the mass of the proton, and exactly
! one twelveth of the ^{12}C isotope. The relation of the atomic mass
! unit and the atomic unit of mass, au_[mass], is:
!
! 1 au_[mass] = 5.485799110e-4 u
!
! The atomic unit of mass is the mass of the electron. Unfortunately, the
! code uses units of mass of (eV/A^2)(h/(2pieV))^2, which are related to
! atomic units through 1cu_[mass] = 7.619963358 au_[mass] . So:
!
! 1u = (1/5.485799110e-4) au_[mass] = (1/5.485799110e-4) * 
!      (1/7.619963358) cu_[mass] = 239.225360 cu_[mass].


module units
use global
use fdf

implicit none

type unit_type
  real(r8) :: factor
  character(len=10) :: abbrev ! common abbreviation of the unit name
  character(len=30) :: name   ! common name
end type unit_type

type unit_system_type
  type(unit_type) :: length
  type(unit_type) :: energy
  type(unit_type) :: time
  type(unit_type) :: velocity
  type(unit_type) :: mass
end type unit_system_type

type(unit_system_type) :: units_inp, units_out

contains

subroutine units_init()
  character(len=10) :: c
  character(len=3) :: cinp, cout

  sub_name = 'units_init'; call push_sub()

  if(fdf_defined("Units")) then
    c = fdf_string("Units", "")
    cinp = c(1:3)
    cout = c(1:3)
  else
    c = fdf_string("UnitsInput", "a.u")
    cinp = c(1:3)
    c = fdf_string("UnitsOutput", "a.u")
    cout = c(1:3)
  end if

  call get_units(units_inp, cinp)
  call get_units(units_out, cout)

  call pop_sub()

contains

  subroutine get_units(u, c)
    type(unit_system_type), intent(out) :: u
    character(len=3) :: c
    
    select case(c)
    case ("a.u")
      call units_atomic(u)
    case ("eVA")
      call units_eV_Ang(u)
    case default
      message(1) = "Invalid unit specification: '"//trim(c)//"'"
      message(2) = "Valid units are: 'a.u', 'eVA'"
      call write_fatal(2)
    end select
  end subroutine get_units
  
end subroutine units_init

! these routines output the unit conversions factors, defined by
! [a.u.] = <input>*u.unit
! <output> = [a.u.]/u.unit

subroutine units_atomic(u)
  type(unit_system_type), intent(out) :: u

  u%length%abbrev = "b"
  u%length%name   = "bohr"
  u%length%factor = 1._r8

  u%energy%abbrev = "H"
  u%energy%name   = "Hartree"
  u%energy%factor = 1._r8

  u%time%abbrev = "hbar/H"
  u%time%name   = "hbar/Hartree"
  u%time%factor = 1._r8/u%energy%factor

  u%velocity%abbrev = "b H (2 pi/h)"
  u%velocity%name   = "bohr times Hartree over h bar"
  u%velocity%factor = 1._r8

  u%mass%abbrev   = "u"
  u%mass%name     = "1/12 of the mass of C^12"
  u%mass%factor   = 1._r8/5.485799110e-4_r8
end subroutine units_atomic

subroutine units_eV_Ang(u)
  type(unit_system_type), intent(out) :: u

  u%length%abbrev = "A"
  u%length%name   = "Angstrom"
  u%length%factor = P_Ang  ! 1 a.u. = 0.529 A

  u%energy%abbrev = "eV"
  u%energy%name   = "electron volt"
  u%energy%factor = 1._r8/(2._r8*P_Ry)   ! 1 a.u. = 27.2 eV

  u%time%abbrev = "hbar/eV"
  u%time%name   = "hbar/electron volt"
  u%time%factor = 1._r8/u%energy%factor

  u%velocity%abbrev = "A eV (2 pi/h)"
  u%velocity%name   = "Angstrom times electron volts over h bar"
  u%velocity%factor = u%length%factor*u%energy%factor

  u%mass%abbrev   = "u"
  u%mass%name     = "1/12 of the mass of C^12"
  u%mass%factor   = 1._r8/5.485799110e-4_r8
end subroutine units_eV_Ang

end module units

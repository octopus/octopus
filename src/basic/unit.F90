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

module unit_m
  use global_m

  implicit none

  private
  public ::            &
    unit_t,            &
    units_to_atomic,   &
    units_from_atomic, &
    units_abbrev,      &
    operator(*),       &
    operator(/),       &
    operator(**),      &
    sqrt

  type unit_t
    real(8)           :: factor
    character(len=12) :: abbrev ! common abbreviation of the unit name
    character(len=50) :: name   ! common name
  end type unit_t

  interface operator (*)
    module procedure units_multiply
  end interface

  interface operator (/)
    module procedure units_divide
  end interface

  interface operator (**)
    module procedure units_pow
  end interface

  interface units_to_atomic
    module procedure dunits_to_atomic, zunits_to_atomic, &
      dunits_to_atomic_4, zunits_to_atomic_4
  end interface

  interface units_from_atomic
    module procedure dunits_from_atomic, zunits_from_atomic, &
      dunits_from_atomic_4, zunits_from_atomic_4
  end interface
  
  interface sqrt
    module procedure units_sqrt
  end interface

contains

  !-----------------------------------------------

  real(8) elemental pure function dunits_to_atomic(this, val) result(res)
    type(unit_t), intent(in) :: this
    real(8),      intent(in) :: val

    res = val*this%factor

  end function dunits_to_atomic
 
  !-----------------------------------------------

  complex(8) elemental pure function zunits_to_atomic(this, val) result(res)
    type(unit_t), intent(in) :: this
    complex(8),   intent(in) :: val

    res = val*this%factor

  end function zunits_to_atomic

  !-----------------------------------------------
  
  real(8) elemental pure function dunits_from_atomic(this, val) result(res)
    type(unit_t), intent(in) :: this
    real(8),      intent(in) :: val

    res = val/this%factor

  end function dunits_from_atomic

  !-----------------------------------------------
  
  complex(8) elemental pure function zunits_from_atomic(this, val) result(res)
    type(unit_t), intent(in) :: this
    complex(8),   intent(in) :: val

    res = val/this%factor

  end function zunits_from_atomic

  !-----------------------------------------------
  ! now the single-precision functions
  !-----------------------------------------------

  real(4) elemental pure function dunits_to_atomic_4(this, val) result(res)
    type(unit_t), intent(in) :: this
    real(4),      intent(in) :: val

    res = val*this%factor

  end function dunits_to_atomic_4
 
  !-----------------------------------------------

  complex(4) elemental pure function zunits_to_atomic_4(this, val) result(res)
    type(unit_t), intent(in) :: this
    complex(4),   intent(in) :: val

    res = val*this%factor

  end function zunits_to_atomic_4

  !-----------------------------------------------
  
  real(4) elemental pure function dunits_from_atomic_4(this, val) result(res)
    type(unit_t), intent(in) :: this
    real(4),      intent(in) :: val

    res = val/this%factor

  end function dunits_from_atomic_4

  !-----------------------------------------------
  
  complex(4) elemental pure function zunits_from_atomic_4(this, val) result(res)
    type(unit_t), intent(in) :: this
    complex(4),   intent(in) :: val

    res = val/this%factor

  end function zunits_from_atomic_4

  !-----------------------------------------------


  character(len=12) pure function units_abbrev(this) result(abbrev)
    type(unit_t), intent(in) :: this
    
    abbrev = this%abbrev
  end function units_abbrev

  !-----------------------------------------------

  type(unit_t) pure function units_multiply(aa, bb) result(cc)
    type(unit_t), intent(in) :: aa
    type(unit_t), intent(in) :: bb

    cc%factor = aa%factor*bb%factor
    cc%abbrev = trim(aa%abbrev)//'*'//trim(bb%abbrev)

  end function units_multiply

  !-----------------------------------------------

  type(unit_t) pure function units_divide(aa, bb) result(cc)
    type(unit_t), intent(in) :: aa
    type(unit_t), intent(in) :: bb

    cc%factor = aa%factor/bb%factor
    cc%abbrev = trim(aa%abbrev)//'/'//trim(bb%abbrev)

  end function units_divide
  !-----------------------------------------------

  type(unit_t) pure function units_pow(aa, nn) result(cc)
    type(unit_t), intent(in) :: aa
    integer,      intent(in) :: nn

    cc%factor = aa%factor**nn

    ! We have to do the conversion by hand. This is ugly, but we
    ! cannot use write here since this function might be called inside
    ! another write (stupid Fortran).

    select case(nn)
    case(-3)
      cc%abbrev = trim(aa%abbrev)//'^-3'
    case(-2)
      cc%abbrev = trim(aa%abbrev)//'^-2'
    case(-1)
      cc%abbrev = trim(aa%abbrev)//'^-1'
    case(0)
      cc%abbrev = '1'
    case(1)
      cc%abbrev = trim(aa%abbrev)
    case(2)
      cc%abbrev = trim(aa%abbrev)//'^2'
    case(3)
      cc%abbrev = trim(aa%abbrev)//'^3'
    case(4)
      cc%abbrev = trim(aa%abbrev)//'^4'
    case(5)
      cc%abbrev = trim(aa%abbrev)//'^5'
    case default
      cc%abbrev = trim(aa%abbrev)//'^n'
    end select

  end function units_pow


  !-----------------------------------------------

  type(unit_t) pure function units_sqrt(aa) result(cc)
    type(unit_t), intent(in) :: aa
    
    cc%factor = sqrt(aa%factor)
    cc%abbrev = 'sqrt('//trim(aa%abbrev)//')'

  end function units_sqrt

end module unit_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

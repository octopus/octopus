!! Copyright (C) 2010 X. Andrade
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
!! $Id: calc_mode.F90 6441 2010-04-03 23:00:15Z xavier $

#include "global.h"

module types_m
  
  implicit none
  
  private

  public ::            &
    type_t,            &
    types_get_size,    &
    operator(==),      &
    operator(/=)

  type type_t
    integer :: itype
  end type type_t

  type(type_t), public :: TYPE_FLOAT   = type_t(1)
  type(type_t), public :: TYPE_CMPLX   = type_t(2)
  type(type_t), public :: TYPE_INTEGER = type_t(3)
  type(type_t), public :: TYPE_BYTE    = type_t(4)

  interface operator(==)
    module procedure types_equal
  end interface operator(==)

  interface operator(/=)
    module procedure types_not_equal
  end interface operator(/=)

#ifdef SINGLE_PRECISION
  integer :: sizes(4) = (/4, 8, 4, 1/)
#else
  integer :: sizes(4) = (/8, 16, 4, 1/)
#endif
  
contains
  
  integer pure function types_get_size(this) result(size)
    type(type_t), intent(in) :: this
    
    size = sizes(this%itype)
  end function types_get_size

  ! -----------------------------------------------------

  logical pure function types_equal(ta, tb) result(equal)
    type(type_t), intent(in) :: ta
    type(type_t), intent(in) :: tb
    
    equal = ta%itype == tb%itype

  end function types_equal
  
  ! -----------------------------------------------------

  logical pure function types_not_equal(ta, tb) result(equal)
    type(type_t), intent(in) :: ta
    type(type_t), intent(in) :: tb
    
    equal = ta%itype /= tb%itype

  end function types_not_equal
end module types_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

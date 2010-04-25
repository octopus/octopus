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
    types_get_size

  integer, public, parameter ::   &
    TYPE_FLOAT            =   1,  &
    TYPE_CMPLX            =   2,  &
    TYPE_INTEGER          =   3
  
#ifdef SINGLE_PRECISION
  integer :: sizes(3) = (/4, 8, 4/)
#else
  integer :: sizes(3) = (/8, 16, 4/)
#endif

  contains

  integer function types_get_size(type) result(size)
    integer, intent(in) :: type
    
    size = sizes(type)
  end function types_get_size

end module types_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

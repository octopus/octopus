!! Copyright (C) 2020 N. Tancogne-Dejean
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

#include "global.h"

module observable_oct_m
  use clock_oct_m
  implicit none

  private
  public ::                   &
     observable_t 

  integer, public, parameter ::         &
    POSITION                     =  1,  &
    VELOCITY                     =  2,  &
    CURRENT                      =  3,  &
    DENSITY                      =  4,  &
    SCALAR_POTENTIAL             =  5,  &
    VECTOR_POTENTIAL             =  6,  &
    E_FIELD                      =  7,  &
    B_FIELD                      =  8,  &
    MAX_OBSERVABLES              =  8


  type observable_t
    private

    type(clock_t), public :: clock
  end type observable_t

contains

end module observable_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

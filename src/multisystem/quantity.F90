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

module quantity_oct_m
  use clock_oct_m
  implicit none

  private
  public ::                   &
     quantity_t

  integer, public, parameter ::         &
    POSITION                     =  1,  &
    VELOCITY                     =  2,  &
    CURRENT                      =  3,  &
    DENSITY                      =  4,  &
    SCALAR_POTENTIAL             =  5,  &
    VECTOR_POTENTIAL             =  6,  &
    E_FIELD                      =  7,  &
    B_FIELD                      =  8,  &
    MASS                         =  9,  &
    CHARGE                       = 10,  &
    PERMITTIVITY                 = 11,  &
    PERMEABILITY                 = 12,  &
    E_CONDUCTIVITY               = 13,  &
    M_CONDUCTIVITY               = 14,  &
    MAX_QUANTITIES               = 14

  character(len=17), public, parameter :: QUANTITY_LABEL(MAX_QUANTITIES) = (/ &
    "position        ", &
    "velocity        ", &
    "current         ", &
    "density         ", &
    "scalar potential", &
    "vector potential", &
    "E field         ", &
    "B field         ", &
    "mass            ", &
    "charge          ", &
    "permittivity    ", &
    "permeability    ", &
    "e_conductivity  ", &
    "m_conductivity  "  &
    /)

  !> Systems can expose quantities that can be used to calculate interactions
  !! with other systems.
  !!
  !! Some quantities are dynamical variables of the system. Such quantities are
  !! usually updated by the propagation algorithm and cannot be calculated
  !! on-demand. Such quantities must be marked as "protected".
  type quantity_t
    private
    type(clock_t), public :: clock               !< Clock storing the time at which the quantity was last updated.
    logical,       public :: required = .false.  !< Should this quantities be calculated?
    logical,       public :: protected = .false. !< Is this quantity protected, i.e., it cannot be updated on-demand
  end type quantity_t

contains

end module quantity_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

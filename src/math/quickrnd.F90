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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#include "global.h"

!> This module is intended to contain "only mathematical" functions
!! and procedures.
module quickrnd_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                     &
    quickrnd
  
contains

  ! ---------------------------------------------------------
  !> a simple congruent random number generator
  subroutine quickrnd(iseed, rnd)
    integer, intent(inout) :: iseed
    FLOAT,   intent(inout) :: rnd

    integer, parameter :: im=6075, ia=106, ic=1283

    ! no PUSH_SUB, called too often

    iseed = mod(iseed*ia + ic, im)
    rnd = real(iseed, REAL_PRECISION)/real(im, REAL_PRECISION)

  end subroutine quickrnd
  
end module quickrnd_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

!! Copyright (C) 2002-2016 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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

! ---------------------------------------------------------
!> a simple congruent random number generator

module quickrnd_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                     &
    quickrnd,                   &
    shiftseed

  interface quickrnd
    module procedure dquickrnd_single, dquickrnd_array, zquickrnd_array
  end interface quickrnd

contains

  subroutine dquickrnd_single(iseed, rnd)
    integer, intent(inout) :: iseed
    FLOAT,   intent(out)   :: rnd

    integer, parameter :: im=6075, ia=106, ic=1283

    ! no PUSH_SUB, called too often

    iseed = mod(iseed*ia + ic, im)
    rnd = real(iseed, REAL_PRECISION)/real(im, REAL_PRECISION)

  end subroutine dquickrnd_single

  ! ---------------------------------------------------------

  subroutine dquickrnd_array(iseed, nn, rnd)
    integer, intent(inout) :: iseed
    integer, intent(in)    :: nn
    FLOAT,   intent(inout) :: rnd(:)

    integer, parameter :: im=6075, ia=106, ic=1283
    integer :: ii
    
    PUSH_SUB(quickrnd_array)

    do ii = 1, nn
      iseed = mod(iseed*ia + ic, im)
      rnd(ii) = real(iseed, REAL_PRECISION)/real(im, REAL_PRECISION)
    end do
    
    POP_SUB(quickrnd_array)
    
  end subroutine dquickrnd_array

  ! ---------------------------------------------------------

  subroutine zquickrnd_array(iseed, nn, rnd)
    integer, intent(inout) :: iseed
    integer, intent(in)    :: nn
    CMPLX,   intent(inout) :: rnd(:)

    integer, parameter :: im=6075, ia=106, ic=1283
    integer :: ii
    
    PUSH_SUB(quickrnd_array)

    do ii = 1, nn
      iseed = mod(iseed*ia + ic, im)
      rnd(ii) = real(iseed, REAL_PRECISION)/real(im, REAL_PRECISION)
      iseed = mod(iseed*ia + ic, im)
      rnd(ii) = rnd(ii) + M_ZI*real(iseed, REAL_PRECISION)/real(im, REAL_PRECISION)
      rnd(ii) = rnd(ii)/sqrt(CNST(2.0))
    end do
    
    POP_SUB(quickrnd_array)
    
  end subroutine zquickrnd_array

  ! ---------------------------------------------------------

  subroutine shiftseed(iseed, n)
    integer, intent(inout) :: iseed
    integer, intent(in)    :: n

    integer, parameter :: im=6075, ia=106, ic=1283
    integer :: ii
    
    PUSH_SUB(shiftseed)

    do ii = 1, n
      iseed = mod(iseed*ia + ic, im) 
    end do    

    POP_SUB(shiftseed)

  end subroutine shiftseed
 
 
end module quickrnd_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

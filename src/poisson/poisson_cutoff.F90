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
!! $Id: math.F90 3142 2007-08-28 09:42:03Z xavier $

#include "global.h"

module poisson_cutoff_m
  use global_m
  use loct_math_m
  use messages_m

  implicit none

  private
  public ::                       &
    poisson_cutoff_3D_0D,         &
    poisson_cutoff_3D_1D,         &
    poisson_cutoff_3D_1D_finite,  &
    poisson_cutoff_3D_2D,         &
    poisson_cutoff_2D_0D,         &
    poisson_cutoff_2D_1D,         &
    poisson_cutoff_1D_0D,         &
    poisson_cutoff_intcoslog


  interface poisson_cutoff_3D_1D_finite
    real(8) function c_poisson_cutoff_3d_1d_finite(gx, gperp, xsize, rsize)
      real(8), intent(in) :: gx, gperp, rsize, xsize
    end function c_poisson_cutoff_3d_1d_finite
    module procedure poisson_cutoff_3d_1d_finite4
  end interface poisson_cutoff_3D_1D_finite

  interface poisson_cutoff_2D_0D
    real(8) function c_poisson_cutoff_2d_0d(x, y)
      real(8), intent(in) :: x, y
    end function c_poisson_cutoff_2d_0d
    module procedure poisson_cutoff_2d_0d4
  end interface poisson_cutoff_2D_0D

  interface poisson_cutoff_2D_1D
    real(8) function c_poisson_cutoff_2d_1d(gy, gx, r_c)
      real(8), intent(in) :: gy, gx, r_c
    end function c_poisson_cutoff_2d_1d
    module procedure poisson_cutoff_2d_1d4
  end interface poisson_cutoff_2D_1D

  interface poisson_cutoff_1D_0D
    real(8) function c_poisson_cutoff_1d_0d(g, a, r_c)
      real(8), intent(in) :: g, a, r_c
    end function c_poisson_cutoff_1d_0d
    module procedure poisson_cutoff_1d_0d4
  end interface poisson_cutoff_1D_0D

  interface poisson_cutoff_intcoslog
    real(8) function intcoslog(mu, gx, gy)
      real(8), intent(in) :: mu, gx, gy
    end function intcoslog
    module procedure poisson_cutoff_intcoslog4
  end interface poisson_cutoff_intcoslog

  
contains


  ! ---------------------------------------------------------
  real(4) function poisson_cutoff_intcoslog4(mu, gx, gy)
    real(4), intent(in) :: mu, gx, gy

    !no PUSH_SUB, called too often

    poisson_cutoff_intcoslog4 = intcoslog(real(mu, 8), real(gx, 8), real(gy, 8))

  end function poisson_cutoff_intcoslog4
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  real(4) function poisson_cutoff_2d_0d4(x, y)
    real(4), intent(in) :: x, y

    real(8) :: res8

    !no PUSH_SUB, called too often

    res8 = c_poisson_cutoff_2d_0d(real(x, 8), real(y, 8))
    poisson_cutoff_2d_0d4 = real(res8, 4)

  end function poisson_cutoff_2d_0d4
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  real(4) function poisson_cutoff_2d_1d4(gy, gx, r_c)
    real(4), intent(in) :: gy, gx, r_c

    real(8) :: res8

    !no PUSH_SUB, called too often

    res8 = c_poisson_cutoff_2d_1d(real(gy, 8), real(gx, 8), real(r_c, 8))
    poisson_cutoff_2d_1d4 = real(res8, 4)

  end function poisson_cutoff_2d_1d4
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  real(4) function poisson_cutoff_1d_0d4(g, a, r_c)
    real(4), intent(in) :: g, a, r_c

    real(8) :: res8

    !no PUSH_SUB, called too often

    res8 = c_poisson_cutoff_1d_0d(real(g, 8), real(a, 8), real(r_c, 8))
    poisson_cutoff_1d_0d4 = real(res8, 4)

  end function poisson_cutoff_1d_0d4
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  real(4) function poisson_cutoff_3d_1d_finite4(gx, gperp, xsize, rsize)
    real(4), intent(in) :: gx, gperp, rsize, xsize

    real(8) :: res8

    !no PUSH_SUB, called too often

    res8 = c_poisson_cutoff_3d_1d_finite(real(gx, 8), real(gperp, 8), real(xsize, 8), real(rsize, 8))
    poisson_cutoff_3d_1d_finite4 = real(res8, 4)

  end function poisson_cutoff_3d_1d_finite4
  ! ---------------------------------------------------------
  

  ! ---------------------------------------------------------
  FLOAT function poisson_cutoff_3D_0D(x, r) result(cutoff)
    FLOAT, intent(in) ::  x, r

    !no PUSH_SUB, called too often

    cutoff = M_ONE - cos(x*r)

  end function poisson_cutoff_3D_0D
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function poisson_cutoff_3D_1D(x, p, rmax) result(cutoff)
    FLOAT, intent(in) ::  x, p, rmax

    integer :: j
    FLOAT :: dr, r, sum

    integer :: nr = CNST(1000)

    !no PUSH_SUB, called too often

    if ( abs(x) < M_EPSILON ) then
      ! Simpson rule for the G_x = 0 contribution -log(r)
      dr = rmax/real(nr)
      sum = M_ZERO
      do j = 1, nr - 1, 2
        r = j*dr
        sum = sum - M_FOUR*r*loct_bessel_j0(p*r)*log(r)
        r = (j+1)*dr
        sum = sum - M_TWO*r*loct_bessel_j0(p*r)*log(r)
      end do
      sum = sum - rmax*loct_bessel_j0(p*rmax)*log(rmax)
      cutoff = (p**2)*M_THIRD*sum*dr
    else
      cutoff = M_ONE + p*rmax*loct_bessel_j1(p*rmax)*loct_bessel_k0(x*rmax) &
        - x*rmax*loct_bessel_j0(p*rmax)*loct_bessel_k1(x*rmax)
    end if

  end function poisson_cutoff_3D_1D


  ! ---------------------------------------------------------
  FLOAT function poisson_cutoff_3D_2D(p, z, r) result(cutoff)
    FLOAT, intent(in) ::  p, z, r

    !no PUSH_SUB, called too often

    if (abs(p) < M_EPSILON) then
      cutoff = M_ONE - cos(z*r) - z*r*sin(z*r)
    else
      cutoff = M_ONE + exp(-p*r)*(z*sin(z*r)/p-cos(z*r))
    end if

  end function poisson_cutoff_3D_2D
  ! ---------------------------------------------------------

end module poisson_cutoff_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

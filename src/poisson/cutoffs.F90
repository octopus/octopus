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

! This module is intended to contain "only mathematical" functions
! and procedures.

module poisson_cutoffs_m
  use global_m
  use loct_math_m
  use messages_m

  implicit none

  private
  public ::                       &
    poisson_cutoff_sphere,        &
    poisson_cutoff_inf_cylinder,  &
    poisson_cutoff_fin_cylinder,  &
    poisson_cutoff_slab,          &
    besselint

  interface poisson_cutoff_fin_cylinder
    real(8) function c_poisson_cutoff_fin_cylinder(gx, gperp, xsize, rsize)
      real(8), intent(in) :: gx, gperp, rsize, xsize
    end function c_poisson_cutoff_fin_cylinder
    module procedure poisson_cutoff_fin_cylinder4
  end interface
  
contains

  real(4) function poisson_cutoff_fin_cylinder4(gx, gperp, xsize, rsize)
    real(4), intent(in) :: gx, gperp, rsize, xsize
    real(8) :: res8
    
    res8 = c_poisson_cutoff_fin_cylinder(real(gx, 8), real(gperp, 8), real(xsize, 8), real(rsize, 8))
    poisson_cutoff_fin_cylinder4 = real(res8, 4)
  end function poisson_cutoff_fin_cylinder4
  
  ! ---------------------------------------------------------
  FLOAT function poisson_cutoff_sphere(x, r) result(cutoff)
    FLOAT, intent(in) ::  x, r

    cutoff = M_ONE - cos(x*r)

  end function poisson_cutoff_sphere


  ! ---------------------------------------------------------
  FLOAT function poisson_cutoff_inf_cylinder(x, p, rmax) result(cutoff)
    FLOAT, intent(in) ::  x, p, rmax

    integer :: j
    FLOAT :: dr, r, sum

    integer :: nr = CNST(1000)

    if ( x == M_ZERO ) then
      ! Simpson rule for the G_x = 0 contribution -log(r)
      dr = rmax/real(nr)
      sum = M_ZERO;
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

  end function poisson_cutoff_inf_cylinder


  ! ---------------------------------------------------------
  FLOAT function poisson_cutoff_slab(p, z, r) result(cutoff)
    FLOAT, intent(in) ::  p, z, r

    if ( p == M_ZERO ) then
      cutoff = M_ONE - cos(z*r) - z*r*sin(z*r)
    else
      cutoff = M_ONE + exp(-p*r)*(z*sin(z*r)/p-cos(z*r))
    end if

  end function poisson_cutoff_slab


  ! ---------------------------------------------------------
  ! F(x) = (1/x) Integrate[ BesselJ[0, r], {0, x, r} ] = 
  !      = HypergeometricPFQ[ {1/2}, {1,3/2}, -x*x/r ] = 
  !      = (1/x) * 2 * sum_{k=0}^{\infty} BesselJ[k, x]
  FLOAT function besselint(x) result(y)
    FLOAT, intent(in) :: x
    integer :: k, nmax
    real(8) :: z, s
    real(8), allocatable :: bess(:)

    real(8), parameter :: large = CNST(1.0e10)

    if(x < CNST(0.2)) then
      y = M_ONE - (M_ONE/CNST(12.0))*x**2
      return
    end if

    nmax = 0

    main_loop: do 
      nmax = nmax + 100
      if(.not.allocated(bess)) ALLOCATE(bess(0:nmax), nmax+1)

      ! We need to do a backwards recursion since otherwise it is unstable.
      bess(0:nmax) = M_ZERO
      bess(nmax) = M_ZERO
      bess(nmax-1) = M_ONE
      s = bess(nmax)
      do k = nmax - 2, 0, -1
         bess(k) = (M_TWO*(k+1)/x)*bess(k+1) - bess(k+2)
         if(bess(k) > large) then
           bess(k:nmax) = bess(k:nmax) / large
           s = s / large
         end if
         if(mod(k,2).eq.0) s = s + bess(k)
      end do
      s = 2*s - bess(0)
      do k = 0, nmax
        bess(k) = bess(k)/s
      end do

      y = CNST(0.0)
      k = 2
      ! In the sum, I use the recursion relation to eliminate half of the terms.
      do 
        if(k + sqrt(CNST(40.0)*k) > nmax) exit ! Beyond this value, bess(k) could be imprecise.
        if(mod(k-2,4).eq.0) then
          z = 2*k*bess(k)/x**2
          y = y + z
          if(abs(z) < CNST(1.0e-9)) exit main_loop
        end if
        k = k + 1
      end do

      deallocate(bess)

    end do main_loop

    if(allocated(bess)) deallocate(bess)
    y = CNST(2.0)*y
  end function besselint


end module poisson_cutoffs_m

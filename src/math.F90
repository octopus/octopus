!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

! This module is intended to contain "only mathematical" functions and        !
!	procedures.                                                           !

module math
  use linalg

  implicit none

contains

! a simple congruent random number generator
subroutine quickrnd(iseed, rnd)
  integer, intent(inout) :: iseed
  real(r8), intent(inout) :: rnd

  integer, parameter :: im=6075, ia=106, ic=1283

  iseed = mod(iseed*ia + ic, im)
  rnd = real(iseed, r8)/real(im, r8)
  
end subroutine quickrnd

! Step function, needed for definition of fermi function.
function stepf(x)
  real(r8), intent(in) ::  x
  real(r8) :: stepf

  if (x.gt.100.0_r8) then
     stepf = 0.0_r8
  elseif (x.lt.-100.0_r8) then
     stepf = 2.0_r8
  else
     stepf = 2.0_r8 / ( 1.0_r8 + exp(x) )
  endif

end function stepf

! This is a Numerical Recipes based subroutine
! computes real spherical harmonics ylm in the direction of vector r:
!    ylm = c * plm( cos(theta) ) * sin(m*phi)   for   m <  0
!    ylm = c * plm( cos(theta) ) * cos(m*phi)   for   m >= 0
! with (theta,phi) the polar angles of r, c a positive normalization
subroutine grylmr(x, y, z, li, mi, ylm, grylm)
  integer, intent(in) :: li, mi
  real(r8), intent(in) :: x, y, z
  real(r8), intent(out) :: ylm, grylm(3)

  integer, parameter :: lmaxd = 20
  real(r8), parameter :: tiny=1.e-30_r8, zero=0.0_r8, half=0.5_r8, one=1.0_r8, &
       two=2.0_r8, three=3.0_r8, six=6.0_r8
  integer :: i, ilm0, l, m, mabs
  integer, save :: lmax = -1

  real(r8) :: cmi, cosm, cosmm1, cosphi, dphase, dplg, fac, &
       fourpi, plgndr, phase, pll, pmm, pmmp1, sinm, &
       sinmm1, sinphi, r2, rsize, Rx, Ry, Rz, xysize
  real(r8), save :: c(0:(lmaxd+1)*(lmaxd+1))
  

! evaluate normalization constants once and for all
  if (li.gt.lmax) then
    fourpi = two**4*datan(one)
    do l = 0, li
      ilm0 = l*l + l
      do m = 0, l
        fac = (2*l+1)/fourpi
        do i = l - m + 1, l + m
          fac = fac/i
        end do
        c(ilm0 + m) = sqrt(fac)
        ! next line because ylm are real combinations of m and -m
        if(m.ne.0) c(ilm0 + m) = c(ilm0 + m)*sqrt(two)
        c(ilm0 - m) = c(ilm0 + m)
      end do
    end do
    lmax = li
  end if

  ! if l=0, no calculations are required
  if (li.eq.0) then
    ylm = c(0)
    grylm(:) = zero
    return
  end if

  ! if r=0, direction is undefined => make ylm=0 except for l=0
  r2 = x**2 + y**2 + z**2
  if(r2.lt.tiny) then
    ylm = zero
    grylm(:) = zero
    return
  endif
  rsize = sqrt(r2)

  Rx = x/rsize
  Ry = y/rsize
  Rz = z/rsize

  ! explicit formulas for l=1 and l=2
  if(li.eq.1) then
    select case(mi)
    case(-1)
      ylm = (-c(1))*Ry
      grylm(1) = c(1)*Rx*Ry/rsize
      grylm(2) = (-c(1))*(one - Ry*Ry)/rsize
      grylm(3) = c(1)*Rz*Ry/rsize 
    case(0)
      ylm = c(2)*Rz
      grylm(1) = (-c(2))*Rx*Rz/rsize
      grylm(2) = (-c(2))*Ry*Rz/rsize
      grylm(3) = c(2)*(one - Rz*Rz)/rsize
    case(1)
      ylm = (-c(3))*Rx
      grylm(1) = (-c(3))*(one - Rx*Rx)/rsize
      grylm(2) = c(3)*Ry*Rx/rsize
      grylm(3) = c(3)*Rz*Rx/rsize
    end select
    return
  end if
         
  if(li.eq.2) then
    select case(mi)
    case(-2)
      ylm = c(4)*six*Rx*Ry
      grylm(1) = (-c(4))*six*(two*Rx*Rx*Ry - Ry)/rsize
      grylm(2) = (-c(4))*six*(two*Ry*Rx*Ry - Rx)/rsize
      grylm(3) = (-c(4))*six*(two*Rz*Rx*Ry)/rsize
    case(-1)
      ylm = (-c(5))*three*Ry*Rz
      grylm(1) = c(5)*three*(two*Rx*Ry*Rz)/rsize
      grylm(2) = c(5)*three*(two*Ry*Ry*Rz - Rz)/rsize
      grylm(3) = c(5)*three*(two*Rz*Ry*Rz - Ry)/rsize
    case(0)
      ylm = c(6)*half*(three*Rz*Rz - one)
      grylm(1) = (-c(6))*three*(Rx*Rz*Rz)/rsize
      grylm(2) = (-c(6))*three*(Ry*Rz*Rz)/rsize
      grylm(3) = (-c(6))*three*(Rz*Rz - one)*Rz/rsize
    case(1)
      ylm = (-c(7))*three*Rx*Rz
      grylm(1) = c(7)*three*(two*Rx*Rx*Rz - Rz)/rsize
      grylm(2) = c(7)*three*(two*Ry*Rx*Rz)/rsize
      grylm(3) = c(7)*three*(two*Rz*Rx*Rz - Rx)/rsize
    case(2)
      ylm = c(8)*three*(Rx*Rx - Ry*Ry)
      grylm(1) = (-c(8))*six*(Rx*Rx - Ry*Ry - one)*Rx/rsize
      grylm(2) = (-c(8))*six*(Rx*Rx - Ry*Ry + one)*Ry/rsize
      grylm(3) = (-c(8))*six*(Rx*Rx - Ry*Ry)*Rz/rsize
    end select
    return
  end if

! general algorithm based on routine plgndr of numerical recipes
  mabs = abs(mi)
  xysize = sqrt(max(Rx*Rx + Ry*Ry, tiny))
  cosphi = Rx/xysize
  sinphi = Ry/xysize
  cosm = one
  sinm = zero
  do m = 1, mabs
    cosmm1 = cosm
    sinmm1 = sinm
    cosm = cosmm1*cosphi - sinmm1*sinphi
    sinm = cosmm1*sinphi + sinmm1*cosphi
  end do
     
  if(mi.lt.0) then
    phase = sinm
    dphase = mabs*cosm
  else
    phase = cosm
    dphase = (-mabs)*sinm
  end if

  pmm = one
  fac = one

  if(mabs.gt.zero) then
    do i = 1, mabs
      pmm = (-pmm)*fac*xysize
      fac = fac + two
    end do
  end if

  if(li.eq.mabs) then
    plgndr = pmm
    dplg = (-li)*Rz*pmm/(xysize**2)
  else
    pmmp1 = Rz*(2*mabs + 1)*pmm
    if(li.eq.mabs + 1) then
      plgndr = pmmp1
      dplg = -((li*Rz*pmmp1 - (mabs + li)*pmm)/(xysize**2))
    else
      do l = mabs + 2, li
        pll = (Rz*(2*l - 1)*pmmp1 - (l + mabs - 1)*pmm)/(l - mabs)
        pmm = pmmp1
        pmmp1 = pll
      end do
      plgndr = pll
      dplg = -((li*Rz*pll - (l + mabs - 1)*pmm)/(xysize**2))
    end if
  end if   

  ilm0 = li*li + li
  cmi = c(ilm0 + mi)
  ylm = cmi*plgndr*phase
  grylm(1) = (-cmi)*dplg*Rx*Rz*phase/rsize     &
       -cmi*plgndr*dphase*Ry/(rsize*xysize**2)
  grylm(2) = (-cmi)*dplg*Ry*Rz*phase/rsize     &
       +cmi*plgndr*dphase*Rx/(rsize*xysize**2)
  grylm(3)= cmi*dplg*(one - Rz*Rz)*phase/rsize
   
  return
end subroutine grylmr

! Compute the Weights for finite-difference calculations:
!
!  N -> highest order fo the derivative to be approximated
!  M -> number of grid points to be used in the approsimation.
!
!  c(j,k,i) -> ith order derivative at kth-order approximation
!              j=0,k: the coefficients acting of each point
subroutine weights(N, M, cc)
  integer, intent(in) :: N, M
  real(r8), intent(out) :: cc(0:M, 0:M, 0:N)

  integer :: i, j, k, mn
  real(r8) :: c1, c2, c3, c4, c5, xi
  real(r8) :: x(0:M)

  ! grid-points for one-side finite-difference formulas on an equi.spaced grid
  ! x(:) = (/(i,i=0,M)/) 

  ! grid-points for centered finite-difference formulas on an equi.spaced grid
  mn = M/2
  x(:) = (/0,(-i,i,i=1,mn)/)

  xi = 0.0_r8  ! point at which the approx. are to be accurate

  cc = 0.0_r8
  cc(0,0,0) = 1.0_r8

  c1 = 1.0_r8
  c4 = x(0) - xi
       
  do j = 1, M
    mn = min(j,N)
    c2 = 1.0_r8
    c5 = c4
    c4 = x(j) - xi
    
    do k = 0, j - 1
      c3 = x(j) - x(k)
      c2 = c2*c3
      
      if (j <= N) cc(k, j - 1, j) = 0.0_r8
      cc(k, j, 0) = c4*cc(k, j - 1, 0)/c3
      
      do i = 1, mn
        cc(k, j, i) = (c4*cc(k, j - 1, i) - i*cc(k, j - 1, i - 1))/c3
      end do
      
      cc(j, j, 0) = -c1*c5*cc(j - 1, j - 1, 0) / c2         
    end do
    
    do i = 1, mn
      cc(j, j, i) = c1*(i*cc(j - 1, j - 1, i - 1) - c5*cc(j - 1, j - 1, i))/c2
    end do

    c1 = c2
  end do

end subroutine weights

#include "undef.F90"
#include "complex.F90"
#include "math_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "math_inc.F90"

end module math

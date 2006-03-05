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
!!
!! $Id$

#include "global.h"

! This module is intended to contain "only mathematical" functions
! and procedures.

module math_m
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_basic_alg_m
  use lib_adv_alg_m

  implicit none

  private
  public ::                     &
    dconjugate_gradients,       &
    zconjugate_gradients,       &
    ddot_product,               &
    zdot_product,               &
    quickrnd,                   &
    stepf,                      &
    grylmr,                     &
    weights,                    &
    cutoff0,                    &
    cutoff1,                    &
    cutoff2,                    &
    besselint,                  &
    dextrapolate, zextrapolate, &
    sort,                       &
    factorial,                  &
    hermite,                    &
    math_divisors,              &
    set_app_threshold,          &
    operator(.app.)


  !------------------------------------------------------------------------------
  ! This is the common interface to a sorting routine.
  ! It performs the shell algorithm, not as fast as the quicksort for large numbers,
  ! but it seems that better for moderate numbers (around 100).
  ! Their possible interfaces are:
  !   subroutine sort(a [, ind] )
  !     FLOAT, intent(inout) :: a(:)
  !     integer, intent(inout), optional :: ind(:)
  !     ! This routine sorts, from smallest to largest, the array a.
  !     ! If the integer array ind is present, it puts in it the indexing
  !     ! of the sorting, so that other arrays can be sorted according to
  !     ! the sorting of a.
  !   end subroutine sort
  !
  !   subroutine sort(a, x)
  !     FLOAT, intent(inout) :: a(:)
  !     FLOAT_OR_COMPLEX, intent(inout) :: x(:, : [, :])
  !     ! This routine sorts, from smallest to largest, the array a.
  !     ! The real or complex array x, which may be two or three dimensional,
  !     ! is sorted according to the ordering of a. The last dimension of x
  !     ! must have the same size as a.
  !   end subroutine sort
  interface sort
    module procedure shellsort, dshellsort1, zshellsort1, dshellsort2, zshellsort2
  end interface
  !------------------------------------------------------------------------------


  ! This common interface applies to the two procedures defined in math_cg_inc.F90
  interface dconjugate_gradients
    module procedure dsym_conjugate_gradients, dbi_conjugate_gradients
  end interface

  interface zconjugate_gradients
    module procedure zsym_conjugate_gradients, zbi_conjugate_gradients
  end interface

  ! Euler McLaurin Integration constants
  FLOAT, public, parameter :: EMcLCoeff(1:5) =  &
    (/                                          &
    CNST( 95.0)/CNST(288.0),                    &
    CNST(   317.0)/CNST(240.0),                 &
    CNST(    23.0)/CNST( 30.0),                 &
    CNST(   793.0)/CNST(720.0),                 &
    CNST(   157.0)/CNST(160.0)                  &
    /)


  ! This operator is .true. if the two operands are approximately
  ! equal (i.e. equal to within APP_THRESHOLD). For arrays, all
  ! elements have to be approximately equal.
  ! The subroutine set_app_threshold changes APP_THRESHOLD.
  FLOAT :: APP_THRESHOLD = CNST(1.0e-10)
  interface operator(.app.)
    module procedure dapproximate_equal
    module procedure dapproximate_equal_1
    module procedure dapproximate_equal_2
    module procedure dapproximate_equal_3
    module procedure zapproximate_equal
    module procedure zapproximate_equal_1
    module procedure zapproximate_equal_2
    module procedure zapproximate_equal_3
  end interface

contains

  ! ---------------------------------------------------------
  recursive function hermite(n, x) result (h)
    integer, intent(in) :: n
    FLOAT,   intent(in) :: x

    FLOAT :: h

    if(n<=0) then
      h = M_ONE
    elseif(n==1) then
      h = M_TWO*x
    else
      h = M_TWO*x*hermite(n-1,x) - M_TWO*(n-1)*hermite(n-2,x)
    end if

  end function hermite


  ! ---------------------------------------------------------
  recursive function factorial (n) RESULT (fac)
    integer, intent(in) :: n
    integer :: fac
    if(n<=1) then 
      fac = 1
    else
      fac = n*factorial(n-1)
    end if
  end function factorial


  ! ---------------------------------------------------------
  ! a simple congruent random number generator
  subroutine quickrnd(iseed, rnd)
    integer, intent(inout) :: iseed
    FLOAT, intent(inout) :: rnd

    integer, parameter :: im=6075, ia=106, ic=1283

    iseed = mod(iseed*ia + ic, im)
    rnd = real(iseed, PRECISION)/real(im, PRECISION)

  end subroutine quickrnd


  ! ---------------------------------------------------------
  ! Step function, needed for definition of fermi function.
  function stepf(x)
    FLOAT, intent(in) ::  x
    FLOAT :: stepf

    if (x.gt.CNST(100.0)) then
      stepf = M_ZERO
    elseif (x.lt.CNST(-100.0)) then
      stepf = M_TWO
    else
      stepf = M_TWO / ( M_ONE + exp(x) )
    end if

  end function stepf


  ! ---------------------------------------------------------
  ! This is a Numerical Recipes based subroutine
  ! computes real spherical harmonics ylm in the direction of vector r:
  !    ylm = c * plm( cos(theta) ) * sin(m*phi)   for   m <  0
  !    ylm = c * plm( cos(theta) ) * cos(m*phi)   for   m >= 0
  ! with (theta,phi) the polar angles of r, c a positive normalization
  subroutine grylmr(x, y, z, li, mi, ylm, grylm)
    integer, intent(in) :: li, mi
    FLOAT, intent(in) :: x, y, z
    FLOAT, intent(out) :: ylm, grylm(3)

    integer, parameter :: lmaxd = 20
    FLOAT,   parameter :: tiny = CNST(1.e-30)
    integer :: i, ilm0, l, m, mabs
    integer, save :: lmax = -1

    FLOAT :: cmi, cosm, cosmm1, cosphi, dphase, dplg, fac, &
      fourpi, plgndr, phase, pll, pmm, pmmp1, sinm, &
      sinmm1, sinphi, r2, rsize, Rx, Ry, Rz, xysize
    FLOAT, save :: c(0:(lmaxd+1)*(lmaxd+1))


    ! evaluate normalization constants once and for all
    if (li.gt.lmax) then
      fourpi = M_FOUR*M_PI
      do l = 0, li
        ilm0 = l*l + l
        do m = 0, l
          fac = (2*l+1)/fourpi
          do i = l - m + 1, l + m
            fac = fac/i
          end do
          c(ilm0 + m) = sqrt(fac)
          ! next line because ylm are real combinations of m and -m
          if(m.ne.0) c(ilm0 + m) = c(ilm0 + m)*sqrt(M_TWO)
          c(ilm0 - m) = c(ilm0 + m)
        end do
      end do
      lmax = li
    end if

    ! if l=0, no calculations are required
    if (li.eq.0) then
      ylm = c(0)
      grylm(:) = M_ZERO
      return
    end if

    ! if r=0, direction is undefined => make ylm=0 except for l=0
    r2 = x**2 + y**2 + z**2
    if(r2.lt.tiny) then
      ylm = M_ZERO
      grylm(:) = M_ZERO
      return
    end if
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
        grylm(2) = (-c(1))*(M_ONE - Ry*Ry)/rsize
        grylm(3) = c(1)*Rz*Ry/rsize
      case(0)
        ylm = c(2)*Rz
        grylm(1) = (-c(2))*Rx*Rz/rsize
        grylm(2) = (-c(2))*Ry*Rz/rsize
        grylm(3) = c(2)*(M_ONE - Rz*Rz)/rsize
      case(1)
        ylm = (-c(3))*Rx
        grylm(1) = (-c(3))*(M_ONE - Rx*Rx)/rsize
        grylm(2) = c(3)*Ry*Rx/rsize
        grylm(3) = c(3)*Rz*Rx/rsize
      end select
      return
    end if

    if(li.eq.2) then
      select case(mi)
      case(-2)
        ylm = c(4)*M_SIX*Rx*Ry
        grylm(1) = (-c(4))*M_SIX*(M_TWO*Rx*Rx*Ry - Ry)/rsize
        grylm(2) = (-c(4))*M_SIX*(M_TWO*Ry*Rx*Ry - Rx)/rsize
        grylm(3) = (-c(4))*M_SIX*(M_TWO*Rz*Rx*Ry)/rsize
      case(-1)
        ylm = (-c(5))*M_THREE*Ry*Rz
        grylm(1) = c(5)*M_THREE*(M_TWO*Rx*Ry*Rz)/rsize
        grylm(2) = c(5)*M_THREE*(M_TWO*Ry*Ry*Rz - Rz)/rsize
        grylm(3) = c(5)*M_THREE*(M_TWO*Rz*Ry*Rz - Ry)/rsize
      case(0)
        ylm = c(6)*M_HALF*(M_THREE*Rz*Rz - M_ONE)
        grylm(1) = (-c(6))*M_THREE*(Rx*Rz*Rz)/rsize
        grylm(2) = (-c(6))*M_THREE*(Ry*Rz*Rz)/rsize
        grylm(3) = (-c(6))*M_THREE*(Rz*Rz - M_ONE)*Rz/rsize
      case(1)
        ylm = (-c(7))*M_THREE*Rx*Rz
        grylm(1) = c(7)*M_THREE*(M_TWO*Rx*Rx*Rz - Rz)/rsize
        grylm(2) = c(7)*M_THREE*(M_TWO*Ry*Rx*Rz)/rsize
        grylm(3) = c(7)*M_THREE*(M_TWO*Rz*Rx*Rz - Rx)/rsize
      case(2)
        ylm = c(8)*M_THREE*(Rx*Rx - Ry*Ry)
        grylm(1) = (-c(8))*M_SIX*(Rx*Rx - Ry*Ry - M_ONE)*Rx/rsize
        grylm(2) = (-c(8))*M_SIX*(Rx*Rx - Ry*Ry + M_ONE)*Ry/rsize
        grylm(3) = (-c(8))*M_SIX*(Rx*Rx - Ry*Ry)*Rz/rsize
      end select
      return
    end if

    ! general algorithm based on routine plgndr of numerical recipes
    mabs = abs(mi)
    xysize = sqrt(max(Rx*Rx + Ry*Ry, tiny))
    cosphi = Rx/xysize
    sinphi = Ry/xysize
    cosm = M_ONE
    sinm = M_ZERO
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

    pmm = M_ONE
    fac = M_ONE

    if(mabs.gt.M_ZERO) then
      do i = 1, mabs
        pmm = (-pmm)*fac*xysize
        fac = fac + M_TWO
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
    grylm(3)= cmi*dplg*(M_ONE - Rz*Rz)*phase/rsize

    return
  end subroutine grylmr


  ! ---------------------------------------------------------
  ! Compute the weights for finite-difference calculations:
  !
  !  N -> highest order fo the derivative to be approximated
  !  M -> number of grid points to be used in the approsimation.
  !
  !  c(j,k,i) -> ith order derivative at kth-order approximation
  !              j=0,k: the coefficients acting of each point
  subroutine weights(N, M, cc)
    integer, intent(in) :: N, M
    FLOAT, intent(out) :: cc(0:M, 0:M, 0:N)

    integer :: i, j, k, mn
    FLOAT :: c1, c2, c3, c4, c5, xi
    FLOAT, allocatable :: x(:)

    ALLOCATE(x(0:M), M+1)
    ! grid-points for one-side finite-difference formulas on an equi.spaced grid
    ! x(:) = (/(i,i=0,M)/)

    ! grid-points for centered finite-difference formulas on an equi.spaced grid
    mn = M/2
    x(:) = (/0,(-i,i,i=1,mn)/)

    xi = M_ZERO  ! point at which the approx. are to be accurate

    cc = M_ZERO
    cc(0,0,0) = M_ONE

    c1 = M_ONE
    c4 = x(0) - xi

    do j = 1, M
      mn = min(j,N)
      c2 = M_ONE
      c5 = c4
      c4 = x(j) - xi

      do k = 0, j - 1
        c3 = x(j) - x(k)
        c2 = c2*c3

        if (j <= N) cc(k, j - 1, j) = M_ZERO
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

    deallocate(x)

  end subroutine weights


  ! ---------------------------------------------------------
  FLOAT function cutoff0(x,r)
    FLOAT, intent(in) ::  x,r

    cutoff0 = M_ONE - cos(x*r)

  end function cutoff0


  ! ---------------------------------------------------------
  FLOAT function cutoff1(x, p, rmax)
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
      cutoff1 = (p**2)*M_THIRD*sum*dr
    else
      cutoff1 = M_ONE + p*rmax*loct_bessel_j1(p*rmax)*loct_bessel_k0(x*rmax) &
        - x*rmax*loct_bessel_j0(p*rmax)*loct_bessel_k1(x*rmax)
    end if

  end function cutoff1


  ! ---------------------------------------------------------
  FLOAT function cutoff2(p, z, r)
    FLOAT, intent(in) ::  p, z, r

    if ( p == M_ZERO ) then
      cutoff2 = M_ONE - cos(z*r) - z*r*sin(z*r)
    else
      cutoff2 = M_ONE + exp(-p*r)*(z*sin(z*r)/p-cos(z*r))
    end if

  end function cutoff2


  ! ---------------------------------------------------------
  FLOAT function besselint(x) result(y)
    FLOAT, intent(in) :: x
    real(8), external :: bessel
    integer :: k
    real(8) :: z
    if(x < CNST(0.2)) then
      y = CNST(2.0)*M_PI - (M_PI/CNST(6.0))*x**2
      return
    end if
    y = CNST(0.0)
    k = 1
    do
      z = loct_bessel(k, x)/x
      y = y + z
      if(abs(z) < CNST(1.0e-9)) exit
      k = k + 2
    end do
    y = CNST(4.0)*M_PI*y
  end function besselint


  ! ---------------------------------------------------------
  FLOAT function ddot_product(a, b) result(r)
    FLOAT, intent(in) :: a(:), b(:)
    r = dot_product(a, b)
  end function ddot_product


  ! ---------------------------------------------------------
  CMPLX function zdot_product(a, b) result(r)
    CMPLX, intent(in) :: a(:), b(:)
    r = sum(conjg(a(:))*b(:))
  end function zdot_product


  ! ---------------------------------------------------------
  subroutine shellsort(a, ind)
    FLOAT, intent(inout) :: a(:)
    integer, intent(inout), optional :: ind(:)

    integer :: i,j,inc,n, indi
    FLOAT   :: v

    n = size(a)

    if(present(ind)) then
      do i = 1, n
        ind(i) = i
      end do
    end if

    inc = 1
    do
      inc=3*inc+1
      if (inc > n) exit
    end do

    do
      inc=inc/3
      do i=inc+1,n
        v=a(i)
        if(present(ind)) indi = ind(i)
        j=i
        do
          if (a(j-inc) <= v) exit
          !if (a(j-inc) >= v) exit
          a(j)=a(j-inc)
          if(present(ind)) ind(j) = ind(j-inc)
          j=j-inc
          if (j <= inc) exit
        end do
        a(j)=v
        if(present(ind)) ind(j) = indi
      end do
      if (inc <= 1) exit
    end do

  end subroutine shellsort


  ! ---------------------------------------------------------
  subroutine math_divisors(n, n_divisors, divisors)
    integer, intent(in)    :: n
    integer, intent(inout) :: n_divisors
    integer, intent(out)   :: divisors(:)

    integer :: i, max_d

    ASSERT(n_divisors > 1)
    max_d = n_divisors

    n_divisors = 1
    divisors(n_divisors) = 1
    do i = 2, n/2
      if(mod(n, i)==0) then
        n_divisors = n_divisors + 1

        if(n_divisors > max_d - 1) then
          message(1) = "Internal error in get_divisors. Please increase n_divisors"
          call write_fatal(1)
        end if

        divisors(n_divisors) = i        
      end if
    end do
    n_divisors = n_divisors + 1
    divisors(n_divisors) = n
  end subroutine math_divisors

  subroutine set_app_threshold(thr)
    FLOAT, intent(in) :: thr
    APP_THRESHOLD = thr
  end subroutine set_app_threshold


#include "undef.F90"
#include "complex.F90"
#include "math_cg_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "math_cg_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "math_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "math_inc.F90"

end module math_m

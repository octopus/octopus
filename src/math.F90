#include "config.h"

! This module is intended to contain "only mathematical" functions and        !
!	procedures.                                                           !

module math

  use global

  implicit none
  
  private
  public :: stepf, derf, derfc, ratint, polint, ylmr, grylmr, sort, &
       gaussj, high_derivative, weights, asinh, quickrnd
  
contains

! as obvious the arcsin hiperbolicus
#include "asinh.F90"

! a simple congruent random number generator
subroutine quickrnd(iseed, rnd)
  integer, intent(inout) :: iseed
  real(r8), intent(inout) :: rnd

  integer, parameter :: im=6075, ia=106, ic=1283

  iseed = mod(iseed*ia + ic, im)
  rnd = real(iseed, r8)/real(im, r8)
  
end subroutine quickrnd

! Step function, needed for definition of fermi function.                     !
function stepf(x)
  real(r8), intent(in) ::  x
  real(r8) :: stepf

! Note:                                                                       ! 
!     complementary error function. ref: fu & ho, prb 28, 5480 (1983)         !
!     stepf=derfc(x)                                                          !
!     improved step function. ref: methfessel & paxton prb40 (15/aug/89)      !
!     parameter (c=0.5641895835_r8)                                           !
!     stepf=derfc(x)-c*x*dexp(-x*x)                                           !

  if (x.gt.100.0_r8) then
     stepf = 0.0_r8
  elseif (x.lt.-100.0_r8) then
     stepf = 2.0_r8
  else
     stepf = 2.0_r8 / ( 1.0_r8 + exp(x) )
  endif

end function stepf

!	complementary error function from "numerical recipes"
!	note: single precision accuracy
function derfc(x)
  real(r8), intent(in) :: x
  real(r8) :: derfc

  real(r8):: t, z

  z = abs(x)
  t = 1.0_r8/(1.0_r8 + 0.5_r8*z)

  derfc=t*exp(-(z*z)-1.26551223_r8+t*(1.00002368_r8+t*(0.37409196_r8+    &
        t*(0.09678418_r8+t*(-0.18628806_r8+t*(0.27886807_r8+             &
        t*(-1.13520398_r8+t*(1.48851587_r8+t*(-0.82215223_r8+            &
        t*.17087277_r8)))))))))

  if (x < 0.0_r8) derfc=2.0_r8-derfc

end function derfc

! Error function (any reasonable compiler will inline this ;)
function derf (x)
  real(r8), intent(in) :: x
  real(r8) :: derf

  derf = 1.0_r8 - derfc(x)

end function derf


! An Adapted Numerical Recipes ratint.f and polint.f  interpolation subroutines  (1998)
subroutine ratint(xa, ya, n, x, y, dy)
  integer, intent(in) :: n
  real(r8), intent(IN) :: xa(n), ya(n)
  real(r8), intent(inout) :: x, y, dy
 
  real(r8), parameter :: tiny = 1.d-20
  real(r8) :: c(n), d(n)
  integer :: i, m, ns 
  real(r8) :: dd, h, hh, t, w

  ns = 1
  hh = abs(x - xa(1))
  do i = 1, n
    h = abs(x - xa(i))
    if(h < tiny)then
      y  = ya(i)
      dy = 0.0_r8
      return
    else if(h < hh) then
      ns = i
      hh = h
    end if
    c(i) = ya(i)
    d(i) = ya(i) + tiny
  end do

  y  = ya(ns)
  ns = ns - 1
  main: do m = 1, n - 1
    do i = 1, n - m
      w = c(i + 1) - d(i)
      h = xa(i + m) - x
      t = (xa(i) - x)*d(i)/h
      dd = t - c(i + 1)
      if(dd == 0.0_r8) exit main
      dd = w/dd
      d(i) = c(i + 1)*dd
      c(i) = t*dd
    end do
    if(2*ns.lt.n - m) then
      dy = c(ns + 1)
    else
      dy = d(ns)
      ns = ns-1
    end if
    y = y + dy
  end do main

  ! as rational interpolation does not converge,we try with a polynomial one
  ! MALM: I think this is wrong, so I comment it out
  !call polint(xa, ya, n, x, y, dy)

end subroutine ratint

! Subroutine polint
subroutine polint(xa, ya, n, x, y, dy)
  integer, intent(in) :: n
  real(r8), intent(IN) ::  xa(n), ya(n)
  real(r8), intent(inout) :: x, y, dy

  integer :: i, m, ns
  real(r8) :: c(n), d(n)
  real(r8) :: den, dif, dift, ho, hp, w

  ns = 1
  dif = abs(x - xa(1))
  do i = 1, n 
    dift = abs(x - xa(i))
    if (dift < dif) then
      ns  = i
      dif = dift
    end if
    c(i) = ya(i)
    d(i) = ya(i)
  end do
  y = ya(ns)
  ns = ns - 1
  do m = 1, n - 1
    do i = 1, n - m
      ho = xa(i) - x
      hp = xa(i + m) - x
      w  = c(i + 1) - d(i)
      den = ho - hp
      if(den == 0.0_r8)stop 'polint: den = 0'
      den = w/den
      d(i) = hp*den
      c(i) = ho*den
    end do
    if (2*ns.lt.n - m)then
      dy = c(ns + 1)
    else
      dy = d(ns)
      ns = ns - 1
    end if
    y = y + dy
  end do

end subroutine polint

! This is a Numerical Recipes based subroutine.
! Computes real spherical harmonics ylm in the direction of vector r:
!    ylm = c * plm( cos(theta) ) * sin(m*phi)   for   m <  0
!    ylm = c * plm( cos(theta) ) * cos(m*phi)   for   m >= 0
! with (theta,phi) the polar angles of r, c a positive normalization
! constant and plm associated legendre polynomials.
subroutine ylmr(x, y, z, li, mi, ylm)
  integer, intent(in) :: li, mi
  real(r8), intent(in) :: x, y, z
  real(r8), intent(out) :: ylm

  integer, parameter :: lmaxd = 20
  real(r8), parameter :: tiny=1.d-30, zero=0.0_r8, half=0.5_r8, one=1.0_r8, &
       two=2.0_r8, three=3.0_r8, six=6.0_r8
  real(r8), save :: c(0:(lmaxd+1)*(lmaxd+1))

  integer :: i, ilm0, l, m, mabs
  integer, save :: lmax = -1
  real(r8) :: cmi, cosm, cosmm1, cosphi, dphase, dplg, fac, &
       fourpi, plgndr, phase, pll, pmm, pmmp1, sinm, &
       sinmm1, sinphi, r2, rsize, Rx, Ry, Rz, xysize

! evaluate normalization constants once and for all                           !
  if(li.gt.lmax) then
    fourpi = two**4*atan(one)
    do l = 0, li
      ilm0 = l*l + l
      do m = 0, l
        fac = (2*l + 1)/fourpi
        do i = l - m + 1, l + m
          fac = fac/i
        end do
        c(ilm0 + m) = sqrt(fac)
        ! next line because ylm are real combinations of m and -m
        if (m.ne.0) c(ilm0+m) = c(ilm0 + m)*sqrt(two)
        c(ilm0 - m) = c(ilm0 + m)
      end do
    end do
    lmax = li
  end if

  ! if l=0, no calculations are required
  if(li.eq.0) then
    ylm = c(0)
    return
  end if

  ! if r=0, direction is undefined => make ylm=0 except for l=0
  r2 = x**2 + y**2 + z**2
  if(r2.lt.tiny) then
    ylm = zero
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
    case(0)
      ylm = c(2)*Rz
    case(1)
      ylm = (-c(3))*Rx
    end select
    return
  end if
         
  if(li.eq.2) then
    select case(mi)
    case(-2)
      ylm = c(4)*six*Rx*Ry
    case(-1)
      ylm = (-c(5))*three*Ry*Rz
    case(0)
      ylm = c(6)*half*(three*Rz*Rz-one)
    case(1)
      ylm = (-c(7))*three*Rx*Rz
    case(2)
      ylm = c(8)*three*(Rx*Rx-Ry*Ry)
    end select
    return
  end if

! general algorithm based on routine plgndr of 'numerical recipes'
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
  enddo
     
  if(mi.lt.0) then
    phase = sinm
  else
    phase = cosm
  endif

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
  else
    pmmp1 = Rz*(2*mabs + 1)*pmm
    if(li.eq.mabs+1) then
      plgndr = pmmp1
    else
      do l = mabs + 2, li
        pll = (Rz*(2*l - 1)*pmmp1 - (l + mabs - 1)*pmm)/(l - mabs)
        pmm = pmmp1
        pmmp1 = pll
      end do
      plgndr = pll
    end if
  end if
        
  ilm0 = li*li + li
  cmi  = c(ilm0 + mi)
  ylm  = cmi*plgndr*phase

  return
end subroutine ylmr

! This is a Numerica Recipes based subroutine
! computes real spherical harmonics ylm in the direction of vector r:
!    ylm = c * plm( cos(theta) ) * sin(m*phi)   for   m <  0
!    ylm = c * plm( cos(theta) ) * cos(m*phi)   for   m >= 0
! with (theta,phi) the polar angles of r, c a positive normalization
subroutine grylmr(x, y, z, li, mi, ylm, grylm)
  integer, intent(in) :: li, mi
  real(r8), intent(in) :: x, y, z
  real(r8), intent(out) :: ylm, grylm(3)

  integer, parameter :: lmaxd = 20
  real(r8), parameter :: tiny=1.d-30, zero=0.0_r8, half=0.5_r8, one=1.0_r8, &
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

! Sorts an arrat by the heapsort method,
! W. H. Preuss et al. (Numerical Recipes)
subroutine sort(n, arrin, indx)
  integer, intent(in) :: n
  real(r8), intent(IN) :: arrin(n)
  integer, intent(out):: indx(n)
  
  real(r8) :: q
  integer :: i, indxt, ir, j, l

  do j = 1, n
    indx(j) = j
  end do
 
  l = n/2 + 1
  ir = n
20 continue
  if (l .gt. 1) then
    l = l - 1
    indxt = indx(l)
    q = arrin(indxt)
  else
    indxt = indx(ir)
    q = arrin(indxt)
    indx(ir) = indx(1)
    ir = ir - 1
    if (ir .eq. 1) then
      indx(1) = indxt
      return
    end if
  end if
  i = l
  j = l + l
30 continue
  if (j .le. ir) then
    if (j .lt. ir) then
      if (arrin(indx(j)) .lt. arrin(indx(j+1))) j = j + 1
    end if
    if (q .lt. arrin(indx(j))) then
      indx(i) = indx(j)
      i = j
      j = j + j
    else
      j = ir + 1
    end if

    go to 30
  end if
  indx(i) = indxt

  go to 20

end subroutine sort

! Numerical Recipes routine. Explained in the book
subroutine gaussj(a, n, np, b, m, mp)
  integer ::  m, mp, n, np
  real(r8) :: a(np, np), b(np, mp)

  integer, parameter :: NMAX=50
  integer :: i, icol, irow, j, k, l, ll, indxc(NMAX), indxr(NMAX), ipiv(NMAX)
  real(r8) :: big, dum, pivinv

  do j = 1, n
    ipiv(j) = 0
  end do
  do i = 1, n
    big = 0._r8
    do j = 1, n
      if(ipiv(j).ne.1)then
        do k = 1, n
          if(ipiv(k).eq.0) then
            if(abs(a(j, k)).ge.big)then
              big = abs(a(j, k))
              irow = j
              icol = k
            end if
          else if(ipiv(k).gt.1) then
            stop 'singular matrix in gaussj'
          end if
        end do
      end if
    end do
    ipiv(icol) = ipiv(icol) + 1
    if(irow.ne.icol) then
      do l = 1, n
        dum = a(irow, l)
        a(irow, l) = a(icol, l)
        a(icol, l) = dum
      end do
      do l = 1, m
        dum = b(irow, l)
        b(irow, l) = b(icol, l)
        b(icol, l) = dum
      end do
    end if
    indxr(i) = irow
    indxc(i) = icol
    if(a(icol, icol).eq.0.) pause 'singular matrix in gaussj'
    pivinv = 1./a(icol, icol)
    a(icol, icol)=1.
    do l = 1, n
      a(icol, l) = a(icol, l)*pivinv
    end do
    do l = 1, m
      b(icol, l) = b(icol, l)*pivinv
    end do
    do ll = 1, n
      if(ll.ne.icol)then
        dum = a(ll, icol)
        a(ll, icol) = 0._r8
        do l = 1, n
          a(ll, l) = a(ll, l) - a(icol, l)*dum
        end do
        do l = 1, m
          b(ll, l) = b(ll, l) - b(icol, l)*dum
        end do
      end if
    end do
  end do

  do l = n, 1, -1
    if(indxr(l).ne.indxc(l))then
      do k = 1, n
        dum = a(k, indxr(l))
        a(k, indxr(l)) = a(k, indxc(l))
        a(k, indxc(l)) = dum
      end do
    end if
  end do

  return
end subroutine gaussj

! computes the coefficients for the high-orde finite differences
! approximation to the kinetic energy operator. (3D-Laplacian)
!
!                      Norder
!                       ---
!        nabla^2 f(i) = \     cN(k)* { f(i+k) + f(i-k) }  + cN(0)*f(i)
!                       /__
!                      k = 1
!
! A general algorithm is implemented based on B. Fornberg and D.M. Sloan.
! The first six-order are tabulated as in ref. J.R. Chelikowsky, N. Troullier,
! K. Wu and Y. Saad, PRB50, 11355 (1994); and references therein)
!
! A. Rubio (September 1999)
subroutine high_derivative(norder, cn)
  integer, intent(in) :: norder           ! high order discretization
  real(r8), intent(out) :: cn(0:norder)
  
  integer :: i, Morder, N_laplac
  real(r8), allocatable :: cc(:,:,:) 


  ! Up to order six, the values are given explicitly
  if (Norder <= 6) then
    select case(Norder)
    case(1) 
      cN(0) = -2.0_r8 
      cN(1) =  1.0_r8       
    case(2) 
      cN(0) = -5.0_r8/2._r8
      cN(1) =  4.0_r8/3.0_r8       
      cN(2) = -1.0_r8/12.0_r8
    case(3) 
      cN(0) = -49.0_r8/18.0_r8
      cN(1) =  3.0_r8/2.0_r8       
      cN(2) = -3.0_r8/20.0_r8
      cN(3) =  1.0_r8/90.0_r8
    case(4) 
      cN(0) = -205.0_r8/72.0_r8
      cN(1) =  8.0_r8/5.0_r8       
      cN(2) = -1.0_r8/5.0_r8
      cN(3) =  8.0_r8/315.0_r8
      cN(4) = -1.0_r8/560.0_r8
    case(5) 
      cN(0) = -5269.0_r8/1800.0_r8
      cN(1) =  5.0_r8/3.0_r8
      cN(2) = -5.0_r8/21.0_r8
      cN(3) =  5.0_r8/126.0_r8
      cN(4) = -5.0_r8/1008.0_r8
      cN(5) =  1.0_r8/3150.0_r8
    case(6) 
      cN(0) = -5369.0_r8/1800._r8
      cN(1) =  12.0_r8/7.0_r8
      cN(2) = -15.0_r8/56.0_r8
      cN(3) =  10.0_r8/189.0_r8
      cN(4) = -1.0_r8/112.0_r8
      cN(5) =  2.0_r8/1925.0_r8
      cN(6) = -1.0_r8/16632.0_r8
    end select


    ! General calculation for centered approximation (the real orde of the
    ! approximation is 2*Norder, as Norder refers to the number of positive
    ! coefficients, in fact the total number of coefficients is 2*Norder +1 in
    ! the finite-difference approach to the Laplacian).
    ! NOTE: the weight subroutine is general for all differential operators!
  else
    Morder = Norder * 2
    N_laplac = 2
    allocate(cc(0:Morder, 0:Morder, 0:N_laplac))

    call weights(N_laplac, Morder, cc)

    do i = 0, Norder
      cN(i) = cc(i*2, Morder, N_laplac)           !for the laplacian
    end do

    deallocate(cc)
  endif

end subroutine high_derivative


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

end module math

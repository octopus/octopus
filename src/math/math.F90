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
!! $Id$

#include "global.h"

! This module is intended to contain "only mathematical" functions
! and procedures.

module math_m
  use global_m
  use lalg_basic_m
  use loct_math_m
  use messages_m
  use loct_m
  use blas_m

  implicit none

  private
  public ::                     &
    dconjugate_gradients,       &
    zconjugate_gradients,       &
    zqmr_sym,                   &
    dparker_traub,              &
    zparker_traub,              &
    dmatrix_newton_raphson,     &
    zmatrix_newton_raphson,     &
    dmatrix_inv_residual,       &
    zmatrix_inv_residual,       &
    ddot_product,               &
    zdot_product,               &
    quickrnd,                   &
    ylmr,                       &
    grylmr,                     &
    weights,                    &
    sort,                       &
    factorial,                  &
    hermite,                    &
    math_divisors,              &
    set_app_threshold,          &
    operator(.app.),            &
    math_xor,                   &
    dcross_product,             &
    zcross_product,             &
    hypot,                      &
    ddelta,                     &
    member,                     &
    make_idx_set,               &
    infinity_norm,              &
    interpolation_coefficients, &
    interpolate


  !------------------------------------------------------------------------------
  ! This is the common interface to a simple-minded polynomical interpolation
  ! procudure (simple use of the classical formula of Lagrange.
  interface interpolate
    module procedure dinterpolate_0, dinterpolate_1, dinterpolate_2
    module procedure zinterpolate_0, zinterpolate_1, zinterpolate_2
  end interface
  !------------------------------------------------------------------------------



  !------------------------------------------------------------------------------
  ! This is the common interface to a sorting routine.
  ! It performs the shell algorithm, not as fast as the quicksort for large numbers,
  ! but it seems that better for moderate numbers (around 100).
  ! Their possible interfaces are:
  !   subroutine sort(a [, ind] )
  !     FLOAT_OR_INTEGER, intent(inout) :: a(:)
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
    module procedure shellsort, dshellsort1, zshellsort1, &
      dshellsort2, zshellsort2, ishellsort
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
    CNST(    95.0)/CNST(288.0),                 &
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

  interface hypot
    real(8) function oct_hypotd(x, y)
      real(8), intent(in) :: x, y
    end function oct_hypotd
    real(4) function oct_hypotf(x, y)
      real(4), intent(in) :: x, y
    end function oct_hypotf
  end interface

  interface infinity_norm
    module procedure dinfinity_norm, zinfinity_norm
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
    rnd = real(iseed, REAL_PRECISION)/real(im, REAL_PRECISION)

  end subroutine quickrnd


  ! ---------------------------------------------------------
  ! Computes spherical harmonics ylm in the direction of vector r
  subroutine ylmr(x, y, z, li, mi, ylm)
    integer, intent(in) :: li, mi
    FLOAT, intent(in) :: x, y, z
    CMPLX, intent(out) :: ylm

    integer :: i
    FLOAT :: dx, dy, dz, r, plm, cosm, sinm, cosmm1, sinmm1, cosphi, sinphi

    ! if l=0, no calculations are required
    if (li == 0) then
      ylm = cmplx(CNST(0.282094791773878), M_ZERO, REAL_PRECISION)
      return
    end if

    r = sqrt(x*x + y*y + z*z)

    ! if r=0, direction is undefined => make ylm=0 except for l=0
    if (r == M_ZERO) then
      ylm = M_z0
      return
    end if

    dx = x/r; dy = y/r; dz = z/r

    ! get the associated Legendre polynomial (including the normalization factor)
    plm = loct_legendre_sphplm(li, abs(mi), dz)

    ! compute sin(|m|*phi) and cos(|m|*phi)
    r = sqrt(dx*dx + dy*dy)
    if (abs(r) < CNST(1e-20)) r = CNST(1e-20)

    cosphi = dx/r; sinphi = dy/r

    cosm = M_ONE; sinm = M_ZERO
    do i = 1, abs(mi)
      cosmm1 = cosm
      sinmm1 = sinm
      cosm = cosmm1*cosphi - sinmm1*sinphi
      sinm = cosmm1*sinphi + sinmm1*cosphi
    end do

    !And now ylm
    ylm = plm*cmplx(cosm, sinm, REAL_PRECISION)

    if (mi < 0) then
      ylm = conjg(ylm)
      do i = 1, abs(mi)
        ylm = -ylm
      end do
    end if

  end subroutine ylmr

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
  !  M -> number of grid points to be used in the approximation.
  !
  !  c(j,k,i) -> ith order derivative at kth-order approximation
  !              j=0,k: the coefficients acting of each point
  !
  !  side -> -1 left sided, +1 right sided, 0 centered (default)
  subroutine weights(N, M, cc, side)
    integer, intent(in) :: N, M
    FLOAT, intent(out) :: cc(0:M, 0:M, 0:N)
    integer, optional, intent(in) :: side

    integer :: i, j, k, mn, side_
    FLOAT :: c1, c2, c3, c4, c5, xi
    FLOAT, allocatable :: x(:)

    ALLOCATE(x(0:M), M+1)

    if (present(side)) then
      side_ = side
    else
      side_ = 0
    end if

    select case(side_)
    case(-1)
      ! grid-points for left-side finite-difference formulas on an equi.spaced grid
      mn = M
      x(:) = (/(-i,i=0,mn)/)
    case(+1)
      ! grid-points for right-side finite-difference formulas on an equi.spaced grid
      mn = M
      x(:) = (/(i,i=0,mn)/)
    case default
      ! grid-points for centered finite-difference formulas on an equi.spaced grid
      mn = M/2
      x(:) = (/0,(-i,i,i=1,mn)/)
    end select



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

    integer :: i,j,inc,n, indi, indj
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

          !workaround to a bug in itanium ifort
          !if(present(ind)) ind(j) = ind(j-inc)
          if(present(ind)) indj = ind(j-inc)
          if(present(ind)) ind(j) = indj

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
  ! Shell sort for integer arrays.
  subroutine ishellsort(a, ind)
    integer, intent(inout)           :: a(:)
    integer, intent(inout), optional :: ind(:)

    integer :: i,j,inc,n, indi, indj
    integer :: v

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
          !workaround to a bug in itanium ifort
          !if(present(ind)) ind(j) = ind(j-inc)
          if(present(ind)) indj = ind(j-inc)
          if(present(ind)) ind(j) = indj

          j=j-inc
          if (j <= inc) exit
        end do
        a(j)=v
        if(present(ind)) ind(j) = indi
      end do
      if (inc <= 1) exit
    end do

  end subroutine ishellsort


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
  
  FLOAT pure function ddelta(i, j)
    integer, intent(in) :: i
    integer, intent(in) :: j
    
    if ( i == j) then 
      ddelta = M_ONE
    else 
      ddelta = M_ZERO
    end if

  end function ddelta

  logical function math_xor(a, b)
    logical, intent(in) :: a
    logical, intent(in) :: b
    
    math_xor = ( a .or. b )
    if ( a .and. b ) math_xor = .false.
    
  end function math_xor


  ! ---------------------------------------------------------
  ! Construct out(1:length) = (/1, ..., n/) if in is not present,
  ! out(1:length) = in otherwise.
  subroutine make_idx_set(n, out, length, in)
    integer,           intent(in)  :: n
    integer,           pointer     :: out(:)
    integer,           intent(out) :: length
    integer, optional, intent(in)  :: in(:)

    integer :: i

    call push_sub('math.idx_set')

    if(present(in)) then
      length = ubound(in, 1)
      ALLOCATE(out(length), length)
      out = in
    else
      length = n
      ALLOCATE(out(length), length)
      do i = 1, length
        out(i) = i
      end do
    end if

    call pop_sub()
  end subroutine make_idx_set


  ! ---------------------------------------------------------
  ! Considers a(1:ubound(a, 1)) as an integer set and checks
  ! if n is a member of it.
  logical function member(n, a)
    integer, intent(in) :: n
    integer, intent(in) :: a(:)

    integer :: i

    call push_sub('states_inc.member')

    member = .false.

    do i = 1, ubound(a, 1)
      if(a(i).eq.n) then
        member = .true.
        exit
      end if
    end do

    call pop_sub()
  end function member

! ---------------------------------------------------------
! QMR (quasi-minimal residual) algorithm for complex symmetric matrices
! algorithm taken from:
! Parallel implementation of efficient preconditioned linear solver for
! grid-based applications in chemical physics. II: QMR linear solver
! Appendix A. Simplified QMR algorithm
subroutine zqmr_sym(np, x, b, op, prec, iter, residue, threshold, showprogress)
  integer, intent(in)             :: np
  CMPLX,  intent(inout)           :: x(:)
  CMPLX,  intent(in)              :: b(:)
  interface
    subroutine op(x, y)
      CMPLX, intent(in)           :: x(:)
      CMPLX, intent(out)          :: y(:)
    end subroutine op
  end interface
  interface
    subroutine prec(x, y)
      CMPLX, intent(in)           :: x(:)
      CMPLX, intent(out)          :: y(:)
    end subroutine prec
  end interface
  integer,          intent(inout) :: iter
  FLOAT, optional, intent(out)    :: residue
  FLOAT, optional,  intent(in)    :: threshold
  logical, optional, intent(in)   :: showprogress

  CMPLX, allocatable  :: r(:), v(:), z(:), q(:), p(:), deltax(:), deltar(:)
  CMPLX               :: eta, delta, epsilon, beta
  FLOAT               :: rho, xsi, gamma, alpha, theta, threshold_, res, oldtheta, oldgamma, oldrho
  integer             :: max_iter, err
  logical             :: showprogress_
  FLOAT               :: log_res, log_thr
  integer             :: ilog_res, ilog_thr

  call push_sub('math_cg_inc.zqmr_sym_x')

  if(present(threshold)) then
    threshold_ = threshold
  else
    threshold_ = CNST(1.0e-6)
  end if

  if(present(showprogress)) then
    showprogress_ = showprogress
  else
    showprogress_ = .false.
  end if

  ALLOCATE( r(np), np)
  ALLOCATE( v(np), np)
  ALLOCATE( z(np), np)
  ALLOCATE( q(np), np)
  ALLOCATE( p(np), np)
  ALLOCATE( deltax(np), np)
  ALLOCATE( deltar(np), np)

  ! use v as temp var
  call lalg_copy(np, b, r)
  call op(x, v)
  call lalg_axpy(np, -M_z1, v, r)
  call lalg_copy(np, r, v)
  rho = lalg_nrm2(np, v)
  call prec(v, z)
  xsi = lalg_nrm2(np, z)
  gamma = M_ONE
  eta = -M_z1
  alpha = M_ONE
  theta = M_ZERO
  
  ! initialize progress bar
  log_thr = -log(threshold_)
  ilog_thr = M_TEN**2*log_thr
  if (showprogress_) call loct_progress_bar(-1, ilog_thr)

  max_iter = iter
  iter     = 0
  err      = 0
  do while(iter < max_iter)
    iter = iter + 1
    if((abs(rho) < M_EPSILON) .or. (abs(xsi) < M_EPSILON)) then
      err = 1
      exit
    end if
    alpha = alpha*xsi/rho
    call lalg_scal(np, M_z1/rho, v)
    call lalg_scal(np, M_z1/xsi, z)
    delta = blas_dotu(np, v(1), 1, z(1), 1)
    if (abs(delta) < M_EPSILON) then
      err = 2
      exit
    end if
    if (iter == 1) then
      call lalg_copy(np, z, q)
    else
      call lalg_scal(np, -rho*delta/epsilon, q)
      call lalg_axpy(np, M_z1, z, q)
    end if
    call op(q, p)
    call lalg_scal(np, alpha, p)
    epsilon = blas_dotu(np, q(1), 1, p(1), 1)
    if (abs(epsilon) < M_EPSILON) then
      err = 3
      exit
    end if
    beta = epsilon/delta
    call lalg_scal(np, -beta, v)
    call lalg_axpy(np, M_z1, p, v)
    oldrho = rho
    rho = lalg_nrm2(np, v)
    call prec(v, z)
    call lalg_scal(np, M_z1/alpha, z)
    xsi = lalg_nrm2(np, z)
    oldtheta = theta
    theta = rho/(gamma*abs(beta))
    oldgamma = gamma
    gamma = M_ONE/sqrt(M_ONE+theta**2)
    if (abs(gamma) < M_EPSILON) then
      err = 4
      exit
    end if
    eta = -eta*oldrho*gamma**2/(beta*oldgamma**2)
    if (iter == 1) then
      call lalg_copy(np, q, deltax)
      call lalg_scal(np, eta*alpha, deltax)
      call lalg_copy(np, p, deltar)
      call lalg_scal(np, eta, deltar)
    else
      call lalg_scal(np, (oldtheta*gamma)**2, deltax)
      call lalg_axpy(np, eta*alpha, q, deltax)
      call lalg_scal(np, (oldtheta*gamma)**2, deltar)
      call lalg_axpy(np, eta, p, deltar)
    end if
    call lalg_axpy(np, M_z1, deltax, x)
    call lalg_axpy(np, -M_z1, deltar, r)
    res = lalg_nrm2(np,r)/lalg_nrm2(np,x)
    log_res = -log(res)
    if (log_res < 0) log_res = 0
    ilog_res = M_TEN**2*log_res
    if (showprogress_)  call loct_progress_bar(ilog_res, ilog_thr)
    if (res < threshold_) exit
    !write(*,*) res
  end do
  if (iter.eq.max_iter) err = 5

  select case(err)
  case(0)
    ! converged
  case(1)
    write(message(1), '(a)') "QMR failure, can't continue: b or P*b is the zero vector!"
    call write_fatal(1)
  case(2)
    write(message(1), '(a)') "QMR failure, can't continue: v^T*z is zero!"
    call write_fatal(1)
  case(3)
    write(message(1), '(a)') "QMR failure, can't continue: q^T*p is zero!"
    call write_fatal(1)
  case(4)
    write(message(1), '(a)') "QMR failure, can't continue: gamma is zero!"
    call write_fatal(1)
  case(5)
    write(message(1), '(a)') "QMR Solver not converged!"
    call write_warning(1)
  end select
  if (showprogress_) write(*,*) ''

  if(present(residue)) residue = res

  deallocate(r, v, z, q, p, deltax, deltar)

  call pop_sub()
end subroutine zqmr_sym

subroutine interpolation_coefficients(nn, xa, xx, cc)
  integer, intent(in)  :: nn    ! the number of points and coefficients
  FLOAT,   intent(in)  :: xa(:) ! the nn points where we know the function
  FLOAT,   intent(in)  :: xx    ! the point where we want the function
  FLOAT,   intent(out) :: cc(:)  ! the coefficients

  integer :: ii, kk

  do ii = 1, nn
    cc(ii) = M_ONE
    do kk = 1, nn
      if(kk .eq. ii) cycle
      cc(ii) = cc(ii)*(xx - xa(kk))/(xa(ii) - xa(kk))
    end do
  end do
  
end subroutine interpolation_coefficients

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

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

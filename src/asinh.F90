! Asinh function taken from netlib
! april 1977 edition.  w. fullerton, c3, los alamos scientific lab.
function asinh (x)
  real(r8), intent(in) :: x
  real(r8) :: asinh

  real(r8), parameter :: asnhcs(20) = (/ &
      -.12820039911738186_r8, -.058811761189951768_r8, &
      .004727465432212481_r8, -.000493836316265361_r8, &
      .000058506207058557_r8, -.000007466998328931_r8, &
      .000001001169358355_r8, -.000000139035438587_r8, &
      .000000019823169483_r8, -.000000002884746841_r8, &
      .000000000426729654_r8, -.000000000063976084_r8, &
      .000000000009699168_r8, -.000000000001484427_r8, &
      .000000000000229037_r8, -.000000000000035588_r8, &
      .000000000000005563_r8, -.000000000000000874_r8, &
      .000000000000000138_r8, -.000000000000000021_r8 /)

  real(r8), parameter :: aln2 = 0.69314718055994530942_r8

! series for asnh       on the interval  0.          to  1.00000d+00
!                                        with weighted error   2.19e-17
!                                         log weighted error  16.66
!                               significant figures required  15.60
!                                    decimal places required  17.31
!

  integer, save :: nterms = 0
  real(r8), save :: xmax = 0._r8, sqeps = 0._r8

  real(r8) :: y

  if (nterms == 0) then
    nterms = inits (asnhcs, 20, 0.1_r8*eps_r8)
    sqeps = sqrt (eps_r8)
    xmax = 1._r8/sqeps
  end if

  y = abs(x)
  if (y.le.1.0) then
    asinh = x
    if (y.gt.sqeps) asinh = x*(1._r8 + csevl (2.*x*x-1._r8, asnhcs, nterms))
    return
  end if

  if (y.lt.xmax) asinh = log (y + sqrt(y**2+1._r8))
  if (y.ge.xmax) asinh = aln2 + log(y)
  asinh = sign (asinh, x)

  return
end function asinh


! april 1977 version.  w. fullerton, c3, los alamos scientific lab.
! evaluate the n-term chebyshev series cs at x.  adapted from
! r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).  also see fox
! and parker, chebyshev polys in numerical analysis, oxford press, p.56.
!
!             input arguments --
! x      value at which the series is to be evaluated.
! cs     array of n terms of a chebyshev series.  in eval-
!        uating cs, only half the first coef is summed.
! n      number of terms in array cs.

function csevl (x, cs, n)
  integer, intent(in) :: n
  real(r8), intent(in) :: x
  real(r8), intent(in) :: cs(n)
  real(r8) :: csevl
  
  integer i, ni
  real(r8) :: b0, b1, b2, twox

  if (n.lt.1) then
    message(1) = "Math::csevl: number of terms le 0"
    call write_fatal(1)
  end if
  if(n.gt.1000) then
    message(1) = "Math::csevl: number of terms gt 1000"
    call write_fatal(1)
  end if

  if (x.lt.(-1.1) .or. x.gt.1.1) then
    message(1) = "Math::csevl:  x outside (-1,+1)"
    call write_fatal(1)
  end if

  b1 = 0._r8
  b0 = 0._r8
  twox = 2.*x
  do i = 1, n
    b2 = b1
    b1 = b0
    ni = n + 1 - i
    b0 = twox*b1 - b2 + cs(ni)
  end do

  csevl = 0.5_r8 * (b0 - b2)

  return
end function csevl


! april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
! initialize the orthogonal series so that inits is the number of terms
! needed to insure the error is no larger than eta.  ordinarily, eta
! will be chosen to be one-tenth machine precision.
!
!             input arguments --
! os     array of nos coefficients in an orthogonal series.
! nos    number of coefficients in os.
! eta    requested accuracy of series.
function inits (os, nos, eta)
  integer, intent(in) :: nos
  real(r8), intent(in) :: os(nos)
  real(r8), intent(in) :: eta
  integer :: inits

  integer :: i, ii
  real(r8) :: err

  if (nos.lt.1) then
    message(1) = "Math::inits: number of terms le 0"
    call write_fatal(1)
  end if

  err = 0.
  do ii=1,nos
    i = nos + 1 - ii
    err = err + abs(os(i))
    if (err.gt.eta) exit
  end do

  if (i.eq.nos) then
    message(1) = "Math::inits: eta may be too small"
    call write_warning(1)
  end if

  inits = i

  return
end function inits

#include "config.h"

module spline
use global
use math, only : ratint

implicit none

type spline_type
  real(r8) :: delt ! distance between two points
  real(r8), pointer :: func(:), der2(:)
end type spline_type

integer, parameter :: ntbmax = 200

contains

subroutine spline_init(spl)
  type(spline_type), intent(inout) :: spl

  allocate(spl%func(ntbmax), spl%der2(ntbmax))
end subroutine spline_init

subroutine spline_end(spl)
  type(spline_type), intent(inout) :: spl

  if(associated(spl%func)) then ! sanity check
    deallocate(spl%func, spl%der2)
    nullify(spl%func, spl%der2)
  end if
end subroutine spline_end

! this function converts ffit from a logarithmic mesh to a linear mesh
! and then fits a spline to it
! parameters are:
!   r_max: maximum value of r for which the function is defined
!   a, b:  parameters for the log grid
!   rofi:  r as a function of index i
!   ffit: function to fit
!   der1, dern: derivatives at points 1 and n
!   spl: spline (output)
subroutine fit_spline(ffit, r_max, a, b, rofi, der1, dern, spl)
  real(r8), intent(in) :: a, b, der1, dern, r_max
  real(r8), intent(IN) :: ffit(:), rofi(:)
  type(spline_type), intent(out) :: spl

  real(r8) :: dy, r
  integer :: itb, nrc, nr, nmin, nmax, niden
  integer, parameter :: npoint = 8
  
  nrc = nint(log(r_max/b + 1.0_r8)/a) + 1
  spl%delt = r_max/dble(ntbmax - 1)

  ! apparently this converts from the logarithm mesh into a linear one
  do itb = 1, ntbmax
    r     = spl%delt*(itb - 1)
    nr    = nint(log(r/b + 1.0_r8)/a) + 1
    nmin  = max(1, nr - npoint)
    nmax  = min(nrc, nr + npoint)
    niden = nmax - nmin + 1
    call ratint(rofi(nmin:nmax), ffit(nmin:nmax), niden, r, spl%func(itb), dy)
  enddo

  call second_der(spl, der1, dern)

  return
end subroutine fit_spline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! An Adapted Numerical Recipes Subroutine spline.f  p.109
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine second_der(spl, yp1, ypn)
  type(spline_type), intent(inout) :: spl
  real(r8), intent(IN) :: yp1, ypn
  
  integer :: i, k
  real(r8) :: p, qn, un, sig
  real(r8), allocatable :: u(:)
  
  allocate(u(ntbmax)) ! auxiliar array
  
  if(yp1.gt..99e30_r8) then
    spl%der2(1) = 0.0_r8
    u(1)  = 0.0_r8
  else
    spl%der2(1) = -0.5_r8
    u(1)  = (3.0d0/spl%delt)*((spl%func(2)-spl%func(1))/spl%delt - yp1)
  end if
  
  do i = 2, ntbmax - 1
    sig   = 0.5_r8
    p     = sig*spl%der2(i - 1) + 2.0_r8
    spl%der2(i) = (sig - 1.0_r8)/p
    u(i)  = (3.0_r8*(spl%func(i+1) + spl%func(i-1) - 2.0_r8*spl%func(i))/(spl%delt**2)  &
         -sig*u(i-1))/p
  end do
  if (ypn.gt..99e30_r8) then
    qn = 0.0d0
    un = 0.0d0
  else
    qn = 0.5d0
    un = (3.0d0/spl%delt)*(ypn - (spl%func(ntbmax) - spl%func(ntbmax-1))/spl%delt)
  end if
  
  spl%der2(ntbmax) = (un - qn*u(ntbmax-1))/(qn*spl%der2(ntbmax-1) + 1.0_r8)
  do k = ntbmax-1, 1, -1
    spl%der2(k) = spl%der2(k)*spl%der2(k+1) + u(k)
  end do
  
  deallocate(u)
  
  return
end subroutine second_der

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! An Adapted Numerical Recipes Subroutine splint.f                            !
! make as a function call.                                                    !
! If the desired value is above the limits of the fitted function             !
! the function returns zero value.         A. Rubio  (1998)                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function splint(spl, x)
  real(r8) :: splint 

  type(spline_type), intent(IN) :: spl
  real(r8), intent(IN) :: x
  
  real(r8) :: a, b
  integer :: nlo, nhi

  if (spl%delt == 0.0_r8) then
    message(1) = 'splint: bad delt input.'
    call write_fatal(1)
  end if

  nlo = int(x/spl%delt) + 1
  nhi = nlo + 1
  if(nhi.gt.ntbmax) then
    splint = 0._r8
    return
  endif

  a = nhi - x/spl%delt - 1._r8
  b = 1._r8 - a
  splint = a*spl%func(nlo) + b*spl%func(nhi) + &
       ((a**3 - a)*spl%der2(nlo) + (b**3 - b)*spl%der2(nhi))*(spl%delt**2)/6._r8

  return
end function splint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return the first derivative of a spline                                     !
! If the desired value is above the limits of the fitted function             !
! the function returns zero value.         M. Marques (2000)                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function spl_der(spl, x)
  real(r8) :: spl_der 

  type(spline_type), intent(IN) :: spl
  real(r8), intent(IN) :: x
  
  real(r8) :: a, b
  integer :: nlo, nhi

  if (spl%delt == 0.0d0) then
    message(1) = 'spl_der: bad delt input.'
    call write_fatal(1)
  end if

  nlo = int(x/spl%delt) + 1
  nhi = nlo + 1
  if(nhi.gt.ntbmax) then
    spl_der = 0._r8
    return
  endif
  
  a = nhi - x/spl%delt - 1._r8
  b = 1._r8 - a
  spl_der = -(spl%func(nlo) - spl%func(nhi))/spl%delt - &
      ((3._r8*a**2 - 1._r8)*spl%der2(nlo) - (3._r8*b**2 - 1._r8)*spl%der2(nhi)) * &
      spl%delt/6._r8

  return
end function spl_der

end module spline

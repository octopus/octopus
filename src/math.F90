#include "config.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module is intended to contain "only mathematical" functions and        !
!	procedures.                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  module math


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Public procedures headers...                                                !
!                                                                             !
!  function stepf(x)                                                          !
!    implicit none                                                            !
!    real(r8), intent(in) :: x                                                !
!    real(r8) :: stepf                                                        !
!  end function stepf                                                         !
!                                                                             !
!  function derfc(x)                                                          !
!    implicit none                                                            !
!    real(r8), intent(in) :: x                                                !
!    real(r8) :: derfc                                                        !
!  end function derfc                                                         !
!                                                                             !
!  function derf(x)                                                           !
!    implicit none                                                            !
!    real(r8), intent(in) :: x                                                !
!    real(r8) :: derf                                                         !
!  end function derf                                                          !
!                                                                             !
!  subroutine polint                                                          !
!    implicit none                                                            !
!    integer, intent(in) :: n                                                 !
!    real(r8), intent(in), dimension(n) ::  xa , ya                           !
!    real(r8), intent(inout) :: x, y, dy                                      !
!  end subroutine polint                                                      !
!                                                                             !
!  subroutine ratint(xa,ya,n,x,y,dy)                                          !
!    implicit  none                                                           !
!    integer, intent(in) :: n                                                 !
!    real(r8), intent(in), dimension(n) :: xa,ya                              !
!    real(r8), intent(inout) :: x,y,dy                                        !
!  end subroutine ratint                                                      !
!                                                                             !
!  subroutine ylmr(x,y,z,li,mi,ylm)                                           !
!    implicit none                                                            !
!    integer, intent(in) :: li, mi                                            !
!    real(r8), intent(in) :: x,y,z                                            !
!    real(r8), intent(out) :: ylm                                             !
!  end subroutine ylmr                                                        !
!                                                                             !
!  subroutine grylmr(x,y,z,li,mi,ylm,grylm)                                   !
!    implicit none                                                            !
!    integer, intent(in) :: li, mi                                            !
!    real(r8), intent(in) :: x,y,z                                            !
!    real(r8), intent(out) :: ylm,grylm(3)                                    !
!  end subroutine grylmr                                                      !
!                                                                             !
!  subroutine sort(n,arrin,indx)                                              !
!    integer, intent(in) :: n                                                 !
!    real(r8), intent(in) :: arrin(n)                                         !
!    integer, intent(out):: indx(n)                                           !
!  end subroutine sort                                                        !
!                                                                             !
!  subroutine gaussj(a,n,np,b,m,mp)                                           !
!    integer, ::  m,mp,n,np                                                   !
!    real(r8), :: a(np,np),b(np,mp)                                           !
!  end subroutine gaussj                                                      !
!                                                                             !
!  subroutine high_derivative(norder,cn)                                      !
!    implicit none                                                            !
!    integer, intent(in) :: norder                                            !
!    real(r8), intent(out) :: cn(0:norder)                                    !
!  end subroutine high_derivative                                             !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
  use global

  implicit none

  private
  public :: stepf,                                                            &
            derf,                                                             &
            derfc,                                                            &
            ratint,                                                           &
            polint,                                                           &
            ylmr,                                                             &
            grylmr,                                                           &
            sort,                                                             &
            gaussj,                                                           &
            high_derivative, &
            weights, &
            asinh

contains

#include "asinh.F90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Step function, needed for definition of fermi function.                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function stepf(x)

  real(r8), intent(in) ::  x
  real(r8) :: stepf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Algorithm.                                                                  !
! Note:                                                                       ! 
!     complementary error function. ref: fu & ho, prb 28, 5480 (1983)         !
!     stepf=derfc(x)                                                          !
!     improved step function. ref: methfessel & paxton prb40 (15/aug/89)      !
!     parameter (c=0.5641895835_r8)                                           !
!     stepf=derfc(x)-c*x*dexp(-x*x)                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (x.gt.100.0_r8) then
     stepf = 0.0_r8
  elseif (x.lt.-100.0_r8) then
     stepf = 2.0_r8
  else
     stepf = 2.0_r8 / ( 1.0_r8 + exp(x) )
  endif


  end function stepf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Complementary step function.                                                !
!	complementary error function from "numerical recipes"                 !
!	note: single precision accuracy                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function derfc (x)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Arguments and local variables                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit none

  real(r8), intent(in) :: x
  real(r8) :: derfc

  real(r8):: t, z


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Algorithm                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  z=abs(x)
  t=1.0_r8/(1.0_r8+0.5_r8*z)

  derfc=t*exp(-(z*z)-1.26551223_r8+t*(1.00002368_r8+t*(0.37409196_r8+    &
        t*(0.09678418_r8+t*(-0.18628806_r8+t*(0.27886807_r8+             &
        t*(-1.13520398_r8+t*(1.48851587_r8+t*(-0.82215223_r8+            &
        t*.17087277_r8)))))))))

  if (x < 0.0_r8) derfc=2.0_r8-derfc

  end function derfc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Error function.                                                             !
!	error function from "numerical recipes"                               !
!	note: single precision accuracy                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function derf (x)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Arguments and local variables                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit none

  real(r8), intent(in) :: x
  real(r8) :: derf

  real(r8) :: t,z

  z=abs(x)
  t=1.0_r8/(1.0_r8+0.5_r8*z)

  derf= t*exp(-(z*z)-1.26551223_r8+t*(1.00002368_r8+t*(0.37409196_r8+   &
        t*(0.09678418_r8+t*(-0.18628806_r8+t*(0.27886807_r8+            &
        t*(-1.13520398_r8+t*(1.48851587_r8+t*(-0.82215223_r8+           &
        t*.17087277_r8)))))))))

  if (x < 0.0_r8) derf=2.0_r8-derf

  derf = 1.0_r8 - derf


  end function derf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!














!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! An Adapted Numerical Recipes ratint.f and polint.f                          !
! interpolation subroutines  (1998)                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

  subroutine ratint(xa,ya,n,x,y,dy)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Local variables                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit  none

  integer, intent(in) :: n
  real(r8), intent(IN), dimension(n) :: xa,ya
  real(r8), intent(inout) :: x,y,dy

 
  real(r8), dimension(n) :: c,d
  real(r8), parameter :: tiny = 1.d-20

  integer :: i, m, ns 
  real(r8) :: dd, h, hh, t, w


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Algorithm                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ns=1
  hh=abs(x-xa(1))
  do i=1,n
     h=dabs(x-xa(i))
     if (h < tiny)then
        y=ya(i)
        dy=0.0d0
        return
     else if (h < hh) then
        ns=i
        hh=h
     endif
     c(i)=ya(i)
     d(i)=ya(i)+tiny
  enddo

  y=ya(ns)
  ns=ns-1
  do m=1,n-1
     do i=1,n-m
        w=c(i+1)-d(i)
        h=xa(i+m)-x
        t=(xa(i)-x)*d(i)/h
        dd=t-c(i+1)
        if(dd == 0.0d0)goto 100
        dd=w/dd
        d(i)=c(i+1)*dd
        c(i)=t*dd
     enddo
     if (2*ns.lt.n-m)then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  enddo

! as rational interpolation does not converge,we try with a polynomial one

100   call polint(xa,ya,n,x,y,dy)

  return

  end subroutine ratint


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine polint                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine polint(xa,ya,n,x,y,dy)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Arguments and local variables                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit none

  integer, intent(in) :: n
  real(r8), intent(IN), dimension(n) ::  xa , ya
  real(r8), intent(inout) :: x, y, dy

  integer :: i, m, ns

  real(r8), dimension(n) :: c,d
  real(r8) :: den, dif, dift, ho, hp, w


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Algorithm                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ns=1
  dif=dabs(x-xa(1))
  do i=1,n 
     dift=dabs(x-xa(i))
     if (dift < dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)
     d(i)=ya(i)
  enddo
  y=ya(ns)
  ns=ns-1
  do m=1,n-1
     do i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den == 0.0d0)stop 'polint: den=0'
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
     enddo
     if (2*ns.lt.n-m)then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  enddo


  return

  end subroutine polint


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is a Numerical Recipes based subroutine.                               !
! Computes real spherical harmonics ylm in the direction of vector r:         !
!    ylm = c * plm( cos(theta) ) * sin(m*phi)   for   m <  0                  !
!    ylm = c * plm( cos(theta) ) * cos(m*phi)   for   m >= 0                  !
! with (theta,phi) the polar angles of r, c a positive normalization          !
! constant and plm associated legendre polynomials.                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ylmr(x,y,z,li,mi,ylm)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Arguments and local variables                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
  implicit none


  integer, intent(in) :: li, mi
  real(r8), intent(in) :: x,y,z
  real(r8), intent(out) :: ylm

  integer, parameter :: lmaxd = 20
  real(r8), parameter :: tiny=1.d-30,                                         &
                         zero=0.0_r8,                                         &
                         half=0.5_r8,                                         &
                         one=1.0_r8,                                          &
                         two=2.0_r8,                                          &
                         three=3.0_r8,                                        &
                         six=6.0_r8
      
  real(r8), save :: c(0:(lmaxd+1)*(lmaxd+1))

  integer :: i, ilm0, l, m, mabs
  integer, save :: lmax = -1
  real(r8) :: cmi, cosm, cosmm1, cosphi, dphase, dplg, fac,                   &
              fourpi, plgndr, phase, pll, pmm, pmmp1, sinm,                   &
              sinmm1, sinphi, r2, rsize, Rx, Ry, Rz, xysize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! evaluate normalization constants once and for all                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (li.gt.lmax) then
     fourpi=two**4*datan(one)
     do l=0,li
        ilm0=l*l+l
        do m=0,l
           fac=(2*l+1)/fourpi
           do i=l-m+1,l+m
                  fac=fac/i
           enddo 
           c(ilm0+m)=sqrt(fac)
! next line because ylm are real combinations of m and -m
           if (m.ne.0) c(ilm0+m)=c(ilm0+m)*sqrt(two)
           c(ilm0-m)=c(ilm0+m)
        enddo 
     enddo 
     lmax=li
  endif

! if l=0, no calculations are required

  if (li.eq.0) then
     ylm=c(0)
     return
  endif


! if r=0, direction is undefined => make ylm=0 except for l=0
  r2=x**2+y**2+z**2
  if (r2.lt.tiny) then
     ylm=zero
     return
  endif
  rsize=sqrt(r2)

  Rx=x/rsize
  Ry=y/rsize
  Rz=z/rsize

! explicit formulas for l=1 and l=2
  if (li.eq.1) then
     if(mi.eq.-1) then 
        ylm=(-c(1))*Ry
     elseif(mi.eq.0) then
        ylm= c(2)*Rz
     elseif(mi.eq.1) then
        ylm=(-c(3))*Rx
     endif
    return
  endif
         
  if (li.eq.2) then
     if(mi.eq.-2)then
        ylm= c(4)*six*Rx*Ry
     elseif(mi.eq.-1)then
        ylm=(-c(5))*three*Ry*Rz
     elseif(mi.eq.0) then
        ylm= c(6)*half*(three*Rz*Rz-one)
     elseif(mi.eq.1) then
        ylm=(-c(7))*three*Rx*Rz
     elseif(mi.eq.2) then
        ylm= c(8)*three*(Rx*Rx-Ry*Ry)
     endif
     return
  endif

! general algorithm based on routine plgndr of 'numerical recipes'

  mabs=abs(mi)
  xysize=sqrt(max(Rx*Rx+Ry*Ry,tiny))
  cosphi=Rx/xysize
  sinphi=Ry/xysize
  cosm=one
  sinm=zero
  do m=1,mabs
     cosmm1=cosm
     sinmm1=sinm
     cosm=cosmm1*cosphi-sinmm1*sinphi
     sinm=cosmm1*sinphi+sinmm1*cosphi
  enddo

     
  if(mi.lt.0) then
    phase=sinm
  else
    phase=cosm
  endif

  pmm=one
  fac=one
  if(mabs.gt.zero) then
    do i=1,mabs
       pmm=(-pmm)*fac*xysize
       fac=fac+two
    enddo 
  endif

  if(li.eq.mabs) then
     plgndr=pmm
  else   
     pmmp1=Rz*(2*mabs+1)*pmm
     if(li.eq.mabs+1) then
        plgndr=pmmp1
     else 
        do l=mabs+2,li
           pll=(Rz*(2*l-1)*pmmp1-(l+mabs-1)*pmm)/(l-mabs)
           pmm=pmmp1
           pmmp1=pll
        enddo 
        plgndr=pll
     endif
  endif         
        

  ilm0=li*li+li
  cmi=c(ilm0+mi)
  ylm=cmi*plgndr*phase



  return

  end subroutine ylmr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is a Numerica Recipes based subroutine                                 !
! computes real spherical harmonics ylm in the direction of vector r:         !
!    ylm = c * plm( cos(theta) ) * sin(m*phi)   for   m <  0                  !
!    ylm = c * plm( cos(theta) ) * cos(m*phi)   for   m >= 0                  !
! with (theta,phi) the polar angles of r, c a positive normalization          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine grylmr(x,y,z,li,mi,ylm,grylm)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Arguments and local variables                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     
  implicit none

  integer, intent(in) :: li, mi
  real(r8), intent(in) :: x,y,z
  real(r8), intent(out) :: ylm,grylm(3)

  integer, parameter :: lmaxd = 20
  real(r8), parameter :: tiny=1.d-30,                                         &
                         zero=0.0_r8,                                         &
                         half=0.5_r8,                                         &
                         one=1.0_r8,                                          &
                         two=2.0_r8,                                          &
                         three=3.0_r8,                                        &
                         six=6.0_r8

  integer :: i, ilm0, l, m, mabs
  integer, save :: lmax = -1

  real(r8) :: cmi, cosm, cosmm1, cosphi, dphase, dplg, fac,                   &
              fourpi, plgndr, phase, pll, pmm, pmmp1, sinm,                   &
              sinmm1, sinphi, r2, rsize, Rx, Ry, Rz, xysize
  real(r8), save :: c(0:(lmaxd+1)*(lmaxd+1))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! evaluate normalization constants once and for all                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (li.gt.lmax) then
     fourpi=two**4*datan(one)
     do l=0,li
        ilm0=l*l+l
        do m=0,l
           fac=(2*l+1)/fourpi
           do i=l-m+1,l+m
              fac=fac/i
           enddo 
           c(ilm0+m)=sqrt(fac)
!          next line because ylm are real combinations of m and -m
           if (m.ne.0) c(ilm0+m)=c(ilm0+m)*sqrt(two)
           c(ilm0-m)=c(ilm0+m)
        enddo 
     enddo 
     lmax=li
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if l=0, no calculations are required                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (li.eq.0) then
     ylm=c(0)
     grylm(1)=0.0
     grylm(2)=0.0
     grylm(3)=0.0
     return
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if r=0, direction is undefined => make ylm=0 except for l=0                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  r2=x**2+y**2+z**2
  if (r2.lt.tiny) then
     ylm=zero
     grylm(1)=zero
     grylm(2)=zero
     grylm(3)=zero
     return
  endif
  rsize=sqrt(r2)

  Rx=x/rsize
  Ry=y/rsize
  Rz=z/rsize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! explicit formulas for l=1 and l=2                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (li.eq.1) then
     if(mi.eq.-1) then 
       ylm=(-c(1))*Ry
       grylm(1)= c(1)*Rx*Ry/rsize
       grylm(2)=(-c(1))*(one-Ry*Ry)/rsize
       grylm(3)= c(1)*Rz*Ry/rsize 
     elseif(mi.eq.0) then
       ylm= c(2)*Rz
       grylm(1)=(-c(2))*Rx*Rz/rsize
       grylm(2)=(-c(2))*Ry*Rz/rsize
       grylm(3)= c(2)*(one-Rz*Rz)/rsize
     elseif(mi.eq.1) then
       ylm=(-c(3))*Rx
       grylm(1)=(-c(3))*(one-Rx*Rx)/rsize
       grylm(2)= c(3)*Ry*Rx/rsize
       grylm(3)= c(3)*Rz*Rx/rsize
     endif
     return
  endif
         
  if (li.eq.2) then
     if(mi.eq.-2)then
       ylm= c(4)*six*Rx*Ry
       grylm(1)=(-c(4))*six*(two*Rx*Rx*Ry-Ry)/rsize
       grylm(2)=(-c(4))*six*(two*Ry*Rx*Ry-Rx)/rsize
       grylm(3)=(-c(4))*six*(two*Rz*Rx*Ry)/rsize
     elseif(mi.eq.-1)then
       ylm=(-c(5))*three*Ry*Rz
       grylm(1)= c(5)*three*(two*Rx*Ry*Rz)/rsize
       grylm(2)= c(5)*three*(two*Ry*Ry*Rz-Rz)/rsize
       grylm(3)= c(5)*three*(two*Rz*Ry*Rz-Ry)/rsize
     elseif(mi.eq.0) then
       ylm= c(6)*half*(three*Rz*Rz-one)
       grylm(1)=(-c(6))*three*(Rx*Rz*Rz)/rsize
       grylm(2)=(-c(6))*three*(Ry*Rz*Rz)/rsize
       grylm(3)=(-c(6))*three*(Rz*Rz-one)*Rz/rsize
     elseif(mi.eq.1) then
       ylm=(-c(7))*three*Rx*Rz
       grylm(1)= c(7)*three*(two*Rx*Rx*Rz-Rz)/rsize
       grylm(2)= c(7)*three*(two*Ry*Rx*Rz)/rsize
       grylm(3)= c(7)*three*(two*Rz*Rx*Rz-Rx)/rsize
     elseif(mi.eq.2) then
       ylm= c(8)*three*(Rx*Rx-Ry*Ry)
       grylm(1)=(-c(8))*six*(Rx*Rx-Ry*Ry-one)*Rx/rsize
       grylm(2)=(-c(8))*six*(Rx*Rx-Ry*Ry+one)*Ry/rsize
       grylm(3)=(-c(8))*six*(Rx*Rx-Ry*Ry)*Rz/rsize
     endif

     return
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! general algorithm based on routine plgndr of 'numerical recipes'            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  mabs=abs(mi)
  xysize=sqrt(max(Rx*Rx+Ry*Ry,tiny))
  cosphi=Rx/xysize
  sinphi=Ry/xysize
  cosm=one
  sinm=zero
  do m=1,mabs
     cosmm1=cosm
     sinmm1=sinm
     cosm=cosmm1*cosphi-sinmm1*sinphi
     sinm=cosmm1*sinphi+sinmm1*cosphi
  enddo

     
  if(mi.lt.0) then
    phase=sinm
    dphase=mabs*cosm
  else
    phase=cosm
    dphase=(-mabs)*sinm
  endif

  pmm=one
  fac=one

  if(mabs.gt.zero) then
    do i=1,mabs
       pmm=(-pmm)*fac*xysize
      fac=fac+two
    enddo
  endif

  if(li.eq.mabs) then
    plgndr=pmm
    dplg=(-li)*Rz*pmm/(xysize**2)
  else   
    pmmp1=Rz*(2*mabs+1)*pmm
    if(li.eq.mabs+1) then
      plgndr=pmmp1
      dplg=-((li*Rz*pmmp1-(mabs+li)*pmm)/(xysize**2))
    else 

      do l=mabs+2,li
         pll=(Rz*(2*l-1)*pmmp1-(l+mabs-1)*pmm)/(l-mabs)
         pmm=pmmp1
         pmmp1=pll
      enddo 
      plgndr=pll
      dplg=-((li*Rz*pll-(l+mabs-1)*pmm)/(xysize**2))
    endif
  endif         
        

  ilm0=li*li+li
  cmi=c(ilm0+mi)
  ylm=cmi*plgndr*phase
  grylm(1)=(-cmi)*dplg*Rx*Rz*phase/rsize                                      &
           -cmi*plgndr*dphase*Ry/(rsize*xysize**2)

  grylm(2)=(-cmi)*dplg*Ry*Rz*phase/rsize                                      &
           +cmi*plgndr*dphase*Rx/(rsize*xysize**2)

  grylm(3)= cmi*dplg*(one-Rz*Rz)*phase/rsize
   

  return
      

  end subroutine grylmr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sorts an arrat by the heapsort method,                                      !
! W. H. Preuss et al. (Numerical Recipes)                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine sort(n,arrin,indx)


  use global


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Arguments and local variables                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  integer, intent(in) :: n
  real(r8), intent(IN) :: arrin(n)
  integer, intent(out):: indx(n)

  real(r8) :: q
  integer :: i, indxt, ir, j, l


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Algorithm                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  do j = 1, n
         indx(j) = j
  enddo
 
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Numerical Recipes routine. Explained in the book                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine gaussj(a,n,np,b,m,mp)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Arguments and variables                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  integer ::  m,mp,n,np
  real(r8) :: a(np,np),b(np,mp)

  integer, parameter :: NMAX=50
  integer :: i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
  real(r8) :: big,dum,pivinv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Algorithm                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  do j=1,n
     ipiv(j)=0
  enddo
  do i=1,n
     big=0.
     do j=1,n
        if(ipiv(j).ne.1)then
        do k=1,n
           if (ipiv(k).eq.0) then
              if (abs(a(j,k)).ge.big)then
                 big=abs(a(j,k))
                 irow=j
                 icol=k
              endif
           else if (ipiv(k).gt.1) then
              stop 'singular matrix in gaussj'
           endif
        enddo
        endif
     enddo
     ipiv(icol)=ipiv(icol)+1
     if (irow.ne.icol) then
        do l=1,n
           dum=a(irow,l)
           a(irow,l)=a(icol,l)
           a(icol,l)=dum
        enddo
        do l=1,m
           dum=b(irow,l)
           b(irow,l)=b(icol,l)
           b(icol,l)=dum
        enddo
     endif
     indxr(i)=irow
     indxc(i)=icol
     if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
     pivinv=1./a(icol,icol)
     a(icol,icol)=1.
     do l=1,n
        a(icol,l)=a(icol,l)*pivinv
     enddo
     do l=1,m
        b(icol,l)=b(icol,l)*pivinv
     enddo
     do ll=1,n
        if(ll.ne.icol)then
        dum=a(ll,icol)
        a(ll,icol)=0.
        do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
        enddo
        do l=1,m
           b(ll,l)=b(ll,l)-b(icol,l)*dum
        enddo
        endif
     enddo

  enddo


  do l=n,1,-1
     if(indxr(l).ne.indxc(l))then
       do k=1,n
          dum=a(k,indxr(l))
          a(k,indxr(l))=a(k,indxc(l))
          a(k,indxc(l))=dum
       enddo
     endif
  enddo


  return

  end subroutine gaussj 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
! computes the coefficients for the high-orde finite differences              !
! approximation to the kinetic energy operator. (3D-Laplacian)                !
!                                                                             !
!                      Norder                                                 !
!                       ---                                                   !
!        nabla^2 f(i) = \     cN(k)* { f(i+k) + f(i-k) }  + cN(0)*f(i)        !
!                       /__                                                   !
!                      k = 1                                                  !
!                                                                             !
! A general algorithm is implemented based on B. Fornberg and D.M. Sloan.     !
! The first six-order are tabulated as in ref. J.R. Chelikowsky, N. Troullier,!
! K. Wu and Y. Saad, PRB50, 11355 (1994); and references therein).            !
!                                                                             !
! A. Rubio (September 1999)                                                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine high_derivative(norder,cn)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Arguments and local variables                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit none

  integer, intent(in) :: norder           ! high order discretization
  real(r8), intent(out) :: cn(0:norder)

  integer :: i, Morder, N_laplac
  real(r8), allocatable :: cc(:,:,:) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Up to order six, the values are given explicitly                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! General calculation for centered approximation (the real orde of the        !
! approximation is 2*Norder, as Norder refers to the number of positive       !
! coefficients, in fact the total number of coefficients is 2*Norder +1 in    !
! the finite-difference approach to the Laplacian).                           !
! NOTE: the weight subroutine is general for all differential operators!      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  else
                          
     Morder = Norder * 2
     N_laplac = 2
     allocate ( cc(0:Morder,0:Morder,0:N_laplac) )

     call weights(N_laplac,Morder,cc)

     do i=0,Norder
        cN(i) = cc(i*2,Morder,N_laplac)           !for the laplacian
     enddo

     deallocate(cc)

  endif

  end subroutine high_derivative

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the Weights for finite-difference calucaltions:                     ! 
!                                                                             !
!  N -> highest order fo the derivative to be approximated                    !
!  M -> number of grid points to be used in the approsimation.                !
!                                                                             !
!  c(j,k,i) -> ith order derivative at kth-order approximation                !
!              j=0,k: the coefficients acting of each point                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine weights (N,M,cc)
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Arguments and local variables                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit none

  integer, intent(in) :: N, M
  real(r8), intent(out) :: cc(0:M,0:M,0:N)

  integer :: i,j,k, mn
  real(r8) :: c1, c2, c3, c4, c5, xi
  real(r8) :: x(0:M)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Algorithm                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!
!   grid-points for one-side finite-difference formulas on an equi.spaced grid
!  x(:) = (/(i,i=0,M)/) 

!   grid-points for centered finite-difference formulas on an equi.spaced grid
  mn = M/2
  x(:) = (/0,(-i,i,i=1,mn)/)

  xi = 0.0_r8               !point at which the approx. are to be accurate

  cc = 0.0_r8
  cc(0,0,0) = 1.0_r8

  c1 = 1.0_r8
  c4 = x(0) - xi
       
  do j = 1 , M

     mn = min(j,N)
     c2 = 1.0_r8
     c5 = c4
     c4 = x(j) - xi
 
     do k = 0, j-1
        c3 = x(j) - x(k)
        c2 = c2*c3
            
        if (j <= N) cc(k,j-1,j) = 0.0_r8
        cc(k,j,0) = c4*cc(k,j-1,0)/c3

        do i = 1, mn
              cc(k,j,i) = ( c4*cc(k,j-1,i) -i*cc(k,j-1,i-1) )/c3
        enddo

        cc(j,j,0) = -c1*c5*cc(j-1,j-1,0) / c2         
     enddo

     do i = 1, mn
        cc(j,j,i) = c1*(i*cc(j-1,j-1,i-1)-c5*cc(j-1,j-1,i))/c2
     enddo

     c1 = c2

  enddo

  end subroutine weights 

end module math


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

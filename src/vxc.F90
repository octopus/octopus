#include "global.h"

module vxc

  use global

  implicit none

  contains
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates the exchange energy and exchange potential for a homogeneous
! electron gas (LDA for exchange).
! It works in both the 3D and 2D cases
! (1D not yet implemented, although one) only has to find out the a_x constant,
! but I do not know what is the value of \int_0^\infty sin**2(x)/x**3 )
!
! Input:
! integer nsp   :: 1 for spin-unpolarized, 2 for spin-polarized.
! integer irel  :: 0 for non-relativistic exchange; relativistic otherwise.
! FLOAT ds(nsp) :: density.
! Output:
! vx(nsp)       :: exchange potential.
! ex            :: exchange energy density.
!
! The basic formulae are (Hartree atomic units are assumed):
!
!    p = ((dim+1)/dim)
!    ex(n) = a_x*n**(1/dim)
!    ex(n,z) = ex(n)*f(z)
!    vx_up(n, z) = ex(n)*( p*f(z) + (df/dz)(z)*(1-z) )
!    vx_do(n, z) = ex(n)*( p*f(z) - (df/dz)(z)*(1+z) )
!    f(z) = (1/2)*( (1+z)**p + (1-z)**p)
!    a_x = -(3/(4*pi))*(3*pi**2)**(1/3) in 3D
!    a_x = -(4/3)*sqrt(2/pi) in 2D
!    a_x = -(1/2) * \int_0^\infty (sin(x))**2/x**3
!
! If irel is not zero, relativistic correction factors have to be applied.
! These are however not implemented in 2D, so nothing is done in that case.
!
! WARNING: Check that the relativistic corrections are OK for the potential
!          in the spin polarized case.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine exchange(nsp, ds, ex, vx, irel)
    integer, intent(in) :: nsp, irel
    FLOAT, intent(in)   :: ds(nsp)
    FLOAT, intent(out)  :: ex, vx(nsp)

    FLOAT, parameter :: c014 = CNST(0.014),                     &
                        opf = CNST(1.5),                        &
                        MINDEN = CNST(1e-15),                   &
                        a_x(3) = (/ -M_ONE,                     &
                                    CNST(-1.06384608107049),    &
                                    CNST(-0.738558766382022) /)
    FLOAT :: d, d1, d2, z, fz, dfdz, beta, alb, sb, rs, p

    p = (conf%dim+M_ONE)/conf%dim

    select case(nsp)
    case(1) ! Spin - unpolarized
       d = max(M_ZERO, ds(1))
       if(d < MINDEN) then
         ex = M_ZERO; vx(1) = M_ZERO
         return
       endif
       ex    = a_x(conf%dim)*d**(M_ONE/conf%dim)
       vx(1) = p*ex

    case(2) ! Spin-polarized
       d1 = max(M_ZERO, ds(1))
       d2 = max(M_ZERO, ds(2))
       d = d1 + d2
       if(d < MINDEN) then
          ex = M_ZERO; vx(1:2) = M_ZERO
          return
       endif
       z = (d1-d2)/(d1+d2)
       fz = M_HALF * ( (1+z)**p+(1-z)**p )
       dfdz = M_HALF * p * ( (1+z)**(p-1) - (1-z)**(p-1) )
       ex = a_x(conf%dim)*d**(M_ONE/conf%dim) ! Unpolarized yet, to be used in the following formulae.
       vx(1) = ex * ( p*fz + dfdz*(1-z) )
       vx(2) = ex * ( p*fz - dfdz*(1+z) )
       ex = ex*fz ! Now it is the correct polarized result.
    case default
       message(1) = 'Wrong "nsp" argument passed to exchange3d'
       call write_fatal(1)
    end select

    ! Relativistic corrections (I believe that they are rather useless).
    if( irel .ne. 0 .and. conf%dim == 3 ) then
      rs = (M_THREE / (M_FOUR*M_PI*d) )**M_THIRD
      beta = c014/rs
      sb = sqrt(1+beta*beta)
      alb = log(beta+sb)
      ex = ex*(M_ONE-opf*((beta*sb-alb)/beta**2)**2)
      vx = vx*(-M_HALF + opf * alb / (beta*sb))
    endif

  end subroutine exchange

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! /* Calculates the LDA corrected exchange suggested in PRA 49, 2421 (1994) */
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine lb94(nspin, dens, gdens, dexdd, r, ip, qtot, modified, beta, threshold)
    integer, intent(in)  :: nspin
    FLOAT, intent(in) :: dens(nspin), gdens(3, nspin)
    FLOAT, intent(out):: dexdd(nspin)!, dexdgd(3, nspin)
    FLOAT, intent(in) :: r, ip, qtot, beta, threshold
    logical, intent(in) :: modified

    integer  :: is
    FLOAT :: alpha, gdm, x, f, gamma

    ! First, get the LDA exchange potential.
    call exchange(nspin, dens, x, dexdd, 0)

    ! If ip is zero, the alpha is set to 0.5. In this way, the boundary 
    ! condition is -1/r. Equivalently, if ip = 1/32, alpha is also 0.5
    if(ip > M_ZERO) then 
      alpha = M_TWO*sqrt(M_TWO*ip)
      gamma = qtot**M_THIRD/(M_TWO*alpha)
    else
      alpha = M_HALF
      gamma = qtot**M_THIRD
    endif

    if(.not.modified) then
      gamma = M_ONE
    endif

    do is = 1, nspin
      gdm   = sqrt( gdens(1,is)**2 + gdens(2,is)**2 + gdens(3,is)**2 )
      !gdm = alpha*dens(is)
      if(dens(is) >= threshold .and. gdm >=threshold) then
        x = gdm / dens(is)**(M_FOUR/M_THREE)
        f = -beta*dens(is)**M_THIRD*&
           x**2/(M_ONE + M_THREE*beta*x*loct_asinh(gamma*x))
        dexdd(is) = dexdd(is) + f !- beta * dens(is)**M_THIRD * f
      elseif(r > M_ZERO .and. dens(is)<= threshold) then
        f = r + (M_THREE/alpha)*log(2*gamma*alpha*qtot**(-M_THIRD))
        !f = f + (qtot*exp(-alpha*r))**M_THIRD/(beta*alpha**2)
        dexdd(is) = dexdd(is) - M_ONE/f
      endif
    enddo

  end subroutine lb94

end module vxc

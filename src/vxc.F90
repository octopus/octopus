#include "global.h"

module vxc

  use global

  implicit none

  private :: change

  contains
   
  subroutine change(nsp, ds, d, z)
    integer, intent(in) :: nsp
    FLOAT, intent(in)   :: ds(nsp)
    FLOAT, intent(out)  :: d, z
    FLOAT :: d1, d2

    ASSERT(nsp==1.or.nsp==2)

    z = M_ZERO
    select case(nsp)
    case(1) ! Spin - unpolarized
      d = max(M_ZERO, ds(1))
    case(2) ! Spin-polarized
      d1 = max(M_ZERO, ds(1))
      d2 = max(M_ZERO, ds(2))
      d = d1 + d2
      if(d > M_ZERO) z = (d1-d2)/d
    end select

  end subroutine change

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
! Finds the exchange and correlation energies at a point, and their           !
! derivatives with respect to density and density gradient, in the            !
! Generalized Gradient Correction approximation.                              !
! Lengths in Bohr, energies in Hartrees                                       !
! Written by L.C.Balbas and J.M.Soler, Dec 96. Version 0.5.                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ggaxc( AUTHOR, IREL, NSPIN, D, GD,                               &
                    EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
  use global


  implicit none
  character(len=*) :: AUTHOR
  integer :: IREL, NSPIN
  FLOAT :: D(NSPIN), DECDD(NSPIN), DECDGD(3,NSPIN),  &
      DEXDD(NSPIN), DEXDGD(3,NSPIN),                 &
      EPSC, EPSX, GD(3,NSPIN)

  IF (NSPIN .GT. 2) STOP 'GGAXC: Not prepared for this NSPIN'

  IF (AUTHOR.EQ.'PBE' .OR. AUTHOR.EQ.'pbe') THEN
    CALL PBEX( IREL, NSPIN, D, GD,   &
        EPSX, DEXDD, DEXDGD)
    CALL PBEC( NSPIN, D, GD,         &
        EPSC, DECDD, DECDGD )
  ELSE
     WRITE(6,*) 'GGAXC: Unknown author ', AUTHOR
     STOP
  ENDIF

  end subroutine ggaxc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation.       !
! Ref: J.P.Perdew, K.Burke & M.Ernzerhof, PRL 77, 3865 (1996)                 !
! Written by L.C.Balbas and J.M.Soler. December 1996. Version 0.5.            !
! NOTE: Split in two: pbec and pbex.                                          !
! ******** INPUT ******************************************************       !
! INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)         !
! INTEGER NSPIN          : Number of spin polarizations (1 or 2)              !
! REAL*8  DENS(NSPIN)    : Total electron density (if NSPIN=1) or             !
!                           spin electron density (if NSPIN=2)                !
! REAL*8  GDENS(3,NSPIN) : Total or spin density gradient                     !
! ******** OUTPUT *****************************************************       !
! REAL*8  EX             : Exchange energy density                            !
! REAL*8  EC             : Correlation energy density                         !
! REAL*8  DEXDD(NSPIN)   : Partial derivative                                 !
!                           d(DensTot*Ex)/dDens(nspin),                       !
!                           where DensTot = Sum_nspin( DENS(nspin) )          !
!                           For a constant density, this is the               !
!                          exchange potential                                 !
! REAL*8  DECDD(NSPIN)   : Partial derivative                                 !
!                           d(DensTot*Ec)/dDens(nspin),                       !
!                           where DensTot = Sum_nspin( DENS(nspin) )          !
!                          For a constant density, this is the                !
!                          correlation potential                              !
! REAL*8  DEXDGD(3,NSPIN): Partial derivative                                 !
!                           d(DensTot*Ex)/d(GradDens(i,nspin))                !
! REAL*8  DECDGD(3,NSPIN): Partial derivative                                 !
!                           d(DensTot*Ec)/d(GradDens(i,nspin))                !
! ********* UNITS ****************************************************        !
! Lengths in Bohr                                                             !
! Densities in electrons per Bohr**3                                          !
! Energies in Hartrees                                                        !
! Gradient vectors in cartesian coordinates                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  subroutine pbex( IREL, NSPIN, DENS, GDENS,                                  &
                    EX, DEXDD, DEXDGD)
  use global

  implicit none 
  integer :: IREL, NSPIN
  FLOAT :: DENS(NSPIN),                                                    &
              DEXDD(NSPIN), DEXDGD(3,NSPIN), GDENS(3,NSPIN)

! Internal variables
  integer :: IS, IX 

  FLOAT :: BETA, D(2), DF1DD, DF1DGD, DFDD, DFDGD, DFXDD(2),       &
       DFXDGD(3,2), DKFDD, DS(2), DSDD, DSDGD, DT, EX, EXUNIF,     &
       F, F1, FX, GAMMA, GD(3,2), GDM(2), GDMS, GDMT, GDS, GDT(3), &
       KAPPA, KFS, MU, S, VXUNIF(2)

! Lower bounds of density and its gradient to avoid divisions by zero
  FLOAT, parameter :: DENMIN = CNST(1.e-12), GDMIN = CNST(1.e-12)

! Fix some more numerical constants
  BETA  = CNST(0.066725)
  GAMMA = (M_ONE - log(M_TWO)) / M_PI**2
  MU    = BETA * M_PI**2 / M_THREE
  KAPPA = CNST(0.8040)

! Translate density and its gradient to new variables
  IF (NSPIN .EQ. 1) THEN
     D(1) = M_HALF * DENS(1)
     D(2) = D(1)
     DT = MAX( DENMIN, DENS(1) )
     DO 10 IX = 1, 3
        GD(IX,1) = M_HALF * GDENS(IX,1)
        GD(IX,2) = GD(IX,1)
        GDT(IX) = GDENS(IX,1)
     10 CONTINUE
  ELSE
     D(1) = DENS(1)
     D(2) = DENS(2)
     DT = MAX( DENMIN, DENS(1)+DENS(2) )
     do IX = 1, 3
       GD(IX,1) = GDENS(IX,1)
       GD(IX,2) = GDENS(IX,2)
       GDT(IX) = GDENS(IX,1) + GDENS(IX,2)
     end do
  ENDIF
  GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
  GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
  GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
  GDMT   = MAX( GDMIN, GDMT )

! Find exchange energy and potential
  FX = 0
  do IS = 1,2
    DS(IS)   = MAX( DENMIN, 2 * D(IS) )
    GDMS = MAX( GDMIN, 2 * GDM(IS) )
    KFS = (3 * M_PI**2 * DS(IS))**M_THIRD
    S = GDMS / (2 * KFS * DS(IS))
    F1 = 1 + MU * S**2 / KAPPA
    F = 1 + KAPPA - KAPPA / F1
    
    !Note nspin=1 in call to exchng...
    call exchange( 1, ds(is), EXUNIF, VXUNIF(IS), irel )
    FX = FX + DS(IS) * EXUNIF * F
    
    DKFDD = M_THIRD * KFS / DS(IS)
    DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
    DF1DD = 2 * (F1-1) * DSDD / S
    DFDD = KAPPA * DF1DD / F1**2
    DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD
    
    do IX = 1, 3
      GDS = M_TWO * GD(IX,IS)
      DSDGD = (S / GDMS) * GDS / GDMS
      DF1DGD = M_TWO * MU * S * DSDGD / KAPPA
      DFDGD = KAPPA * DF1DGD / F1**2
      DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
    end do
  end do
  FX = M_HALF * FX / DT

! Set output arguments
  EX = FX
  do IS = 1,NSPIN
    DEXDD(IS) = DFXDD(IS)
    do IX = 1,3
      DEXDGD(IX,IS) = DFXDGD(IX,IS)
    end do
  end do
  
  end subroutine pbex

  subroutine pbec(  NSPIN, DENS, GDENS,                                 &
                    EC, DECDD, DECDGD )
  implicit none 
  integer ::  NSPIN
  FLOAT :: DENS(NSPIN), DECDD(NSPIN), DECDGD(3,NSPIN),                     &
              GDENS(3,NSPIN)

! Internal variables
  integer :: IS, IX 

  FLOAT :: A, BETA, D(2), DADD, DECUDD, DF1DD, DF2DD, DF3DD,  &
       DF4DD, DF3DGD, DF4DGD, DFCDD(2), DFCDGD(3,2), DHDD,    &
       DHDGD, DKFDD, DKSDD, DPDD, DPDZ, DRSDD, DT, DTDD,      &
       DTDGD, DZDD(2), EC, ECUNIF, F1, F2, F3, F4, FC, GAMMA, & 
       GD(3,2), GDM(2), GDMT, GDT(3), H, KAPPA, KF, KS, MU,   &
       PHI, RS, T, VCUNIF(2), ZETA

! Lower bounds of density and its gradient to avoid divisions by zero

! Lower bounds of density and its gradient to avoid divisions by zero
  FLOAT, parameter :: DENMIN = CNST(1.e-12), GDMIN = CNST(1.e-12)

! Fix some more numerical constants
  BETA  = CNST(0.066725)
  GAMMA = (M_ONE - log(M_TWO)) / M_PI**2
  MU    = BETA * M_PI**2 / M_THREE
  KAPPA = CNST(0.8040)

! Translate density and its gradient to new variables
  IF (NSPIN .EQ. 1) THEN
     D(1) = M_HALF * DENS(1)
     D(2) = D(1)
     DT = MAX( DENMIN, DENS(1) )
     DO 10 IX = 1, 3
        GD(IX,1) = M_HALF * GDENS(IX,1)
        GD(IX,2) = GD(IX,1)
        GDT(IX) = GDENS(IX,1)
     10 CONTINUE
  ELSE
     D(1) = DENS(1)
     D(2) = DENS(2)
     DT = MAX( DENMIN, DENS(1)+DENS(2) )
     do IX = 1, 3
       GD(IX,1) = GDENS(IX,1)
       GD(IX,2) = GDENS(IX,2)
       GDT(IX) = GDENS(IX,1) + GDENS(IX,2)
     end do
  ENDIF
  GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
  GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
  GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
  GDMT   = MAX( GDMIN, GDMT )

! Find local correlation energy and potential
  CALL correlation_lda_3D_pw92( 2, D, ECUNIF, VCUNIF )

! Find total correlation energy
  RS = (M_THREE / (4*M_PI*DT) )**M_THIRD
  KF = (M_THREE * M_PI**2 * DT)**M_THIRD
  KS = sqrt( M_FOUR * KF / M_PI )
  ZETA = ( D(1) - D(2) ) / DT
  ZETA = MAX( -M_ONE + DENMIN, ZETA )
  ZETA = MIN(  M_ONE - DENMIN, ZETA )
  PHI = M_HALF * ( (M_ONE + ZETA)**M_TWOTHIRD + (M_ONE - ZETA)**M_TWOTHIRD )
  T = GDMT / (M_TWO * PHI * KS * DT)
  F1 = ECUNIF / GAMMA / PHI**3
  F2 = EXP(-F1)
  A = BETA / GAMMA / (F2-1)
  F3 = T**2 + A * T**4
  F4 = BETA/GAMMA * F3 / (M_ONE + A*F3)
  H = GAMMA * PHI**3 * LOG( M_ONE + F4 )
  FC = ECUNIF + H


! Find correlation energy derivatives
  DRSDD = - (M_THIRD * RS / DT)
  DKFDD =    M_THIRD * KF / DT
  DKSDD = M_HALF * KS * DKFDD / KF
  DZDD(1) =   1 / DT - ZETA / DT
  DZDD(2) = - (1 / DT) - ZETA / DT
  DPDZ = M_HALF * M_TWOTHIRD * ( 1/(1+ZETA)**M_THIRD - 1/(1-ZETA)**M_THIRD )
  do IS = 1,2
    DECUDD = ( VCUNIF(IS) - ECUNIF ) / DT
    DPDD = DPDZ * DZDD(IS)
    DTDD = (- T) * ( DPDD/PHI + DKSDD/KS + 1/DT )
    DF1DD = F1 * ( DECUDD/ECUNIF - 3*DPDD/PHI )
    DF2DD = (- F2) * DF1DD
    DADD = (- A) * DF2DD / (F2-1)
    DF3DD = (2*T + 4*A*T**3) * DTDD + DADD * T**4
    DF4DD = F4 * ( DF3DD/F3 - (DADD*F3+A*DF3DD)/(1+A*F3) )
    DHDD = 3 * H * DPDD / PHI
    DHDD = DHDD + GAMMA * PHI**3 * DF4DD / (1+F4)
    DFCDD(IS) = VCUNIF(IS) + H + DT * DHDD

    do IX = 1,3
      DTDGD = (T / GDMT) * GDT(IX) / GDMT
      DF3DGD = DTDGD * ( 2 * T + 4 * A * T**3 )
      DF4DGD = F4 * DF3DGD * ( 1/F3 - A/(1+A*F3) ) 
      DHDGD = GAMMA * PHI**3 * DF4DGD / (1+F4)
      DFCDGD(IX,IS) = DT * DHDGD
    end do
  end do

! Set output arguments
  EC = FC
  do IS = 1,NSPIN
    DECDD(IS) = DFCDD(IS)
    do IX = 1,3
      DECDGD(IX,IS) = DFCDGD(IX,IS)
    end do
  end do
  
  end subroutine pbec

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates the exchange energy and correlation potential for a homogeneous
! electron gas (LDA for exchange).
! Ref: J.P.Perdew & Y.Wang, Phys. Rev. B 45, 13244 (1992).
!
! Input:
! integer nsp   :: 1 for spin-unpolarized, 2 for spin-polarized.
! FLOAT ds(nsp) :: density.
! Output:
! vc(nsp)       :: correlation potential.
! ec            :: correlation energy density.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine correlation_lda_3D_pw92(nsp, ds, ec, vc)
    implicit none
    integer, intent(in) :: nsp
    FLOAT, intent(in)   :: ds(nsp)
    FLOAT, intent(out)  :: ec, vc(nsp)

    FLOAT :: d, z, ec0, ec1, alphac, rs, fz, fpz, &
             dec0drs, dec1drs, dalphacdrs, decdrs, decdz
    FLOAT, parameter :: &
       a(1:3)       = (/ CNST(0.031091), CNST(0.015545), CNST(0.016887) /), &
       alpha(1:3)   = (/ CNST(0.21370),  CNST(0.20548),  CNST(0.11125) /), &
       beta(1:4, 1:3) = reshape( (/ CNST(7.5957),   CNST(3.5876),   CNST(1.6382),  CNST(0.49294) ,  &
                                    CNST(14.1189),  CNST(6.1977),   CNST(3.3662),  CNST(0.62517) ,  &
                                    CNST(10.357),   CNST(3.6231),   CNST(0.88026), CNST(0.49671) /),&
                                 (/4, 3/) )
    FLOAT, parameter :: MINDEN = CNST(1e-12), &
                        fzconst = CNST(1.92366105093154), & ! = 1/(2**(4/3)-2)
                        fz20 = CNST(1.709921) ! (d^2f/dz^2)(z=0)

    select case(nsp)
    case(1) ! Spin - unpolarized
       d = max(MINDEN, ds(1))
       z   = M_ZERO
       fz  = M_ZERO 
       fpz = M_ZERO
    case(2) ! Spin-polarized
       d  = max(MINDEN, ds(1)+ds(2))
       z = (ds(1)-ds(2))/d
       fz  = fzconst*((M_ONE+z)**(M_FOUR/M_THREE) + (M_ONE-z)**(M_FOUR/M_THREE) - M_TWO)
       fpz = fzconst*(M_FOUR/M_THREE)*((M_ONE+z)**(M_ONE/M_THREE) - (M_ONE-z)**(M_ONE/M_THREE))
    case default
       message(1) = 'Wrong "nsp" argument passed to exchange3d'
       call write_fatal(1)
    end select

    rs = (M_THREE / (M_FOUR*M_PI*d) )**M_THIRD

    ec0        =  g(1)
    ec1        =  g(2)
    alphac     = -g(3)
    dec0drs    =  dgdrs(1)
    dec1drs    =  dgdrs(2)
    dalphacdrs = -dgdrs(3)

    ec = ec0 + alphac*fz*(M_ONE-z**4)/fz20 + (ec1-ec0)*fz*z**4
    decdrs = dec0drs*(M_ONE - fz*z**4) + dec1drs*fz*z**4 + dalphacdrs*fz*(M_ONE-z**4)/fz20
    decdz  = M_FOUR*z**3*fz*(ec1-ec0-alphac/fz20) + &
             fpz*(z**4*ec1 - z**4*ec0 + (M_ONE-z**4)*alphac/fz20)
    vc(1) = ec - (rs/M_THREE)*decdrs - (z-M_ONE)*decdz
    if(nsp==2) then
       vc(2) = ec - (rs/M_THREE)*decdrs - (z+M_ONE)*decdz
    endif

    contains

    ! Function g defined Eq. 10 of the original paper.
    FLOAT function g(k)
      integer, intent(in) :: k
      g = - M_TWO * a(k) * ( M_ONE + alpha(k) * rs) * log( M_ONE + &
          M_ONE/(M_TWO*a(k)*(beta(1, k)*sqrt(rs) + beta(2, k)*rs + beta(3, k)*sqrt(rs)**3 + beta(4, k)*rs**2)) )
    end function g
    ! The derivative of previous function with respect to rs (Eq. A5)
    FLOAT function dgdrs(k)
      integer, intent(in) :: k
      FLOAT :: q0, q1, q1p, b
      b = (beta(1, k)*sqrt(rs) + beta(2, k)*rs + beta(3, k)*sqrt(rs)**3 + beta(4, k)*rs**2)
      q0 = -M_TWO * a(k) * ( M_ONE + alpha(k) * rs)
      q1 = M_TWO * a(k) * b
      q1p = a(k) * (beta(1, k)/sqrt(rs) + M_TWO*beta(2, k) + M_THREE*beta(3, k)*sqrt(rs) + M_FOUR*beta(4, k)*rs)
      dgdrs = -M_TWO * a(k) * alpha(k) * log(M_ONE + M_ONE/q1) - (q0*q1p)/(q1**2+q1)
    end function dgdrs
  end subroutine correlation_lda_3D_pw92

end module vxc

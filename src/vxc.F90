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
    select case(nsp)
    case(1) ! Spin - unpolarized
       d = max(M_ZERO, ds(1))
       z = M_ZERO
    case(2) ! Spin-polarized
       d1 = max(M_ZERO, ds(1))
       d2 = max(M_ZERO, ds(2))
       d = d1 + d2
       z = (d1-d2)/(d1+d2)
    case default
       message(1) = 'Wrong "nsp" argument passed to change'
       call write_fatal(1)
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
! Finds the exchange and correlation energies and potentials, in the          !
! Local (spin) Density Approximation.                                         !
! Written by L.C.Balbas and J.M.Soler, Dec 96.                                !
! Non-collinear spin added by J.M.Soler, May 98                               !
! *********** INPUT ************************************************          !
! CHARACTER*(*) AUTHOR : Parametrization desired:                             !
!     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)                !
!           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992)                !
!                     Uppercase is optional                                   !
! INTEGER IREL     : Relativistic exchange? (0=>no, 1=>yes)                   !
! INTEGER NSPIN    : NSPIN=1 => unpolarized; NSPIN=2 => polarized;            !
!                    NSPIN=4 => non-collinear polarization                    !
! REAL*8  D(NSPIN) : Local (spin) density. For non-collinear                  !
!                    polarization, the density matrix is given by:            !
!                    D(1)=D11, D(2)=D22, D(3)=Real(D12), D(4)=Im(D12)         !
! *********** OUTPUT ***********************************************          !
! REAL*8 EPSX, EPSC : Exchange and correlation energy densities               !
! REAL*8 VX(NSPIN), VC(NSPIN) : Exchange and correlation potentials,          !
!                               defined as dExc/dD(nspin)                     !
! *********** UNITS ************************************************          !
! Lengths in Bohr, energies in Hartrees                                       !
! ******************************************************************          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ldaxc( AUTHOR, IREL, NSPIN, D, EPSX, EPSC, VX, VC )

  use global

  IMPLICIT NONE
  CHARACTER(len=*) :: AUTHOR
  INTEGER :: IREL, NSPIN
  FLOAT :: D(NSPIN), EPSC, EPSX, VX(NSPIN), VC(NSPIN)

  INTEGER :: IS, NS
  FLOAT :: DD(2), DPOL, DTOT, VCD(2), VPOL, VXD(2)

  FLOAT, parameter :: TINY = CNST(1.e-12)

  IF (NSPIN .EQ. 4) THEN
!       Find eigenvalues of density matrix (up and down densities
!       along the spin direction)
!       Note: D(1)=D11, D(2)=D22, D(3)=Real(D12), D(4)=Im(D12)
     NS = 2
     DTOT = D(1) + D(2)
     DPOL = SQRT( (D(1)-D(2))**2 + M_FOUR*(D(3)**2+D(4)**2) )
     DD(1) = M_HALF * ( DTOT + DPOL )
     DD(2) = M_HALF * ( DTOT - DPOL )
  ELSE
     NS = NSPIN
     DO 10 IS = 1,NSPIN
        DD(IS) = D(IS)
     10 CONTINUE
  ENDIF

  IF ( AUTHOR.EQ.'CA' .OR. AUTHOR.EQ.'ca' .OR.                                &
       AUTHOR.EQ.'PZ' .OR. AUTHOR.EQ.'pz') THEN
     call exchange(ns, dd, epsx, vxd, irel)
     CALL correlation_lda_3D_pz   (      NS, DD, EPSC, VCD)
  ELSEIF ( AUTHOR.EQ.'PW92' .OR. AUTHOR.EQ.'pw92' ) THEN
     call exchange(ns, dd, epsx, vxd, irel)
     CALL correlation_lda_3D_pw92 ( NS, DD, EPSC, VCD )
  ELSE
     WRITE(6,*) 'LDAXC: Unknown author ', AUTHOR
     STOP
  ENDIF

  IF (NSPIN .EQ. 4) THEN
!       Find dE/dD(nspin) = dE/dDup * dDup/dD(nspin) +
!                           dE/dDdown * dDown/dD(nspin)
     VPOL  = (VXD(1)-VXD(2)) * (D(1)-D(2)) / (DPOL+TINY)
     VX(1) = M_HALF * ( VXD(1) + VXD(2) + VPOL )
     VX(2) = M_HALF * ( VXD(1) + VXD(2) - VPOL )
     VX(3) = (VXD(1)-VXD(2)) * D(3) / (DPOL+TINY)
     VX(4) = (VXD(1)-VXD(2)) * D(4) / (DPOL+TINY)
     VPOL  = (VCD(1)-VCD(2)) * (D(1)-D(2)) / (DPOL+TINY)
     VC(1) = M_HALF * ( VCD(1) + VCD(2) + VPOL )
     VC(2) = M_HALF * ( VCD(1) + VCD(2) - VPOL )
     VC(3) = (VCD(1)-VCD(2)) * D(3) / (DPOL+TINY)
     VC(4) = (VCD(1)-VCD(2)) * D(4) / (DPOL+TINY)
  ELSE
     DO 20 IS = 1,NSPIN
        VX(IS) = VXD(IS)
        VC(IS) = VCD(IS)
     20 CONTINUE
  ENDIF

  end subroutine ldaxc

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
! Correlation energy per particle and potentials for a homogeneous electron
! gas in 2D, as parametrized by Attacalite et al.
! Refs: [1] C. Attacalite et al, Phys. Rev. Lett. 88, 256601 (2002) for 2D.
!       [2] C. Attacalite, PhD thesis. 
  subroutine correlation_lda_2D_atta(nsp, ds, ec, vc)
    implicit none
    integer, intent(in) :: nsp
    FLOAT, intent(in)   :: ds(nsp)
    FLOAT, intent(out)  :: ec, vc(nsp)

    FLOAT, parameter :: a(0:2) = (/ CNST(-0.1925),   CNST(0.117331),    CNST(0.0234188) /), &
                        b(0:2) = (/ CNST(0.0863136), CNST(-3.394e-2),   CNST(-0.037093) /), &
                        c(0:2) = (/ CNST(0.0572384), CNST(-7.66765e-3), CNST(0.0163618) /), &
                        e(0:2) = (/ CNST(1.0022),    CNST(0.4133),      CNST(1.424301) /), &
                        f(0:2) = (/ CNST(-0.02069),  CNST(0.0),         CNST(0.0) /), &
                        g(0:2) = (/ CNST(0.33997),   CNST(6.68467e-2),  CNST(0.0) /), &
                        h(0:2) = (/ CNST(1.747e-2),  CNST(7.799e-4),    CNST(1.163099) /), &
                        d(0:2) = -a(0:2)*h(0:2)
    FLOAT, parameter :: MINDEN = CNST(1e-15), &
                        beta = CNST(1.3386)
    FLOAT :: rs, dens, d1, d2, zeta, ex, ex0, vx(2), dd(2), ex6, calf, calfp, ax, decdrs, decdz

    ax = -M_FOUR/(M_THREE*M_PI*sqrt(M_TWO))

    ! Get the trace and the polarization of the density
    call change(nsp, ds, dens, zeta)
    ! Wigner radius
    rs = sqrt(M_ONE / (M_PI*dens))

    ! If the density is too small, return zero  
    if(dens < MINDEN) then
       ec = M_ZERO; vc = M_ZERO
       return
    endif

    ! In unpolarized cases, expressions are fairly simple.
    if(nsp==1) then
      ec = alpha(0)
      vc(1) = dalphadrs(0)
      return
    endif

    ! In the case of spin-polarized calculations, all this is necessary...
    ! Unpolarized exchange energy...
    ex0 = ((-M_FOUR*sqrt(M_TWO))/(M_THREE*M_PI*rs))
    ! Polarized exchange energy
    ex  = M_HALF*((M_ONE+zeta)**(M_HALF*M_THREE)+(M_ONE-zeta)**(M_HALF*M_THREE))*ex0
    ! Taylor expansion of ex, in zeta, beyond fourth order.
    ex6 = ex - (M_ONE  + (M_THREE/CNST(8.0))*zeta**2 + (M_THREE/CNST(128.0))*zeta**4)*ex0
    ! Correlation energy (Eq.[1]3)
    ec = alpha(0) + alpha(1)*zeta**2 + alpha(2)*zeta**4 + (exp(-beta*rs)-M_ONE)*ex6
    ! Function calf (Eq.[2]4.10)
    calf = (M_ONE+zeta)**(M_THREE/M_TWO)+(M_ONE-zeta)**(M_THREE/M_TWO) - &
           (M_TWO + (M_THREE/M_FOUR)*zeta**2 + (M_THREE/CNST(64.0))*zeta**4)
    ! Function calfp (Eq.[2]C8)
    calfp = (M_THREE/M_TWO)*(sqrt(M_ONE+zeta)-sqrt(M_ONE-zeta)) &
            - (M_THREE/M_TWO)*zeta - (M_THREE/CNST(16.0))*zeta**3
    ! Derivative of the correlation energy with respect to rs (Eq.[2]C2)
    decdrs = ax*calf*(M_ONE-exp(-beta*rs)*(M_ONE+beta*rs))/rs**2 + &
             dalphadrs(0) + dalphadrs(1)*zeta**2 + dalphadrs(2)*zeta**4
    ! Derivative of the correlation energy with respect to zeta (Eq.[2]C7)
    decdz  = ax*(exp(-beta*rs)-M_ONE)*calfp/rs + M_TWO*alpha(1)*zeta + M_FOUR*alpha(2)*rs**3
    ! And finally, the potentials (Eq.[2]C1)
    vc(1) = ec - M_HALF*rs*decdrs - (zeta-M_ONE)*decdz
    vc(1) = ec - M_HALF*rs*decdrs - (zeta+M_ONE)*decdz

    contains
    FLOAT function alpha(i) ! Eq.[1]4
      integer, intent(in) :: i
      alpha = a(i) + (b(i)*rs + c(i)*rs**2 + d(i)*rs**3) * &
              log(M_ONE + M_ONE/(e(i)*rs + f(i)*sqrt(rs)**3 + g(i)*rs**2 + h(i)*rs**3) )
    end function alpha
    FLOAT function dalphadrs(i) ! Eq.[2]C3
      integer, intent(in) :: i
      FLOAT :: efe, efep, lg, x
      efe  = e(i)*rs + f(i)*sqrt(rs)**3 + g(i)*rs**2 + h(i)*rs**3 ! Eq.[2]C5
      efep = e(i) + M_HALF*M_THREE*f(i)*sqrt(rs) + M_TWO*g(i)*rs + M_THREE*h(i)*rs**2 ! Eq. [2]C6
      lg = log(M_ONE + M_ONE/efe)
      x  = ((b(i)+c(i)*rs**2+d(i)*rs**3)*efep)/(efe**2+efe)
      dalphadrs = (b(i) + M_TWO*c(i)*rs + M_THREE*d(i)*rs**2)*lg - x ! Eq.[2]C3
    end function dalphadrs
  end subroutine correlation_lda_2D_atta

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Correlation energy per-particle and potential of a HEG as parameterized by
! Perdew & Zunger (Phys. Rev. B 23, 5048).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define eclow(i)  (gamma(i)/(M_ONE + beta1(i)*sqrs + beta2(i)*rs))
#define echigh(i) (a(i)*rslog + b(i) + c(i)*rs*rslog + d(i)*rs)
#define potlow(i) ((M_ONE + (CNST(7.0)/CNST(6.0))*beta1(i)*sqrs + (M_FOUR/M_THREE)*beta2(i)*rs) / \
                  (M_ONE + beta1(i)*sqrs + beta2(i)*rs))
#define pothigh(i) (a(i)*rslog + (b(i)-M_THIRD*a(i)) + M_TWO*c(i)*rs*rslog/M_THREE + \
                   M_THIRD*(M_TWO*d(i)-c(i))*rs)
  subroutine correlation_lda_3d_pz(nsp, ds, ec, vc)
    integer, intent(in)   :: nsp
    FLOAT, intent(in)  :: ds(nsp)
    FLOAT, intent(out) :: ec, vc(nsp)

    FLOAT, parameter :: TFTM = CNST(0.519842099789746380 ), &
                        FTRD = M_FOUR*M_THIRD
    FLOAT, parameter :: gamma(1:2) = (/ -CNST(0.1423), -CNST(0.0843) /), &
                        beta1(1:2) = (/  CNST(1.0529),  CNST(1.3981) /), &
                        beta2(1:2) = (/  CNST(0.3334),  CNST(0.2611) /), &
                        a(1:2)     = (/  CNST(0.0311),  CNST(0.015555) /), &
                        b(1:2)     = (/ -CNST(0.048),  -CNST(0.0269) /), &
                        c(1:2)     = (/  CNST(0.0020191519406228),  CNST(0.00069255121311694) /), &
                        d(1:2)     = (/ -CNST(0.0116320663789130), -CNST(0.00480126353790614) /)
    FLOAT :: dens, z, rs, sqrs, rslog, fz, fzp, ecp, vcp, ecu, vcu, x

    call change(nsp, ds, dens, z)
    if(dens.le.M_ZERO) then
       ec = M_ZERO
       vc = M_ZERO
       return
    endif

    rs = (M_THREE / (M_FOUR*M_PI*dens) )**M_THIRD
    sqrs = sqrt(rs)
    rslog = log(rs)

    if(rs>M_ONE) then
       ec = eclow(1)
       vc(1) = ec*potlow(1)
    else
       ec = echigh(1)
       vc(1) = pothigh(1)
    endif

    if(nsp==1) return

    fz = ((M_ONE+z)**FTRD+(M_ONE-z)**FTRD-M_TWO)/TFTM
    fzp = FTRD*((M_ONE+z)**M_THIRD-(M_ONE-z)**M_THIRD)/TFTM 
    ecu = ec
    vcu = vc(1)
    if(rs>M_ONE) then
       ecp = eclow(2)
       vcp = ecp*potlow(2)
    else
       ecp = echigh(2)
       vcp = pothigh(2)
    endif
    ec = ecu + fz*(ecp-ecu)
    x = vcu + fz*(vcp-vcu) - z*(ecp-ecu)*fzp
    vc(1) = x + (ecp-ecu)*fzp
    vc(2) = x - (ecp-ecu)*fzp

  end subroutine correlation_lda_3d_pz
#undef eclow
#undef echigh
#undef potlow
#undef pothigh

end module vxc

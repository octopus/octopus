#include "global.h"

module vxc

  use global

  implicit none

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Finds local exchange energy density and potential                          !
!  Adapted by J.M.Soler from routine velect of Froyens                        !
!    pseudopotential generation program. Madrid, Jan 97. Version 0.5.         !
! **** Input ******************************************************           !
! INTEGER IREL    : relativistic-exchange switch (0=no, 1=yes)                !
! INTEGER NSP     : spin-polarizations (1=>unpolarized, 2=>polarized)         !
! REAL*8  DS(NSP) : total (nsp=1) or spin (nsp=2) electron density            !
! **** Output *****************************************************           !
! REAL*8  EX      : exchange energy density                                   !
! REAL*8  VX(NSP) : (spin-dependent) exchange potential                       !
! **** Units ******************************************************           !
! Densities in electrons/Bohr**3                                              !
! Energies in Hartrees                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine exchng( irel, nsp, ds, ex, vx )
    integer, intent(in)   :: irel, nsp
    FLOAT, intent(in)  :: ds(nsp)
    FLOAT, intent(out) :: vx(nsp), ex

    FLOAT, parameter :: opf = CNST(1.5), c014 = CNST(0.014),  &
        ftrd = M_FOUR/M_THREE, alp  = M_TWO*M_THIRD
    FLOAT :: d1, d2, d, z, fz, fzp, rs, vxp, exp, vxf, exf,   &
                beta, sb, alb, tftm, a0

    a0   = (M_FOUR/(9*M_PI))**M_THIRD
    tftm = M_TWO**FTRD - M_TWO

    IF (NSP .EQ. 2) THEN
      D1 = MAX(DS(1), M_ZERO)
      D2 = MAX(DS(2), M_ZERO)
      D = D1 + D2
      IF (D .LE. M_ZERO) THEN
        EX    = M_ZERO
        VX(1) = M_ZERO
        VX(2) = M_ZERO
        RETURN
      ENDIF
      Z = (D1 - D2) / D
      FZ = ((1+Z)**FTRD+(1-Z)**FTRD-2)/TFTM
      FZP = FTRD*((1+Z)**M_THIRD-(1-Z)**M_THIRD)/TFTM 
    ELSE
      D = DS(1)
      IF (D .LE. M_ZERO) THEN
        EX = M_ZERO
        VX(1) = M_ZERO
        RETURN
      ENDIF
      Z = M_ZERO
      FZ = M_ZERO
      FZP = M_ZERO
    ENDIF
    RS = (3 / (4*M_PI*D) )**M_THIRD
    VXP = -(3*ALP/(2*M_PI*A0*RS))
    EXP = 3*VXP/4
    IF (IREL .EQ. 1) THEN
      BETA = C014/RS
      SB = SQRT(1+BETA*BETA)
      ALB = LOG(BETA+SB)
      VXP = VXP * (-M_HALF + OPF * ALB / (BETA*SB))
      EXP = EXP * (M_ONE-OPF*((BETA*SB-ALB)/BETA**2)**2) 
    ENDIF
    VXF = 2**M_THIRD*VXP
    EXF = 2**M_THIRD*EXP
    IF (NSP .EQ. 2) THEN
      VX(1) = VXP + FZ*(VXF-VXP) + (1-Z)*FZP*(EXF-EXP)
      VX(2) = VXP + FZ*(VXF-VXP) - (1+Z)*FZP*(EXF-EXP)
      EX    = EXP + FZ*(EXF-EXP)
    ELSE
      VX(1) = VXP
      EX    = EXP
    ENDIF
    
  end subroutine exchng


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
     CALL EXCHNG(IREL, NS, DD, EPSX, VXD)
     CALL PZC   (      NS, DD, EPSC, VCD)
  ELSEIF ( AUTHOR.EQ.'PW92' .OR. AUTHOR.EQ.'pw92' ) THEN
     CALL EXCHNG(IREL, NS, DD, EPSX, VXD)
     CALL PW92C (      NS, DD, EPSC, VCD )
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

  FLOAT :: A, BETA, D(2), DADD, DECUDD,                                       & 
              DF1DD, DF2DD, DF3DD, DF4DD, DF1DGD, DF3DGD, DF4DGD,             &
              DFCDD(2), DFCDGD(3,2), DFDD, DFDGD, DFXDD(2), DFXDGD(3,2),      &
              DHDD, DHDGD, DKFDD, DKSDD, DPDD, DPDZ, DRSDD,                   &
              DS(2), DSDD, DSDGD, DT, DTDD, DTDGD, DZDD(2),                   &
              EC, ECUNIF, EX, EXUNIF,                                         &
              F, F1, F2, F3, F4, FC, FX,                                      &
              GAMMA, GD(3,2), GDM(2), GDMS, GDMT, GDS, GDT(3),                &
              H, KAPPA, KF, KFS, KS, MU, PHI, RS, S,                          &
              T, VCUNIF(2), VXUNIF(2), ZETA

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
    CALL EXCHNG( IREL, 1, DS(IS), EXUNIF, VXUNIF(IS) )
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

  FLOAT :: A, BETA, D(2), DADD, DECUDD,                                       & 
              DF1DD, DF2DD, DF3DD, DF4DD, DF1DGD, DF3DGD, DF4DGD,             &
              DFCDD(2), DFCDGD(3,2), DFDD, DFDGD, DFXDD(2), DFXDGD(3,2),      &
              DHDD, DHDGD, DKFDD, DKSDD, DPDD, DPDZ, DRSDD,                   &
              DS(2), DSDD, DSDGD, DT, DTDD, DTDGD, DZDD(2),                   &
              EC, ECUNIF, EX, EXUNIF,                                         &
              F, F1, F2, F3, F4, FC, FX,                                      &
              GAMMA, GD(3,2), GDM(2), GDMS, GDMT, GDS, GDT(3),                &
              H, KAPPA, KF, KFS, KS, MU, PHI, RS, S,                          &
              T, VCUNIF(2), VXUNIF(2), ZETA

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
  CALL PW92C( 2, D, ECUNIF, VCUNIF )

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
    FLOAT :: alpha, gdm, gd(3, nspin), x, f, gamma

    ! First, get the LDA exchange potential.
    call exchng(0, nspin, dens, x, dexdd)

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
                x**2/(M_ONE + M_THREE*beta*x*oct_asinh(gamma*x))
           dexdd(is) = dexdd(is) + f !- beta * dens(is)**M_THIRD * f
      elseif(r > M_ZERO .and. dens(is)<= threshold) then
           f = r + (M_THREE/alpha)*log(2*gamma*alpha*qtot**(-M_THIRD))
           !f = f + (qtot*exp(-alpha*r))**M_THIRD/(beta*alpha**2)
           dexdd(is) = dexdd(is) - M_ONE/f
      endif
    enddo

  end subroutine lb94

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Implements the Perdew-Wang 92 local correlation (beyond RPA).               !
  ! Ref: J.P.Perdew & Y.Wang, PRB, 45, 13244 (1992)                             !
  ! Written by L.C.Balbas and J.M.Soler. Dec 96.  Version 0.5.                  !
  ! ********* INPUT ****************************************************        !
  ! INTEGER NSPIN       : Number of spin polarizations (1 or 2)                 !
  ! REAL*8  DENS(NSPIN) : Local (spin) density                                  !
  ! ********* OUTPUT ***************************************************        !
  ! REAL*8  EC        : Correlation energy density                              !
  ! REAL*8  VC(NSPIN) : Correlation (spin) potential                            !
  ! ********* UNITS ****************************************************        ! 
  ! Densities in electrons per Bohr**3                                          ! 
  ! Energies in Hartrees                                                        !
  ! ********* ROUTINES CALLED ******************************************        !
  ! None                                                                        !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pw92c( NSPIN, DENS, EC, VC )
  implicit none 
  integer :: NSPIN
  FLOAT :: DENS(NSPIN), EC, VC(NSPIN)

  ! Internal variable declarations
  integer :: IG
  FLOAT :: A(0:2), ALPHA1(0:2), B, BETA(0:2,4), C,     &
       DBDRS, DECDD(2), DECDRS, DECDZ, DFDZ,           &
       DGDRS(0:2), DCDRS, DRSDD, DTOT, DZDD(2),        &
       F, FPP0, G(0:2), P(0:2), RS, ZETA

  ! Fix lower bound of density to avoid division by zero
  FLOAT, parameter :: DENMIN = CNST(1e-12), &
      FOUTHD = M_FOUR/M_THREE, THRHLF = M_THREE/M_TWO

  ! Parameters from Table I of Perdew & Wang, PRB, 45, 13244 (92)
  DATA P      / CNST(1.00    ), CNST( 1.00    ), CNST( 1.00    )  /
  DATA A      / CNST(0.031091), CNST( 0.015545), CNST( 0.016887)  /
  DATA ALPHA1 / CNST(0.21370 ), CNST( 0.20548 ), CNST( 0.11125 )  /
  DATA BETA   / CNST(7.5957  ), CNST(14.1189  ), CNST(10.357   ), &
                CNST(3.5876  ), CNST( 6.1977  ), CNST( 3.6231  ), &
                CNST(1.6382  ), CNST( 3.3662  ), CNST( 0.88026 ), &
                CNST(0.49294 ), CNST( 0.62517 ), CNST( 0.49671 )  /    

  ! Find rs and zeta
  IF (NSPIN .EQ. 1) THEN
     DTOT = MAX( DENMIN, DENS(1) )
     ZETA = 0
     RS = ( 3 / (4*M_PI*DTOT) )**M_THIRD
  !      Find derivatives dRs/dDens and dZeta/dDens
     DRSDD = (- RS) / DTOT / 3
     DZDD(1) = 0
  ELSE
     DTOT = MAX( DENMIN, DENS(1)+DENS(2) )
     ZETA = ( DENS(1) - DENS(2) ) / DTOT
     RS = ( 3 / (4*M_PI*DTOT) )**M_THIRD
     DRSDD = (- RS) / DTOT / 3
     DZDD(1) =   1 / DTOT - ZETA / DTOT
     DZDD(2) = - (1 / DTOT) - (ZETA / DTOT)
  ENDIF

  ! Find eps_c(rs,0)=G(0), eps_c(rs,1)=G(1) and -alpha_c(rs)=G(2)
  ! using eq.(10) of cited reference (Perdew & Wang, PRB, 45, 13244 (92))

  DO 20 IG = 0,2
     B = BETA(IG,1) * RS**M_HALF   +                                            &
         BETA(IG,2) * RS         +                                            &
         BETA(IG,3) * RS**THRHLF +                                            &
         BETA(IG,4) * RS**(P(IG)+1)
     DBDRS = BETA(IG,1) * M_HALF      / RS**M_HALF +                              &
         BETA(IG,2)                         +                                 &
         BETA(IG,3) * THRHLF    * RS**M_HALF +                                  &
         BETA(IG,4) * (P(IG)+1) * RS**P(IG)
     C = 1 + 1 / (2 * A(IG) * B)
     DCDRS = - ( (C-1) * DBDRS / B )
     G(IG) = (- 2) * A(IG) * ( 1 + ALPHA1(IG)*RS ) * LOG(C)
     DGDRS(IG) = (- 2) *A(IG) * ( ALPHA1(IG) * LOG(C) +                       &
                                (1+ALPHA1(IG)*RS) * DCDRS / C )
  20 CONTINUE

  ! Find f''(0) and f(zeta) from eq.(9)
  C = 1 / (2**FOUTHD - 2)
  FPP0 = 8 * C / 9
  F = ( (1+ZETA)**FOUTHD + (1-ZETA)**FOUTHD - 2 ) * C
  DFDZ = FOUTHD * ( (1+ZETA)**M_THIRD - (1-ZETA)**M_THIRD ) * C

  ! Find eps_c(rs,zeta) from eq.(8)
  EC = G(0) - G(2) * F / FPP0 * (1-ZETA**4) +                                 &
       (G(1)-G(0)) * F * ZETA**4
  DECDRS = DGDRS(0) - DGDRS(2) * F / FPP0 * (1-ZETA**4) +                     &
           (DGDRS(1)-DGDRS(0)) * F * ZETA**4
  DECDZ = (- G(2)) / FPP0 * ( DFDZ*(1-ZETA**4) - F*4*ZETA**3 ) +              &
          (G(1)-G(0)) * ( DFDZ*ZETA**4 + F*4*ZETA**3 )
      
  ! Find correlation potential
  IF (NSPIN .EQ. 1) THEN
     DECDD(1) = DECDRS * DRSDD
     VC(1) = EC + DTOT * DECDD(1)
  ELSE
     DECDD(1) = DECDRS * DRSDD + DECDZ * DZDD(1)
     DECDD(2) = DECDRS * DRSDD + DECDZ * DZDD(2)
     VC(1) = EC + DTOT * DECDD(1)
     VC(2) = EC + DTOT * DECDD(2)
  ENDIF

  end subroutine pw92c 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Perdew-Zunger parameterization of Ceperley-Alder                           !
  !  correlation. Ref: Perdew & Zunger, Phys. Rev. B 23 5075 (1981).            !
  !  Adapted by J.M.Soler from routine velect of Froyen s                       !
  !    pseudopotential generation program.                                      !
  ! **** Input *****************************************************            !
  ! INTEGER NSP     : spin-polarizations (1=>unpolarized, 2=>polarized)         !
  ! REAL*8  DS(NSP) : total (nsp=1) or spin (nsp=2) electron density            !
  ! **** Output *****************************************************           !
  ! REAL*8  EC      : correlation energy density                                !
  ! REAL*8  VC(NSP) : (spin-dependent) correlation potential                    !
  ! **** Units *******************************************************          !
  ! Densities in electrons/Bohr**3                                              !
  ! Energies in Hartrees                                                        !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pzc( NSP, DS, EC, VC )
    integer, intent(in)   :: nsp
    FLOAT, intent(in)  :: ds(nsp)
    FLOAT, intent(out) :: ec, vc(nsp)

    FLOAT, parameter :: &
        FIVE  = CNST(0.50    ), OPF   = CNST(1.50),    PNN  = CNST(.99    ), &
        PTHREE= CNST( 0.30   ), PSEVF = CNST(0.750),   C0504= CNST(0.0504 ), &
        C0254 = CNST( 0.02540), C014  = CNST(0.0140),  C0406= CNST(0.04060), &
        C15P9 = CNST(15.90   ), C0666 = CNST(0.06660), C11P4= CNST(11.40  ), &
        C045  = CNST( 0.0450 ), C7P8  = CNST(7.80),    C88  = CNST(0.880  ), C20P59 = CNST(20.5920), &
        C3P52 = CNST( 3.520  ), C0311 = CNST(0.03110), C0014= CNST(0.00140), &
        C0538 = CNST( 0.05380), C0096 = CNST(0.00960), C096 = CNST(0.0960 ), &
        C0622 = CNST( 0.06220), C004  = CNST(0.0040),  C0232= CNST(0.02320), &
        C1686 = CNST( 0.16860), C1P398= CNST(1.39810), C2611= CNST(0.26110), &
        C2846 = CNST( 0.28460), C1P053= CNST(1.05290), C3334= CNST(0.33340)

    ! Ceperly-Alder 'ca' constants. Internal energies in Rydbergs.
    FLOAT, parameter :: &
        CON1 = M_ONE/M_SIX,           CON2 = CNST(0.0080 )/M_THREE, &
        CON3 = CNST(0.35020)/M_THREE, CON4 = CNST(0.05040)/M_THREE, &
        CON5 = CNST(0.00280)/M_THREE, CON6 = CNST(0.19250)/M_THREE, &
        CON7 = CNST(0.02060)/M_THREE, CON8 = CNST(9.78670)/M_SIX,   &
        CON9 = CNST(1.0444 )/M_THREE, CON10= CNST(7.37030)/M_SIX,   &
        CON11= CNST(1.33360)/M_THREE

    ! X-alpha parameter:
    FLOAT, parameter :: ALP = M_TWOTHIRD

    ! Other variables converted into parameters by J.M.Soler
    FLOAT, parameter :: TRD  = M_THIRD, FTRD = M_FOUR / M_THREE, &
        TFTM = CNST(0.519842099789746380 ),  A0 = CNST(0.521061761197848080), &
        CRS  = CNST(0.6203504908994000870), CXP = -M_THREE*ALP/(M_PI*A0),     &
        CXF  = CNST(1.259921049894873190 )

    FLOAT :: d1, d2, d, z, fz, fzp, rs, sqrs, te, be, ecp, vcp, ecf, vcf, rslog
    integer :: isp

    !      Find density and polarization
    IF (NSP .EQ. 2) THEN
     D1 = MAX(DS(1), M_ZERO)
     D2 = MAX(DS(2), M_ZERO)
     D = D1 + D2
     IF (D .LE. M_ZERO) THEN
        EC = M_ZERO
        VC(1) = M_ZERO
        VC(2) = M_ZERO
        RETURN
     ENDIF
     Z = (D1 - D2) / D
     FZ = ((1+Z)**FTRD+(1-Z)**FTRD-2)/TFTM
     FZP = FTRD*((1+Z)**TRD-(1-Z)**TRD)/TFTM 
    ELSE
     D = DS(1)
     IF (D .LE. M_ZERO) THEN
        EC = M_ZERO
        VC(1) = M_ZERO
        RETURN
     ENDIF
     Z = M_ZERO
     FZ = M_ZERO
     FZP = M_ZERO
    ENDIF
    RS = CRS / D**TRD

    !      Correlation 
    IF (RS .GT. M_ONE) THEN  
      SQRS=SQRT(RS)
      TE = M_ONE+CON10*SQRS+CON11*RS
      BE = M_ONE+C1P053*SQRS+C3334*RS
      ECP = -(C2846/BE)
      VCP = ECP*TE/BE
      TE = M_ONE+CON8*SQRS+CON9*RS
      BE = M_ONE+C1P398*SQRS+C2611*RS
      ECF = -(C1686/BE)
      VCF = ECF*TE/BE
    ELSE
      RSLOG=LOG(RS)
      ECP=(C0622+C004*RS)*RSLOG-C096-C0232*RS
      VCP=(C0622+CON2*RS)*RSLOG-CON3-CON4*RS
      ECF=(C0311+C0014*RS)*RSLOG-C0538-C0096*RS
      VCF=(C0311+CON5*RS)*RSLOG-CON6-CON7*RS
    ENDIF

    !      Find up and down potentials
    IF (NSP .EQ. 2) THEN
      EC    = ECP + FZ*(ECF-ECP)
      VC(1) = VCP + FZ*(VCF-VCP) + (1-Z)*FZP*(ECF-ECP)
      VC(2) = VCP + FZ*(VCF-VCP) - (1+Z)*FZP*(ECF-ECP)
    ELSE
      EC    = ECP
      VC(1) = VCP
    ENDIF

    !      Change from Rydbergs to Hartrees
    EC = M_HALF * EC
    do ISP = 1,NSP
       VC(ISP) = M_HALF * VC(ISP)
    end do

  end subroutine pzc 

end module vxc

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

! Return LDA exchange energy and potential
! density, potential and energy in a.u.
subroutine xc_x_lda(rel, nspin, rho, vx, ex)
  logical, intent(in) :: rel
  integer, intent(in) :: nspin
  real(r8), intent(IN) :: rho(1:nspin)
  real(r8), intent(inout) :: vx(1:nspin), ex

  real(r8), parameter :: &
      HALF=.5_r8, OPF=1.5_r8, C014=0.014_r8, &
      TRD = ONE/3.0_r8, FTRD = 4.0_r8*TRD, ALP = 2.0_r8 * TRD, &
      TFTM = 0.519842099789746380_r8, A0   = 0.521061761197848080_r8, &
      CRS  = 0.6203504908994000870_r8, CXP  = (- 3.0_r8) * ALP / (M_PI*A0), &
      CXF  = 1.259921049894873190_r8

  real(r8) :: d1, d2, d, z, fz, fzp, rs, vxp, exp, beta, sb, alb, vxf, exf, bxc, &
              bin, vin, m(3)

  if(.not.fzeta(nspin, rho, d, z, fz, fzp)) then
    vx = ZERO; ex = ZERO; return
  endif

  rs = CRS / d**TRD
  vxp = CXP / rs
  exp = 0.750_r8*vxp
  if(rel) then
    beta = C014/rs
    sb = sqrt(1 + beta*beta)
    alb = log(beta + sb)
    vxp = vxp * (-HALF + OPF * alb / (beta*sb))
    exp = exp * (ONE - OPF*((beta*sb-alb)/beta**2)**2) 
  endif
  vxf = CXF * vxp
  exf = CXF * exp

  vin = vxp + fz*(vxf-vxp)-z*fzp*(exf-exp)
  bin = fzp*(exf-exp)
  if(nspin>2) m(1:3) = rho(2:4)
  call xc_matrix(nspin, vin, bin, m, vx)
  ex = exp + fz*(exf - exp)

  !Change from Rydbergs to Hartrees
  ex    = HALF * ex
  vx(:) = HALF * vx(:) 
  return
end subroutine xc_x_lda


! LSD Perdew & Zunger, PRB 23, 5048 (1981)
subroutine xc_c_pz(nspin, rho, vc, ec)
  integer, intent(in) :: nspin
  real(r8), intent(IN) :: rho(nspin)
  real(r8), intent(inout) :: vc(nspin), ec

!      X-alpha parameter:
  real(r8), parameter :: ALP = 2.0_r8 / 3.0_r8

  real(r8), parameter :: TRD  = 1.0_r8 / 3.0_r8,        &
       FTRD = 4.0_r8 / 3.0_r8, TFTM = 0.519842099789746380_r8,          &
       A0   = 0.521061761197848080_r8, CRS = 0.6203504908994000870_r8

  real(r8), parameter :: &
      ZERO=0.0_r8,      HALF = 0.5_r8, ONE = 1.0_r8, &
      C0311=0.03110_r8, C0014=0.00140_r8, &
      C0538=0.05380_r8, C0096=0.00960_r8, C096=0.0960_r8, &
      C0622=0.06220_r8, C004=0.0040_r8, C0232=0.02320_r8, &
      C1686=0.16860_r8, C1P398=1.39810_r8, C2611=0.26110_r8, &
      C2846=0.28460_r8, C1P053=1.05290_r8, C3334=0.33340_r8

!    Ceperly-Alder 'ca' constants. Internal energies in Rydbergs.
  real(r8), parameter :: &
      CON2=0.0080_r8/3, CON3=0.35020_r8/3, &
      CON4=0.05040_r8/3, CON5=0.00280_r8/3, CON6=0.19250_r8/3, &
      CON7=0.02060_r8/3, CON8=9.78670_r8/6, CON9=1.0444_r8/3, &
      CON10=7.37030_r8/6, CON11=1.33360_r8/3

  real(r8) :: d, d1, d2, z, fz, fzp, rs, sqrs, te, be, ecp, vcp, &
      ecf, vcf, rslog, bcc, bin, vin, m(3)

  if(.not.fzeta(nspin, rho, d, z, fz, fzp)) then
    vc = ZERO; ec = ZERO; return
  endif

  rs = CRS / d**TRD
  if (rs .gt. ONE) then  
    sqrs = sqrt(rs)
    te   = ONE + CON10*sqrs  + CON11*rs
    be   = ONE + C1P053*sqrs + C3334*rs
    ecp  = -(C2846/be)
    vcp  = ecp*te/be
    te   = ONE + CON8*sqrs   + CON9*rs
    be   = ONE + C1P398*sqrs + C2611*rs
    ecf  = -(C1686/be)
    vcf  = ecf*te/be
  else
    rslog = log(rs)
    ecp   = (C0622+C004 *rs)*rslog - C096  - C0232*rs
    vcp   = (C0622+CON2 *rs)*rslog - CON3  - CON4*rs
    ecf   = (C0311+C0014*rs)*rslog - C0538 - C0096*rs
    vcf   = (C0311+CON5 *rs)*rslog - CON6  - CON7*rs
  endif

  vin = vcp + fz*(vcf-vcp)-z*fzp*(ecf-ecp)
  bin = fzp*(ecf-ecp)
  if(nspin>2) m(1:3) = rho(2:4)
  call xc_matrix(nspin, vin, bin, m, vc)
  ec = ecp + fz*(ecf - ecp)

  !Change from Rydbergs to Hartrees (Yeah, the Good Old Rydbergs...)
  ec = HALF * ec
  vc = HALF * vc
  return
end subroutine xc_c_pz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Implements the Perdew-Wang 92 local correlation (beyond RPA).
! Ref: J.P.Perdew & Y.Wang, PRB, 45, 13244 (1992)
! Written by L.C.Balbas and J.M.Soler. Dec 96.  Version 0.5.
!
! Modified slightly later.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine xc_c_pw92( NSPIN, DENS, VC, EC)
  integer :: NSPIN
  real(r8) :: DENS(NSPIN), EC, VC(NSPIN)
  ! Internal variable declarations
  integer :: IG
  real(r8) :: A(0:2), ALPHA1(0:2), B, BETA(0:2,4), C,         &
      DBDRS, DECDD(2), DECDRS, DECDZ, DENMIN, DFDZ,           &
      DGDRS(0:2), DCDRS, DRSDD, DTOT, DZDD(2),                &
      F, FPP0, FOUTHD, G(0:2), HALF,                          &
      P(0:2), RS, THD, THRHLF, ZETA, bin, vin, m(3)
  ! Fix lower bound of density to avoid division by zero
  PARAMETER ( DENMIN = 1.E-12 )
  ! Fix some numerical constants
  PARAMETER ( FOUTHD=4.0_r8/3.0_r8, HALF=0.50_r8,             &
      THD=1.0_r8/3.0_r8, THRHLF=1.50_r8 )
  ! Parameters from Table I of Perdew & Wang, PRB, 45, 13244 (92)
  DATA P      / 1.00_r8,     1.00_r8,     1.00_r8     /
  DATA A      / 0.031091_r8, 0.015545_r8, 0.016887_r8 /
  DATA ALPHA1 / 0.21370_r8,  0.20548_r8,  0.11125_r8  /
  DATA BETA   / 7.5957_r8,  14.1189_r8,  10.357_r8,           &
      3.5876_r8,   6.1977_r8,   3.6231_r8,                    &
      1.6382_r8,   3.3662_r8,   0.88026_r8,                   &
      0.49294_r8,  0.62517_r8,  0.49671_r8 /    

  ! Find rs and zeta
  IF (NSPIN .EQ. 1) THEN
    DTOT = MAX( DENMIN, DENS(1) )
    ZETA = 0
    RS = ( 3 / (4*M_PI*DTOT) )**THD
    ! Find derivatives dRs/dDens and dZeta/dDens
    DRSDD = (- RS) / DTOT / 3
    DZDD(1) = 0
  ELSE
    DTOT = MAX(DENMIN, DENS(1))
    IF (NSPIN .EQ. 2) THEN
       ZETA = DENS(2) / DTOT
    ELSE
       ZETA = sqrt(DENS(2)**2+DENS(3)**2+DENS(4)**2)/DTOT
    ENDIF
    RS = ( 3 / (4*M_PI*DTOT) )**THD
    DRSDD = (- RS) / DTOT / 3
    DZDD(1) =   1 / DTOT - ZETA / DTOT
    DZDD(2) = - (1 / DTOT) - (ZETA / DTOT)
  ENDIF

  ! Find eps_c(rs,0)=G(0), eps_c(rs,1)=G(1) and -alpha_c(rs)=G(2)
  ! using eq.(10) of cited reference (Perdew & Wang, PRB, 45, 13244 (92))
  DO IG = 0,2
    B = BETA(IG,1) * RS**HALF   +                                            &
        BETA(IG,2) * RS         +                                            &
        BETA(IG,3) * RS**THRHLF +                                            &
        BETA(IG,4) * RS**(P(IG)+1)
    DBDRS = BETA(IG,1) * HALF      / RS**HALF +                              &
        BETA(IG,2)                         +                                 &
        BETA(IG,3) * THRHLF    * RS**HALF +                                  &
        BETA(IG,4) * (P(IG)+1) * RS**P(IG)
    C = 1 + 1 / (2 * A(IG) * B)
    DCDRS = - ( (C-1) * DBDRS / B )
    G(IG) = (- 2) * A(IG) * ( 1 + ALPHA1(IG)*RS ) * LOG(C)
    DGDRS(IG) = (- 2) *A(IG) * ( ALPHA1(IG) * LOG(C) +                       &
        (1+ALPHA1(IG)*RS) * DCDRS / C )
  end DO
  
  ! Find f''(0) and f(zeta) from eq.(9)
  C = 1 / (2**FOUTHD - 2)
  FPP0 = 8 * C / 9
  F = ( (1+ZETA)**FOUTHD + (1-ZETA)**FOUTHD - 2 ) * C
  DFDZ = FOUTHD * ( (1+ZETA)**THD - (1-ZETA)**THD ) * C

  ! Find eps_c(rs,zeta) from eq.(8)
  EC = G(0) - G(2) * F / FPP0 * (1-ZETA**4) +                                 &
      (G(1)-G(0)) * F * ZETA**4
  DECDRS = DGDRS(0) - DGDRS(2) * F / FPP0 * (1-ZETA**4) +                     &
      (DGDRS(1)-DGDRS(0)) * F * ZETA**4
  DECDZ = (- G(2)) / FPP0 * ( DFDZ*(1-ZETA**4) - F*4*ZETA**3 ) +              &
      (G(1)-G(0)) * ( DFDZ*ZETA**4 + F*4*ZETA**3 )

  vin = ec + dtot*decdrs*drsdd - zeta*decdz
  bin = decdz
  if(nspin>2) m(1:3) = dens(2:4)
  call xc_matrix(nspin, vin, bin, m, vc)
  
  return
end subroutine xc_c_pw92

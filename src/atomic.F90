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

#include "global.h"

module atomic
  use vxc
  use logrid

  implicit none

  private
  public :: atomxc, atomhxc, vhrtre, egofv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Next stuff is for the "valence_conf" data type, made to handle atomic configurations.
  public :: valconf, &
            valconf_copy, &
            write_valconf, &
            read_valconf, &
            valconf_null, &
            VALCONF_STRING_LENGTH
  character(len=1), parameter :: spec_notation(0:3) = (/ 's', 'p', 'd', 'f' /)
  integer, parameter :: VALCONF_STRING_LENGTH = 80

  type valconf
    integer           :: z
    character(len=2)  :: symbol
    integer           :: type     ! 0 for the most normal valence configuration, 1 for semicore.
    integer           :: p        ! number of orbitals.
    integer           :: n(6)     ! n quantum number
    integer           :: l(6)     ! l quantum number
    real(r8)          :: occ(6,2) ! occupations of each level
  end type
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutines to write and read valence configurations.
  subroutine valconf_null(c)
    type(valconf) :: c
    c%z = 0; c%symbol = ""; c%type = 0; c%p = 0; c%n = 0; c%l = 0; c%occ = M_ZERO
  end subroutine valconf_null

  subroutine valconf_copy(cout, cin)
    type(valconf), intent(out) :: cout
    type(valconf), intent(in)  :: cin
    cout%z      = cin%z
    cout%symbol = cin%symbol
    cout%type   = cin%type
    cout%p      = cin%p
    cout%n      = cin%n
    cout%l      = cin%l
    cout%occ    = cin%occ
  end subroutine valconf_copy

  subroutine write_valconf(c, s)
    type(valconf), intent(in) :: c
    character(len=VALCONF_STRING_LENGTH) :: s
    integer :: j
    write(s,'(i2,1x,a2,i1,1x,i1,a1,6(i1,a1,f6.3,a1))') c%z, c%symbol, c%type, c%p, ':',&
         (c%n(j),spec_notation(c%l(j)),c%occ(j,1),',',j=1,c%p)
  end subroutine write_valconf

  subroutine read_valconf(s, c)
    character(len=VALCONF_STRING_LENGTH), intent(in) :: s
    type(valconf), intent(out) :: c
    integer :: j
    character(len=1) :: lvalues(1:6)
    read (s,'(i2,1x,a2,i1,1x,i1,1x,6(i1,a1,f6.3,1x))') c%z, c%symbol, c%type, c%p,&
         (c%n(j),lvalues(j),c%occ(j,1),j=1,c%p)
    do j = 1, c%p
       select case(lvalues(j))
       case('s'); c%l(j) = 0
       case('p'); c%l(j) = 1
       case('d'); c%l(j) = 2
       case('f'); c%l(j) = 3
       case default; stop 'Error'
       end select
    enddo
  end subroutine read_valconf

  subroutine atomhxc(functl, irel, g, nspin, dens, v, extra)
    character(len=*), intent(in)   :: functl
    integer, intent(in)            :: irel, nspin
    type(logrid_type), intent(in)  :: g
    real(r8), intent(in)           :: dens(g%nrval, nspin)
    real(r8), intent(out)          :: v(g%nrval, nspin)
    real(r8), intent(in), optional :: extra(g%nrval)

    character(len=5) :: xcfunc, xcauth
    integer :: is, ir
    real(r8), allocatable :: xc(:, :), ve(:, :), rho(:, :)
    real(r8) :: r2, ex, ec, dx, dc

    call push_sub('atomhxc')

    allocate(ve(g%nrval, nspin), xc(g%nrval, nspin), rho(g%nrval, nspin))
             ve = M_ZERO; xc = M_ZERO; rho = M_ZERO

    ! To calculate the Hartree term, we put all the density in one variable.
    do is = 1, nspin
       rho(:, 1) = rho(:, 1) + dens(:, is)
    enddo
    ve = M_ZERO
    do is = 1, nspin
       call vhrtre(rho(:, 1), ve(:, is), g%rofi, g%drdi, g%s, g%nrval, g%a)
    enddo
    ve(1, 1:nspin) = ve(2, 1:nspin)

    select case(functl)
    case('LDA')
      xcfunc = 'LDA'; xcauth = 'PZ'
    case('GGA')
      xcfunc = 'GGA'; xcauth = 'PBE'
    case default
      message(1) = 'Internal Error'
      call write_fatal(1)
    end select

    rho = dens
    do is = 1, nspin
       if(present(extra)) rho(:, is) = rho(:, is) + extra(:)/nspin
       rho(2:g%nrval, is) = rho(2:g%nrval, is)/(M_FOUR*M_PI*g%rofi(2:g%nrval)**2)
    enddo
    r2 = g%rofi(2)/(g%rofi(3)-g%rofi(2))
    rho(1, 1:nspin) = rho(2, 1:nspin) - (rho(3, 1:nspin)-rho(2, 1:nspin))*r2

    do is = 1, nspin
       do ir = 1, g%nrval
          if(rho(ir, is) < 1.0e-15_r8) rho(ir, is) = 0.0_r8
       enddo
    enddo

    call atomxc(xcfunc, xcauth, irel, g%nrval, g%nrval, g%rofi, &
                nspin, rho, ex, ec, dx, dc, xc)

    v = ve + xc

    deallocate(ve, xc, rho)
    call pop_sub()
  end subroutine atomhxc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finds total exchange-correlation energy and potential for a                 !
! spherical electron density distribution.                                    !
! This version implements the Local (spin) Density Approximation and          !
! the Generalized-Gradient-Aproximation with the explicit mesh                !
! functional method of White & Bird, PRB 50, 4954 (1994).                     !
! Gradients are defined by numerical derivatives, using 2*NN+1 mesh         ! 
!   points, where NN is a parameter defined below                             !
! Coded by L.C.Balbas and J.M.Soler. December 1996. Version 0.5.              !
!                                                                             !
! CHARACTER*(*) FUNCTL : Functional to be used:                               !
!              LDA or LSD => Local (spin) Density Approximation           !
!                     GGA => Generalized Gradient Corrections             !
!                                Uppercase is optional                        !
! CHARACTER*(*) AUTHOR : Parametrization desired:                             !
!     CA or PZ => LSD Perdew & Zunger, PRB 23, 5075 (1981)                !
!         PW92 => LSD Perdew & Wang, PRB, 45, 13244 (1992). This is       !
!                     the local density limit of the next:                    !
!          PBE => GGA Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)      !
!                     Uppercase is optional                                   !
! INTEGER IREL         : Relativistic exchange? (0=>no, 1=>yes)               !
! INTEGER NR           : Number of radial mesh points                         !
! INTEGER MAXR         : Physical first dimension of RMESH, DENS and V_XC     !
! REAL*8  RMESH(MAXR)  : Radial mesh points                                   !
! INTEGER NSPIN        : NSPIN=1 => unpolarized; NSPIN=2 => polarized         !
! REAL*8  DENS(MAXR,NSPIN) : Total (NSPIN=1) or spin (NSPIN=2) electron       !
!                            density at mesh points                           !
! ************************* OUTPUT **********************************         !
! REAL*8  EX              : Total exchange energy                             !
! REAL*8  EC              : Total correlation energy                          !
! REAL*8  DX              : IntegralOf( rho * (eps_x - v_x) )                 !
! REAL*8  DC              : IntegralOf( rho * (eps_c - v_c) )                 !
! REAL*8  V_XC(MAXR,NSPIN): (Spin) exch-corr potential                        !
! ************************ UNITS ************************************         !
! Distances in atomic units (Bohr).                                           !
! Densities in atomic units (electrons/Bohr**3)                               !
! Energy unit depending of parameter EUNIT below                              !
! ********* ROUTINES CALLED *****************************************         !
! GGAXC, LDAXC                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine atomxc( FUNCTL, AUTHOR, IREL,                                    &
                     NR, MAXR, RMESH, NSPIN, DENS,                            &
                     EX, EC, DX, DC, V_XC )
  use global

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Argument types and dimensions                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none

  character(len=*) :: FUNCTL, AUTHOR
  integer :: IREL, MAXR, NR, NSPIN
  real(r8) :: DENS(MAXR,NSPIN), RMESH(MAXR), V_XC(MAXR,NSPIN)
  real(r8) :: DC, DX, EC, EX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal parameters                                                         !
! NN     : order of the numerical derivatives: the number of radial           !
!          points used is 2*NN+1                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, parameter :: NN = 5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fix energy unit:  EUNIT=1.0 => Hartrees,                                    !
!                   EUNIT=0.5 => Rydbergs,                                    !
!                   EUNIT=0.03674903 => eV                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(r8), parameter :: EUNIT = 0.5_r8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DVMIN is added to differential of volume to avoid division by zero          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(r8), parameter :: DVMIN = 1.0E-12

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Local variables and arrays                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical :: GGA 
  integer :: IN, IN1, IN2, IR, IS, JN 
  real(r8) :: AUX(MAXR), D(NSPIN), DECDD(NSPIN), DECDGD(3,NSPIN),             & 
              DEXDD(NSPIN), DEXDGD(3,NSPIN),                                  &
              DGDM(-NN:NN), DGIDFJ(-NN:NN), DRDM, DVOL,                       & 
              EPSC, EPSX, F1, F2, GD(3,NSPIN), PI
  !external :: GGAXC, LDAXC 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set GGA switch                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF ( FUNCTL.EQ.'LDA' .OR. FUNCTL.EQ.'lda' .OR.                              &
       FUNCTL.EQ.'LSD' .OR. FUNCTL.EQ.'lsd' ) THEN
     GGA = .FALSE.
  ELSEIF ( FUNCTL.EQ.'GGA' .OR. FUNCTL.EQ.'gga') THEN
     GGA = .TRUE.
  ELSE
     WRITE(6,*) 'atomxc: Unknown functional ', FUNCTL
     STOP
  ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize output                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  EX = 0
  EC = 0
  DX = 0
  DC = 0
  DO 20 IS = 1,NSPIN
     DO 10 IR = 1,NR
        V_XC(IR,IS) = 0
     10   CONTINUE
  20 CONTINUE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get number pi                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PI = 4 * ATAN(1.0_r8)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loop on mesh points                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO 140 IR = 1,NR

!       Find interval of neighbour points to calculate derivatives
     IN1 = MAX(  1, IR-NN ) - IR
     IN2 = MIN( NR, IR+NN ) - IR

!    Find weights of numerical derivation from Lagrange
!    interpolation formula
     DO 50 IN = IN1,IN2
        IF (IN .EQ. 0) THEN
           DGDM(IN) = 0
           DO 30 JN = IN1,IN2
              IF (JN.NE.0) DGDM(IN) = DGDM(IN) + 1.0_r8 / (0 - JN)
           30 CONTINUE
        ELSE
           F1 = 1
           F2 = 1
           DO 40 JN = IN1,IN2
              IF (JN.NE.IN .AND. JN.NE.0) F1 = F1 * (0  - JN)
              IF (JN.NE.IN)               F2 = F2 * (IN - JN)
           40 CONTINUE
           DGDM(IN) = F1 / F2
        ENDIF
     50 CONTINUE

!       Find dr/dmesh
     DRDM = 0
     DO 60 IN = IN1,IN2
        DRDM = DRDM + RMESH(IR+IN) * DGDM(IN)
     60 CONTINUE

!       Find differential of volume. Use trapezoidal integration rule
     DVOL = 4 * PI * RMESH(IR)**2 * DRDM
!    DVMIN is a small number added to avoid a division by zero
     DVOL = DVOL + DVMIN
     IF (IR.EQ.1 .OR. IR.EQ.NR) DVOL = DVOL / 2
     IF (GGA) AUX(IR) = DVOL

!       Find the weights for the derivative d(gradF(i))/d(F(j)), of
!       the gradient at point i with respect to the value at point j
     IF (GGA) THEN
        DO 80 IN = IN1,IN2
           DGIDFJ(IN) = DGDM(IN) / DRDM
        80 CONTINUE
     ENDIF

!    Find density and gradient of density at this point
     DO 90 IS = 1,NSPIN
        D(IS) = DENS(IR,IS)
     90 CONTINUE
     IF (GGA) THEN
        DO 110 IS = 1,NSPIN
           GD(1,IS) = 0
           GD(2,IS) = 0
           GD(3,IS) = 0
           DO 100 IN = IN1,IN2
              GD(3,IS) = GD(3,IS) + DGIDFJ(IN) * DENS(IR+IN,IS)
           100 CONTINUE
        110 CONTINUE
     ENDIF

!    Find exchange and correlation energy densities and their 
!    derivatives with respect to density and density gradient
     IF (GGA) THEN
        CALL GGAXC( AUTHOR, IREL, NSPIN, D, GD,                               &
                    EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
     ELSE
        CALL LDAXC( AUTHOR, IREL, NSPIN, D, EPSX, EPSC, DEXDD, DECDD )
     ENDIF

!       Add contributions to exchange-correlation energy and its
!       derivatives with respect to density at all points
     DO 130 IS = 1,NSPIN
        EX = EX + DVOL * D(IS) * EPSX
        EC = EC + DVOL * D(IS) * EPSC
        DX = DX + DVOL * D(IS) * (EPSX - DEXDD(IS))
        DC = DC + DVOL * D(IS) * (EPSC - DECDD(IS))
        IF (GGA) THEN
            V_XC(IR,IS) = V_XC(IR,IS) + DVOL * ( DEXDD(IS) + DECDD(IS) )
            DO 120 IN = IN1,IN2
               DX= DX - DVOL * DENS(IR+IN,IS) * DEXDGD(3,IS) * DGIDFJ(IN)
               DC= DC - DVOL * DENS(IR+IN,IS) * DECDGD(3,IS) * DGIDFJ(IN)
               V_XC(IR+IN,IS) = V_XC(IR+IN,IS) + DVOL *                         &
                  (DEXDGD(3,IS) + DECDGD(3,IS)) * DGIDFJ(IN)
            120 CONTINUE
        ELSE
            V_XC(IR,IS) = DEXDD(IS) + DECDD(IS)
        ENDIF
     130 CONTINUE
  140 CONTINUE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Divide by volume element to obtain the potential (per electron)             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (GGA) THEN
     DO 160 IS = 1,NSPIN
        DO 150 IR = 1,NR
           DVOL = AUX(IR)
           V_XC(IR,IS) = V_XC(IR,IS) / DVOL
        150 CONTINUE
     160 CONTINUE
  ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Divide by energy unit                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  EX = EX / EUNIT
  EC = EC / EUNIT
  DX = DX / EUNIT
  DC = DC / EUNIT
  DO 180 IS = 1,NSPIN
     DO 170 IR = 1,NR
        V_XC(IR,IS) = V_XC(IR,IS) / EUNIT
     170 CONTINUE
  180 CONTINUE

  end subroutine atomxc



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   VHRTRE CONSTRUCTS THE ELECTROSTATIC POTENTIAL DUE TO A SUPPLIED           !
!   ELECTRON DENSITY.  THE NUMEROV METHOD IS USED TO INTEGRATE                !
!   POISSONS EQN.                                                             !
!                                                                             !
!   DESCRIPTION OF ARGUMENTS:                                                 !
!      RHO....4*PI*R**2 * THE ELECTRON DENSITY FOR WHICH WE CALCULATING       !
!             THE ELECTROSTATIC POTENTIAL                                     !
!      V......THE ELECTROSTATIC POTENTIAL DUE TO THE ELECTRON DENSITY         !
!             RHO.  THE CONSTANTS OF INTEGRATION ARE FIXED SO THAT THE        !
!             POTENTIAL TENDS TO A CONSTANT AT THE ORIGIN AND TO              !
!             2*Q/R AT R=R(NR), WHERE Q IS THE INTEGRATED CHARGE              !
!             CONTAINED IN RHO(R)                                             !
!      R......THE RADIAL MESH R(I) = B*(EXP(A(I-1))-1)                        !
!      NR.....THE NUMBER OF RADIAL MESH POINTS                                !
!      DRDI...DR(I)/DI                                                        !
!     SRDRDI.SQRT(DR/DI)                                                      !
!      A......THE PARAMETER APPEARING IN R(I) = B*(EXP(A(I-1))-1)             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vhrtre(rho, v, r, drdi, srdrdi, nr, a)
  real(r8), intent(IN), dimension(*) :: rho, r, drdi, srdrdi
  real(r8), intent(in) :: a
  real(r8), intent(out), dimension(*) :: v
  integer, intent(in) :: nr

  real(r8) :: x,y,q,a2by4,ybyq,qbyy,qpartc,v0,qt,dz,t,beta,dv
  integer :: nrm1,nrm2,ir

  NRM1 = NR - 1
  NRM2 = NR - 2
  A2BY4 = A*A/4.0_r8
  YBYQ = 1.0_r8 - A*A/48.0_r8
  QBYY = 1.0_r8/YBYQ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SIMPSONS RULE IS USED TO PERFORM TWO INTEGRALS OVER THE ELECTRON          ! 
!  DENSITY.  THE TOTAL CHARGE QT IS USED TO FIX THE POTENTIAL AT R=R(NR)      !
!  AND V0 (THE INTEGRAL OF THE ELECTRON DENSITY DIVIDED BY R) FIXES           !
!  THE ELECTROSTATIC POTENTIAL AT THE ORIGIN                                  !
!  Modified to be consistent with pseudopotential generation (use the         !
!  trapeziodal rule for integration). A. Rubio. Jan. 2000                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  V0 = 0.0_r8
  QT = 0.0_r8
  do IR = 2, NRM1, 2
    DZ = DRDI(IR)*RHO(IR)
    QT = QT + DZ
    V0 = V0 + DZ/R(IR)
  end do

  do IR=3, NRM2, 2
    DZ = DRDI(IR)*RHO(IR)
    QT = QT + DZ
    V0 = V0 + DZ/R(IR)
  end do
  DZ = DRDI(NR)*RHO(NR)

  QT = (QT + QT + DZ)*M_HALF
  V0 = (V0 + V0 + DZ/R(NR))
  V(1) = 2.0_r8*V0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  THE ELECTROSTATIC POTENTIAL AT R=0 IS SET EQUAL TO                         !
!                       THE AVERAGE VALUE OF RHO(R)/R                         !
!  BEGIN CONSTRUCTION OF THE POTENTIAL AT FINITE                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IR   = 2
  T    = SRDRDI(IR)/R(IR)
  BETA = DRDI(IR)*T*RHO(IR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  THE NEXT 4 STATEMENTS INDICATE THAT WE FIRST FIND THE PARTICULAR           !
!  SOLUTION TO THE INHOMOGENEOUS EQN. FOR WHICH Q(2)=0, WE THEN               !
!  ADD TO THIS PARTICULAR SOLUTION A SOLUTION OF THE HOMOGENEOUS EQN.         !
!  (A CONSTANT IN V OR A Q PROPORTIONAL TO R)                                 !
!  WHICH WHEN DIVIDED BY R IN GOING FROM Q TO V GIVES                         !
!  THE POTENTIAL THE DESIRED COULOMB TAIL OUTSIDE THE ELECTRON DENSITY.       !
!  THE SIGNIFICANCE OF THE SOLUTION VANISHING AT THE SECOND RADIAL            !
!  MESH POINT IS THAT, SINCE ALL REGULAR SOLUTIONS OF THE EQUATION            !
!  FOR Q=R*V VANISH AT THE ORIGIN, THE KNOWLEDGE OF THE SOLUTION              !
!  VALUE AT THE SECOND MESH POINT PROVIDES THE TWO SOLUTION VALUES            !
!  REQUIRED TO START THE NUMEROV PROCEDURE.                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  X = 0.0_r8
  Y = 0.0_r8
  Q = (Y - BETA/12.0_r8)*QBYY
  V(IR) = 2.0_r8*T*Q

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  BEGINNING OF THE NUMEROV ALGORITHM                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

3 X  = X + A2BY4*Q - BETA
  Y  = Y + X
  IR = IR + 1
  T  = SRDRDI(IR)/R(IR)
  BETA = T*DRDI(IR)*RHO(IR)
  Q = (Y-BETA/12.0_r8)*QBYY
  V(IR) = 2.0_r8*T*Q
  IF(IR.LT.NR) GO TO 3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  END OF THE NUMEROV ALGORITHM                                               !
!                                                                             !
!  WE HAVE NOW FOUND A PARTICULAR SOLUTION TO THE INHOMOGENEOUS EQN.          !
!  FOR WHICH Q(R) AT THE SECOND RADIAL MESH POINT EQUALS ZERO.                !
!  NOTE THAT ALL REGULAR SOLUTIONS TO THE EQUATION FOR Q=R*V                  !
!  VANISH AT THE ORIGIN.                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  QPARTC = R(NR)*V(NR)/2.0_r8
  DZ = QT - QPARTC
  DV = 2.0_r8*DZ/R(NR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  THE LOOP FOLLOWING ADDS THE CONSTANT SOLUTION OF THE HOMOGENEOUSi          !
!  EQN TO THE PARTICULAR SOLUTION OF THE INHOMOGENEOUS EQN.                   !
!  NOTE THAT V(1) IS CONSTRUCTED INDEPENDENTLY                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do IR = 2, NR
    V(IR) = V(IR) + DV
  end do

  return
  end subroutine vhrtre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  egofv determines the eigenenergy and wavefunction corresponding            !
!  to a particular l, principal quantum number and boundary condition.        !
!                                                                             !
!  two fundamental techniques are used to locate the solution:                ! 
!      1) node counting and bisection                                         !
!       2) variational estimate based on a slope discontinuity in psi         !
!  the arguments are defined as follows:                                      !
!       h,s: g = (h-e*s)*g                                                   !
!       nr: maximum allowed number of radial points                           !
!       e: e(i) is the i-th energy found                                      !
!       ne: number of energies found                                          !
!       l: the angular momentum                                               !
!       ncor: the number of lower-energy state                                !
!                                                                             ! 
!  the individual energies are resolved by performing a fixed number          !
!  of bisections after a given eigenvalue has been isolated                   ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine egofv(h,s,n,e,g,y,l,z,a,b,rmax,nprin,nnode,dr)

  implicit real(r8) (a-h,o-z)
  integer :: i,n,l,nprin,nnode,ncor,n1,n2,niter,nt

  real(r8),dimension(*) ::h,s,g,y

  real(r8), parameter :: tol = 1.0e-5_r8

  ncor=nprin-l-1
  n1=nnode
  n2=nnode-1
  e1=e
  e2=e


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the labels 1 and 2 refer to the bisection process, defining the            !
!  range in which the desired solution is located.  the initial               !
!  settings of n1, n2, e1 and e2 are not consistent with the bisection        !
!  algorithm; they are set to consistent values when the desired              !
!  energy interval has been located.                                          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  del = M_HALF
  de  = M_ZERO
  niter = 0
1 niter = niter + 1
  if(niter.gt.40) go to 3
  et = e + de
! the following line is the fundamental "bisection"
  e = 0.5*(e1+e2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the following concatenation of logical ors ensures that node               !
!  counting is used unless the previous integration of the radial             !
!  eq produced both the correct number of nodes and a sensible                !
!  prediction for the energy.                                                 !
!                                                                             !
!     sensible means that et must be greater than e1 and less than e2         !
!     correct number of nodes means that nt = nnode or nnode-1.               !
!                                                                             !
!     leaving e set to its bisection value, and transfering to                !
!     the call to yofe means that we are performing bisection,                !
!     whereas setting e to et is use of the variational estimate.             !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if(et.le.e1 .or. et.ge.e2 .or.                                              &
     nt.lt.nnode-1 .or. nt.gt.nnode) go to 2
  e=et
  if(dabs(de).lt.tol) go to 6
2 call yofe(e,de,dr,rmax,h,s,y,n,l,ncor,nt,z,a,b)
!     write(6,101) l,dr,n1,nt,nnode,n2,e1,e,e2,de
!101  format('  l     dr     n1  nt   n  n2       e1           e',
!    1       '          e2          de'/i3,d10.3,4i4,4f12.5)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  yofe integrates the schro eq.; now the bisection logic                     !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if(nt.ge.nnode) go to 5
!  too few nodes; set e1 and n1
  e1=e
  n1=nt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  at this point, we have just set the bottom of the bisection range;         !
!  if the top is also set, we procede.  if the top of the range has not       !
!  been set, it means that we have yet to find an e greater than the          !
!  desired energy.  the upper end of the range is extended.                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if(n2.ge.nnode) go to 1
  del=del*M_TWO
  e2=e1+del
  go to 1
!  too many nodes; set e2 and n2
5 e2=e
  n2=nt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  at this point, we have just set the top of the bisection range;            !
!  if the top is also set, we procede.  if the top of the range has           !
!  not been set, it means that we have yet to find an e less than the         !
!  desired energy.  the lower end of the range is extended.                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if(n1.lt.nnode) go to 1
  del=del*M_TWO
  e1=e2-del
  go to 1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the numerov method uses a transformation of the radial wave fcn.           !
!  that we call "y".  having located the eigenenergy, we transform            !
!  y to "g", from which the density is easily constructed.                    !
!  finally, the call to "nrmlzg" normalizes g to one electron.                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


6 g(1) = M_ZERO
  do 7 i=2,n
  t=h(i)-e*s(i)
  g(i)=y(i)/(M_ONE-t/12._r8)
7 continue
  call nrmlzg(g,s,n)
  return
3 write(6,4) z,l,nnode,e,de
4 format(' egofv: too many iterations; execution stopping'/                   &
            ' z=',f3.0,'  l=',i2,'  nnode=',i2,'  e=',f12.5,                  &
            '  de=',f12.5)
  stop 8


  end subroutine egofv


  subroutine yofe(e,de,dr,rmax,h,s,y,nmax,l,ncor,nnode,z,a,b)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   yofe integrates the radial schrodinger eqn using the numerov              !
!   method.                                                                   !
!                                                                             !
!   the arguments are defined as follows:                                     !
!       e is the old energy(overwritten) by the new energy                    !
!      de is the e change predicted to elim the kink in psi                   !
!       dr is the log deriv (the boundary condition)                          !
!       gpp = (h-es)g (all diagonal in i (radius) )                            !
!       y is the numerov independent variable y = g - gpp/12                   !
!       n is the number of radial mesh points                                 !
!       l is the angular momentum                                             !
!       ncor is the number of states of lower energy                          !
!       nnode is 1 + the number of interior nodes in psi                      !
!       z is the atomic number                                                !
!       a and b specify the radial mesh r(i)=(exp(a*(i-1))-1)*b               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit real(r8) (a-h,o-z)
  real(r8) :: h(*),s(*),y(*)
  integer :: nmax,l,ncor,nnode,n,knk,nndin,i

  zdr = z*a*b
  n=nmax
8 if( h(n)-e*s(n) .lt. M_ONE ) go to 9
  y(n)=M_ZERO
  n=n-1
  go to 8
9 continue
  call bcorgn(e,h,s,l,zdr,y2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  bcorgn computes y2, which embodies the boundary condition                  !
!  satisfied by the radial wave function at the origin                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  knk=n
  call numout(e,h,s,y,ncor,knk,nnode,y2,g,gsg,x)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the outward integration is now complete                                    !
!                                                                             !
!     we first decide if the kinetic energy is sufficiently non               !
!     negative to permit use of the numerov eq at rmax.  if                   !
!     it is not, then zero-value boundary condition is used                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  yn = M_ZERO
  if(n.lt.nmax .or. dabs(dr).gt.1.e3_r8) go to 7
  call bcrmax(e,dr,rmax,h,s,n,yn,a,b)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  bcrmax computes yn, which embodies the boundary condition                  !
!  satisfied by the radial wave function at rmax                              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


7 call numin(e,h,s,y,n,nndin,yn,gin,gsgin,xin,knk)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  numin performs the inward integration by the numerov method                !
!                                                                             !
!  the energy increment is now evaluated from the kink in psi                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ratio = g/gin
  xin=xin*ratio
  gsg=gsg+gsgin*ratio*ratio
  t=h(knk)-e*s(knk)
  de=g*(x+xin+t*g)/gsg
  nnode=nnode+nndin
  if(de.lt.M_ZERO) nnode=nnode+1
  do 6 i=knk,n
  y(i) = y(i)*ratio
6 continue
  return


  end subroutine yofe




  subroutine nrmlzg(g,s,n)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   nrmlzg normalizes the supplied radial wave function                       !
!                                                                             !
!   the arguments are defined as follows:                                     !
!       g is the radial wave function appropriate to the numerov              !
!             representation of the radial schrodinger equation               !
!             that is, the radial fcn r(r) = (drdi)**1/2 g(i) / r(i)          !
!       gpp = (h-es)g (all diagonal in i (radius) )                           !
!       n1 is the number of radial mesh points corresponding to               !
!             the portion of the radial mesh on which the norm                !
!             is defined                                                      !
!       n2 is the number of radial mesh points corresponding to               !
!             the portion of the radial mesh on which the wave                !
!             function is defined.  for the intended use of this              !
!             routine, n1 = nrval and n2 = nrcor                              !
!       a and b are the radial mesh parameters                                !
!             r(i) = ( exp(a*(i-1)) - 1 ) * b                                 !
!             (dr/di = a*b at the origin)                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit real(r8) (a-h,o-z)
  real(r8) s(*),g(*),norm
  integer :: n,nm1,nm2,i


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  determine the norm of g(i) using simpsons rule                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (mod(n,2).ne.1) write(6,*) ' nrmlzg: n should be odd. n =',n
  norm = M_ZERO
  nm1 = n - 1
  do 2 i = 2,nm1,2
  norm=norm + g(i)*s(i)*g(i)
2 continue
  norm = norm * M_TWO
  nm2  = n - 2
  do 3 i = 3,nm2,2
  norm=norm + g(i)*s(i)*g(i)
3 continue
  norm = norm * M_TWO
  nm2  = n - 2
  do 4 i = 1,n,nm1
  norm=norm + g(i)*s(i)*g(i)
4 continue
  norm = norm/M_THREE
  srnrm = dsqrt(norm)
  do 5 i=1,n
  g(i) = g(i)/srnrm
5 continue
  return


  end subroutine nrmlzg


  subroutine bcorgn(e,h,s,l,zdr,y2)

  implicit real(r8) (a-h,o-z)
  real(r8) h(*),s(*)
  integer :: l


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   the quantity called d(i) in the program is actually the inverse           !
!   of the diagonal of the tri-diagonal numerov matrix                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  t2=h(2)-e*s(2)
  d2=-((24._r8+10._r8*t2)/(12._r8-t2))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  the following section deals with the fact that the independent             !
!  variable "y" in the numerov equation is not zero at the origin             !
!  for l less than 2                                                          !
!  the l=0 solution g vanishes, but the first and second                      !
!  derivatives are finite, making the numerov variable y finite               !
!  the l=1 solution g vanishes, and gprime also vanishes, but                 !
!  the second derivative gpp is finite making y finite.  for l > 1,           !
!  g and its first two derivatives vanish, making y zero.                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if(l.ge.2) goto 3
  if(l.gt.0) goto 1
  c0=zdr/6._r8
  c0=c0/(M_ONE-0.75*zdr)
  go to 2
1 c0=M_ONE/12._r8
  c0=(-c0)*8._r8/M_THREE
2 c1=c0*(12._r8+13._r8*t2)/(12._r8-t2)
  t3=h(3)-e*s(3)
  c2=(-M_HALF)*c0*(24._r8-t3)/(12._r8-t3)
  d2=(d2-c1)/(M_ONE-c2)
3 y2=(-M_ONE)/d2
  return

  end subroutine bcorgn



  subroutine bcrmax(e,dr,rmax,h,s,n,yn,a,b)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 22.7.85                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit real(r8) (a-h,o-z)
  real(r8)  h(*),s(*),                                                          &
   e,dr,rmax,yn,a,b,tnm1,tn,tnp1,beta,dg,c1,c2,c3,dn
  integer :: n

!     write(6,*) 'bcrmax:',dr
  tnm1=h(n-1)-e*s(n-1)
  tn  =h(n  )-e*s(n  )
  tnp1=h(n+1)-e*s(n+1)
  beta=M_ONE+b/rmax
  dg=a*beta*(dr+M_ONE-M_HALF/beta)


  c2=24._r8*dg/(12._r8-tn)
  dn=-((24._r8+10._r8*tn)/(12._r8-tn))

  c1= (M_ONE-tnm1/6._r8)/(M_ONE-tnm1/12._r8)
  c3=-((M_ONE-tnp1/6._r8)/(M_ONE-tnp1/12._r8))
  yn=-((M_ONE-c1/c3)/(dn-c2/c3))


  return

  end subroutine bcrmax



  subroutine numin(e,h,s,y,n,nnode,yn,g,gsg,x,knk)

  implicit real(r8) (a-h,o-z)
  integer :: i,n,nnode,knk
  real(r8) :: h(n),s(n),y(n)

  y(n)=yn
  t=h(n)-e*s(n)
  g=y(n)/(M_ONE-t/12._r8)
  gsg=g*s(n)*g
  i=n-1
  y(i)=M_ONE
  t=h(i)-e*s(i)
  g=y(i)/(M_ONE-t/12._r8)
  gsg=gsg+g*s(i)*g
  x=y(i)-y(n)
  nnode=0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  begin the inward integrationby the numerov method                          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


1 x=x+t*g
  i=i-1
  y(i)=y(i+1)+x
  if( y(i)*y(i+1) .lt. M_ZERO) nnode=nnode+1
  t=h(i)-e*s(i)
  g=y(i)/(M_ONE-t/12._r8)
  gsg=gsg+g*s(i)*g
  if(i.gt.knk) go to 1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the last statement defines the kink radius as the point where              !
!  psi first turns downward.  this usually means at the outermost             !
!  maximum                                                                    !
!                                                                             !
!  the inward integration is now complete                                     !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  knk=i
  return

  end subroutine numin


  subroutine numout(e,h,s,y,ncor,knk,nnode,y2,g,gsg,x)

  implicit real(r8) (a-h,o-z)
  integer :: ncor,nnode,knk,i,nm4
  real(r8) :: h(knk),s(knk),y(knk)

  y(1)=M_ZERO
  y(2)=y2
  t=h(2)-e*s(2)
  g=y(2)/(M_ONE-t/12._r8)
  gsg=g*s(2)*g
  y(3)=M_ONE
  t=h(3)-e*s(3)
  g=y(3)/(M_ONE-t/12._r8)
  gsg=gsg+g*s(3)*g
  x=y(3)-y(2)
  i=3
  nnode=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  begin the outward integrationby the numerov method                         !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nm4=knk-4
1 xl=x
  x=x+t*g
  i=i+1
  y(i)=y(i-1)+x
!     write(6,300) i,y(i),x,t,h(i),s(i)
!300  format(i5,5d14.5)
  if( y(i)*y(i-1) .lt. M_ZERO) nnode=nnode+1
  t=h(i)-e*s(i)
  g=y(i)/(M_ONE-t/12._r8)
  gsg=gsg+g*s(i)*g
  if(i.eq.nm4) go to 2
  if(nnode.lt.ncor) go to 1
  if(xl*x.gt.M_ZERO) go to 1
2 knk=i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the outward integration is now complete                                    !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  return
  end subroutine numout

  end module atomic

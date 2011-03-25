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
#define R_TOPREC(x) real(x, REAL_PRECISION)

module atomic_m
  use global_m
  use logrid_m
  use messages_m
  use periodic_table_m
  use profiling_m
  use xc_f90_lib_m ! this is the double-precision version

  implicit none

  private
  public ::        &
    atomxc,        &
    atomhxc,       &
    vhrtre, egofv

  ! Following stuff is for the "val_conf_t" data type, made to handle atomic configurations.
  public ::        &
    valconf_t,     &
    valconf_copy,  &
    write_valconf, &
    read_valconf,  &
    valconf_null,  &
    VALCONF_STRING_LENGTH

  character(len=1), parameter :: &
    spec_notation(0:3) = (/ 's', 'p', 'd', 'f' /)

  integer, parameter :: VALCONF_STRING_LENGTH = 80

  type valconf_t
    integer           :: z
    character(len=3)  :: symbol
    integer           :: type     ! 0 for the most normal valence configuration, 1 for semicore.
    integer           :: p        ! number of orbitals.
    integer           :: n(6)     ! n quantum number
    integer           :: l(6)     ! l quantum number
    FLOAT             :: occ(6,2) ! occupations of each level
  end type valconf_t


contains

  ! ---------------------------------------------------------
  ! Subroutines to write and read valence configurations.
  subroutine valconf_null(c)
    type(valconf_t), intent(out) :: c

    PUSH_SUB(valconf_null)
    c%z = 0
    c%symbol = ""
    c%type = 0
    c%p = 0
    c%n = 0
    c%l = 0
    c%occ = M_ZERO

    POP_SUB(valconf_null)
  end subroutine valconf_null


  ! ---------------------------------------------------------
  subroutine valconf_copy(cout, cin)
    type(valconf_t), intent(out) :: cout
    type(valconf_t), intent(in)  :: cin

    PUSH_SUB(valconf_copy)

    cout%z      = cin%z
    cout%symbol = cin%symbol
    cout%type   = cin%type
    cout%p      = cin%p
    cout%n      = cin%n
    cout%l      = cin%l
    cout%occ    = cin%occ

    POP_SUB(valconf_copy)
  end subroutine valconf_copy


  ! ---------------------------------------------------------
  subroutine write_valconf(c, s)
    type(valconf_t), intent(in) :: c
    character(len=VALCONF_STRING_LENGTH), intent(out) :: s

    integer :: j

    PUSH_SUB(write_valconf)

    write(s,'(i2,1x,a2,i1,1x,i1,a1,6(i1,a1,f6.3,a1))') c%z, c%symbol, c%type, c%p, ':',&
         (c%n(j),spec_notation(c%l(j)),c%occ(j,1),',',j=1,c%p)

    POP_SUB(write_valconf)
  end subroutine write_valconf


  ! ---------------------------------------------------------
  subroutine read_valconf(s, c)
    character(len=VALCONF_STRING_LENGTH), intent(in) :: s
    type(valconf_t), intent(out) :: c

    integer :: j
    character(len=1) :: lvalues(1:6)

    PUSH_SUB(read_valconf)

    read (s,'(i2,1x,a2,i1,1x,i1,1x,6(i1,a1,f6.3,1x))') c%z, c%symbol, c%type, c%p,&
         (c%n(j),lvalues(j),c%occ(j,1),j=1,c%p)
    do j = 1, c%p
       select case(lvalues(j))
       case('s'); c%l(j) = 0
       case('p'); c%l(j) = 1
       case('d'); c%l(j) = 2
       case('f'); c%l(j) = 3
       case default
          message(1) = 'Error: read_valconf.'
          call messages_fatal(1)
       end select
    end do

    POP_SUB(read_valconf)
  end subroutine read_valconf


  ! ---------------------------------------------------------
  subroutine atomhxc(functl, g, nspin, dens, v, extra)
    character(len=*),  intent(in)  :: functl
    type(logrid_t),    intent(in)  :: g
    integer,           intent(in)  :: nspin
    FLOAT,             intent(in)  :: dens(g%nrval, nspin)
    FLOAT,             intent(out) :: v(g%nrval, nspin)
    FLOAT,             intent(in), optional :: extra(g%nrval)

    character(len=5) :: xcfunc, xcauth
    integer :: is, ir
    FLOAT, allocatable :: xc(:, :), ve(:, :), rho(:, :)
    FLOAT :: r2, ex, ec, dx, dc

    PUSH_SUB(atomhxc)

    SAFE_ALLOCATE( ve(1:g%nrval, 1:nspin))
    SAFE_ALLOCATE( xc(1:g%nrval, 1:nspin))
    SAFE_ALLOCATE(rho(1:g%nrval, 1:nspin))
    ve = M_ZERO
    xc = M_ZERO
    rho = M_ZERO

    ! To calculate the Hartree term, we put all the density in one variable.
    do is = 1, nspin
      rho(:, 1) = rho(:, 1) + dens(:, is)
    end do
    ve = M_ZERO
    do is = 1, nspin
       call vhrtre(rho(:, 1), ve(:, is), g%rofi, g%drdi, g%s, g%nrval, g%a)
    end do
    ve(1, 1:nspin) = ve(2, 1:nspin)

    select case(functl)
    case('LDA')
      xcfunc = 'LDA'
      xcauth = 'PZ'
    case('GGA')
      xcfunc = 'GGA'
      xcauth = 'PBE'
    case default
      message(1) = 'Internal Error in atomhxc: unknown functl'
      call messages_fatal(1)
    end select

    rho = dens
    do is = 1, nspin
      if(present(extra)) rho(:, is) = rho(:, is) + extra(:)/nspin
      rho(2:g%nrval, is) = rho(2:g%nrval, is)/(M_FOUR*M_PI*g%rofi(2:g%nrval)**2)
    end do
    r2 = g%rofi(2)/(g%rofi(3)-g%rofi(2))
    rho(1, 1:nspin) = rho(2, 1:nspin) - (rho(3, 1:nspin)-rho(2, 1:nspin))*r2

    do is = 1, nspin
      do ir = 1, g%nrval
        if(rho(ir, is) < M_EPSILON) rho(ir, is) = M_ZERO
      end do
    end do

    call atomxc(xcfunc, xcauth, g%nrval, g%nrval, g%rofi, nspin, rho, ex, ec, dx, dc, xc)
    v = ve + xc

    SAFE_DEALLOCATE_A(ve)
    SAFE_DEALLOCATE_A(xc)
    SAFE_DEALLOCATE_A(rho)

    POP_SUB(atomhxc)
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
  subroutine atomxc( FUNCTL, AUTHOR, NR, MAXR, RMESH, NSPIN, DENS, EX, EC, DX, DC, V_XC )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Argument types and dimensions                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=*), intent(in) :: FUNCTL, AUTHOR
    integer, intent(in) :: NR, MAXR
    FLOAT, intent(in) :: RMESH(MAXR)
    integer, intent(in) :: NSPIN
    FLOAT, intent(in) :: DENS(MAXR,NSPIN)
    FLOAT, intent(out) :: DC, DX, EC, EX, V_XC(MAXR,NSPIN)

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

    FLOAT, parameter :: EUNIT = M_HALF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DVMIN is added to differential of volume to avoid division by zero          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    FLOAT, parameter :: DVMIN = CNST(1.0E-12)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Local variables and arrays                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    logical :: GGA
    integer :: IN, IN1, IN2, IR, IS, JN
    REAL_DOUBLE :: AUX(MAXR), D(NSPIN),                              &
      DGDM(-NN:NN), DGIDFJ(-NN:NN), DRDM, DVOL,                      &
      F1, F2, GD(3,NSPIN), sigma(3), vxsigma(3), vcsigma(3),         &
      EPSX, EPSC,                                                    &
      DEXDD(NSPIN), DEXDGD(3,NSPIN), DECDD(NSPIN), DECDGD(3,NSPIN)
    type(xc_f90_pointer_t) :: x_conf, c_conf
    type(xc_f90_pointer_t) :: x_info, c_info

    PUSH_SUB(atomxc)

    ! sanity check
    ASSERT(NSPIN==1.or.NSPIN==2)

    ! Set GGA switch
    IF ( FUNCTL.EQ.'LDA' .OR. FUNCTL.EQ.'lda' .OR.                              &
      FUNCTL.EQ.'LSD' .OR. FUNCTL.EQ.'lsd' ) THEN
      GGA = .FALSE.
    ELSEIF ( FUNCTL.EQ.'GGA' .OR. FUNCTL.EQ.'gga') THEN
      GGA = .TRUE.
    ELSE
      GGA = .FALSE.
      write(message(1),'(a,a)') 'Error: atomxc: Unknown functional ', FUNCTL
      call messages_fatal(1)
    ENDIF
    
    ! initialize xc functional
    if(GGA) then
      call xc_f90_func_init(x_conf, x_info, XC_GGA_X_PBE, NSPIN)
      call xc_f90_func_init(c_conf, c_info, XC_GGA_C_PBE, NSPIN)
    else
      call xc_f90_func_init(x_conf, x_info, XC_LDA_X, NSPIN)
      if(AUTHOR.EQ.'CA' .OR. AUTHOR.EQ.'ca' .OR.                                &
        AUTHOR.EQ.'PZ' .OR. AUTHOR.EQ.'pz') THEN
        call xc_f90_func_init(c_conf, c_info, XC_LDA_C_PZ, NSPIN)
      else IF ( AUTHOR.EQ.'PW92' .OR. AUTHOR.EQ.'pw92' ) THEN
        call xc_f90_func_init(c_conf, c_info, XC_LDA_C_PW, NSPIN)
      else
        write(message(1),'(a,a)') 'Error: LDAXC: Unknown author ', AUTHOR
        call messages_fatal(1)
      end if
    end if
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize output                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    EX = M_ZERO
    EC = M_ZERO
    DX = M_ZERO
    DC = M_ZERO
    do IS = 1,NSPIN
      do IR = 1,NR
        V_XC(IR,IS) = M_ZERO
      end do
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loop on mesh points                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do IR = 1,NR
    
      ! Find interval of neighbour points to calculate derivatives
      IN1 = MAX(  1, IR-NN ) - IR
      IN2 = MIN( NR, IR+NN ) - IR
      
      ! Find weights of numerical derivative from Lagrange
      ! interpolation formula
      do IN = IN1,IN2
        if (IN .eq. 0) then
          DGDM(IN) = 0
          do JN = IN1,IN2
            if (JN.ne.0) DGDM(IN) = DGDM(IN) + M_ONE / (0 - JN)
          end do
        else
          F1 = 1
          F2 = 1
          do JN = IN1,IN2
            if (JN.ne.IN .and. JN.ne.0) F1 = F1 * (0  - JN)
            if (JN.ne.IN)               F2 = F2 * (IN - JN)
          end do
          DGDM(IN) = F1 / F2
        end if
      end do
      
      ! Find dr/dmesh
      DRDM = 0
      do IN = IN1,IN2
        DRDM = DRDM + RMESH(IR+IN) * DGDM(IN)
      end do
      
      ! Find differential of volume. Use trapezoidal integration rule
      DVOL = 4 * M_PI * RMESH(IR)**2 * DRDM
      
      ! DVMIN is a small number added to avoid a division by zero
      DVOL = DVOL + DVMIN
      if (IR.eq.1 .or. IR.eq.NR) DVOL = DVOL / 2
      if (GGA) AUX(IR) = DVOL
      
      ! Find the weights for the derivative d(gradF(i))/d(F(j)), of
      ! the gradient at point i with respect to the value at point j
      if (GGA) then
        do IN = IN1,IN2
          DGIDFJ(IN) = DGDM(IN) / DRDM
        end do
      end if
      
      ! Find density and gradient of density at this point
      do IS = 1,NSPIN
        D(IS) = DENS(IR,IS)
      end do
      if (GGA) then
        do IS = 1,NSPIN
          GD(:,IS) = M_ZERO
          do IN = IN1,IN2
            GD(3,IS) = GD(3,IS) + DGIDFJ(IN) * DENS(IR+IN,IS)
          end do
        end do
      else
        GD = M_ZERO
      end if
      
      ! The xc_f90_gga routines need as input these combinations of 
      ! the gradient of the densities.
      sigma(1) = gd(3, 1)*gd(3, 1)
      if(NSPIN > 1) then
        sigma(2) = gd(3, 1) * gd(3, 2)
        sigma(3) = gd(3, 2) * gd(3, 2)
      end if
      
      ! Find exchange and correlation energy densities and their
      ! derivatives with respect to density and density gradient
      if (GGA) then
        call xc_f90_gga_exc_vxc(x_conf, 1, D(1), sigma(1), EPSX, DEXDD(1), vxsigma(1))
        call xc_f90_gga_exc_vxc(c_conf, 1, D(1), sigma(1), EPSC, DECDD(1), vcsigma(1))
      else
        call xc_f90_lda_exc_vxc(x_conf, 1, D(1), EPSX, DEXDD(1))
        call xc_f90_lda_exc_vxc(c_conf, 1, D(1), EPSC, DECDD(1))
      end if
      
      ! The derivatives of the exchange and correlation energies per particle with
      ! respect to the density gradient have to be recovered from the derivatives
      ! with respect to the sigma values defined above.
      if(gga) then
        dexdgd(:, 1) = M_TWO * vxsigma(1)*gd(:, 1)
        decdgd(:, 1) = M_TWO * vcsigma(1)*gd(:, 1)
        if(NSPIN .eq. 2) then
          dexdgd(:, 1) = dexdgd(:, 1) + vxsigma(2)*gd(:, 2)
          dexdgd(:, 2) = M_TWO * vxsigma(3)*gd(:, 2) + vxsigma(2)*gd(:, 1)
          decdgd(:, 1) = decdgd(:, 1) + vcsigma(2)*gd(:, 2)
          decdgd(:, 2) = M_TWO * vcsigma(3)*gd(:, 2) + vcsigma(2)*gd(:, 1)
        end if
      end if
      
      ! Add contributions to exchange-correlation energy and its
      ! derivatives with respect to density at all points
      do is = 1, NSPIN
        EX = EX + DVOL * D(IS) * EPSX
        EC = EC + DVOL * D(IS) * EPSC
        DX = DX + DVOL * D(IS) * (EPSX - DEXDD(IS))
        DC = DC + DVOL * D(IS) * (EPSC - DECDD(IS))
        if (GGA) then
          V_XC(IR,IS) = V_XC(IR,IS) + DVOL * ( DEXDD(IS) + DECDD(IS) )
          do IN = IN1,IN2
            DX = DX - DVOL * DENS(IR+IN,IS) * DEXDGD(3,IS) * DGIDFJ(IN)
            DC = DC - DVOL * DENS(IR+IN,IS) * DECDGD(3,IS) * DGIDFJ(IN)
            V_XC(IR+IN,IS) = V_XC(IR+IN,IS) + DVOL *                         &
              (DEXDGD(3,IS) + DECDGD(3,IS)) * DGIDFJ(IN)
          end do
        else
          V_XC(IR,IS) = DEXDD(IS) + DECDD(IS)
        end if
      end do
      
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Divide by volume element to obtain the potential (per electron)             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if(GGA) then
      do IS = 1,NSPIN
        do IR = 1,NR
          DVOL = AUX(IR)
          V_XC(IR,IS) = V_XC(IR,IS) / DVOL
        end do
      end do
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Divide by energy unit                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    EX = EX / EUNIT
    EC = EC / EUNIT
    DX = DX / EUNIT
    DC = DC / EUNIT
    do IS = 1,NSPIN
      do IR = 1,NR
        V_XC(IR,IS) = V_XC(IR,IS) / EUNIT
      end do
    end do

    call xc_f90_func_end(x_conf)
    call xc_f90_func_end(c_conf)
    
    POP_SUB(atomxc)
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
    FLOAT,   intent(in)  :: rho(:)
    FLOAT,   intent(out) :: v(:)
    FLOAT,   intent(in)  :: r(:), drdi(:), srdrdi(:)
    integer, intent(in)  :: nr
    FLOAT,   intent(in)  :: a
    
    FLOAT :: x,y,q,a2by4,ybyq,qbyy,qpartc,v0,qt,dz,t,beta,dv
    integer :: nrm1,nrm2,ir
    
    PUSH_SUB(vhrtre)
    
    NRM1 = NR - 1
    NRM2 = NR - 2
    A2BY4 = A*A/M_FOUR
    YBYQ = M_ONE - A*A/CNST(48.0)
    QBYY = M_ONE/YBYQ
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SIMPSON`S RULE IS USED TO PERFORM TWO INTEGRALS OVER THE ELECTRON          !
!  DENSITY.  THE TOTAL CHARGE QT IS USED TO FIX THE POTENTIAL AT R=R(NR)      !
!  AND V0 (THE INTEGRAL OF THE ELECTRON DENSITY DIVIDED BY R) FIXES           !
!  THE ELECTROSTATIC POTENTIAL AT THE ORIGIN                                  !
!  Modified to be consistent with pseudopotential generation (use the         !
!  trapezoidal rule for integration). A. Rubio. Jan. 2000                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    V0 = M_ZERO
    QT = M_ZERO
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
    V(1) = M_TWO*V0

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

    X = M_ZERO
    Y = M_ZERO
    Q = (Y - BETA/CNST(12.0))*QBYY
    V(IR) = M_TWO*T*Q

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  BEGINNING OF THE NUMEROV ALGORITHM                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do IR = 2, NR
      X  = X + A2BY4*Q - BETA
      Y  = Y + X
      T  = SRDRDI(IR)/R(IR)
      BETA = T*DRDI(IR)*RHO(IR)
      Q = (Y-BETA/CNST(12.0))*QBYY
      V(IR) = M_TWO*T*Q
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  END OF THE NUMEROV ALGORITHM                                               !
!                                                                             !
!  WE HAVE NOW FOUND A PARTICULAR SOLUTION TO THE INHOMOGENEOUS EQN.          !
!  FOR WHICH Q(R) AT THE SECOND RADIAL MESH POINT EQUALS ZERO.                !
!  NOTE THAT ALL REGULAR SOLUTIONS TO THE EQUATION FOR Q=R*V                  !
!  VANISH AT THE ORIGIN.                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    QPARTC = R(NR)*V(NR)/M_TWO
    DZ = QT - QPARTC
    DV = M_TWO*DZ/R(NR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  THE LOOP FOLLOWING ADDS THE CONSTANT SOLUTION OF THE HOMOGENEOUS           !
!  EQN TO THE PARTICULAR SOLUTION OF THE INHOMOGENEOUS EQN.                   !
!  NOTE THAT V(1) IS CONSTRUCTED INDEPENDENTLY                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do IR = 2, NR
      V(IR) = V(IR) + DV
    end do

    POP_SUB(vhrtre)
  end subroutine vhrtre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FROM NOW, ON, ALL SUBROUTINES ARE IN DOUBLE PRECISION, NO MATTER IF THE CODE
! IS COMPILED WITH SINGLE PRECISION OR NOT. APPARENTLY DOUBLE PRECISION IS
! NEEDED FOR THESE PROCEDURES.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  egofv determines the eigenenergy and wavefunction corresponding            !
!  to a particular l, principal quantum number and boundary condition.        !
!                                                                             !
!  two fundamental techniques are used to locate the solution:                !
!      1) node counting and bisection                                         !
!      2) variational estimate based on a slope discontinuity in psi          !
!  the arguments are defined as follows:                                      !
!       h,s: g = (h-e*s)*g                                                    !
!       nr: maximum allowed number of radial points                           !
!       e: e is the energy found                                              !
!       ne: number of energies found                                          !
!       l: the angular momentum                                               !
!       ncor: the number of lower-energy state                                !
!                                                                             !
!  the individual energies are resolved by performing a fixed number          !
!  of bisections after a given eigenvalue has been isolated                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine egofv(h,s,n,e,g,y,l,z,a,b,rmax,nprin,nnode,dr,ierr)
    REAL_DOUBLE, dimension(n), intent(in)    :: h, s
    integer,                   intent(in)    :: n
    REAL_DOUBLE,               intent(inout) :: e
    REAL_DOUBLE, dimension(n), intent(out)   :: g, y
    integer,                   intent(in)    :: l
    REAL_DOUBLE,               intent(in)    :: z, a, b, rmax
    integer,                   intent(in)    :: nprin, nnode
    REAL_DOUBLE,               intent(in)    :: dr
    integer,                   intent(out)   :: ierr
    
    integer :: i,ncor,n1,n2,niter,nt
    REAL_DOUBLE               :: e1, e2, del, de, et, t
    REAL_DOUBLE, parameter    :: tol = CNST(1.0e-5)
    
    PUSH_SUB(egofv)
    
    ncor=nprin-l-1
    n1=nnode
    n2=nnode-1
    e1=e
    e2=e
    ierr = 1
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the labels 1 and 2 refer to the bisection process, defining the            !
!  range in which the desired solution is located.  the initial               !
!  settings of n1, n2, e1 and e2 are not consistent with the bisection        !
!  algorithm; they are set to consistent values when the desired              !
!  energy interval has been located.                                          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nt = 0
    del = M_HALF
    de  = M_ZERO

    do niter = 0, 40

      et = e + de
      ! the following line is the fundamental "bisection"
      e = M_HALF*(e1+e2)

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


      if(.not. (et.le.e1 .or. et.ge.e2 .or. nt.lt.nnode-1 .or. nt.gt.nnode)) then
        e=et
        if(abs(de).lt.tol) then
          ierr = 0
          exit
        endif
      endif
      call yofe(e,de,dr,rmax,h,s,y,n,l,ncor,nt,z,a,b)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  yofe integrates the schro eq.; now the bisection logic                     !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(nt .lt. nnode) then
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

        if(n2.ge.nnode) cycle
        del=del*M_TWO
        e2=e1+del
        cycle
      endif

      !  too many nodes; set e2 and n2
      e2=e
      n2=nt
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  at this point, we have just set the top of the bisection range;            !
!  if the top is also set, we procede.  if the top of the range has           !
!  not been set, it means that we have yet to find an e less than the         !
!  desired energy.  the lower end of the range is extended.                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(n1 .lt. nnode) cycle
      del=del*M_TWO
      e1=e2-del

    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the Numerov method uses a transformation of the radial wavefunction.       !
!  that we call "y".  having located the eigenenergy, we transform            !
!  y to "g", from which the density is easily constructed.                    !
!  finally, the call to "nrmlzg" normalizes g to one electron.                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(ierr .ne. 1) then
      g(1) = M_ZERO
      do i = 2, n
        t=h(i)-e*s(i)
        g(i)=y(i)/(M_ONE-t/CNST(12.))
      enddo
      call nrmlzg(g,s,n)
    endif
    
    POP_SUB(egofv)
  end subroutine egofv


  subroutine yofe(e,de,dr,rmax,h,s,y,nmax,l,ncor,nnode,z,a,b)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   yofe integrates the radial Schroedinger eqn using the Numerov             !
!   method.                                                                   !
!                                                                             !
!   the arguments are defined as follows:                                     !
!       e is the old energy (overwritten) by the new energy                   !
!       de is the e change predicted to elim the kink in psi                  !
!       dr is the log deriv (the boundary condition)                          !
!       gpp = (h-es)g (all diagonal in i (radius) )                           !
!       y is the numerov independent variable y = g - gpp/12                  !
!       n is the number of radial mesh points                                 !
!       l is the angular momentum                                             !
!       ncor is the number of states of lower energy                          !
!       nnode is 1 + the number of interior nodes in psi                      !
!       z is the atomic number                                                !
!       a and b specify the radial mesh r(i)=(exp(a*(i-1))-1)*b               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    REAL_DOUBLE, intent(inout) :: e
    REAL_DOUBLE, intent(out)   :: de
    REAL_DOUBLE, intent(in)    :: dr, rmax
    REAL_DOUBLE, intent(in)    :: h(nmax), s(nmax)
    REAL_DOUBLE, intent(out)   :: y(nmax)
    integer,     intent(in)    :: nmax, l
    integer,     intent(in)    :: ncor
    integer,     intent(out)   :: nnode
    REAL_DOUBLE, intent(in)    :: z, a, b

    REAL_DOUBLE :: y2, g, gsg, x, gin, gsgin, xin
    integer :: n,knk,nndin,i
    REAL_DOUBLE :: zdr, yn, ratio, t
    
    ! No PUSH SUB, called too often.
    
    zdr = z*a*b
    
    do n = nmax, 1, -1
      if( h(n)-e*s(n) .lt. M_ONE ) then
        knk = n
        exit
      endif
      y(n)=M_ZERO
    enddo
    
    call bcorgn(e,h,s,l,zdr,y2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  bcorgn computes y2, which embodies the boundary condition                  !
!  satisfied by the radial wavefunction at the origin                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    if(.not. (n.lt.nmax .or. abs(dr).gt.CNST(1.e3))) then
      call bcrmax(e,dr,rmax,h,s,n,yn,a,b)
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  bcrmax computes yn, which embodies the boundary condition                  !
!  satisfied by the radial wavefunction at rmax                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call numin(e,h,s,y,n,nndin,yn,gin,gsgin,xin,knk)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  numin performs the inward integration by the Numerov method                !
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
    
    do i=knk,n
      y(i) = y(i)*ratio
    enddo

  end subroutine yofe


  subroutine nrmlzg(g,s,n)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   nrmlzg normalizes the supplied radial wavefunction                        !
!                                                                             !
!   the arguments are defined as follows:                                     !
!       g is the radial wavefunction appropriate to the Numerov               !
!             representation of the radial Schroedinger equation              !
!             that is, the radial fcn r(r) = (drdi)**1/2 g(i) / r(i)          !
!       gpp = (h-es)g (all diagonal in i (radius) )                           !
!       n1 is the number of radial mesh points corresponding to               !
!             the portion of the radial mesh on which the norm                !
!             is defined                                                      !
!       n2 is the number of radial mesh points corresponding to               !
!             the portion of the radial mesh on which the wavefunction        !
!             is defined. For the intended use of this                        !
!             routine, n1 = nrval and n2 = nrcor                              !
!       a and b are the radial mesh parameters                                !
!             r(i) = ( exp(a*(i-1)) - 1 ) * b                                 !
!             (dr/di = a*b at the origin)                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    REAL_DOUBLE, intent(inout) :: g(:)
    REAL_DOUBLE, intent(in)    :: s(:)
    integer,     intent(in)    :: n

    integer :: nm1,nm2,i
    REAL_DOUBLE :: norm, srnrm
    
    PUSH_SUB(nrmlzg)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  determine the norm of g(i) using Simpson`s rule                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (mod(n,2).ne.1) then
      write(message(1),'(a,i6)') ' nrmlzg: n should be odd. n =', n
      call messages_warning(1)
    endif

    norm = M_ZERO
    nm1 = n - 1
    do i = 2,nm1,2
      norm=norm + g(i)*s(i)*g(i)
    enddo
    norm = norm * M_TWO
    nm2  = n - 2
    do i = 3,nm2,2
      norm=norm + g(i)*s(i)*g(i)
    enddo
    norm = norm * M_TWO
    nm2  = n - 2
    do i = 1,n,nm1
      norm=norm + g(i)*s(i)*g(i)
    enddo
    norm = norm/M_THREE
    srnrm = sqrt(norm)
    do i=1,n
      g(i) = g(i)/srnrm
    enddo
    
    POP_SUB(nrmlzg)
  end subroutine nrmlzg


  subroutine bcorgn(e,h,s,l,zdr,y2)
    REAL_DOUBLE, intent(in) :: e
    REAL_DOUBLE, intent(in) :: h(:), s(:)
    integer,     intent(in) :: l
    REAL_DOUBLE, intent(in) :: zdr
    REAL_DOUBLE, intent(out) :: y2

    REAL_DOUBLE :: t2, t3, d2, c0, c1, c2
    
    ! no PUSH SUB, called too often

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   the quantity called d(i) in the program is actually the inverse           !
!   of the diagonal of the tri-diagonal Numerov matrix                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    t2=h(2)-e*s(2)
    d2=-((CNST(24.)+M_TEN*t2)/(CNST(12.)-t2))
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  The following section deals with the fact that the independent             !
!  variable "y" in the Numerov equation is not zero at the origin             !
!  for l less than 2.                                                         !
!  The l=0 solution g vanishes, but the first and second                      !
!  derivatives are finite, making the Numerov variable y finite.              !
!  The l=1 solution g vanishes, and gprime also vanishes, but                 !
!  The second derivative gpp is finite making y finite. For l > 1,            !
!  g and its first two derivatives vanish, making y zero.                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(l .lt. 2) then
      if(l .le. 0) then
        c0=zdr/M_SIX
        c0=c0/(M_ONE-CNST(0.75)*zdr)
      else
        c0=M_ONE/CNST(12.)
        c0=(-c0)*M_EIGHT/M_THREE
      endif

      c1=c0*(CNST(12.)+CNST(13.)*t2)/(CNST(12.)-t2)
      t3=h(3)-e*s(3)
      c2=(-M_HALF)*c0*(CNST(24.)-t3)/(CNST(12.)-t3)
      d2=(d2-c1)/(M_ONE-c2)
    endif
    y2=(-M_ONE)/d2
    
    if(l .lt. 2) then
      if(l .le. 0) then
        c0=zdr/M_SIX
        c0=c0/(M_ONE-CNST(0.75)*zdr)
      else
        c0=M_ONE/CNST(12.)
        c0=(-c0)*M_EIGHT/M_THREE
      endif
      c1=c0*(CNST(12.)+CNST(13.)*t2)/(CNST(12.)-t2)
      t3=h(3)-e*s(3)
      c2=(-M_HALF)*c0*(CNST(24.)-t3)/(CNST(12.)-t3)
      d2=(d2-c1)/(M_ONE-c2)
    endif
    y2=(-M_ONE)/d2
    
  end subroutine bcorgn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 22.7.85                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine bcrmax(e, dr, rmax, h, s, n, yn, a, b)
    REAL_DOUBLE, intent(in)  :: e, dr, rmax
    REAL_DOUBLE, intent(in)  :: h(:), s(:)
    integer,     intent(in)  :: n
    REAL_DOUBLE, intent(out) :: yn
    REAL_DOUBLE, intent(in)  :: a, b

    REAL_DOUBLE :: tnm1, tn, tnp1, beta, dg, c1, c2, c3, dn
    
    PUSH_SUB(bcrmax)
    
    tnm1=h(n-1)-e*s(n-1)
    tn  =h(n  )-e*s(n  )
    tnp1=h(n+1)-e*s(n+1)
    beta=M_ONE+b/rmax
    dg=a*beta*(dr+M_ONE-M_HALF/beta)
    
    c2=CNST(24.)*dg/(CNST(12.)-tn)
    dn=-((CNST(24.)+CNST(10.)*tn)/(CNST(12.)-tn))
    
    c1= (M_ONE-tnm1/M_SIX)/(M_ONE-tnm1/CNST(12.))
    c3=-((M_ONE-tnp1/M_SIX)/(M_ONE-tnp1/CNST(12.)))
    yn=-((M_ONE-c1/c3)/(dn-c2/c3))
    
    POP_SUB(bcrmax)
  end subroutine bcrmax


  subroutine numin(e,h,s,y,n,nnode,yn,g,gsg,x,knk)
    REAL_DOUBLE, intent(in)  :: e, h(n), s(n)
    REAL_DOUBLE, intent(out) :: y(n)
    integer,     intent(in)  :: n
    integer,     intent(out) :: nnode
    REAL_DOUBLE, intent(in)  :: yn
    REAL_DOUBLE, intent(out) :: g, gsg, x
    integer,     intent(in)  :: knk

    integer :: i
    REAL_DOUBLE :: t
    
    ! no PUSH SUB, called too often
   
    y(n)=yn
    t=h(n)-e*s(n)
    g=y(n)/(M_ONE-t/CNST(12.))
    gsg=g*s(n)*g
    i=n-1
    y(i)=M_ONE
    t=h(i)-e*s(i)
    g=y(i)/(M_ONE-t/CNST(12.))
    gsg=gsg+g*s(i)*g
    x=y(i)-y(n)
    nnode=0
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  begin the inward integration by the Numerov method                         !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    do i = n - 2, knk, -1
      x=x+t*g
      y(i)=y(i+1)+x
      if( y(i)*y(i+1) .lt. M_ZERO) nnode=nnode+1
      t=h(i)-e*s(i)
      g=y(i)/(M_ONE-t/CNST(12.))
      gsg=gsg+g*s(i)*g
    enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  The kink radius is defined as the point where                              !
!  psi first turns downward.  this usually means at the outermost             !
!  maximum                                                                    !
!                                                                             !
!  the inward integration is now complete                                     !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine numin


  subroutine numout(e, h, s, y, ncor, knk, nnode, y2, g, gsg, x)
    REAL_DOUBLE, intent(in)    :: e, h(:), s(:)
    REAL_DOUBLE, intent(out)   :: y(:)
    integer,     intent(in)    :: ncor
    integer,     intent(inout) :: knk
    integer,     intent(out)   :: nnode
    REAL_DOUBLE, intent(in)    :: y2
    REAL_DOUBLE, intent(out)   :: g, gsg, x

    integer :: i, nm4
    REAL_DOUBLE :: t, xl
    
    ! no PUSH SUB, called too often
    
    y(1)=M_ZERO
    y(2)=y2
    t=h(2)-e*s(2)
    g=y(2)/(M_ONE-t/CNST(12.))
    gsg=g*s(2)*g
    y(3)=M_ONE
    t=h(3)-e*s(3)
    g=y(3)/(M_ONE-t/CNST(12.))
    gsg=gsg+g*s(3)*g
    x=y(3)-y(2)
    nnode=0
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  begin the outward integration by the Numerov method                        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nm4=knk-4
    knk = nm4
    do i = 4, nm4
      xl=x
      x=x+t*g
      y(i)=y(i-1)+x
      if( y(i)*y(i-1) .lt. M_ZERO) nnode=nnode+1
      t=h(i)-e*s(i)
      g=y(i)/(M_ONE-t/CNST(12.))
      gsg=gsg+g*s(i)*g
      if(.not. (nnode.lt.ncor .or. xl*x.gt.M_ZERO)) then
        knk = i
        exit
      endif
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the outward integration is now complete                                    !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine numout

end module atomic_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

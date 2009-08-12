#include <global.h>

module poisson_sete_m

  implicit none

  private

  public ::             &
    poisson_sete_init,  &
    pois,               &
    cbsurf_end

  public ::             &
    rho_nuc,            &
    esurf


  FLOAT, DIMENSION(:,:,:), ALLOCATABLE :: RHO, VH

  FLOAT, DIMENSION(:), ALLOCATABLE :: AW,AWP,Q2,X2,QS
  INTEGER, DIMENSION(:), ALLOCATABLE :: IAD,JAD,IY,JY,KY
  FLOAT :: TOL=0.01,HARTREE=2.0*13.60569193,BOHR=0.52917720859 ! Bohr radius in nm
  INTEGER :: IPOISSON_SETE_ON,NELT,NTOT,LENW,LENIW,IDEV,NGATES, &
    ISYM=0,ITOL=2,ITMAX=201,ITERMIN=5,ITER,IERR,IUNIT=0, &
    NXBOT,NYBOT,NZBOT,NXL,NYL,NXTOT,NYTOT,NZTOT,MD,PI,PCONST
  FLOAT, DIMENSION(:,:,:), ALLOCATABLE :: VBOUND,DIELECTRIC, &
    rhotest,VH_BIG 
  FLOAT, DIMENSION(:), ALLOCATABLE :: XG,YG,ZG, &
    DXG,DYG,DZ,RWORK,DXL,DYL,VT,VTV,ADIAG
  FLOAT :: XWIDTH,YWIDTH,ZWIDTH,DIELECTRIC0
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IPIO
  INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK,IDIAG
  INTEGER :: NX2,NY2,NZ2
  FLOAT, DIMENSION(:,:,:,:), ALLOCATABLE :: SIG
  FLOAT :: ESURF,CHARGE_TOP,CHARGE_BOT,CHARGE_TOT,BORDER,&
    CS1,CS2,CS3,CHARGE_SURF
  FLOAT, DIMENSION(:), ALLOCATABLE :: rho_nuc(:)!, v_nuc(:)

contains

  subroutine poisson_sete_init(nx, ny, nz, xl, yl, zl)
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    FLOAT,   intent(in) :: xl
    FLOAT,   intent(in) :: yl
    FLOAT,   intent(in) :: zl

    CHARACTER*40 :: CDUM,FIL1
    FLOAT :: BOHRNM, ANGSNM
    INTEGER :: J1,K1, M, I,J,K
    FLOAT :: ZCEN, XCEN, YCEN, VHMIN, VHMAX

    PI=DACOS(-1.D0)
    PCONST=4*PI ! need 4 pi for Hartrees I think (??)
    BOHRNM=BOHR/10.0
    ANGSNM=10.0/BOHR
    open(56,file='poisq',status='unknown')
    FIL1='test1.pois'
    OPEN(unit=57,FILE=FIL1,STATUS='UNKNOWN') ! status info file
    ! READ FILE
    READ(57,'(A40)')CDUM
    READ(57,'(A40)')CDUM
    READ(57,'(A40)')CDUM
    READ(57,*)IDEV, IPOISSON_SETE_ON ! for future use - device configuration
    write(*,*)"ISETE is ON", IPOISSON_SETE_ON
    IF(IDEV.EQ.1)THEN
      ! physical device characteristics
      READ(57,*)ZWIDTH;ZWIDTH=ZWIDTH*ANGSNM  ! distance between plates (nm)

      READ(57,*) ZCEN
      zcen = zcen*angsnm

      READ(57,*)NGATES;
      allocate(VTV(NGATES)); allocate(VT(NGATES))
      READ(57,*)VTV(1),VTV(2) ! voltage on plates in volts and nm?
      VT=VTV/HARTREE
      READ(57,*)DIELECTRIC0    ! dielectric constant of region betw plates
      DIELECTRIC0=DIELECTRIC0/PCONST
      ! mesh related
      READ(57,*)XWIDTH,YWIDTH ! bounding box around octopus box
      XWIDTH=XWIDTH*ANGSNM
      YWIDTH=YWIDTH*ANGSNM
      READ(57,*)XCEN,YCEN
      XCEN=XCEN*ANGSNM
      YCEN=YCEN*ANGSNM
      READ(57,*)NX2,NY2,NZ2   ! additional mesh points in bd box
      READ(57,*)NXL,NYL       ! padding points (to sim zone edge)
      allocate(DXL(NXL)); allocate(DYL(NYL))
      READ(57,*)(DXL(I),I=1,NXL)
      READ(57,*)(DYL(J),J=1,NYL)
    ELSE
      WRITE(6,*)' IDEV not supported '
      RETURN
    ENDIF

    call xyzgrid(nx, ny, nz, xl, yl, zl, xcen, ycen, zcen)  ! establish x-y-z grids

    CALL CBSURF !Assign boundaries 
    NELT=32+20*((NXTOT-2)+(NYTOT-2)+(NZTOT-2))+ &
      12*((NXTOT-2)*(NYTOT-2)+(NXTOT-2)*(NZTOT-2)+(NYTOT-2)*(NZTOT-2)) + &
      7*(NXTOT-2)*(NYTOT-2)*(NZTOT-2)
    LENW=NELT+9*NXTOT*NYTOT*NZTOT
    LENIW=NELT+5*NXTOT*NYTOT*NZTOT+12
    allocate(q2(ntot))
    allocate(AW(NELT))
    allocate(QS(NTOT))
    allocate(ADIAG(NTOT))
    allocate(IDIAG(NTOT))
    allocate(IAD(NELT))
    allocate(JAD(NELT))
    allocate(X2(NTOT)) 
    X2=0.0
    CALL POISSONM
    QS=Q2 ! store Dirichlet BC info in QS
    deallocate(Q2)

  end subroutine poisson_sete_init

  !---------------------------------------

  subroutine pois(icase, rho, vh, nx, ny, nz, xl, yl, zl, icalc)
    integer, intent(in)    :: icase
    FLOAT,   intent(in)    :: rho(:, :, :)
    FLOAT,   intent(inout) :: vh(:, :, :)
    integer, intent(in)    :: nx
    integer, intent(in)    :: ny
    integer, intent(in)    :: nz
    FLOAT,   intent(in)    :: xl
    FLOAT,   intent(in)    :: yl
    FLOAT,   intent(in)    :: zl
    integer, intent(in)    :: icalc

    ! Input: NX, NY, NZ -> No. of octopus grid points in each dimension
    ! Input: XL, YL, ZL -> Size of the octopus box

    CHARACTER*40 :: CDUM,FIL1
    FLOAT        :: bohrnm, angsnm, err
    integer      :: j1, k1, m, i, j, k, idum,i1
    FLOAT :: ZCEN, XCEN, YCEN, VHMIN, VHMAX
    PI=DACOS(-1.D0)
    PCONST=4*PI ! need 4 pi for Hartrees I think (??)
    BOHRNM=BOHR/10.0
    ANGSNM=10.0/BOHR

    write(*,*) "LENW LENIW IELT", LENW, LENIW, NELT
    allocate(Q2(NTOT))
    allocate(rhotest(nxtot,nytot,nztot))
    rhotest=0
    do i= 1,ntot
      Q2(i)=QS(i) ! recall stored BC info
    enddo
    DO M = 1,NTOT
      I=IY(M); J=JY(M); K=KY(M)
      IF(IPIO(I,J,K).EQ.0)THEN
        IF((I.GT.NXL+NXBOT.AND.I.LE.NXL+NXBOT+NX).AND. &
          (J.GT.NYL+NYBOT.AND.J.LE.NYL+NYBOT+NY).AND. &
          (K.GT.NZBOT.AND.K.LE.NZBOT+NZ))THEN
          I1=I-(NXL+NXBOT)
          J1=J-(NYL+NYBOT)
          K1=K-(NZBOT)
          IF(I1.GT.0.AND.J1.GT.0.AND.K1.GT.0)THEN ! this is redundant
            Q2(M)=Q2(M)+DYG(J)*DXG(I)*RHO(I1,J1,K1)
            rhotest(i,j,k)=rho(i1,j1,k1)
          ENDIF
        ENDIF
      ELSE
        Q2(M)=0.0
      ENDIF
    ENDDO
    write(*,*) size(Q2), size(ADIAG)
    do i=1,NTOT
      Q2(i)=Q2(i)/ADIAG(i)!! Laplacian scaled with diagonal elements
    enddo

    allocate(rwork(lenw))
    allocate(iwork(leniw))

    call dslucs(NTOT,Q2,X2,NELT,IAD,JAD,AW,ISYM,ITOL, &
      TOL,ITMAX,ITER,ITERMIN,ERR,IERR,IUNIT,RWORK,LENW, &
      IWORK,LENIW)

    leniw = iwork(9)
    lenw  = iwork(10) ! new work array lengths

    deallocate(rwork)
    deallocate(iwork)

    DO M=1,NTOT ! store X2(M) i,j,k wise
      I=IY(M); J=JY(M); K=KY(M)
      VH_BIG(I,J,K)=X2(M)
    ENDDO
    DO M = 1,NTOT
      I=IY(M); J=JY(M); K=KY(M)
      IF(IPIO(I,J,K).EQ.0)THEN
        IF((I.GT.NXL+NXBOT.AND.I.LE.NXL+NXBOT+NX).AND. &
          (J.GT.NYL+NYBOT.AND.J.LE.NYL+NYBOT+NY).AND. &
          (K.GT.NZBOT.AND.K.LE.NZBOT+NZ))THEN
          I1=I-(NXL+NXBOT)
          J1=J-(NYL+NYBOT)
          K1=K-(NZBOT)
          VH(I1,J1,K1)=X2(M)
        ENDIF
      ENDIF
    ENDDO
    call EGATE
    VHMIN=MINVAL(VH)
    VHMAX=MAXVAL(VH)

    deallocate(Q2)
    deallocate(rhotest)
    ! deallocate(VH_BIG);

  end subroutine pois


  subroutine poissonm
    ! version g created 03/25/08 - modified to include 
    !                              DIELECTRIC(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1)
    !           based on ~/LEVEL1/photovolt/poissonpv_v2.f90

    ! creates Laplacian.
    ! inputs:  NTOT,MD,NXTOT,NYTOT,NZTOT,IPIO(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1),
    !          VBOUND(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1),Q2(NTOT),
    !          DXG(NXTOT),DYG(NYTOT),DZ(NZTOT),DIELELCTRIC(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1)
    !
    ! outputs: AW(NELT),ADIAG(NTOT),Q2,
    !          NELT,IDIAG(NTOT),IAD(NELT),JAD(NELT)
    !
    ! local:   JOFF(7),AV(7),DCONST,D1,D2,DEL,
    !          I,J,K,IEC,IAV,N

    !    Note: original treatment of variable grid spacing documented
    !          in RIKEN I pg 50.

    FLOAT   :: AV(7)
    INTEGER :: JOFF(7)
    FLOAT   :: d1, d2, dconst, del
    integer :: i, j, k, ia, iav, icol, idebug, iec, n, ii
    !
    idebug=0
    if(idebug.eq.1)then
      open(57,file='aw.dat',status='unknown')
    endif
    !
    JOFF(1)=-MD ! used to compute row indices JAD
    JOFF(2)=-NZTOT; JOFF(3)=-1; JOFF(4)=0; JOFF(5)=1; JOFF(6)=NZTOT; JOFF(7)=MD
    I=1; J=1; K=1
    IEC=1 ! counter for number of non-zero entries in AW.

    q2(1:ntot) = M_ZERO

    DO IA=1,NTOT
      AV=0.D0
      IF(IPIO(I,J,K).EQ.0)THEN ! the point itself must be
        !                                   in the grid
        DO N=1,6 ! for all neighbors of this IA
          IF(N.GT.3)THEN
            IAV=N+1
          ELSE
            IAV=N
          ENDIF
          ! set grid intervals based on which direction
          ! this neighbor is located in. (+-I,+-J,+-K).
          SELECT CASE(IAV)
          CASE(1) ! I-1,J,K
            D1=DZ(K); D2=DYG(J)
            SELECT CASE(IPIO(I-1,J,K))
            CASE(1) ! Dirichlet boundary point
              DEL=0.5D0*DXG(I)*DZ(K)
              DCONST=DIELECTRIC(I,J,K)
              AV(IAV)=0.D0; AV(4)=AV(4)+DCONST*D1*D2/DEL
              Q2(IA)=Q2(IA)+VBOUND(I-1,J,K)*DCONST*D1*D2/DEL
            CASE(2) ! Neumann boundary point
              AV(IAV)=0.D0; DEL=0.D0
            CASE DEFAULT ! interior point
              DEL=0.5D0*(DXG(I)+DXG(I-1))*DZ(K)
              DCONST=0.5D0*(DIELECTRIC(I-1,J,K)+DIELECTRIC(I,J,K))
              AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
            END SELECT
          CASE(2) ! I,J-1,K
            D1=DZ(K); D2=DXG(I)
            SELECT CASE(IPIO(I,J-1,K))
            CASE(1) ! Dirichlet
              DEL=0.5D0*DYG(J)*DZ(K)
              DCONST=DIELECTRIC(I,J,K)
              AV(IAV)=0.D0; AV(4)=AV(4)+DCONST*D1*D2/DEL
              Q2(IA)=Q2(IA)+VBOUND(I,J-1,K)*DCONST*D1*D2/DEL
            CASE(2) ! Neumann
              AV(IAV)=0.D0; DEL=0.D0
            CASE DEFAULT! interior point
              DEL=0.5D0*(DYG(J)+DYG(J-1))*DZ(K)
              DCONST=0.5D0*(DIELECTRIC(I,J,K)+DIELECTRIC(I,J-1,K))
              AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
            END SELECT
          CASE(3) ! I,J,K-1
            D1=DXG(I); D2=DYG(J)
            SELECT CASE(IPIO(I,J,K-1))
            CASE(1) ! Dirichlet
              DEL=0.5D0*DZ(K)*DZ(K)
              DCONST=DIELECTRIC(I,J,K)
              AV(IAV)=0.D0; AV(4)=AV(4)+DCONST*D1*D2/DEL
              Q2(IA)=Q2(IA)+VBOUND(I,J,K-1)*DCONST*D1*D2/DEL
            CASE(2) ! Neumann
              AV(IAV)=0.D0; DEL=0.D0
            CASE DEFAULT! interior point
              DEL=0.5D0*(DZ(K)+DZ(K-1))*DZ(K)
              DCONST=0.5D0*(DIELECTRIC(I,J,K)+DIELECTRIC(I,J,K-1))
              AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
            END SELECT
          CASE(5) ! I,J,K+1
            D1=DXG(I); D2=DYG(J)
            SELECT CASE(IPIO(I,J,K+1))
            CASE(1) ! Dirichlet
              DEL=0.5D0*DZ(K)*DZ(K)
              DCONST=DIELECTRIC(I,J,K)
              AV(IAV)=0.D0; AV(4)=AV(4)+DCONST*D1*D2/DEL
              Q2(IA)=Q2(IA)+VBOUND(I,J,K+1)*DCONST*D1*D2/DEL
            CASE(2) ! Neumann
              AV(IAV)=0.D0; DEL=0.D0
            CASE DEFAULT! interior point
              DEL=0.5D0*(DZ(K)+DZ(K+1))*DZ(K)
              DCONST=0.5D0*(DIELECTRIC(I,J,K)+DIELECTRIC(I,J,K+1))
              AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
            END SELECT
          CASE(6) ! I,J+1,K
            D1=DZ(K); D2=DXG(I)
            SELECT CASE(IPIO(I,J+1,K))
            CASE(1) ! Dirichlet
              DEL=0.5D0*DYG(J)*DZ(K)
              DCONST=DIELECTRIC(I,J,K)
              AV(IAV)=0.D0; AV(4)=AV(4)+DCONST*D1*D2/DEL
              Q2(IA)=Q2(IA)+VBOUND(I,J+1,K)*DCONST*D1*D2/DEL
            CASE(2) ! Neumann
              AV(IAV)=0.D0; DEL=0.D0
            CASE DEFAULT
              DEL=0.5D0*(DYG(J)+DYG(J+1))*DZ(K)
              DCONST=0.5D0*(DIELECTRIC(I,J,K)+DIELECTRIC(I,J+1,K))
              AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
            END SELECT
          CASE(7) ! I+1,J,K
            D1=DZ(K); D2=DYG(J)
            SELECT CASE(IPIO(I+1,J,K))
            CASE(1) ! Dirichlet
              DEL=0.5D0*DXG(I)*DZ(K)
              DCONST=DIELECTRIC(I,J,K)
              AV(IAV)=0.D0; AV(4)=AV(4)+DCONST*D1*D2/DEL
              Q2(IA)=Q2(IA)+VBOUND(I+1,J,K)*DCONST*D1*D2/DEL
            CASE(2) ! Neumann
              AV(IAV)=0.D0; DEL=0.D0
            CASE DEFAULT
              DEL=0.5D0*(DXG(I)+DXG(I+1))*DZ(K)
              DCONST=0.5D0*(DIELECTRIC(I,J,K)+DIELECTRIC(I+1,J,K))
              AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
            END SELECT
          END SELECT
        ENDDO
      ELSE
        AV(4)=1.D0/PCONST ! for points that are out of the grid.
      ENDIF
      !                                 DSLUCS matrix and normalizations
      ADIAG(IA)=AV(4)
      IF(ADIAG(IA).EQ.0.D0)THEN
        WRITE(56,*)' DIAG ELEM ZERO ',IA
        STOP ' GOTTA QUIT '
      ENDIF

      if(idebug.eq.1)then
        write(57,101)iec,(av(ii),ii=1,7),q2(ia)
101     format(1x,i8,8(1x,e10.4))
      endif

      DO ICOL=1,7  !  THIS IS SLAP TRIAD FORMAT
        IF(AV(ICOL).NE.0.D0)THEN
          AW(IEC)=AV(ICOL)/ADIAG(IA)
          IAD(IEC)=IA
          JAD(IEC)=IA+JOFF(ICOL)
          IF(ICOL.EQ.4)THEN
            IDIAG(IA)=IEC
          ENDIF
          IEC=IEC+1
        ENDIF
      ENDDO
      !                                     increment I,J,K
      K=K+1
      IF(K.EQ.NZTOT+1)THEN
        K=1
        J=J+1
        IF(J.EQ.NYTOT+1)THEN
          J=1
          I=I+1
        ENDIF
      ENDIF


    ENDDO
    if(idebug.eq.1)then
      close(57)
    endif

    IEC=IEC-1 ! overcounted non-zero elements by one.
    write(6,*)' nelt, iec ',nelt,iec
    NELT=IEC

    ! in case NELT in static differs from that counted here
!!$        allocate(AWP(NELT))
!!$        DO I=1,NELT
!!$           AWP(I)=AW(I)
!!$        ENDDO
!!$        deallocate(AW)
!!$        allocate(AW(NELT))
!!$        AW=AWP
!!$        deallocate(AWP)

    if(idebug.eq.1)then
      close(57)
    endif

  end subroutine poissonm


  subroutine cap

    integer :: i, j, k

    !  9/24/08 - THIS VERSION OF CAP is very generic. It is
    !  for surfER. Assumes just IPIO(0:NX+1,0:NY+1,0:NZ+1), 
    !  VBOUND(0:NX+1,0:NY+1,0:NZ+1),DIELECTRIC(NX,NY,NZ).
    !  Outputs a single 3D array SIG(NX,NY,NZ)
    !
    !  for calculating capacitances
    !
    !  SIG(I,J,K) is capacitance
    !
    !  input: NX,NY,NZ
    !         IPIO,VBOUND
    !         ZG,DXG,DYG,DZ
    !         VH_BIG,DIELECTRIC,PCONST
    !
    !  output: SIG(I,J,K,6) 6 faces of grid cube 
    !          used in ETOT
    !  1 I-1, 2 I+1, 3 J-1, 4 J+1, 5 K-1, 6 K+1
    !
    !  local: BORDER,PCONST
    !
    !  a homogeeneous grid approximation is made: derivative
    !  at metal-semi border calculated discretely with 1/DXG (1/DYG,1/DZ)
    !  rather than properly averaged spacings.
    write(*,*) "In capser.F90"
    allocate(SIG(NXTOT,NYTOT,NZTOT,6))
    write(*,*) "Just called SIG"
    SIG=0.0
    DO I=1,NXTOT 
      DO J=1,NYTOT
        DO K=1,NZTOT
          IF(IPIO(I,J,K).EQ.0)THEN ! I,J,K is in Poisson grid
            IF(IPIO(I-1,J,K).EQ.1)THEN
              BORDER=2*(VH_BIG(I,J,K)-VBOUND(I-1,J,K))/DXG(I)
              SIG(I,J,K,1)=SIG(I,J,K,1)-BORDER*DYG(J)*DZ(K)* &
                DIELECTRIC(I,J,K)
            ENDIF
            IF(IPIO(I+1,J,K).EQ.1)THEN
              BORDER=2*(VH_BIG(I,J,K)-VBOUND(I+1,J,K))/DXG(I)
              SIG(I,J,K,2)=SIG(I,J,K,2)-BORDER*DYG(J)*DZ(K)* &
                DIELECTRIC(I,J,K)
            ENDIF
            IF(IPIO(I,J-1,K).EQ.1)THEN
              BORDER=2*(VH_BIG(I,J,K)-VBOUND(I,J-1,K))/DYG(J)
              SIG(I,J,K,3)=SIG(I,J,K,3)-BORDER*DXG(I)*DZ(K)* &
                DIELECTRIC(I,J,K)
            ENDIF
            IF(IPIO(I,J+1,K).EQ.1)THEN
              BORDER=2*(VH_BIG(I,J,K)-VBOUND(I,J+1,K))/DYG(J)
              SIG(I,J,K,4)=SIG(I,J,K,4)-BORDER*DXG(I)*DZ(K)* &
                DIELECTRIC(I,J,K)
            ENDIF
            IF(IPIO(I,J,K-1).EQ.1)THEN
              BORDER=2*(VH_BIG(I,J,K)-VBOUND(I,J,K-1))/DZ(K)
              SIG(I,J,K,5)=SIG(I,J,K,5)-BORDER*DXG(I)*DYG(J)* &
                DIELECTRIC(I,J,K)
            ENDIF
            IF(IPIO(I,J,K+1).EQ.1)THEN
              BORDER=2*(VH_BIG(I,J,K)-VBOUND(I,J,K+1))/DZ(K)
              SIG(I,J,K,6)=SIG(I,J,K,6)-BORDER*DXG(I)*DYG(J)* &
                DIELECTRIC(I,J,K)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  end subroutine cap

  subroutine cbsurf

    integer :: k

    !  IFILL=1 parallel plate
    !
    ! input:  IFILL,IFREE_SURFACE,NX,NY,NXL,NYL
    !         VFILL,XG(NX),YG(NY),XNG(:),YNG(:),VT(NGATES)
    !         XMID,YMID
    ! output: VSURF,VBOUND(0:NX+1,0:NY+1,0:NZ+1)
    !         IPATTERN,IPIO


    ! This is a combination of cbsurfser and spaceser located on
    ! /home/stopa/surfER/build
    write(*,*) "Allocating data"
    allocate(VH_BIG(NXTOT,NYTOT,NZTOT))
    allocate(IPIO(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1))
    allocate(VBOUND(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1))
    allocate(DIELECTRIC(NXTOT,NYTOT,NZTOT))
    IPIO=0; VBOUND=0.0
    IF(IDEV.EQ.1)THEN ! parallel plates
      DIELECTRIC=DIELECTRIC0
      IPIO(:,:,0)=1
      IPIO(:,:,NZTOT+1)=1
      VBOUND(:,:,0)=VT(1)
      VBOUND(:,:,NZTOT+1)=VT(2)
      DO K=1,NZTOT
        IF(K.NE.0.AND.K.NE.NZTOT+1)THEN
          ! lateral BC's
          IPIO(0,:,K)=2
          IPIO(NXTOT+1,:,K)=2
          IPIO(:,0,K)=2
          IPIO(:,NYTOT+1,K)=2
        ENDIF
      ENDDO
    ENDIF

  end subroutine cbsurf

  subroutine cbsurf_end

    deallocate(IPIO)
    deallocate(VBOUND)
    deallocate(VT);deallocate(VTV);deallocate(QS)
    deallocate(IAD);deallocate(JAD);deallocate(ADIAG)
    deallocate(IDIAG);deallocate(X2)
    deallocate(DXL);deallocate(DYL)
    !	deallocate(VH_BIG)
  end subroutine cbsurf_end


  subroutine egate

    integer :: i, j, k

    !
    !     to do energy calculations previously performed in setr
    !  compute capacitances of six surfaces
    call CAP
    ESURF=0.0
    CHARGE_SURF=0.0
    CHARGE_TOP=0.0

    CHARGE_BOT=0.0
    CS1=0.0; CS2=0.0; CS3=0.0
    !write(*,*) "Begin cycle ", NXTOT, NYTOT, NZTOT
    DO I=1,NXTOT
      DO J=1,NYTOT
        DO K=1,NZTOT ! N.B. [SIG]=charge (not charge/Area)
          ESURF=ESURF+VBOUND(I-1,J,K)*SIG(I,J,K,1)
          ESURF=ESURF+VBOUND(I+1,J,K)*SIG(I,J,K,2)
          ESURF=ESURF+VBOUND(I,J-1,K)*SIG(I,J,K,3)
          ESURF=ESURF+VBOUND(I,J+1,K)*SIG(I,J,K,4)
          ESURF=ESURF+VBOUND(I,J,K-1)*SIG(I,J,K,5)
          ESURF=ESURF+VBOUND(I,J,K+1)*SIG(I,J,K,6)
          CS1=CS1+SIG(I,J,K,1)+SIG(I,J,K,2)
          CS2=CS2+SIG(I,J,K,3)+SIG(I,J,K,4)
          CS3=CS3+SIG(I,J,K,5)+SIG(I,J,K,6)
        ENDDO
        CHARGE_TOP=CHARGE_TOP+SIG(I,J,1,5)+SIG(I,J,1,6)
        CHARGE_BOT=CHARGE_BOT+SIG(I,J,NZTOT,5)+SIG(I,J,NZTOT,6)
      ENDDO
    ENDDO
    deallocate(SIG)

    ESURF=0.5*ESURF
    CHARGE_SURF=CS1+CS2+CS3

    WRITE(520,*)'Esurf', ESURF*27.2,' CS1 ',CS1,' CS2 ',CS2,' CS3 ',CS3
    RETURN
  end subroutine egate

  subroutine xyzgrid(nx, ny, nz, xl, yl, zl, xcen, ycen, zcen)
    ! Input: NX, NY, NZ -> No. of octopus grid points in each dimension
    ! Input: XL, YL, ZL -> Size of the octopus box

    integer, intent(in) :: nx, ny, nz
    FLOAT,   intent(in) :: xl, yl, zl
    FLOAT,   intent(in) :: xcen, ycen, zcen

    FLOAT   :: dx1b, dx1t, dx_oct, xbot, xtest, xtop, xwidth_big
    FLOAT   :: dy1b, dy1t, dy_oct, ybot, ytop, ywidth_big
    FLOAT   :: dz1b, dz1t, dz_oct, zbot, ztest, ztop
    integer :: i, j, k, m
    integer :: nxtop, nytop, nztop
    integer :: nxtot_0, nytot_0


    ! inputs:  NZ,NZ2,NX,NX2,NY,NY2
    !          XCEN,YCEN,ZCEN
    !          XWIDTH,YWIDTH,ZWIDTH
    !          XL,YL,ZL
    !
    ! outputs: ZG,XG,YG        
    !write(62,*)"XL,YL,ZL ",XL,YL,ZL
    ! z-mesh -- differs from x,y mesh because no border points
    !           hence NZTOT_0 = NZTOT (ie there IS no NZTOT_0)
    nztot = nz + nz2

    allocate(zg(nztot))
    allocate(dz(nztot))

    !top and bottom z coordinates.
    !Depend on POISSON_SETE (ZWIDTH) 
    !and Octopus (ZL)parameters.
    ztop = CNST(0.5)*zwidth - zcen - CNST(0.5)*zl
    zbot = CNST(0.5)*zwidth + zcen - CNST(0.5)*zl

    nztop = int(nz2*ztop/(ztop + zbot))
    nzbot = nz2 - nztop

    DZ1B=ZBOT/TOFLOAT(NZBOT)
    DO K=1,NZBOT
      DZ(K)=DZ1B
    ENDDO
    DZ_OCT=ZL/TOFLOAT(NZ)

    do k = nzbot + 1, nzbot + nz
      dz(k) = dz_oct
    end do

    DZ1T=ZTOP/TOFLOAT(NZTOP)
    DO K=NZBOT+NZ+1,NZTOT
      DZ(K)=DZ1T
    ENDDO

    ZG(1)=-0.5*ZWIDTH+0.5*DZ(1)
    DO K=2,NZTOT-1
      ZG(K)=ZG(K-1)+0.5*(DZ(K-1)+DZ(K))
    ENDDO
    ZG(NZTOT)=0.5*ZWIDTH-0.5*DZ(NZTOT)

    ZTEST=ZG(NZTOT)
    write(6,*)' nzbot ',nzbot
    write(56,*)' zg nztot ',nztot
    write(56,101)(zg(k),k=1,nztot)
101 format(10(1x,e10.4))
    write(56,*)' dz '
    write(56,101)(dz(k),k=1,nztot)

    ! x-mesh
    NXTOT_0=NX+NX2
    NXTOT=NXTOT_0+2*NXL
    allocate(xg(nxtot))
    allocate(DXG(NXTOT))
    DO I=1,NXL  !  border points
      DXG(I)=DXL(I)
    ENDDO
    DO I=NXL+NXTOT_0+1,NXTOT
      DXG(I)=DXL(NXTOT-I+1)
    ENDDO
    ! outside Octopus mesh
    XTOP=0.5*XWIDTH-XCEN-0.5*XL
    XBOT=0.5*XWIDTH+XCEN-0.5*XL
    NXTOP=INT(NX2*XTOP/(XTOP+XBOT))
    NXBOT=NX2-NXTOP
    DX1B=XBOT/TOFLOAT(NXBOT) ! "bottom" region
    DO I=NXL+1,NXL+NXBOT
      DXG(I)=DX1B
    ENDDO
    DX_OCT=XL/TOFLOAT(NX)  ! Octopus region
    DO I=NXL+NXBOT+1,NXL+NXBOT+NX
      DXG(I)=DX_OCT
    ENDDO
    DX1T=XTOP/TOFLOAT(NXTOP)  ! "top" region
    DO I=NXL+NXBOT+NX+1,NXL+NXBOT+NX+NXTOP
      DXG(I)=DX1T
    ENDDO

    XWIDTH_BIG=SUM(DXG)
    XG(1)=-0.5*XWIDTH_BIG+0.5*DXG(1)
    DO I=2,NXTOT-1
      XG(I)=XG(I-1)+0.5*(DXG(I-1)+DXG(I))
    ENDDO
    XG(NXTOT)=0.5*XWIDTH_BIG-0.5*DXG(NXTOT)

    XTEST=XG(NXTOT)
    write(56,*)' xg nxtot ',nxtot
    write(56,101)(xg(k),k=1,nxtot)
    write(56,*)' dxg '
    write(56,101)(dxg(k),k=1,nxtot)

    ! Y-mesh
    NYTOT_0=NY+NY2  ! excluding border points
    NYTOT=NYTOT_0+2*NYL  ! all Poisson Y mesh point
    allocate(yg(nytot))
    allocate(dyg(nytot))
    DO I=1,NYL  !  border points
      DYG(I)=DYL(I)
    ENDDO
    DO I=NYL+NYTOT_0+1,NYTOT
      DYG(I)=DYL(NYTOT-I+1)
    ENDDO
    ! outside Octopus mesh
    YTOP=0.5*YWIDTH-YCEN-0.5*YL
    YBOT=0.5*YWIDTH+YCEN-0.5*YL
    NYTOP=INT(NY2*YTOP/(YTOP+YBOT))
    NYBOT=NY2-NYTOP
    DY1B=YBOT/TOFLOAT(NYBOT) ! "bottom" region
    DO I=NYL+1,NYL+NYBOT
      DYG(I)=DY1B
    ENDDO
    DY_OCT=YL/TOFLOAT(NY)  ! Octopus region
    DO I=NYL+NYBOT+1,NYL+NYBOT+NY
      DYG(I)=DY_OCT
    ENDDO
    DY1T=YTOP/TOFLOAT(NYTOP)  ! "top" region
    DO I=NYL+NYBOT+NY+1,NYL+NYBOT+NY+NYTOP
      DYG(I)=DY1T
    ENDDO

    YWIDTH_BIG=SUM(DYG)
    YG(1)=-0.5*YWIDTH_BIG+0.5*DYG(1)
    DO I=2,NYTOT-1
      YG(I)=YG(I-1)+0.5*(DYG(I-1)+DYG(I))
    ENDDO
    YG(NYTOT)=0.5*YWIDTH_BIG-0.5*DYG(NYTOT)


    NTOT=NXTOT*NYTOT*NZTOT
    MD=NZTOT*NYTOT

    allocate(iy(ntot))
    allocate(jy(ntot))
    allocate(ky(ntot))
    m=0
    DO I=1,NXTOT
      DO J=1,NYTOT
        DO K=1,NZTOT
          M=M+1
          IY(M)=I
          JY(M)=J
          KY(M)=K
        ENDDO
      ENDDO
    ENDDO

!!$        dxg=1.0; dyg=1.0; dz=1.0
!!$        zg(1)=dz(1)*0.5
!!$        do i=2,nztot-1
!!$           zg(i)=zg(i-1)+0.5*(dz(i)+dz(i-1))
!!$        enddo
!!$        zg(nztot)=zg(nztot-1)+0.5*dz(nztot)
!!$
!!$        xg(1)=dxg(1)*0.5
!!$        do i=2,nxtot-1
!!$           xg(i)=xg(i-1)+0.5*(dxg(i)+dxg(i-1))
!!$        enddo
!!$        xg(nxtot)=xg(nxtot-1)+0.5*dxg(nxtot)
!!$
!!$        yg(1)=dyg(1)*0.5
!!$        do i=2,nytot-1
!!$           yg(i)=yg(i-1)+0.5*(dyg(i)+dyg(i-1))
!!$        enddo
!!$        yg(nytot)=yg(nytot-1)+0.5*dyg(nytot)
!!$           

    write(56,*)' yg '
    write(56,101)(yg(k),k=1,nytot)
    write(56,*)' dyg '
    write(56,101)(dyg(k),k=1,nytot)

  end subroutine xyzgrid


end module poisson_sete_m


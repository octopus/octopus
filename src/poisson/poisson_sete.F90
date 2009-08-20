!! Copyright (C) 2009 R. Olivares, M. Stopa, X. Andrade
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
!! $Id: poisson_multigrid.F90 5486 2009-05-24 14:09:52Z xavier $

#include <global.h>

module poisson_sete_m
  use global_m
  use math_m
  use messages_m
  use profiling_m

  implicit none

  private

  public ::              &
    poisson_sete_t,      &
    poisson_sete_init,   &
    poisson_sete_solve,  &
    poisson_sete_end,    &
    poisson_sete_energy

  public ::              &
    rho_nuc,             &
    count_atoms

  type poisson_sete_t
    integer          :: isete_on
    integer, pointer :: ipio(:,:,:)
    FLOAT,   pointer :: adiag(:)
    FLOAT,   pointer :: qs(:)
    FLOAT,   pointer :: vbound(:,:,:)
    integer, pointer :: iad(:)
    integer, pointer :: jad(:)
    FLOAT,   pointer :: aw(:)
  end type poisson_sete_t

  integer, allocatable :: iy(:), jy(:), ky(:)
  integer :: nelt, ntot, lenw, leniw, idev, isym = 0, itol = 2, itmax = 201, itermin = 5, iter, ierr, iunit = 0, &
    nxbot, nybot, nzbot, nxl, nyl, nxtot, nytot, nztot
  FLOAT, allocatable :: dielectric(:,:,:), vh_big(:,:,:)
  FLOAT, allocatable :: xg(:), yg(:), zg(:), dxg(:), dyg(:), dz(:), dxl(:), dyl(:)
  FLOAT :: xwidth, ywidth, zwidth, pconst
  FLOAT :: esurf, charge_top, charge_bot, charge_tot, border, charge_surf
  FLOAT, allocatable :: rho_nuc(:)
  integer :: noatoms, count_atoms

  FLOAT, parameter :: tol = 0.01, hartree = CNST(2.0*13.60569193), BOHR = CNST(0.52917720859)

contains

  subroutine poisson_sete_init(this, nx, ny, nz, xl, yl, zl, number_atoms)
    type(poisson_sete_t), intent(out)   :: this
    integer,              intent(in)    :: nx
    integer,              intent(in)    :: ny
    integer,              intent(in)    :: nz
    FLOAT,                intent(in)    :: xl
    FLOAT,                intent(in)    :: yl
    FLOAT,                intent(in)    :: zl
    integer,              intent(in)    :: number_atoms

    character(len=40) :: cdum, fil1
    FLOAT :: bohrnm, angsnm
    integer :: i, j, ngates
    FLOAT :: zcen, xcen, ycen, dielectric0
    FLOAT, allocatable :: vtv(:), vt(:), q2(:)
    integer :: nx2, ny2, nz2, md

    pconst = CNST(4.0)*M_PI ! need 4 pi for Hartrees I think (??)
    bohrnm = BOHR/CNST(10.0)
    angsnm = CNST(10.0)/BOHR
    noatoms = number_atoms
    count_atoms = 0

    open(56, file = 'poisq', status = 'unknown')
    FIL1='test1.pois'
    OPEN(unit=57, FILE = FIL1, STATUS = 'UNKNOWN') ! status info file

    ! READ FILE
    read(57, '(a40)') cdum
    read(57, '(a40)') cdum
    read(57, '(a40)') cdum
    read(57,*) idev, this%isete_on ! for future use - device configuration
    write(*,*) "ISETE is ON", this%isete_on

    if(idev == 1) then

      ! physical device characteristics
      read(57,*) zwidth

      zwidth = zwidth*angsnm  ! distance between plates (nm)

      read(57,*) zcen
      zcen = zcen*angsnm

      read(57,*) ngates

      SAFE_ALLOCATE(VTV(1:NGATES))
      SAFE_ALLOCATE(VT(1:NGATES))

      READ(57,*)VTV(1),VTV(2) ! voltage on plates in volts and nm?
      VT=VTV/HARTREE

      SAFE_DEALLOCATE_A(VTV)

      READ(57,*) dielectric0    ! dielectric constant of region betw plates
      dielectric0 = dielectric0/pconst

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
    end if

    call xyzgrid(this, nx, ny, nz, xl, yl, zl, xcen, ycen, zcen, nx2, ny2, nz2)  ! establish x-y-z grids

    call cbsurf(this, vt) !Assign boundaries 

    NELT=32+20*((NXTOT-2)+(NYTOT-2)+(NZTOT-2))+ &
      12*((NXTOT-2)*(NYTOT-2)+(NXTOT-2)*(NZTOT-2)+(NYTOT-2)*(NZTOT-2)) + &
      7*(NXTOT-2)*(NYTOT-2)*(NZTOT-2)
    LENW=NELT+9*NXTOT*NYTOT*NZTOT
    LENIW=NELT+5*NXTOT*NYTOT*NZTOT+12

    SAFE_ALLOCATE(q2(1:ntot))
    SAFE_ALLOCATE(this%aw(1:nelt))
    SAFE_ALLOCATE(this%qs(1:ntot))
    SAFE_ALLOCATE(this%adiag(1:ntot))
    SAFE_ALLOCATE(this%iad(1:nelt))
    SAFE_ALLOCATE(this%jad(1:nelt))

    call poissonm(this)

    this%qs = q2 ! store Dirichlet BC info in this%QS

    SAFE_DEALLOCATE_A(q2)

  contains

    subroutine poissonm(this)
      type(poisson_sete_t), intent(in)    :: this

      ! version g created 03/25/08 - modified to include 
      !                              dielectric(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1)
      !           based on ~/LEVEL1/photovolt/poissonpv_v2.f90

      ! creates Laplacian.
      ! inputs:  NTOT,MD,NXTOT,NYTOT,NZTOT,this%ipio(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1),
      !          this%VBOUND(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1),Q2(NTOT),
      !          DXG(NXTOT),DYG(NYTOT),DZ(NZTOT),DIELELCTRIC(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1)
      !
      ! outputs: this%AW(NELT),this%ADIAG(NTOT),Q2,
      !          NELT, IAD(NELT),this%JAD(NELT)
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
      if(idebug == 1)then
        open(57,file='this%aw.dat',status='unknown')
      endif
      !
      JOFF(1)=-MD ! used to compute row indices this%JAD
      JOFF(2)=-NZTOT; JOFF(3)=-1; JOFF(4)=0; JOFF(5)=1; JOFF(6)=NZTOT; JOFF(7)=MD
      I=1; J=1; K=1
      IEC=1 ! counter for number of non-zero entries in this%AW.

      q2(1:ntot) = M_ZERO

      DO IA=1,NTOT
        AV=M_ZERO
        IF(this%ipio(I,J,K) == 0)THEN ! the point itself must be
          !                                   in the grid
          DO N=1,6 ! for all neighbors of this IA
            IF(N.GT.3)THEN
              IAV=N+1
            ELSE
              IAV=N
            end if
            ! set grid intervals based on which direction
            ! this neighbor is located in. (+-I,+-J,+-K).
            SELECT CASE(IAV)
            CASE(1) ! I-1,J,K
              D1=DZ(K); D2=DYG(J)
              SELECT CASE(this%ipio(I-1,J,K))
              CASE(1) ! Dirichlet boundary point
                DEL=CNST(0.5)*DXG(I)*DZ(K)
                DCONST=dielectric(I,J,K)
                AV(IAV)=M_ZERO; AV(4)=AV(4)+DCONST*D1*D2/DEL
                Q2(IA)=Q2(IA)+this%VBOUND(I-1,J,K)*DCONST*D1*D2/DEL
              CASE(2) ! Neumann boundary point
                AV(IAV)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT ! interior point
                DEL=CNST(0.5)*(DXG(I)+DXG(I-1))*DZ(K)
                DCONST=CNST(0.5)*(dielectric(I-1,J,K)+dielectric(I,J,K))
                AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
              END SELECT
            CASE(2) ! I,J-1,K
              D1=DZ(K); D2=DXG(I)
              SELECT CASE(this%ipio(I,J-1,K))
              CASE(1) ! Dirichlet
                DEL=CNST(0.5)*DYG(J)*DZ(K)
                DCONST=dielectric(I,J,K)
                AV(IAV)=M_ZERO; AV(4)=AV(4)+DCONST*D1*D2/DEL
                Q2(IA)=Q2(IA)+this%VBOUND(I,J-1,K)*DCONST*D1*D2/DEL
              CASE(2) ! Neumann
                AV(IAV)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT! interior point
                DEL=CNST(0.5)*(DYG(J)+DYG(J-1))*DZ(K)
                DCONST=CNST(0.5)*(dielectric(I,J,K)+dielectric(I,J-1,K))
                AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
              END SELECT
            CASE(3) ! I,J,K-1
              D1=DXG(I); D2=DYG(J)
              SELECT CASE(this%ipio(I,J,K-1))
              CASE(1) ! Dirichlet
                DEL=CNST(0.5)*DZ(K)*DZ(K)
                DCONST=dielectric(I,J,K)
                AV(IAV)=M_ZERO; AV(4)=AV(4)+DCONST*D1*D2/DEL
                Q2(IA)=Q2(IA)+this%VBOUND(I,J,K-1)*DCONST*D1*D2/DEL
              CASE(2) ! Neumann
                AV(IAV)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT! interior point
                DEL=CNST(0.5)*(DZ(K)+DZ(K-1))*DZ(K)
                DCONST=CNST(0.5)*(dielectric(I,J,K)+dielectric(I,J,K-1))
                AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
              END SELECT
            CASE(5) ! I,J,K+1
              D1=DXG(I); D2=DYG(J)
              SELECT CASE(this%ipio(I,J,K+1))
              CASE(1) ! Dirichlet
                DEL=CNST(0.5)*DZ(K)*DZ(K)
                DCONST=dielectric(I,J,K)
                AV(IAV)=M_ZERO; AV(4)=AV(4)+DCONST*D1*D2/DEL
                Q2(IA)=Q2(IA)+this%VBOUND(I,J,K+1)*DCONST*D1*D2/DEL
              CASE(2) ! Neumann
                AV(IAV)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT! interior point
                DEL=CNST(0.5)*(DZ(K)+DZ(K+1))*DZ(K)
                DCONST=CNST(0.5)*(dielectric(I,J,K)+dielectric(I,J,K+1))
                AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
              END SELECT
            CASE(6) ! I,J+1,K
              D1=DZ(K); D2=DXG(I)
              SELECT CASE(this%ipio(I,J+1,K))
              CASE(1) ! Dirichlet
                DEL=CNST(0.5)*DYG(J)*DZ(K)
                DCONST=dielectric(I,J,K)
                AV(IAV)=M_ZERO; AV(4)=AV(4)+DCONST*D1*D2/DEL
                Q2(IA)=Q2(IA)+this%VBOUND(I,J+1,K)*DCONST*D1*D2/DEL
              CASE(2) ! Neumann
                AV(IAV)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT
                DEL=CNST(0.5)*(DYG(J)+DYG(J+1))*DZ(K)
                DCONST=CNST(0.5)*(dielectric(I,J,K)+dielectric(I,J+1,K))
                AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
              END SELECT
            CASE(7) ! I+1,J,K
              D1=DZ(K); D2=DYG(J)
              SELECT CASE(this%ipio(I+1,J,K))
              CASE(1) ! Dirichlet
                DEL=CNST(0.5)*DXG(I)*DZ(K)
                DCONST=dielectric(I,J,K)
                AV(IAV)=M_ZERO; AV(4)=AV(4)+DCONST*D1*D2/DEL
                Q2(IA)=Q2(IA)+this%VBOUND(I+1,J,K)*DCONST*D1*D2/DEL
              CASE(2) ! Neumann
                AV(IAV)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT
                DEL=CNST(0.5)*(DXG(I)+DXG(I+1))*DZ(K)
                DCONST=CNST(0.5)*(dielectric(I,J,K)+dielectric(I+1,J,K))
                AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
              END SELECT
            END SELECT
          end do
        ELSE
          write(57,*) "At else", I, J, K
          AV(4) = M_ONE/PCONST ! for points that are out of the grid.
        end if
        !                                 DSLUCS matrix and normalizations
        this%ADIAG(IA)=AV(4)
        IF(this%ADIAG(IA) == M_ZERO)THEN
          WRITE(*,*)' DIAG ELEM ZERO ',IA
          STOP ' GOTTA QUIT '
        end if

        if(idebug == 1) then
          write(57, '(1x,i8,8(1x,e10.4))') iec, (av(ii),ii=1,7), q2(ia)
        endif

        DO ICOL=1,7  !  this IS SLAP TRIAD FORMAT
          IF(AV(ICOL) /= M_ZERO)THEN
            this%AW(IEC)=AV(ICOL)/this%ADIAG(IA)
            this%IAD(IEC)=IA
            this%JAD(IEC)=IA+JOFF(ICOL)
            IEC=IEC+1
          end if
        end do
        !                                     increment I,J,K
        K=K+1
        IF(K == NZTOT+1)THEN
          K=1
          J=J+1
          IF(J == NYTOT+1)THEN
            J=1
            I=I+1
          end if
        end if


      end do
      if(idebug == 1)then
        close(57)
      endif

      IEC=IEC-1 ! overcounted non-zero elements by one.
      write(6,*)' nelt, iec ',nelt,iec
      NELT=IEC

      ! in case NELT in static differs from that counted here
!!$        allocate(this%AWP(NELT))
!!$        DO I=1,NELT
!!$           this%AWP(I)=this%AW(I)
!!$        end do
!!$        deallocate(this%AW)
!!$        allocate(this%AW(NELT))
!!$        this%AW=this%AWP
!!$        deallocate(this%AWP)

      if(idebug == 1)then
        close(57)
      endif
    end subroutine poissonm

    !---------------------------------------------------------

    subroutine cbsurf(this, vt)
      type(poisson_sete_t), intent(in)    :: this
      FLOAT,                intent(in)    :: vt(:)

      integer :: k

      !  IFILL=1 parallel plate
      !
      ! input:  IFILL,IFREE_SURFACE,NX,NY,NXL,NYL
      !         VFILL,XG(NX),YG(NY),XNG(:),YNG(:),VT(NGATES)
      !         XMID,YMID
      ! output: VSURF,this%VBOUND(0:NX+1,0:NY+1,0:NZ+1)
      !         IPATTERN,this%ipio


      ! This is a combination of cbsurfser and spaceser located on
      ! /home/stopa/surfER/build
      allocate(this%ipio(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1))
      allocate(VH_BIG(1:NXTOT,1:NYTOT,1:NZTOT))
      vh_big(1:NXTOT,1:NYTOT,1:NZTOT)= M_ZERO
      allocate(this%VBOUND(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1))
      allocate(dielectric(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1))
      this%ipio=0; this%VBOUND= M_ZERO
      IF(IDEV == 1)THEN ! parallel plates
        dielectric=dielectric0
        dielectric(:,:,0)=6.5/pconst  !For platinum
        dielectric(:,:,NZTOT+1)=6.5/pconst !For platinum
        this%ipio(:,:,0)=1
        this%ipio(:,:,NZTOT+1)=1
        this%VBOUND(:,:,0)=VT(1)
        this%VBOUND(:,:,NZTOT+1)=VT(2)
        DO K=1,NZTOT
          IF(K /= 0 .AND. K /= NZTOT+1)THEN !Redundant?
            ! lateral BC`s
            this%ipio(0,:,K)=2
            this%ipio(NXTOT+1,:,K)=2
            this%ipio(:,0,K)=2
            this%ipio(:,NYTOT+1,K)=2
          end if
        end do
      end if

    end subroutine cbsurf

    !---------------------------------------------------------

    subroutine xyzgrid(this, nx, ny, nz, xl, yl, zl, xcen, ycen, zcen, nx2, ny2, nz2)
      type(poisson_sete_t), intent(in)    :: this
      integer,              intent(in)    :: nx
      integer,              intent(in)    :: ny
      integer,              intent(in)    :: nz
      FLOAT,                intent(in)    :: xl
      FLOAT,                intent(in)    :: yl
      FLOAT,                intent(in)    :: zl
      FLOAT,                intent(in)    :: xcen
      FLOAT,                intent(in)    :: ycen
      FLOAT,                intent(in)    :: zcen
      integer,              intent(in)    :: nx2
      integer,              intent(in)    :: ny2
      integer,              intent(in)    :: nz2

      ! Input: NX, NY, NZ -> No. of octopus grid points in each dimension
      ! Input: XL, YL, ZL -> Size of the octopus box

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
      end do
      DZ_OCT=ZL/TOFLOAT(NZ)

      do k = nzbot + 1, nzbot + nz
        dz(k) = dz_oct
      end do

      DZ1T=ZTOP/TOFLOAT(NZTOP)
      DO K=NZBOT+NZ+1,NZTOT
        DZ(K)=DZ1T
      end do

      ZG(1)=-CNST(0.5)*ZWIDTH+CNST(0.5)*DZ(1)
      DO K=2,NZTOT-1
        ZG(K)=ZG(K-1)+CNST(0.5)*(DZ(K-1)+DZ(K))
      end do
      ZG(NZTOT)=CNST(0.5)*ZWIDTH-CNST(0.5)*DZ(NZTOT)

      ZTEST=ZG(NZTOT)
      !    write(6,*)' nzbot ',nzbot
      !    write(6,*)' zg nztot ',nztot
      !    write(6,101)(zg(k),k=1,nztot)
      ! 101 format(10(1x,e10.4))
      !    write(6,*)' dz '
      !    write(6,101)(dz(k),k=1,nztot)

      ! x-mesh
      NXTOT_0=NX+NX2
      NXTOT=NXTOT_0+2*NXL
      allocate(xg(nxtot))
      allocate(DXG(NXTOT))
      DO I=1,NXL  !  border points
        DXG(I)=DXL(I)
      end do
      DO I=NXL+NXTOT_0+1,NXTOT
        DXG(I)=DXL(NXTOT-I+1)
      end do
      ! outside Octopus mesh
      XTOP=CNST(0.5)*XWIDTH-XCEN-CNST(0.5)*XL
      XBOT=CNST(0.5)*XWIDTH+XCEN-CNST(0.5)*XL
      NXTOP=INT(NX2*XTOP/(XTOP+XBOT))
      NXBOT=NX2-NXTOP
      DX1B=XBOT/TOFLOAT(NXBOT) ! "bottom" region
      DO I=NXL+1,NXL+NXBOT
        DXG(I)=DX1B
      end do
      DX_OCT=XL/TOFLOAT(NX)  ! Octopus region
      DO I=NXL+NXBOT+1,NXL+NXBOT+NX
        DXG(I)=DX_OCT
      end do
      DX1T=XTOP/TOFLOAT(NXTOP)  ! "top" region
      DO I=NXL+NXBOT+NX+1,NXL+NXBOT+NX+NXTOP
        DXG(I)=DX1T
      end do

      XWIDTH_BIG=SUM(DXG)
      XG(1)=-CNST(0.5)*XWIDTH_BIG+CNST(0.5)*DXG(1)
      DO I=2,NXTOT-1
        XG(I)=XG(I-1)+CNST(0.5)*(DXG(I-1)+DXG(I))
      end do
      XG(NXTOT)=CNST(0.5)*XWIDTH_BIG-CNST(0.5)*DXG(NXTOT)

      XTEST=XG(NXTOT)
      !    write(6,*)' xg nxtot ',nxtot
      !    write(6,101)(xg(k),k=1,nxtot)
      !    write(6,*)' dxg '
      !    write(6,101)(dxg(k),k=1,nxtot)
      ! Y-mesh
      NYTOT_0=NY+NY2  ! excluding border points
      NYTOT=NYTOT_0+2*NYL  ! all Poisson Y mesh point
      allocate(yg(nytot))
      allocate(dyg(nytot))
      DO I=1,NYL  !  border points
        DYG(I)=DYL(I)
      end do
      DO I=NYL+NYTOT_0+1,NYTOT
        DYG(I)=DYL(NYTOT-I+1)
      end do
      ! outside Octopus mesh
      YTOP=CNST(0.5)*YWIDTH-YCEN-CNST(0.5)*YL
      YBOT=CNST(0.5)*YWIDTH+YCEN-CNST(0.5)*YL
      NYTOP=INT(NY2*YTOP/(YTOP+YBOT))
      NYBOT=NY2-NYTOP
      DY1B=YBOT/TOFLOAT(NYBOT) ! "bottom" region
      DO I=NYL+1,NYL+NYBOT
        DYG(I)=DY1B
      end do
      DY_OCT=YL/TOFLOAT(NY)  ! Octopus region
      DO I=NYL+NYBOT+1,NYL+NYBOT+NY
        DYG(I)=DY_OCT
      end do
      DY1T=YTOP/TOFLOAT(NYTOP)  ! "top" region
      DO I=NYL+NYBOT+NY+1,NYL+NYBOT+NY+NYTOP
        DYG(I)=DY1T
      end do

      YWIDTH_BIG=SUM(DYG)
      YG(1)=-CNST(0.5)*YWIDTH_BIG+CNST(0.5)*DYG(1)
      DO I=2,NYTOT-1
        YG(I)=YG(I-1)+CNST(0.5)*(DYG(I-1)+DYG(I))
      end do
      YG(NYTOT)=CNST(0.5)*YWIDTH_BIG-CNST(0.5)*DYG(NYTOT)


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
          end do
        end do
      end do

!!$        dxg=1.0; dyg=1.0; dz=1.0
!!$        zg(1)=dz(1)*CNST(0.5)
!!$        do i=2,nztot-1
!!$           zg(i)=zg(i-1)+CNST(0.5)*(dz(i)+dz(i-1))
!!$        enddo
!!$        zg(nztot)=zg(nztot-1)+CNST(0.5)*dz(nztot)
!!$
!!$        xg(1)=dxg(1)*CNST(0.5)
!!$        do i=2,nxtot-1
!!$           xg(i)=xg(i-1)+CNST(0.5)*(dxg(i)+dxg(i-1))
!!$        enddo
!!$        xg(nxtot)=xg(nxtot-1)+CNST(0.5)*dxg(nxtot)
!!$
!!$        yg(1)=dyg(1)*CNST(0.5)
!!$        do i=2,nytot-1
!!$           yg(i)=yg(i-1)+CNST(0.5)*(dyg(i)+dyg(i-1))
!!$        enddo
!!$        yg(nytot)=yg(nytot-1)+CNST(0.5)*dyg(nytot)
!!$           

      !    write(6,*)' yg '
      !    write(6,101)(yg(k),k=1,nytot)
      !    write(6,*)' dyg '
      !    write(6,101)(dyg(k),k=1,nytot)
    end subroutine xyzgrid

  end subroutine poisson_sete_init

  !---------------------------------------

  subroutine poisson_sete_solve(this, icase, rho, vh, nx, ny, nz, xl, yl, zl, icalc)
    type(poisson_sete_t), intent(in)    :: this
    integer,              intent(in)    :: icase
    FLOAT,                intent(in)    :: rho(:, :, :)
    FLOAT,                intent(inout) :: vh(:, :, :)
    integer,              intent(in)    :: nx
    integer,              intent(in)    :: ny
    integer,              intent(in)    :: nz
    FLOAT,                intent(in)    :: xl
    FLOAT,                intent(in)    :: yl
    FLOAT,                intent(in)    :: zl
    integer,              intent(in)    :: icalc

    ! Input: NX, NY, NZ -> No. of octopus grid points in each dimension
    ! Input: XL, YL, ZL -> Size of the octopus box

    FLOAT        :: bohrnm, angsnm, err
    integer      :: j1, k1, m, i, j, k, i1,idum
    FLOAT :: VHMIN, VHMAX

    integer, allocatable :: iwork(:)
    FLOAT,   allocatable :: rhotest(:,:,:), rwork(:), q2(:), x2(:)

    PCONST=CNST(4.0)*M_PI ! need 4 pi for Hartrees I think (??)
    BOHRNM=BOHR/CNST(10.0)
    ANGSNM=CNST(10.0)/BOHR

    SAFE_ALLOCATE(q2(1:NTOT))
    SAFE_ALLOCATE(rhotest(1:nxtot, 1:nytot, 1:nztot))

    write(*,*) " Count atoms", count_atoms

    if (count_atoms <= noatoms) then
      idum = count_atoms

      SAFE_ALLOCATE(X2(1:NTOT))
      
      do i = 1, ntot
        call quickrnd(IDUM,X2(I))
        if (x2(i) > M_ZERO) then
            x2(i) = - x2(i)
        endif
      enddo

    else if (count_atoms == noatoms+1) then

      SAFE_ALLOCATE(X2(1:NTOT))

      x2(1:NTOT) = M_ZERO

      idum = count_atoms + 1

      do i = 1, ntot
        x2(i) = M_ZERO

        call quickrnd(IDUM,X2(I))

        if (x2(i) < M_ZERO) then
            x2(i) = -x2(i)
        endif

      enddo

    endif

    rhotest = M_ZERO

    do i = 1, ntot
      q2(i) = this%qs(i) ! recall stored BC info
    enddo

    DO M = 1,NTOT
      I=IY(M); J=JY(M); K=KY(M)
      IF(this%ipio(I,J,K) == 0)THEN
        IF((I.GT.NXL+NXBOT.AND.I.LE.NXL+NXBOT+NX).AND. &
          (J.GT.NYL+NYBOT.AND.J.LE.NYL+NYBOT+NY).AND. &
          (K.GT.NZBOT.AND.K.LE.NZBOT+NZ))THEN
          I1=I-(NXL+NXBOT)
          J1=J-(NYL+NYBOT)
          K1=K-(NZBOT)
          IF(I1.GT.0.AND.J1.GT.0.AND.K1.GT.0)THEN ! this is redundant
            Q2(M)=Q2(M)+DYG(J)*DXG(I)*RHO(I1,J1,K1)
            rhotest(i,j,k)=rho(i1,j1,k1)
          end if
        end if
      ELSE
        Q2(M)= M_ZERO
      end if
    end do
    do i=1,NTOT
      Q2(i)=Q2(i)/this%ADIAG(i)!! Laplacian scaled with diagonal elements
    enddo

    allocate(rwork(lenw))
    allocate(iwork(leniw))

    call dslucs(NTOT,Q2,X2,NELT,this%IAD,this%JAD,this%AW,ISYM,ITOL, &
      TOL,ITMAX,ITER,ITERMIN,ERR,IERR,IUNIT,RWORK,LENW, &
      IWORK,LENIW)

    leniw = iwork(9)
    lenw  = iwork(10) ! new work array lengths

    deallocate(rwork)
    deallocate(iwork)

    DO M=1,NTOT ! store X2(M) i,j,k wise
      I=IY(M); J=JY(M); K=KY(M)
      VH_BIG(I,J,K)=X2(M)
    end do
    DO M = 1,NTOT
      I=IY(M); J=JY(M); K=KY(M)
      IF(this%ipio(I,J,K) == 0)THEN
        IF((I.GT.NXL+NXBOT.AND.I.LE.NXL+NXBOT+NX).AND. &
          (J.GT.NYL+NYBOT.AND.J.LE.NYL+NYBOT+NY).AND. &
          (K.GT.NZBOT.AND.K.LE.NZBOT+NZ))THEN
          I1=I-(NXL+NXBOT)
          J1=J-(NYL+NYBOT)
          K1=K-(NZBOT)
          VH(I1,J1,K1)=X2(M)
        end if
      end if
    end do

    call egate(this)

    write(358,*) "#x,z, vh_big(x,11,k)"
    do I = 1, NXTOT
      do K = 1, NZTOT
        write(358,*) I,K,VH_BIG(I,11,K)
      end do
    end do

    VHMIN=MINVAL(VH)
    VHMAX=MAXVAL(VH)

    deallocate(Q2)
    deallocate(rhotest)
    if (count_atoms <= noatoms ) then
      deallocate(x2)
    endif
    !deallocate(VH_BIG);

  end subroutine poisson_sete_solve

  !------------------------------------------------------------

  subroutine poisson_sete_end(this)
    type(poisson_sete_t), intent(inout) :: this

    SAFE_DEALLOCATE_P(this%ipio)
    SAFE_DEALLOCATE_P(this%vbound)
    SAFE_DEALLOCATE_P(this%qs)
    SAFE_DEALLOCATE_P(this%iad)
    SAFE_DEALLOCATE_P(this%JAD)

    deallocate(this%adiag)
!    if (count_atoms > noatoms ) then
!     deallocate(X2)
!    endif
    deallocate(DXL);deallocate(DYL)
    deallocate(dielectric)
    deallocate(zg)
    deallocate(xg)
    deallocate(yg)
    deallocate(dz)
    deallocate(dyg)
    deallocate(dxg)
    deallocate(ky)
    deallocate(iy)
    deallocate(jy)
    deallocate(VH_BIG)
    deallocate(this%aw)
  end subroutine poisson_sete_end

  !--------------------------------------------

  subroutine egate(this)
    type(poisson_sete_t), intent(in)    :: this

    integer :: i, j, k
    FLOAT, allocatable :: sig(:, :, :, :)
    FLOAT :: cs1, cs2, cs3

    !
    !     to do energy calculations previously performed in setr
    !  compute capacitances of six surfaces
    
    SAFE_ALLOCATE(sig(1:nxtot, 1:nytot, 1:nztot, 1:6))

    call cap(this, sig)

    ESURF = M_ZERO
    CHARGE_SURF = M_ZERO
    CHARGE_TOP = M_ZERO

    CHARGE_BOT = M_ZERO
    CS1 = M_ZERO
    CS2 = M_ZERO
    CS3 = M_ZERO
    !write(*,*) "Begin cycle ", NXTOT, NYTOT, NZTOT
    DO I=1,NXTOT
      DO J=1,NYTOT
        DO K=1,NZTOT ! N.B. [SIG]=charge (not charge/Area)
          ESURF=ESURF+this%VBOUND(I-1,J,K)*SIG(I,J,K,1)
          ESURF=ESURF+this%VBOUND(I+1,J,K)*SIG(I,J,K,2)
          ESURF=ESURF+this%VBOUND(I,J-1,K)*SIG(I,J,K,3)
          ESURF=ESURF+this%VBOUND(I,J+1,K)*SIG(I,J,K,4)
          ESURF=ESURF+this%VBOUND(I,J,K-1)*SIG(I,J,K,5)
          ESURF=ESURF+this%VBOUND(I,J,K+1)*SIG(I,J,K,6)
          CS1=CS1+SIG(I,J,K,1)+SIG(I,J,K,2)
          CS2=CS2+SIG(I,J,K,3)+SIG(I,J,K,4)
          CS3=CS3+SIG(I,J,K,5)+SIG(I,J,K,6)
        end do
        CHARGE_TOP=CHARGE_TOP+SIG(I,J,1,5)+SIG(I,J,1,6)
        CHARGE_BOT=CHARGE_BOT+SIG(I,J,NZTOT,5)+SIG(I,J,NZTOT,6)
      end do
    end do

    SAFE_DEALLOCATE_A(sig)

    ESURF=CNST(0.5)*ESURF
    CHARGE_SURF=CS1+CS2+CS3

    WRITE(520,*)'Esurf', ESURF*27.2
    write(521,*)' CS1 ',CS1,' CS2 ',CS2,' CS3 ',CS3

  contains

    subroutine cap(this, sig)
      type(poisson_sete_t), intent(in)    :: this
      FLOAT,                intent(out)   :: sig(:, :, :, :)

      integer :: i, j, k

      !  9/24/08 - this VERSION OF CAP is very generic. It is
      !  for surfER. Assumes just this%ipio(0:NX+1,0:NY+1,0:NZ+1), 
      !  this%VBOUND(0:NX+1,0:NY+1,0:NZ+1),dielectric(NX,NY,NZ).
      !  Outputs a single 3D array SIG(NX,NY,NZ)
      !
      !  for calculating capacitances
      !
      !  SIG(I,J,K) is capacitance
      !
      !  input: NX,NY,NZ
      !         this%ipio,this%VBOUND
      !         ZG,DXG,DYG,DZ
      !         VH_BIG,dielectric,PCONST
      !
      !  output: SIG(I,J,K,6) 6 faces of grid cube 
      !          used in ETOT
      !  1 I-1, 2 I+1, 3 J-1, 4 J+1, 5 K-1, 6 K+1
      !
      !  local: BORDER,PCONST
      !
      !  a homogeneous grid approximation is made: derivative
      !  at metal-semi border calculated discretely with 1/DXG (1/DYG,1/DZ)
      !  rather than properly averaged spacings.

      SIG= M_ZERO

      DO I=1,NXTOT 
        DO J=1,NYTOT
          DO K=1,NZTOT
            IF(this%ipio(I,J,K) == 0)THEN ! I,J,K is in Poisson grid
              IF(this%ipio(I-1,J,K) == 1)THEN
                BORDER=2*(VH_BIG(I,J,K)-this%VBOUND(I-1,J,K))/DXG(I)
                SIG(I,J,K,1)=SIG(I,J,K,1)-BORDER*DYG(J)*DZ(K)* &
                  dielectric(I,J,K)
              end if
              IF(this%ipio(I+1,J,K) == 1)THEN
                BORDER=2*(VH_BIG(I,J,K)-this%VBOUND(I+1,J,K))/DXG(I)
                SIG(I,J,K,2)=SIG(I,J,K,2)-BORDER*DYG(J)*DZ(K)* &
                  dielectric(I,J,K)
              end if
              IF(this%ipio(I,J-1,K) == 1)THEN
                BORDER=2*(VH_BIG(I,J,K)-this%VBOUND(I,J-1,K))/DYG(J)
                SIG(I,J,K,3)=SIG(I,J,K,3)-BORDER*DXG(I)*DZ(K)* &
                  dielectric(I,J,K)
              end if
              IF(this%ipio(I,J+1,K) == 1)THEN
                BORDER=2*(VH_BIG(I,J,K)-this%VBOUND(I,J+1,K))/DYG(J)
                SIG(I,J,K,4)=SIG(I,J,K,4)-BORDER*DXG(I)*DZ(K)* &
                  dielectric(I,J,K)
              end if
              IF(this%ipio(I,J,K-1) == 1)THEN
                BORDER=2*(VH_BIG(I,J,K)-this%VBOUND(I,J,K-1))/DZ(K)
                SIG(I,J,K,5)=SIG(I,J,K,5)-BORDER*DXG(I)*DYG(J)* &
                  dielectric(I,J,K)
              end if
              IF(this%ipio(I,J,K+1) == 1)THEN
                BORDER=2*(VH_BIG(I,J,K)-this%VBOUND(I,J,K+1))/DZ(K)
                SIG(I,J,K,6)=SIG(I,J,K,6)-BORDER*DXG(I)*DYG(J)* &
                  dielectric(I,J,K)
              end if
            end if
          end do
        end do
      end do

    end subroutine cap


  end subroutine egate

  FLOAT function poisson_sete_energy(this) result(energy)
    type(poisson_sete_t), intent(in)    :: this

    energy = esurf
  end function poisson_sete_energy

end module poisson_sete_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

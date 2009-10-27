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
    integer, pointer :: ipio(:,:,:)
    FLOAT,   pointer :: adiag(:)
    FLOAT,   pointer :: qs(:)
    FLOAT,   pointer :: vbound(:,:,:)
    integer, pointer :: iad(:)
    integer, pointer :: jad(:)
    FLOAT,   pointer :: aw(:)
    FLOAT,   pointer :: dielectric(:,:,:)
    FLOAT,   pointer :: vh_big(:,:,:)
    integer, pointer :: iy(:)
    integer, pointer :: jy(:)
    integer, pointer :: ky(:)
    FLOAT            :: esurf
    integer          :: nxtot
    integer          :: nytot
    integer          :: nztot
    integer          :: nxbot
    integer          :: nybot
    integer          :: nzbot
    integer          :: noatoms
    integer          :: nelt
    integer          :: ntot
    integer          :: counter 
    FLOAT            :: tot_nuc_charge_energy
  end type poisson_sete_t

  integer :: lenw, leniw, idev, isym = 0, itol = 2, itmax = 2001, itermin = 5, iter, ierr, iunit = 0, nxl, nyl
  FLOAT, allocatable :: xg(:), yg(:), zg(:), dxg(:), dyg(:), dz(:), dxl(:), dyl(:)

  ! these variables must be local for the moment
  FLOAT, allocatable :: rho_nuc(:)
  integer :: count_atoms

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

    integer           :: isete_on
    character(len=40) :: cdum, fil1
    FLOAT :: bohrnm, angsnm
    integer :: i, j, ngates
    FLOAT :: zcen, xcen, ycen, dielectric0
    FLOAT, allocatable :: vtv(:), vt(:), q2(:)
    integer :: nx2, ny2, nz2, md
    FLOAT :: xwidth, ywidth, zwidth, pconst

    pconst = CNST(4.0)*M_PI ! need 4 pi for Hartrees I think (??)
    bohrnm = BOHR/CNST(10.0)
    angsnm = CNST(10.0)/BOHR
    this%noatoms = number_atoms
    count_atoms = 0
    this%counter = 0
    this%tot_nuc_charge_energy=0.0

    open(56, file = 'poisq', status = 'unknown')
    FIL1='test1.pois'
    OPEN(unit=57, FILE = FIL1, STATUS = 'UNKNOWN') ! status info file

    ! READ FILE
    read(57, '(a40)') cdum
    read(57, '(a40)') cdum
    read(57, '(a40)') cdum
    read(57,*) idev, isete_on ! for future use - device configuration
    write(*,*) "ISETE is ON", isete_on

    if(idev == 1) then

      ! physical device characteristics
      read(57,*) zwidth

      zwidth = zwidth*angsnm  ! distance between plates (nm)

      read(57,*) zcen
      zcen = zcen*angsnm

      read(57,*) ngates

      SAFE_ALLOCATE(VTV(1:NGATES))
      SAFE_ALLOCATE(VT(1:NGATES))

      READ(57,*)VTV(1),VTV(2) ! voltage on plates 
      VT=VTV/HARTREE

      SAFE_DEALLOCATE_A(VTV)

      READ(57,*) dielectric0    ! dielectric constant of region betw plates
      dielectric0 = dielectric0/pconst

      ! mesh related
      READ(57,*)XWIDTH,YWIDTH ! bounding box around octopus box
      xwidth=xwidth*angsnm
      ywidth=ywidth*angsnm
      READ(57,*)xcen,ycen
      xcen=xcen*angsnm
      ycen=ycen*angsnm
      READ(57,*)nx2,ny2,nz2   ! additional mesh points in bd box
      READ(57,*)nxl,nyl       ! padding points (to sim zone edge)
      allocate(dxl(nxl)); allocate(dyl(nyl))
      READ(57,*)(dxl(i),i=1,nxl)
      READ(57,*)(dyl(j),j=1,nyl)
    ELSE
      WRITE(6,*)' IDEV not supported '
      RETURN
    end if

    call xyzgrid(this, nx, ny, nz, xl, yl, zl, xcen, ycen, zcen, nx2, ny2, nz2)  ! establish x-y-z grids

    call cbsurf(this, vt) !Assign boundaries 

    this%nelt=32+20*((this%nxtot-2)+(this%nytot-2)+(this%nztot-2))+ &
      12*((this%nxtot-2)*(this%nytot-2)+(this%nxtot-2)*(this%nztot-2)+(this%nytot-2)*(this%nztot-2)) + &
      7*(this%nxtot-2)*(this%nytot-2)*(this%nztot-2)
    lenw=this%nelt+9*this%nxtot*this%nytot*this%nztot
    leniw=this%nelt+5*this%nxtot*this%nytot*this%nztot+12

    SAFE_ALLOCATE(q2(1:this%ntot))
    SAFE_ALLOCATE(this%aw(1:this%nelt))
    SAFE_ALLOCATE(this%qs(1:this%ntot))
    SAFE_ALLOCATE(this%adiag(1:this%ntot))
    SAFE_ALLOCATE(this%iad(1:this%nelt))
    SAFE_ALLOCATE(this%jad(1:this%nelt))

    call poissonm(this)

    this%qs = q2 ! store Dirichlet BC info in this%QS

    SAFE_DEALLOCATE_A(q2)

  contains

    subroutine poissonm(this)
      type(poisson_sete_t), intent(inout) :: this

      ! version g created 03/25/08 - modified to include 
      !                              this%dielectric(0:THIS%NXTOT+1,0:THIS%NYTOT+1,0:THIS%NZTOT+1)
      !           based on ~/LEVEL1/photovolt/poissonpv_v2.f90

      ! creates Laplacian.
      ! inputs:  THIS%NTOT,MD,THIS%NXTOT,THIS%NYTOT,THIS%NZTOT,this%ipio(0:THIS%NXTOT+1,0:THIS%NYTOT+1,0:THIS%NZTOT+1),
      !          this%VBOUND(0:THIS%NXTOT+1,0:THIS%NYTOT+1,0:THIS%NZTOT+1),Q2(THIS%NTOT),
      !          DXG(THIS%NXTOT),DYG(THIS%NYTOT),DZ(THIS%NZTOT),DIELELCTRIC(0:THIS%NXTOT+1,0:THIS%NYTOT+1,0:THIS%NZTOT+1)
      !
      ! outputs: this%AW(THIS%NELT),this%ADIAG(THIS%NTOT),Q2,
      !          THIS%NELT, IAD(THIS%NELT),this%JAD(THIS%NELT)
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
      JOFF(2)=-THIS%NZTOT; JOFF(3)=-1; JOFF(4)=0; JOFF(5)=1; JOFF(6)=THIS%NZTOT; JOFF(7)=MD
      I=1; J=1; K=1
      IEC=1 ! counter for number of non-zero entries in this%AW.

      q2(1:this%ntot) = M_ZERO

      DO IA=1,THIS%NTOT
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
                DCONST=this%dielectric(I,J,K)
                AV(IAV)=M_ZERO; AV(4)=AV(4)+DCONST*D1*D2/DEL
                Q2(IA)=Q2(IA)+this%VBOUND(I-1,J,K)*DCONST*D1*D2/DEL
              CASE(2) ! Neumann boundary point
                AV(IAV)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT ! interior point
                DEL=CNST(0.5)*(DXG(I)+DXG(I-1))*DZ(K)
                DCONST=CNST(0.5)*(this%dielectric(I-1,J,K)+this%dielectric(I,J,K))
                AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
              END SELECT
            CASE(2) ! I,J-1,K
              D1=DZ(K); D2=DXG(I)
              SELECT CASE(this%ipio(I,J-1,K))
              CASE(1) ! Dirichlet
                DEL=CNST(0.5)*DYG(J)*DZ(K)
                DCONST=this%dielectric(I,J,K)
                AV(IAV)=M_ZERO; AV(4)=AV(4)+DCONST*D1*D2/DEL
                Q2(IA)=Q2(IA)+this%VBOUND(I,J-1,K)*DCONST*D1*D2/DEL
              CASE(2) ! Neumann
                AV(IAV)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT! interior point
                DEL=CNST(0.5)*(DYG(J)+DYG(J-1))*DZ(K)
                DCONST=CNST(0.5)*(this%dielectric(I,J,K)+this%dielectric(I,J-1,K))
                AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
              END SELECT
            CASE(3) ! I,J,K-1
              D1=DXG(I); D2=DYG(J)
              SELECT CASE(this%ipio(I,J,K-1))
              CASE(1) ! Dirichlet
                DEL=CNST(0.5)*DZ(K)*DZ(K)
                DCONST=this%dielectric(I,J,K)
                AV(IAV)=M_ZERO; AV(4)=AV(4)+DCONST*D1*D2/DEL
                Q2(IA)=Q2(IA)+this%VBOUND(I,J,K-1)*DCONST*D1*D2/DEL
              CASE(2) ! Neumann
                AV(IAV)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT! interior point
                DEL=CNST(0.5)*(DZ(K)+DZ(K-1))*DZ(K)
                DCONST=CNST(0.5)*(this%dielectric(I,J,K)+this%dielectric(I,J,K-1))
                AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
              END SELECT
            CASE(5) ! I,J,K+1
              D1=DXG(I); D2=DYG(J)
              SELECT CASE(this%ipio(I,J,K+1))
              CASE(1) ! Dirichlet
                DEL=CNST(0.5)*DZ(K)*DZ(K)
                DCONST=this%dielectric(I,J,K)
                AV(IAV)=M_ZERO; AV(4)=AV(4)+DCONST*D1*D2/DEL
                Q2(IA)=Q2(IA)+this%VBOUND(I,J,K+1)*DCONST*D1*D2/DEL
              CASE(2) ! Neumann
                AV(IAV)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT! interior point
                DEL=CNST(0.5)*(DZ(K)+DZ(K+1))*DZ(K)
                DCONST=CNST(0.5)*(this%dielectric(I,J,K)+this%dielectric(I,J,K+1))
                AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
              END SELECT
            CASE(6) ! I,J+1,K
              D1=DZ(K); D2=DXG(I)
              SELECT CASE(this%ipio(I,J+1,K))
              CASE(1) ! Dirichlet
                DEL=CNST(0.5)*DYG(J)*DZ(K)
                DCONST=this%dielectric(I,J,K)
                AV(IAV)=M_ZERO; AV(4)=AV(4)+DCONST*D1*D2/DEL
                Q2(IA)=Q2(IA)+this%VBOUND(I,J+1,K)*DCONST*D1*D2/DEL
              CASE(2) ! Neumann
                AV(IAV)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT
                DEL=CNST(0.5)*(DYG(J)+DYG(J+1))*DZ(K)
                DCONST=CNST(0.5)*(this%dielectric(I,J,K)+this%dielectric(I,J+1,K))
                AV(IAV)=-DCONST*D1*D2/DEL; AV(4)=AV(4)+DCONST*D1*D2/DEL
              END SELECT
            CASE(7) ! I+1,J,K
              D1=DZ(K); D2=DYG(J)
              SELECT CASE(this%ipio(I+1,J,K))
              CASE(1) ! Dirichlet
                DEL=CNST(0.5)*DXG(I)*DZ(K)
                DCONST=this%dielectric(I,J,K)
                AV(IAV)=M_ZERO; AV(4)=AV(4)+DCONST*D1*D2/DEL
                Q2(IA)=Q2(IA)+this%VBOUND(I+1,J,K)*DCONST*D1*D2/DEL
              CASE(2) ! Neumann
                AV(IAV)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT
                DEL=CNST(0.5)*(DXG(I)+DXG(I+1))*DZ(K)
                DCONST=CNST(0.5)*(this%dielectric(I,J,K)+this%dielectric(I+1,J,K))
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
        IF(K == THIS%NZTOT+1)THEN
          K=1
          J=J+1
          IF(J == THIS%NYTOT+1)THEN
            J=1
            I=I+1
          end if
        end if


      end do
      if(idebug == 1)then
        close(57)
      endif

      IEC=IEC-1 ! overcounted non-zero elements by one.
      write(6,*)' this%nelt, iec ',this%nelt,iec
      this%nelt = iec

      ! in case THIS%NELT in static differs from that counted here
!!$        allocate(this%AWP(THIS%NELT))
!!$        DO I=1,THIS%NELT
!!$           this%AWP(I)=this%AW(I)
!!$        end do
!!$        deallocate(this%AW)
!!$        allocate(this%AW(THIS%NELT))
!!$        this%AW=this%AWP
!!$        deallocate(this%AWP)

      if(idebug == 1)then
        close(57)
      endif
    end subroutine poissonm

    !---------------------------------------------------------

    subroutine cbsurf(this, vt)
      type(poisson_sete_t), intent(inout) :: this
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
      allocate(this%ipio(0:THIS%NXTOT+1,0:THIS%NYTOT+1,0:THIS%NZTOT+1))
      allocate(THIS%VH_BIG(1:THIS%NXTOT,1:THIS%NYTOT,1:THIS%NZTOT))
      this%vh_big(1:THIS%NXTOT,1:THIS%NYTOT,1:THIS%NZTOT)= M_ZERO
      allocate(this%VBOUND(0:THIS%NXTOT+1,0:THIS%NYTOT+1,0:THIS%NZTOT+1))
      allocate(this%dielectric(0:THIS%NXTOT+1,0:THIS%NYTOT+1,0:THIS%NZTOT+1))
      this%ipio=0; this%VBOUND= M_ZERO
      IF(IDEV == 1)THEN ! parallel plates
        this%dielectric=dielectric0
        this%ipio(:,:,0)=1
        this%ipio(:,:,THIS%NZTOT+1)=1
        this%VBOUND(:,:,0)=VT(1)
        this%VBOUND(:,:,THIS%NZTOT+1)=VT(2)
        DO K=1,THIS%NZTOT
          IF(K /= 0 .AND. K /= THIS%NZTOT+1)THEN !Redundant?
            ! lateral BC`s
            this%ipio(0,:,K)=2
            this%ipio(THIS%NXTOT+1,:,K)=2
            this%ipio(:,0,K)=2
            this%ipio(:,THIS%NYTOT+1,K)=2
          end if
        end do
      end if

    end subroutine cbsurf

    !---------------------------------------------------------

    subroutine xyzgrid(this, nx, ny, nz, xl, yl, zl, xcen, ycen, zcen, nx2, ny2, nz2)
      type(poisson_sete_t), intent(inout) :: this
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
      !           hence NZTOT_0 = THIS%NZTOT (ie there IS no NZTOT_0)
      this%nztot = nz + nz2

      allocate(zg(this%nztot))
      allocate(dz(this%nztot))

      !top and bottom z coordinates.
      !Depend on POISSON_SETE (ZWIDTH) 
      !and Octopus (ZL)parameters.
      ztop = CNST(0.5)*zwidth - zcen - CNST(0.5)*zl
      zbot = CNST(0.5)*zwidth + zcen - CNST(0.5)*zl

      nztop = int(nz2*ztop/(ztop + zbot))
      this%nzbot = nz2 - nztop
      dz1b=zbot/TOFLOAT(this%nzbot)
      DO K=1,this%nzbot
        dz(k)=dz1b
      end do
      write(67,*) zl
      write(67,*) nz, (nz-1)/2
      dz_oct=(zl/TOFLOAT((nz-1)))

      do k = this%nzbot + 1, this%nzbot + nz
        dz(k) = dz_oct
        write(67,*) dz(k)
      end do

      dz1t=ztop/TOFLOAT(nztop)
      DO k=this%nzbot+nz+1,this%nztot
        dz(k)=dz1t
        write(67,*) dz(k)
      end do

      zg(1)=-CNST(0.5)*zwidth+CNST(0.5)*dz(1)
      write(67,*) zg(1)
      DO k=2,this%nztot-1
        zg(k)=zg(k-1)+CNST(0.5)*(dz(k-1)+dz(k))
        write(67,*) zg(k)
      end do
      zg(this%nztot)=CNST(0.5)*zwidth-CNST(0.5)*dz(this%nztot)
      write(67,*) zg(this%nztot)

      ZTEST=ZG(THIS%NZTOT)
          write(67,*)' this%nzbot ',this%nzbot
          write(67,*)' zg this%nztot ',this%nztot
          write(67,101)(zg(k),k=1,this%nztot)
       101 format(10(1x,e10.4))
          write(67,*)' dz '
          write(67,101)(dz(k),k=1,this%nztot)

      ! x-mesh
      NXTOT_0=NX+NX2
      THIS%NXTOT=NXTOT_0+2*NXL
      allocate(xg(this%nxtot))
      allocate(DXG(THIS%NXTOT))
      DO I=1,NXL  !  border points
        DXG(I)=DXL(I)
      end do
      DO I=NXL+NXTOT_0+1,THIS%NXTOT
        DXG(I)=DXL(THIS%NXTOT-I+1)
      end do
      ! outside Octopus mesh
      XTOP=CNST(0.5)*XWIDTH-XCEN-CNST(0.5)*XL
      XBOT=CNST(0.5)*XWIDTH+XCEN-CNST(0.5)*XL
      NXTOP=INT(NX2*XTOP/(XTOP+XBOT))
      THIS%NXBOT=NX2-NXTOP
      DX1B=XBOT/TOFLOAT(THIS%NXBOT) ! "bottom" region
      DO I=NXL+1,NXL+THIS%NXBOT
        DXG(I)=DX1B
      end do
      DX_OCT=XL/TOFLOAT(NX)  ! Octopus region
      DO I=NXL+THIS%NXBOT+1,NXL+THIS%NXBOT+NX
        DXG(I)=DX_OCT
      end do
      DX1T=XTOP/TOFLOAT(NXTOP)  ! "top" region
      DO I=NXL+THIS%NXBOT+NX+1,NXL+THIS%NXBOT+NX+NXTOP
        DXG(I)=DX1T
      end do

      xwidth_big=SUM(dxg)
      xg(1)=-CNST(0.5)*xwidth_big+CNST(0.5)*dxg(1)
      DO I=2,THIS%NXTOT-1
        XG(I)=XG(I-1)+CNST(0.5)*(DXG(I-1)+DXG(I))
      end do
      XG(THIS%NXTOT)=CNST(0.5)*XWIDTH_BIG-CNST(0.5)*DXG(THIS%NXTOT)

      XTEST=XG(THIS%NXTOT)
          write(57,*)' xg this%nxtot ',this%nxtot
          write(57,101)(xg(k),k=1,this%nxtot)
          write(57,*)' dxg '
          write(57,101)(dxg(k),k=1,this%nxtot)
      ! Y-mesh
      NYTOT_0=NY+NY2  ! excluding border points
      THIS%NYTOT=NYTOT_0+2*NYL  ! all Poisson Y mesh point
      allocate(yg(this%nytot))
      allocate(dyg(this%nytot))
      DO I=1,NYL  !  border points
        DYG(I)=DYL(I)
      end do
      DO I=NYL+NYTOT_0+1,THIS%NYTOT
        DYG(I)=DYL(THIS%NYTOT-I+1)
      end do
      ! outside Octopus mesh
      YTOP=CNST(0.5)*YWIDTH-YCEN-CNST(0.5)*YL
      YBOT=CNST(0.5)*YWIDTH+YCEN-CNST(0.5)*YL
      NYTOP=INT(NY2*YTOP/(YTOP+YBOT))
      THIS%NYBOT=NY2-NYTOP
      DY1B=YBOT/TOFLOAT(THIS%NYBOT) ! "bottom" region
      DO I=NYL+1,NYL+THIS%NYBOT
        DYG(I)=DY1B
      end do
      DY_OCT=YL/TOFLOAT(NY)  ! Octopus region
      DO I=NYL+THIS%NYBOT+1,NYL+THIS%NYBOT+NY
        DYG(I)=DY_OCT
      end do
      DY1T=YTOP/TOFLOAT(NYTOP)  ! "top" region
      DO I=NYL+THIS%NYBOT+NY+1,NYL+THIS%NYBOT+NY+NYTOP
        DYG(I)=DY1T
      end do

      YWIDTH_BIG=SUM(DYG)
      YG(1)=-CNST(0.5)*YWIDTH_BIG+CNST(0.5)*DYG(1)
      DO I=2,THIS%NYTOT-1
        YG(I)=YG(I-1)+CNST(0.5)*(DYG(I-1)+DYG(I))
      end do
      YG(THIS%NYTOT)=CNST(0.5)*YWIDTH_BIG-CNST(0.5)*DYG(THIS%NYTOT)
          write(57,*)' yg this%nytot ',this%nytot
          write(57,101)(yg(k),k=1,this%nytot)
          write(57,*)' dyg '
          write(57,101)(dyg(k),k=1,this%nytot)


      THIS%NTOT=THIS%NXTOT*THIS%NYTOT*THIS%NZTOT
      MD=THIS%NZTOT*THIS%NYTOT

      allocate(this%iy(this%ntot))
      allocate(this%jy(this%ntot))
      allocate(this%ky(this%ntot))
      m=0
      DO I=1,THIS%NXTOT
        DO J=1,THIS%NYTOT
          DO K=1,THIS%NZTOT
            M=M+1
            THIS%IY(M)=I
            THIS%JY(M)=J
            THIS%KY(M)=K
          end do
        end do
      end do

!!$        dxg=1.0; dyg=1.0; dz=1.0
!!$        zg(1)=dz(1)*CNST(0.5)
!!$        do i=2,this%nztot-1
!!$           zg(i)=zg(i-1)+CNST(0.5)*(dz(i)+dz(i-1))
!!$        enddo
!!$        zg(this%nztot)=zg(this%nztot-1)+CNST(0.5)*dz(this%nztot)
!!$
!!$        xg(1)=dxg(1)*CNST(0.5)
!!$        do i=2,this%nxtot-1
!!$           xg(i)=xg(i-1)+CNST(0.5)*(dxg(i)+dxg(i-1))
!!$        enddo
!!$        xg(this%nxtot)=xg(this%nxtot-1)+CNST(0.5)*dxg(this%nxtot)
!!$
!!$        yg(1)=dyg(1)*CNST(0.5)
!!$        do i=2,this%nytot-1
!!$           yg(i)=yg(i-1)+CNST(0.5)*(dyg(i)+dyg(i-1))
!!$        enddo
!!$        yg(this%nytot)=yg(this%nytot-1)+CNST(0.5)*dyg(this%nytot)
!!$           

      !    write(6,*)' yg '
      !    write(6,101)(yg(k),k=1,this%nytot)
      !    write(6,*)' dyg '
      !    write(6,101)(dyg(k),k=1,this%nytot)
    end subroutine xyzgrid

  end subroutine poisson_sete_init

  !---------------------------------------

  subroutine poisson_sete_solve(this, icase, rho, vh, nx, ny, nz, xl, yl, zl, icalc)
    type(poisson_sete_t), intent(inout) :: this
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
    FLOAT :: pconst

    integer, allocatable :: iwork(:)
    FLOAT,   allocatable :: rhotest(:,:,:), rwork(:), q2(:), x2(:)

    PCONST=CNST(4.0)*M_PI ! need 4 pi for Hartrees I think (??)
    BOHRNM=BOHR/CNST(10.0)
    ANGSNM=CNST(10.0)/BOHR

    SAFE_ALLOCATE(q2(1:THIS%NTOT))
    SAFE_ALLOCATE(rhotest(1:this%nxtot, 1:this%nytot, 1:this%nztot))

    write(*,*) " Count atoms", count_atoms

    if (count_atoms <= this%noatoms) then
      idum = count_atoms

      SAFE_ALLOCATE(X2(1:THIS%NTOT))
      
      do i = 1, this%ntot
        call quickrnd(IDUM,X2(I))
        if (x2(i) > M_ZERO) then
            x2(i) = - x2(i)
        endif
      enddo

    else if (count_atoms == this%noatoms+1) then

      SAFE_ALLOCATE(X2(1:THIS%NTOT))

      x2(1:THIS%NTOT) = M_ZERO

      idum = count_atoms + 1

      do i = 1, this%ntot
        x2(i) = M_ZERO

        call quickrnd(idum, x2(i))

        if (x2(i) < M_ZERO) then
            x2(i) = -x2(i)
        endif

      enddo

    endif

    rhotest = M_ZERO

    do i = 1, this%ntot
      q2(i) = this%qs(i) ! recall stored BC info
    enddo

    DO M = 1,THIS%NTOT
      I=THIS%IY(M); J=THIS%JY(M); K=THIS%KY(M)
      IF(this%ipio(I,J,K) == 0)THEN
        IF((I.GT.NXL+THIS%NXBOT.AND.I.LE.NXL+THIS%NXBOT+NX).AND. &
          (J.GT.NYL+THIS%NYBOT.AND.J.LE.NYL+THIS%NYBOT+NY).AND. &
          (K.GT.THIS%NZBOT.AND.K.LE.THIS%NZBOT+NZ))THEN
          I1=I-(NXL+THIS%NXBOT)
          J1=J-(NYL+THIS%NYBOT)
          K1=K-(THIS%NZBOT)
          IF(I1.GT.0.AND.J1.GT.0.AND.K1.GT.0)THEN ! this is redundant
            Q2(M)=Q2(M)+DYG(J)*DXG(I)*RHO(I1,J1,K1)
            rhotest(i,j,k)=rho(i1,j1,k1)
          end if
        end if
      ELSE
        Q2(M)= M_ZERO
      end if
    end do
    do i=1,THIS%NTOT
      Q2(i)=Q2(i)/this%ADIAG(i)!! Laplacian scaled with diagonal elements
    enddo

    allocate(rwork(lenw))
    allocate(iwork(leniw))

    call dslucs(THIS%NTOT,Q2,X2,THIS%NELT,this%IAD,this%JAD,this%AW,ISYM,ITOL, &
      TOL,ITMAX,ITER,ITERMIN,ERR,IERR,IUNIT,RWORK,LENW, &
      IWORK,LENIW)
    if (ierr /= 0) then
	    write(*,*) "DSLUCS error number:", ierr, err
	    STOP
    endif
    leniw = iwork(9)
    lenw  = iwork(10) ! new work array lengths

    deallocate(rwork)
    deallocate(iwork)

    do m = 1, this%ntot ! store x2(m) i,j,k wise
      i = this%iy(m)
      j = this%jy(m)
      k = this%ky(m)
      this%vh_big(i, j, k) = x2(m)
    end do
    DO m = 1, this%ntot
      I=THIS%IY(M); J=THIS%JY(M); K=THIS%KY(M)
      IF(this%ipio(I,J,K) == 0)THEN
        IF((I.GT.NXL+THIS%NXBOT.AND.I.LE.NXL+THIS%NXBOT+NX).AND. &
          (J.GT.NYL+THIS%NYBOT.AND.J.LE.NYL+THIS%NYBOT+NY).AND. &
          (K.GT.THIS%NZBOT.AND.K.LE.THIS%NZBOT+NZ))THEN
          I1=I-(NXL+THIS%NXBOT)
          J1=J-(NYL+THIS%NYBOT)
          K1=K-(THIS%NZBOT)
          VH(I1,J1,K1)=X2(M)
        end if
      end if
    end do

    call egate(this)

!    write(358,*) "#x,z, this%vh_big(x,11,k)"
    do I = 1, THIS%NXTOT
      do K = 1, THIS%NZTOT
 !       write(358,*) I,K,THIS%VH_BIG(I,11,K)
      end do
    end do

    VHMIN=MINVAL(VH)
    VHMAX=MAXVAL(VH)

    deallocate(Q2)
    deallocate(rhotest)
    if (count_atoms <= this%noatoms ) then
      deallocate(x2)
    endif
    !deallocate(THIS%VH_BIG);

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
!    if (count_atoms > this%noatoms ) then
!     deallocate(X2)
!    endif
    deallocate(DXL);deallocate(DYL)
    deallocate(this%dielectric)
    deallocate(zg)
    deallocate(xg)
    deallocate(yg)
    deallocate(dz)
    deallocate(dyg)
    deallocate(dxg)
    deallocate(this%iy)
    deallocate(this%jy)
    deallocate(this%ky)
    deallocate(this%vh_big)
    deallocate(this%aw)
  end subroutine poisson_sete_end

  !--------------------------------------------

  subroutine egate(this)
    type(poisson_sete_t), intent(inout) :: this

    integer :: i, j, k
    FLOAT, allocatable :: sig(:, :, :, :) 
    FLOAT :: cs1, cs2, cs3,temporary
    FLOAT :: charge_top, charge_surf, charge_bot

    !
    !     to do energy calculations previously performed in setr
    !  compute capacitances of six surfaces
    
    SAFE_ALLOCATE(sig(1:this%nxtot, 1:this%nytot, 1:this%nztot, 1:6))

    call cap(this, sig)

    THIS%ESURF = M_ZERO
    CHARGE_SURF = M_ZERO
    CHARGE_TOP = M_ZERO

    CHARGE_BOT = M_ZERO
    CS1 = M_ZERO
    CS2 = M_ZERO
    CS3 = M_ZERO
    !write(*,*) "Begin cycle ", THIS%NXTOT, THIS%NYTOT, THIS%NZTOT
    DO I=1,THIS%NXTOT
      DO J=1,THIS%NYTOT
        DO K=1,THIS%NZTOT ! N.B. [SIG]=charge (not charge/Area)
          THIS%ESURF=THIS%ESURF+this%VBOUND(I-1,J,K)*SIG(I,J,K,1)
          THIS%ESURF=THIS%ESURF+this%VBOUND(I+1,J,K)*SIG(I,J,K,2)
          THIS%ESURF=THIS%ESURF+this%VBOUND(I,J-1,K)*SIG(I,J,K,3)
          THIS%ESURF=THIS%ESURF+this%VBOUND(I,J+1,K)*SIG(I,J,K,4)
          THIS%ESURF=THIS%ESURF+this%VBOUND(I,J,K-1)*SIG(I,J,K,5)
          THIS%ESURF=THIS%ESURF+this%VBOUND(I,J,K+1)*SIG(I,J,K,6)
          CS1=CS1+SIG(I,J,K,1)+SIG(I,J,K,2)
          CS2=CS2+SIG(I,J,K,3)+SIG(I,J,K,4)
          CS3=CS3+SIG(I,J,K,5)+SIG(I,J,K,6)
        end do
        CHARGE_TOP=CHARGE_TOP+SIG(I,J,1,5)+SIG(I,J,1,6)
        CHARGE_BOT=CHARGE_BOT+SIG(I,J,THIS%NZTOT,5)+SIG(I,J,THIS%NZTOT,6)
      end do
    end do

    SAFE_DEALLOCATE_A(sig)

      THIS%ESURF=CNST(0.5)*THIS%ESURF
      CHARGE_SURF=CS1+CS2+CS3
!    WRITE(520,*)THIS%ESURF*hartree*CNST(0.5)
!    write(521,*)' CS1 ',CS1,' CS2 ',CS2,' CS3 ',CS3
    write(*,*) "COUNTING ATOMS",this%counter+1
    if (this%counter<this%noatoms) then
      this%counter=this%counter+1
      this%tot_nuc_charge_energy=this%tot_nuc_charge_energy+THIS%ESURF
!      write(523,*)this%tot_nuc_charge_energy,THIS%ESURF
    else if (this%counter == this%noatoms) then
      temporary=THIS%ESURF
      THIS%ESURF=THIS%ESURF+this%tot_nuc_charge_energy 
!      WRITE(522,*)this%counter,this%tot_nuc_charge_energy*hartree,temporary*hartree,THIS%ESURF*hartree
      endif 
  contains

    subroutine cap(this, sig)
      type(poisson_sete_t), intent(in)    :: this
      FLOAT,                intent(out)   :: sig(:, :, :, :)

      integer :: i, j, k
      FLOAT   :: border

      !  9/24/08 - this VERSION OF CAP is very generic. It is
      !  for surfER. Assumes just this%ipio(0:NX+1,0:NY+1,0:NZ+1), 
      !  this%VBOUND(0:NX+1,0:NY+1,0:NZ+1),this%dielectric(NX,NY,NZ).
      !  Outputs a single 3D array SIG(NX,NY,NZ)
      !
      !  for calculating capacitances
      !
      !  SIG(I,J,K) is capacitance
      !
      !  input: NX,NY,NZ
      !         this%ipio,this%VBOUND
      !         ZG,DXG,DYG,DZ
      !         THIS%VH_BIG,this%dielectric,PCONST
      !
      !  output: SIG(I,J,K,6) 6 faces of grid cube 
      !          used in ETOT
      !  1 I-1, 2 I+1, 3 J-1, 4 J+1, 5 K-1, 6 K+1
      !
      !  local: BORDER
      !
      !  a homogeneous grid approximation is made: derivative
      !  at metal-semi border calculated discretely with 1/DXG (1/DYG,1/DZ)
      !  rather than properly averaged spacings.

      SIG= M_ZERO

      DO I=1,THIS%NXTOT 
        DO J=1,THIS%NYTOT
          DO K=1,THIS%NZTOT
            IF(this%ipio(I,J,K) == 0)THEN ! I,J,K is in Poisson grid
              IF(this%ipio(I-1,J,K) == 1)THEN
                BORDER=2*(THIS%VH_BIG(I,J,K)-this%VBOUND(I-1,J,K))/DXG(I)
                SIG(I,J,K,1)=SIG(I,J,K,1)-BORDER*DYG(J)*DZ(K)* &
                  this%dielectric(I,J,K)
              end if
              IF(this%ipio(I+1,J,K) == 1)THEN
                BORDER=2*(THIS%VH_BIG(I,J,K)-this%VBOUND(I+1,J,K))/DXG(I)
                SIG(I,J,K,2)=SIG(I,J,K,2)-BORDER*DYG(J)*DZ(K)* &
                  this%dielectric(I,J,K)
              end if
              IF(this%ipio(I,J-1,K) == 1)THEN
                BORDER=2*(THIS%VH_BIG(I,J,K)-this%VBOUND(I,J-1,K))/DYG(J)
                SIG(I,J,K,3)=SIG(I,J,K,3)-BORDER*DXG(I)*DZ(K)* &
                  this%dielectric(I,J,K)
              end if
              IF(this%ipio(I,J+1,K) == 1)THEN
                BORDER=2*(THIS%VH_BIG(I,J,K)-this%VBOUND(I,J+1,K))/DYG(J)
                SIG(I,J,K,4)=SIG(I,J,K,4)-BORDER*DXG(I)*DZ(K)* &
                  this%dielectric(I,J,K)
              end if
              IF(this%ipio(I,J,K-1) == 1)THEN
                BORDER=2*(THIS%VH_BIG(I,J,K)-this%VBOUND(I,J,K-1))/DZ(K)
                SIG(I,J,K,5)=SIG(I,J,K,5)-BORDER*DXG(I)*DYG(J)* &
                  this%dielectric(I,J,K)
              end if
              IF(this%ipio(I,J,K+1) == 1)THEN
                BORDER=2*(THIS%VH_BIG(I,J,K)-this%VBOUND(I,J,K+1))/DZ(K)
                SIG(I,J,K,6)=SIG(I,J,K,6)-BORDER*DXG(I)*DYG(J)* &
                  this%dielectric(I,J,K)
              end if
            end if
          end do
        end do
      end do

    end subroutine cap


  end subroutine egate

  FLOAT function poisson_sete_energy(this) result(energy)
    type(poisson_sete_t), intent(in)    :: this

    energy = this%esurf
  end function poisson_sete_energy

end module poisson_sete_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

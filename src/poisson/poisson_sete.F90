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
    v_es,                &
    v_es3,               &
    count_atoms,         &
    calc_gate_energy,    &
    es_energy

  type poisson_sete_t
    integer, pointer :: ipio(:,:,:)
    FLOAT,   pointer :: adiag(:)
    FLOAT,   pointer :: qs(:)
    FLOAT,   pointer :: vbound(:,:,:)
    FLOAT,   pointer :: vbound_lap(:,:,:)
    FLOAT,   pointer :: top(:,:)
    FLOAT,   pointer :: bottom(:,:)
    integer, pointer :: iad(:)
    integer, pointer :: jad(:)
    FLOAT,   pointer :: aw(:)
    FLOAT,   pointer :: dielectric(:,:,:)
    FLOAT,   pointer :: vh_big(:,:,:)
    FLOAT,   pointer :: vh_lap(:,:,:)
    integer, pointer :: iy(:)
    integer, pointer :: jy(:)
    integer, pointer :: ky(:)
    FLOAT            :: esurf
    FLOAT            :: tol
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
    integer          :: nx2
    integer          :: ny2
    integer          :: nz2
    integer          :: md
    FLOAT            :: dx
    FLOAT            :: dy
    FLOAT            :: dz
    FLOAT            :: zwidth
    integer          :: po
    integer          :: free
    FLOAT,   pointer :: vt(:)
    FLOAT,   pointer :: nucp(:,:,:)
    FLOAT,   pointer :: v_lap(:)
  end type poisson_sete_t

  integer :: lenw, leniw, idev, isym = 0, itol = 2, itmax = 2001, itermin = 5, iter, ierr, iunit = 0 
  integer :: nxl, nyl, test=0, nuc_or_elec=0, add_bias=0, calc_egate=1, egate_just_field=0, nuc_egate=1
  FLOAT, allocatable :: xg(:), yg(:), zg(:), dxg(:), dyg(:), dz(:), dxl(:), dyl(:),x2(:)

  ! these variables must be local for the moment
  FLOAT, allocatable :: rho_nuc(:), v_es3(:,:,:), v_es(:)
  integer :: count_atoms, calc_gate_energy
  FLOAT :: es_energy

  FLOAT, parameter :: hartree = CNST(2.0*13.60569193), BOHR = CNST(0.52917720859) 
  

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
    FLOAT, allocatable :: vtv(:), q2(:)
    FLOAT :: xwidth, ywidth,  pconst

    PUSH_SUB(poisson_sete_init)

    pconst = CNST(4.0)*M_PI ! need 4 pi for Hartrees I think (??)
    bohrnm = BOHR*CNST(0.1) 
    angsnm = CNST(10.0)/BOHR
    this%tol = 0.0001
    this%noatoms = number_atoms
    count_atoms = 0
    calc_gate_energy = 0
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
      read(57,*) this%zwidth

      this%zwidth = this%zwidth*angsnm  ! distance between plates (nm)

      read(57,*) zcen
      zcen = zcen*angsnm

      read(57,*) ngates

      write(68,*) "Allocating v_es3", nx, ny, nz
      SAFE_ALLOCATE(v_es3(1:nx, 1:ny, 1:nz))
      v_es3(1:nx, 1:ny, 1:nz)=M_ZERO
      write(68,*) " Allocated"

      SAFE_ALLOCATE(vtv(1:ngates))
      SAFE_ALLOCATE(this%vt(1:ngates))

      read(57,*)vtv(1),vtv(2) ! voltage on plates 
      !this%vt=vtv/(hartree*(this%noatoms+1))
      this%vt=vtv/(hartree)!*(this%noatoms+1))

      SAFE_DEALLOCATE_A(vtv)

      READ(57,*) dielectric0    ! dielectric constant of region betw plates
      dielectric0 = dielectric0/pconst

      ! mesh related
      READ(57,*)xwidth,ywidth ! bounding box around octopus box
      xwidth=xwidth*angsnm
      ywidth=ywidth*angsnm
      READ(57,*)xcen,ycen
      xcen=xcen*angsnm
      ycen=ycen*angsnm
      READ(57,*)this%nx2,this%ny2,this%nz2   ! additional mesh points in bd box
      READ(57,*)this%dx,this%dy,this%dz      ! Differences, experimental
      this%dx=this%dx*angsnm
      this%dy=this%dy*angsnm
      this%dz=this%dz*angsnm
      READ(57,*)nxl,nyl       ! padding points (to sim zone edge)
      allocate(dxl(nxl)); 
      allocate(dyl(nyl))
      READ(57,*)(dxl(i),i=1,nxl)
      READ(57,*)(dyl(j),j=1,nyl)
      read(57,*)this%po
      read(57,*)this%free
    ELSE
      WRITE(6,*)' IDEV not supported '
      RETURN
    end if

    call xyzgrid(this, nx, ny, nz, xl, yl, zl, xcen, ycen, zcen, xwidth, ywidth)!, this%zwidth)  ! establish x-y-z grids

    call cbsurf(this) !Assign boundaries 

    this%nelt=32+20*((this%nxtot-2)+(this%nytot-2)+(this%nztot-2))+ &
      12*((this%nxtot-2)*(this%nytot-2)+(this%nxtot-2)*(this%nztot-2)+(this%nytot-2)*(this%nztot-2)) + &
      7*(this%nxtot-2)*(this%nytot-2)*(this%nztot-2)
    if (this%nelt<=0) then
	  write(6,*) "this%nelt", this%nelt
          stop ' GOTTA QUIT '
  end if

    lenw=this%nelt+9*this%nxtot*this%nytot*this%nztot
    leniw=this%nelt+5*this%nxtot*this%nytot*this%nztot+12

    SAFE_ALLOCATE(q2(1:this%ntot))
    SAFE_ALLOCATE(this%aw(1:this%nelt))
    SAFE_ALLOCATE(this%qs(1:this%ntot))
    SAFE_ALLOCATE(this%adiag(1:this%ntot))
    SAFE_ALLOCATE(this%iad(1:this%nelt))
    SAFE_ALLOCATE(this%jad(1:this%nelt))
    
    SAFE_ALLOCATE(this%v_lap(1:this%ntot))
    this%v_lap(1:this%ntot)=M_ZERO

    call poissonm(this)

    this%qs = q2 ! store Dirichlet BC info in this%QS

    SAFE_DEALLOCATE_A(q2)
    POP_SUB(poisson_sete_init)

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

      FLOAT   :: av(7)
      INTEGER :: joff(7)
      FLOAT   :: d1, d2, dconst, del, err
      integer :: i, j, k, ia, iav, icol, idebug, iec, n, ii, m
      FLOAT,   allocatable :: q2_lap(:), rwork(:)
      integer, allocatable :: iwork(:) 
      integer :: i1, j1, k1
      !

      PUSH_SUB(poisson_sete_init.poissonm)
      SAFE_ALLOCATE(q2_lap(1:this%ntot))

      idebug=0
      if(idebug == 1)then
        open(57,file='this%aw.dat',status='unknown')
      endif
      !
      joff(1)=-this%md ! used to compute row indices this%JAD
      joff(2)=-this%nztot 
      joff(3)=-1 
      joff(4)=0 
      joff(5)=1 
      joff(6)=this%nztot 
      joff(7)=this%md
      I=1; J=1; K=1
      IEC=1 ! counter for number of non-zero entries in this%AW.

      q2(1:this%ntot) = M_ZERO
      q2_lap(1:this%ntot) = M_ZERO

      DO ia=1,this%ntot
        av=M_ZERO
        IF(this%ipio(i,j,k) == 0)THEN ! the point itself must be
          !                                   in the grid
          DO n=1,6 ! for all neighbors of this IA
            IF(n.GT.3)THEN
              iav=n+1
            ELSE
              iav=n
            end if
            ! set grid intervals based on which direction
            ! this neighbor is located in. (+-I,+-J,+-K).
            SELECT CASE(iav)
            CASE(1) ! I-1,J,K
              d1=dz(k); d2=dyg(j)
              SELECT CASE(this%ipio(I-1,J,K))
              CASE(1) ! Dirichlet boundary point
                del=M_HALF*dxg(i)*dz(k)
                dconst=this%dielectric(I,J,K)
                av(iav)=M_ZERO 
                av(4)=av(4)+dconst*d1*d2/del
                q2(ia)=q2(iA)+this%vbound(i-1,j,k)*dconst*d1*d2/del
                q2_lap(ia)=q2_lap(iA)+this%vbound_lap(i-1,j,k)*dconst*d1*d2/del
              CASE(2) ! Neumann boundary point
                av(iav)=M_ZERO
                del=M_ZERO
              CASE DEFAULT ! interior point
                del=M_HALF*(dxg(i)+dxg(i-1))*dz(k)
                dconst=M_HALF*(this%dielectric(i-1,j,k)+this%dielectric(I,J,K))
                av(iav)=-dconst*d1*d2/del 
                av(4)=av(4)+dconst*d1*d2/del
              END SELECT
            CASE(2) ! I,J-1,K
              d1=dz(k)
              d2=dxg(i)
              SELECT CASE(this%ipio(I,J-1,K))
              CASE(1) ! Dirichlet
                del=M_HALF*dyg(j)*dz(k)
                dconst=this%dielectric(I,J,K)
                av(iav)=M_ZERO 
                av(4)=av(4)+dconst*d1*d2/del
                q2(ia)=q2(ia)+this%vbound(i,j-1,k)*dconst*d1*d2/del
                q2_lap(ia)=q2_lap(ia)+this%vbound_lap(i,j-1,k)*dconst*d1*d2/del
              CASE(2) ! Neumann
                av(iav)=M_ZERO 
		del=M_ZERO
              CASE DEFAULT! interior point
                del=M_HALF*(dyg(j)+dyg(j-1))*dz(k)
                dconst=M_HALF*(this%dielectric(i ,j ,k)+this%dielectric(i,j-1,k))
                av(iav)=-dconst*d1*d2/del 
                av(4)=av(4)+dconst*d1*d2/del
              END SELECT
            CASE(3) ! I,J,K-1
              D1=DXG(I); D2=DYG(J)
              SELECT CASE(this%ipio(I,J,K-1))
              CASE(1) ! Dirichlet
                DEL=M_HALF*DZ(K)*DZ(K)
                DCONST=this%dielectric(I,J,K)
                AV(IAV)=M_ZERO 
		av(4)=av(4)+DCONST*D1*D2/DEL
                q2(ia)=q2(ia)+this%vbound(i,j,k-1)*dconst*d1*d2/del
                q2_lap(ia)=q2_lap(ia)+this%vbound_lap(i,j,k-1)*dconst*d1*d2/del
              CASE(2) ! Neumann
                av(iav)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT! interior point
                DEL=M_HALF*(DZ(K)+DZ(K-1))*DZ(K)
                DCONST=M_HALF*(this%dielectric(I,J,K)+this%dielectric(I,J,K-1))
                av(iav)=-dconst*d1*d2/del 
		av(4)=av(4)+dconst*d1*d2/del
              END SELECT
            CASE(5) ! I,J,K+1
              D1=DXG(I); D2=DYG(J)
              SELECT CASE(this%ipio(I,J,K+1))
              CASE(1) ! Dirichlet
                DEL=M_HALF*DZ(K)*DZ(K)
                DCONST=this%dielectric(I,J,K)
                AV(IAV)=M_ZERO 
		av(4)=av(4)+DCONST*D1*D2/DEL
                q2(ia)=q2(ia)+THIS%vbound(i,j,k+1)*dconst*d1*d2/del
                q2_lap(ia)=q2_lap(ia)+THIS%vbound_lap(i,j,k+1)*dconst*d1*d2/del
              CASE(2) ! Neumann
                AV(IAV)=M_ZERO; DEL=M_ZERO
              CASE DEFAULT! interior point
                del=M_HALF*(dz(k)+dz(k+1))*dz(k)
                dconst=M_HALF*(this%dielectric(i,j,k)+this%dielectric(i,j,k+1))
                av(iav)=-dconst*d1*d2/del 
                av(4)=av(4)+dconst*d1*d2/del
              END SELECT
            CASE(6) ! I,J+1,K
              D1=DZ(K); D2=DXG(I)
              SELECT CASE(this%ipio(I,J+1,K))
              CASE(1) ! Dirichlet
                DEL=M_HALF*DYG(J)*DZ(K)
                DCONST=this%dielectric(I,J,K)
                AV(IAV)=M_ZERO 
		AV(4)=AV(4)+DCONST*D1*D2/DEL
                Q2(IA)=Q2(IA)+this%VBOUND(I,J+1,K)*DCONST*D1*D2/DEL
                q2_lap(IA)=q2_lap(IA)+this%vbound_lap(I,J+1,K)*dconst*d1*d2/del
              CASE(2) ! Neumann
                AV(IAV)=M_ZERO 
		DEL=M_ZERO
              CASE DEFAULT
                DEL=M_HALF*(DYG(J)+DYG(J+1))*DZ(K)
                DCONST=M_HALF*(this%dielectric(I,J,K)+this%dielectric(I,J+1,K))
                AV(IAV)=-DCONST*D1*D2/DEL 
		AV(4)=AV(4)+DCONST*D1*D2/DEL
              END SELECT
            CASE(7) ! I+1,J,K
              D1=DZ(K); D2=DYG(J)
              SELECT CASE(this%ipio(I+1,J,K))
              CASE(1) ! Dirichlet
                del=M_HALF*dxg(i)*dz(k)
                dconst=this%dielectric(i,j,k)
                av(iav)=M_ZERO 
		av(4)=av(4)+dconst*d1*d2/del
                q2(ia)=q2(ia)+this%vbound(i+1,j,k)*dconst*d1*d2/del
                q2_lap(ia)=q2_lap(ia)+this%vbound_lap(i+1,j,k)*dconst*d1*d2/del
              CASE(2) ! Neumann
                av(iav)=M_ZERO 
		del=M_ZERO
              CASE DEFAULT
                del=M_HALF*(dxg(i)+dxg(i+1))*dz(k)
                dconst=M_HALF*(this%dielectric(I,J,K)+this%dielectric(i+1,j,k))
                av(iav)=-dconst*d1*d2/del
                av(4)=av(4)+dconst*d1*d2/del
              END SELECT
            END SELECT
          end do
        ELSE
          av(4) = M_ONE/pconst ! for points that are out of the grid.
        end if
        !                                 DSLUCS matrix and normalizations
        this%adiag(ia)=av(4)
        IF(this%adiag(ia) == M_ZERO)THEN
          WRITE(*,*)' DIAG ELEM ZERO ',ia
          STOP ' GOTTA QUIT '
        end if

        if(idebug == 1) then
          write(57, '(1x,i8,8(1x,e11.4))') iec, (av(ii),ii=1,7), q2(ia)
        endif

        DO ICOL=1,7  !  this IS SLAP TRIAD FORMAT
          IF(AV(ICOL) /= M_ZERO)THEN
            this%AW(IEC)=AV(ICOL)/this%adiag(ia)
            this%IAD(IEC)=IA
            this%JAD(IEC)=IA+JOFF(ICOL)
            IEC=IEC+1
          end if
        end do
        !increment I,J,K
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

      IEC=IEC-1 ! overcounted non-zero elements by one.
      !write(6,*)' this%nelt, iec ',this%nelt,iec
      this%nelt = iec
      write(6,*) "this%nelt", this%nelt


      ! in case THIS%NELT in static differs from that counted here
!!$        allocate(this%AWP(THIS%NELT))
!!$        DO I=1,THIS%NELT
!!$           this%AWP(I)=this%AW(I)
!!$        end do
!!$        deallocate(this%AW)
!!$        allocate(this%AW(THIS%NELT))
!!$        this%AW=this%AWP
!!$        deallocate(this%AWP)

     
    do i=1,this%ntot
      q2_lap(i)=q2_lap(i)/this%adiag(i)!! Laplacian scaled with diagonal elements
    enddo
    allocate(rwork(lenw))
    allocate(iwork(leniw))

	this%tol=0.000001
    call dslucs(this%ntot,q2_lap,this%v_lap,this%nelt,this%iad,this%jad,this%aw, &
    isym,itol, this%tol,itmax,iter,itermin,err,ierr,iunit,rwork,lenw, &
      iwork,leniw)
    do m=1, this%ntot
	      i = this%iy(m)
	      j = this%jy(m)
	      k = this%ky(m)
              IF((I.GT.NXL+THIS%NXBOT.AND.I.LE.NXL+THIS%NXBOT+NX).AND. &
                (J.GT.NYL+THIS%NYBOT.AND.J.LE.NYL+THIS%NYBOT+NY).AND. &
                (K.GT.THIS%NZBOT.AND.K.LE.THIS%NZBOT+NZ))THEN
                   I1=I-(NXL+THIS%NXBOT)
                   J1=J-(NYL+THIS%NYBOT)
                   K1=K-(THIS%NZBOT)
		   v_es3(i1, j1, k1) = this%v_lap(m)
		   !if(i1==1.and.j1==1) then
		   !write(225,*) k1, zg(k1), v_es3(i1,j1,k1)*27.2
	           !endif
      endif
    enddo


    deallocate(rwork)
    deallocate(iwork)
    SAFE_DEALLOCATE_P(this%vbound)
    SAFE_DEALLOCATE_A(q2_lap)
      if(idebug == 1)then
        close(57)
      endif

    POP_SUB(poisson_sete_init.poissonm)
    end subroutine poissonm

    !---------------------------------------------------------

    subroutine cbsurf(this)
      type(poisson_sete_t), intent(inout) :: this

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

      PUSH_SUB(poisson_sete_init.cbsurf)

      allocate(this%ipio(0:THIS%NXTOT+1,0:THIS%NYTOT+1,0:THIS%NZTOT+1))
      allocate(this%vh_big(1:this%nxtot,1:this%nytot,1:this%nztot))
      this%vh_big(1:this%nxtot,1:this%nytot,1:this%nztot)= M_ZERO
      if (egate_just_field==0) then
      allocate(this%nucp(1:this%nxtot,1:this%nytot,1:this%nztot))
      this%nucp(1:this%nxtot,1:this%nytot,1:this%nztot)= M_ZERO
      endif
      allocate(this%vbound(0:THIS%NXTOT+1,0:THIS%NYTOT+1,0:THIS%NZTOT+1))
      allocate(this%vbound_lap(0:THIS%NXTOT+1,0:THIS%NYTOT+1,0:THIS%NZTOT+1))
      allocate(this%dielectric(0:THIS%NXTOT+1,0:THIS%NYTOT+1,0:THIS%NZTOT+1))
      this%ipio=0; 
      this%vbound= M_ZERO
      this%vbound_lap= M_ZERO
      
      if (this%po==1) then
      SAFE_ALLOCATE(this%top(1:this%nxtot,1:this%nytot))
      SAFE_ALLOCATE(this%bottom(1:this%nxtot,1:this%nytot))
      do i=1,this%nxtot
      do j=1,this%nytot
      this%top(i,j)=M_ZERO
      this%bottom(i,j)=M_ZERO
      enddo
      enddo
      endif
      !write(68,*) "IDEV," , idev
      IF(IDEV == 1)THEN ! parallel plates
        this%dielectric=dielectric0
        this%ipio(:,:,0)=1
        this%ipio(:,:,THIS%NZTOT+1)=1

        this%vbound(:,:,0)=0.0
        this%vbound(:,:,this%nztot+1)=0.0
        this%vbound_lap(:,:,0)=this%vt(1)
        this%vbound_lap(:,:,this%nztot+1)=this%vt(2)

        DO k=1,this%nztot
          IF(k /= 0 .AND. k /= this%nztot+1)THEN !Redundant?
            ! lateral BC`s
            this%ipio(0,:,K)=2
            this%ipio(THIS%NXTOT+1,:,K)=2
            this%ipio(:,0,K)=2
            this%ipio(:,THIS%NYTOT+1,K)=2
          end if
        end do
      else if(idev == 0) then ! one plate only
	write(*,*) "There are still some parameters that need to be"
	write(*,*) "correctly accounted for."
	write(*,*) "Exiting."
	stop
        this%dielectric=dielectric0
        this%ipio(:,:,0)=1

        this%vbound(:,:,0)=0.0
        this%vbound_lap(:,:,0)=this%vt(1)

        DO k=1,this%nztot
          IF(k /= 0 )THEN !Redundant?
            ! lateral BC`s
            this%ipio(0,:,K)=2
            this%ipio(THIS%NXTOT+1,:,K)=2
            this%ipio(:,0,K)=2
            this%ipio(:,THIS%NYTOT+1,K)=2
          end if
        end do
      end if

      POP_SUB(poisson_sete_init.cbsurf)
    end subroutine cbsurf

    !---------------------------------------------------------

    subroutine xyzgrid(this, nx, ny, nz, xl, yl, zl, xcen, ycen, zcen, xwidth, ywidth)
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
      FLOAT,                intent(in)    :: xwidth
      FLOAT,                intent(in)    :: ywidth

      ! Input: NX, NY, NZ -> No. of octopus grid points in each dimension
      ! Input: XL, YL, ZL -> Size of the octopus box
      ! Input : zwidth is now this%zwidth

      FLOAT   :: dx1b, dx1t, dx_oct, xbot, xtest, xtop, xwidth_big
      FLOAT   :: dy1b, dy1t, dy_oct, ybot, ytop, ywidth_big
      FLOAT   :: dz1b, dz1t, dz_oct, zbot, ztest, ztop
      integer :: i, j, k, m
      integer :: nxtop, nytop, nztop
      integer :: nxtot_0, nytot_0

      ! inputs:  NZ,NZ2,NX,NX2,NY,NY2
      !          XCEN,YCEN,ZCEN
      !          XWIDTH,YWIDTH,this%zwidth
      !          XL,YL,ZL
      !
      ! outputs: ZG,XG,YG        
      !write(62,*)"XL,YL,ZL ",XL,YL,ZL
      ! z-mesh -- differs from x,y mesh because no border points
      !           hence NZTOT_0 = THIS%NZTOT (ie there IS no NZTOT_0)


      !top and bottom z coordinates.
      !Depend on POISSON_SETE (this%zwidth) 
      !and Octopus (ZL)parameters.

      PUSH_SUB(poisson_sete_init.xyzgrid)

      dz_oct=(zl/TOFLOAT((nz-1))) 
      ztop = M_HALF*this%zwidth - zcen - ((nz-1)/2)*dz_oct -M_HALF*dz_oct
      zbot = M_HALF*this%zwidth + zcen - ((nz-1)/2)*dz_oct -M_HALF*dz_oct
      if (ztop<M_ZERO.and.zbot<M_ZERO) then
	      write(*,*) "ztop and zbot are less than zero", ztop, zbot
	      stop
      endif
      
      if (ztop<=M_ZERO) then
	      ztop=M_ZERO
	      nztop=0
	      dz1t=M_ZERO
      else
	      dz1t=this%dz
	      write(70,*) 'dz1t', dz1t
	      nztop=ztop/dz1t
	      write(70,*) 'nztop', nztop
	      if (mod(ztop,this%dz)/=M_ZERO) then
		      nztop=nztop+1
		      write(70,*) 'nztop', nztop
		      dz1t=ztop/nztop
		      write(70,*) 'dz1t', dz1t
	      endif
      endif
      if (zbot<=M_ZERO) then
	      zbot=M_ZERO
	      this%nzbot=0
	      dz1b=M_ZERO
      else
	      dz1b=this%dz
	      this%nzbot = zbot/dz1b
	      if (mod(zbot,this%dz)/=M_ZERO) then
		      this%nzbot=this%nzbot+1
		      dz1b=zbot/this%nzbot
	      endif
      endif

      this%nz2=this%nzbot+nztop
      this%nztot = nz + this%nz2
      allocate(zg(this%nztot))
      allocate(dz(this%nztot))
      DO k=1,this%nzbot
        dz(k)=dz1b
      end do
      !DZ_OCT=ZL/TOFLOAT(NZ)

      do k = this%nzbot + 1, this%nzbot + nz
        dz(k) = dz_oct
      end do

      DO K=THIS%NZBOT+NZ+1,THIS%NZTOT
        DZ(K)=DZ1T
      end do
      write(70,*)"this%zwidth, ",this%zwidth 
      write(70,*)"zbot ztop",zbot,ztop ! Differences, experimental
      write(70,*)"dz1b, dz1t, dz_oct",dz1b,dz1t, dz_oct  ! Differences, experimental
      write(70,*)"nztot, nzbot nztop", this%nztot, this%nzbot, nztop

      ZG(1)=-M_HALF*this%zwidth+M_HALF*DZ(1)
      DO K=2,THIS%NZTOT-1 
      ZG(K)=ZG(K-1)+M_HALF*(DZ(K-1)+DZ(K))
      end do
      ZG(THIS%NZTOT)=M_HALF*this%zwidth-M_HALF*DZ(THIS%NZTOT)

      do i=1,this%nztot
      write(70,*) "i,zg, dz",i, zg(i), dz(i)
      enddo
      ztest=zg(this%nztot)
      ! x-mesh
      dx_oct=xl/TOFLOAT(nx-1)  ! Octopus region
      xtop=M_HALF*xwidth-xcen-((nx-1)/2)*dx_oct-M_HALF*dx_oct
      xbot=M_HALF*xwidth+xcen-((nx-1)/2)*dx_oct-M_HALF*dx_oct
      if (xtop<M_ZERO.and.xbot<M_ZERO) then
	      write(*,*) "ztop and zbot are less than zero", ztop, zbot
	      stop
      endif
      if (xtop<=M_ZERO) then
	      xtop=M_ZERO
	      nxtop=0
	      dx1t=M_ZERO
      else
	      dx1t=this%dx
	      write(70,*) 'dx1t', dx1t
	      nxtop=xtop/dx1t
	      write(70,*) 'nxtop', nxtop
	      if (mod(xtop,this%dx)/=M_ZERO) then
		      nxtop=nxtop+1
		      write(70,*) 'nxtop', nxtop
		      dx1t=xtop/nxtop
		      write(70,*) 'dx1t', dx1t
	      endif
      endif
      if (xbot<=M_ZERO) then
	      xbot=M_ZERO
	      this%nxbot=0
	      dx1b=M_ZERO
      else
	      dx1b=this%dx
	      this%nxbot = xbot/dx1b
	      if (mod(xbot,this%dx)/=M_ZERO) then
		      this%nxbot=this%nxbot+1
		      dx1b=xbot/this%nxbot
	      endif
      endif
      this%nx2=this%nxbot+nxtop
      nxtot_0=nx+this%nx2
      this%nxtot=nxtot_0+2*nxl
      allocate(xg(this%nxtot))
      allocate(dxg(this%nxtot))

      DO I=1,NXL  !  border points
	DXG(I)=DXL(I)
      end do
      DO I=NXL+NXTOT_0+1,THIS%NXTOT
	DXG(I)=DXL(THIS%NXTOT-I+1)
      enddo
      NXTOP=INT(this%NX2*XTOP/(XTOP+XBOT))
      THIS%NXBOT=this%NX2-NXTOP
      DX1B=XBOT/TOFLOAT(THIS%NXBOT) ! "bottom" region
      DO I=NXL+1,NXL+THIS%NXBOT
        DXG(I)=DX1B
      end do
      DX1T=XTOP/TOFLOAT(NXTOP)  ! "top" region
      DO I=NXL+THIS%NXBOT+NX+1,NXL+THIS%NXBOT+NX+NXTOP
        DXG(I)=DX1T
      end do
      write(70,*)"xwidth ",xwidth ! Differences, experimental
      write(70,*)"xbot xtop",xbot,xtop ! Differences, experimental
      write(70,*)"dx1b, dx1t, dx_oct",dx1b,dx1t, dx_oct  ! Differences, experimental
      write(70,*)"nxtot, nxbot, nxtop nx2",this%nxtot,this%nxbot, nxtop, this%nx2
      DO I=NXL+THIS%NXBOT+1,NXL+THIS%NXBOT+NX
        DXG(I)=DX_OCT
      end do
      XWIDTH_BIG=SUM(DXG)
      XG(1)=-M_HALF*XWIDTH_BIG+M_HALF*DXG(1)
      DO I=2,THIS%NXTOT-1
        XG(I)=XG(I-1)+M_HALF*(DXG(I-1)+DXG(I))
      end do
      XG(THIS%NXTOT)=M_HALF*XWIDTH_BIG-M_HALF*DXG(THIS%NXTOT)
      do i=1,this%nxtot
      write(70,*) "i, xg, dx",i, xg(i), dxg(i)
      enddo
          
      XTEST=XG(THIS%NXTOT)
      ! Y-mesh
      dy_oct=yl/TOFLOAT(ny-1)  ! Octopus region
      ytop=M_HALF*ywidth-ycen-((ny-1)/2)*dy_oct-M_HALF*dy_oct
      ybot=M_HALF*ywidth+ycen-((ny-1)/2)*dy_oct-M_HALF*dy_oct
      if (ytop<M_ZERO.and.ybot<M_ZERO) then
	      write(*,*) "ztop and zbot are less than zero", ztop, zbot
	      stop
      endif
      if (ytop<=M_ZERO) then
	      ytop=M_ZERO
	      nytop=0
	      dy1t=M_ZERO
      else
	      dy1t=this%dy
	      write(70,*) 'dy1t', dy1t
	      nytop=ytop/dy1t
	      write(70,*) 'nytop', nytop
	      if (mod(ytop,this%dy)/=M_ZERO) then
		      nytop=nytop+1
		      write(70,*) 'nytop', nytop
		      dy1t=ytop/nytop
		      write(70,*) 'dy1t', dy1t
	      endif
      endif
      if (ybot<=M_ZERO) then
	      ybot=M_ZERO
	      this%nybot=0
	      dy1b=M_ZERO
      else
	      dy1b=this%dy
	      this%nybot = ybot/dy1b
	      if (mod(ybot,this%dy)/=M_ZERO) then
		      this%nybot=this%nybot+1
		      dy1b=ybot/this%nybot
	      endif
      endif
      this%ny2=this%nybot+nytop
      nytot_0=ny+THIS%ny2  ! excluding border points
      this%nytot=nytot_0+2*nyl  ! all Poisson Y mesh point
      allocate(yg(this%nytot))
      allocate(dyg(this%nytot))
      DO I=1,NYL  !  border points
        DYG(I)=DYL(I)
      end do
      do i=nyl+nytot_0+1,this%nytot
        DYG(I)=DYL(THIS%NYTOT-I+1)
      end do
      NYTOP=INT(this%NY2*YTOP/(YTOP+YBOT))
      THIS%NYBOT=this%NY2-NYTOP
      DY1B=YBOT/TOFLOAT(THIS%NYBOT) ! "bottom" region
      DO I=NYL+1,NYL+THIS%NYBOT
        DYG(I)=DY1B
      end do
      DO I=NYL+THIS%NYBOT+1,NYL+THIS%NYBOT+NY
        DYG(I)=DY_OCT
      end do
      DY1T=YTOP/TOFLOAT(NYTOP)  ! "top" region
      DO I=NYL+THIS%NYBOT+NY+1,NYL+THIS%NYBOT+NY+NYTOP
        DYG(I)=DY1T
      end do

      YWIDTH_BIG=SUM(DYG)
      YG(1)=-M_HALF*YWIDTH_BIG+M_HALF*DYG(1)
      DO I=2,THIS%NYTOT-1
        YG(I)=YG(I-1)+M_HALF*(DYG(I-1)+DYG(I))
      end do
      YG(THIS%NYTOT)=M_HALF*YWIDTH_BIG-M_HALF*DYG(THIS%NYTOT)

      write(70,*)"ywidth", ywidth ! Differences, experimental
      write(70,*)"ybot ytop", ybot,ytop ! Differences, experimental
      write(70,*)"dy1b, dy1t, dy_oct",dy1b,dy1t, dy_oct  ! Differences, experimental
      write(70,*)"nytot, nybot, nytop ny2",this%nytot,this%nybot, nytop, this%ny2
      do i=1,this%nytot
      write(70,*) "i, yg, dy",i, yg(i), dyg(i)
      enddo

      this%ntot=this%nxtot*this%nytot*this%nztot
      this%md=this%nztot*this%nytot

      allocate(this%iy(this%ntot))
      allocate(this%jy(this%ntot))
      allocate(this%ky(this%ntot))
      m=0
      DO i=1,this%nxtot
        DO j=1,this%nytot
          DO k=1,this%nztot
            m=m+1
            this%iy(m)=i
            this%jy(m)=j
            this%ky(m)=k
          end do
        end do
      end do

      POP_SUB(poisson_sete_init.xyzgrid)
    end subroutine xyzgrid

  end subroutine poisson_sete_init

  !---------------------------------------

  subroutine poisson_sete_solve(this, rho, vh, nx, ny, nz)
    type(poisson_sete_t), intent(inout) :: this
    FLOAT,                intent(in)    :: rho(:, :, :)
    FLOAT,                intent(inout) :: vh(:, :, :)
    integer,              intent(in)    :: nx
    integer,              intent(in)    :: ny
    integer,              intent(in)    :: nz

    ! Input: NX, NY, NZ -> No. of octopus grid points in each dimension

    FLOAT        :: bohrnm, angsnm, err
    integer      :: j1, k1, m, i, j, k, i1,idum
    FLOAT :: VHMIN, VHMAX
    FLOAT :: pconst

    integer, allocatable :: iwork(:)
    FLOAT,   allocatable :: rhotest(:,:,:), rwork(:), q2(:)!, q2_lap(:)!, x2(:)

    PUSH_SUB(poisson_sete_solve)

    PCONST=CNST(4.0)*M_PI ! need 4 pi for Hartrees I think (??)
    BOHRNM=BOHR/CNST(10.0)
    ANGSNM=CNST(10.0)/BOHR

    SAFE_ALLOCATE(q2(1:THIS%NTOT))
!    SAFE_ALLOCATE(rhotest(1:this%nxtot, 1:this%nytot, 1:this%nztot))


    if (rho(nx/2,ny/2,nz/2)<0.and.nuc_or_elec==1) then
	!When the electron calculation ends and starts calling the 
	!nuclear energy term.
	nuc_or_elec=0
	calc_egate=0
        add_bias=0
	this%tol=0.000001
        SAFE_DEALLOCATE_A(X2)
    endif

    if (rho(nx/2,ny/2,nz/2)<M_ZERO.and.nuc_or_elec==0) then
	    !Nuclear potential calculation
      nuc_or_elec=0
	this%tol=0.000001
      idum = count_atoms
      add_bias=0
      !To ensure that it does not add the laplacian solver to v0
      !if (calc_gate_energy==1) then
!	      add_bias=1
!      endif
      SAFE_ALLOCATE(X2(1:THIS%NTOT))
      do i = 1, this%ntot
        call quickrnd(idum,x2(i))
        if (x2(i) > M_ZERO) then
            x2(i) = - x2(i)
        endif
      idum = count_atoms
      enddo

    else if (rho(nx/2,ny/2,nz/2)>0.and.nuc_or_elec==0) then
	    !Electronic potential calculation
      nuc_or_elec=1
      calc_egate=1
      write(200,*) "Electronic"
      add_bias=1
      this%tol=0.0001
      SAFE_ALLOCATE(x2(1:this%ntot))
      x2(1:this%ntot) = M_ZERO
      idum = count_atoms + 1
      do i = 1, this%ntot
        x2(i) = M_ZERO
        call quickrnd(idum, x2(i))
        if (x2(i) < M_ZERO) then
            x2(i) = -x2(i)
        endif
     enddo
    endif

!    do i=1,this%nxtot
!     do j=1, this%nytot
!      do k=1, this%nztot
!        rhotest(i,j,k) = M_ZERO
!      enddo
!     enddo
!    enddo

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
          if(i1.gt.0.and.j1.gt.0.and.k1.gt.0)then ! this is redundant
            q2(m)=q2(m)+dyg(j)*dxg(i)*rho(i1,j1,k1)
!            rhotest(i,j,k)=rho(i1,j1,k1)
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
      this%TOL,ITMAX,ITER,ITERMIN,ERR,IERR,IUNIT,RWORK,LENW, &
      IWORK,LENIW)
    if (ierr /= 0) then
	    write(*,*) "DSLUCS error number:", ierr, err
	    stop
    endif
    leniw = iwork(9)
    lenw  = iwork(10) ! new work array lengths


    deallocate(rwork)
    deallocate(iwork)

    if ((add_bias==1).or.(add_bias==2)) then
	    do m=1,this%ntot
	    i=this%iy(m) 
	    j=this%jy(m)
 	    k=this%ky(m)
	    x2(m)=x2(m)+this%v_lap(m)
	 enddo
 endif

    !Only calculate once.
    !Test 1: Just use the voltage bias
    if (egate_just_field==1) then
    if (rho(nx/2,ny/2,nz/2)>0.and.calc_egate==1) then
	    !field=(this%vt(2)-this%vt(1))/(this%zwidth)!zg(this%nztot)-zg(1))
	    do m = 1, this%ntot ! store x2(m) i,j,k wise
	      i = this%iy(m)
	      j = this%jy(m)
	      k = this%ky(m)
              this%vh_big(i, j, k) = this%v_lap(m)
	      enddo
    endif 
    else
      if (calc_egate==1) then
	      if (rho(nx/2,ny/2,nz/2)<0.and.nuc_egate==1) then
		      do m = 1, this%ntot ! store x2(m) i,j,k wise
		      i = this%iy(m)
		      j = this%jy(m)
		      k = this%ky(m)
		      this%nucp(i,j,k) = this%nucp(i,j,k)+x2(m)
		      enddo
		if (idum==this%noatoms.and.nuc_egate==1) then
		    nuc_egate=0
		    endif
               else if (rho(nx/2,ny/2,nz/2)>0) then
		      !The electronic potential should already have the field in.
		      do m = 1, this%ntot ! store x2(m) i,j,k wise
		      i = this%iy(m)
		      j = this%jy(m)
		      k = this%ky(m)
                      this%vh_big(i, j, k) = x2(m) + this%nucp(i, j, k)
		      enddo
               endif
        endif
    endif



    do m = 1, this%ntot
      I=THIS%IY(M); J=THIS%JY(M); K=THIS%KY(M)
      if(this%ipio(I,J,K) == 0) then
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


    if (rho(nx/2,ny/2,nz/2)>0.and.calc_egate==1) then
    !if (rho(nx/2,ny/2,nz/2)>0) then
    !if (calc_egate==1) then
    call egate(this)
    endif
    vhmin = minval(vh)
    vhmax = maxval(vh)

    deallocate(Q2)
!    deallocate(rhotest)
    if (rho(nx/2,ny/2,nz/2)<0.and.nuc_or_elec==0) then
          SAFE_DEALLOCATE_A(x2)
    endif

    !deallocate(THIS%VH_BIG);

    POP_SUB(poisson_sete_solve)
  end subroutine poisson_sete_solve

  !------------------------------------------------------------

  subroutine poisson_sete_end(this)
    type(poisson_sete_t), intent(inout) :: this

    PUSH_SUB(poisson_sete_end)

    SAFE_DEALLOCATE_P(this%ipio)
    SAFE_DEALLOCATE_P(this%qs)
    SAFE_DEALLOCATE_P(this%vbound_lap)
    SAFE_DEALLOCATE_P(this%iad)
    SAFE_DEALLOCATE_P(this%JAD)

    deallocate(this%adiag)
    if (nuc_or_elec==1) then
            nuc_or_elec=0
            SAFE_DEALLOCATE_A(X2)
	    write(77,*) "deAllocating x2 poisson_sete_end"
    endif
    deallocate(DXL)
    deallocate(DYL)
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
    deallocate(this%v_lap)
    deallocate(this%aw)
    if (this%po==1) then
	    deallocate(this%top)
	    deallocate(this%bottom)
	    deallocate(this%vt)
    endif
    if (egate_just_field==0) then
	    deallocate(this%nucp)
    endif

    POP_SUB(poisson_sete_end)
  end subroutine poisson_sete_end

  !--------------------------------------------
subroutine egate(this)
type(poisson_sete_t), intent(inout) :: this

    integer :: i, j, k
    FLOAT, allocatable :: sig(:, :, :, :) 
    FLOAT :: cs1, cs2, cs3,temporary
    FLOAT :: charge_top, charge_surf, charge_bot

    PUSH_SUB(egate)

    !
    !     to do energy calculations previously performed in setr
    !  compute capacitances of six surfaces
    
    SAFE_ALLOCATE(sig(1:this%nxtot, 1:this%nytot, 1:this%nztot, 1:6))


    !this%vbound(:,:,0)=this%vt(1)
    !this%vbound(:,:,this%nztot+1)=this%vt(2)
    call cap(this, sig)

    this%esurf = M_ZERO
    charge_surf = M_ZERO
    charge_top = M_ZERO
    charge_bot = M_ZERO
    cs1 = M_ZERO
    cs2 = M_ZERO
    cs3 = M_ZERO
    !write(*,*) "Begin cycle ", THIS%NXTOT, THIS%NYTOT, THIS%NZTOT
    do i=1,this%nxtot
      do j=1,this%nytot
        do k=1,this%nztot ! n.b. [sig]=CHARGE (NOT CHARGE/aREA)
          this%esurf=this%esurf+this%vbound_lap(i-1,j,k)*sig(i,j,k,1)
          this%esurf=this%esurf+this%vbound_lap(i+1,j,k)*sig(i,j,k,2)
          this%esurf=this%esurf+this%vbound_lap(i,j-1,k)*sig(i,j,k,3)
          this%esurf=this%esurf+this%vbound_lap(i,j+1,k)*sig(i,j,k,4)
          this%esurf=this%esurf+this%vbound_lap(i,j,k-1)*sig(i,j,k,5)
          this%esurf=this%esurf+this%vbound_lap(i,j,k+1)*sig(i,j,k,6)
          cs1=cs1+sig(i,j,k,1)+sig(i,j,k,2)
          cs2=cs2+sig(i,j,k,3)+sig(i,j,k,4)
          cs3=cs3+sig(i,j,k,5)+sig(i,j,k,6)
        END DO
        charge_top=charge_top+sig(i,j,1,5)+sig(i,j,1,6)
        charge_bot=charge_bot+sig(i,j,this%nztot,5)+sig(i,j,this%nztot,6)
      end do
    end do
    !this%vbound(:,:,0)=0.0!this%vt(1)
    !this%vbound(:,:,this%nztot+1)=0.0!this%vt(2)

    
    this%esurf = M_HALF*this%esurf
    charge_surf = cs1 + cs2 + cs3

    !    WRITE(520,*)THIS%ESURF*hartree*M_HALF
        write(521,*)' CS1 ',CS1,' CS2 ',CS2,' CS3 ',CS3
	write(521,*) this%counter, 'Charge top', charge_top, sig(this%nxtot/2,this%nytot/2,1,5)
	write(521,*) this%counter, 'Charge bot', charge_bot, sig(this%nxtot/2,this%nytot/2,this%nztot,6)


!    if (this%counter<this%noatoms) then
!      this%counter=this%counter+1
!      this%tot_nuc_charge_energy=this%tot_nuc_charge_energy+this%esurf
!      write(523,*) "nuc",this%tot_nuc_charge_energy*CNST(2.0*13.60569193), this%esurf
!      if (this%po==1) then
!!       open(unit=550+this%counter,action="write",status='replace')
!!       write(550+this%counter,*) "#Gate 1, x, y, Total Charge"
!      do i=1,this%nxtot
!        do j=1,this%nytot
!	 this%top(i,j)=this%top(i,j)+sig(i,j,1,5)+sig(i,j,1,6)
!	 this%bottom(i,j)=this%bottom(i,j)+sig(i,j,this%nztot,5)+sig(i,j,this%nztot,6)
!           write(550+this%counter,50) xg(i), yg(j), sig(i,j,1,5)
!	enddo
!           write(550+this%counter,30)
!      enddo
!      close(550+this%counter)
!      endif
!    else if ((this%counter == this%noatoms).and.(sig(this%nxtot/2,this%nytot/2,1,5)<0)) then
!      temporary=this%esurf
!      this%esurf=this%esurf+this%tot_nuc_charge_energy 

      write(523,*)this%tot_nuc_charge_energy*CNST(2.0*13.60569193), this%esurf
      if (this%po==1) then
       open(unit=524, file="V1.charge",action="write",status='replace')
       open(unit=525, file="V2.charge",action="write",status='replace')
       write(524,*) "#Gate 1, x, y, Total Charge, Electron Charge, Nuclear Charge", this%counter
       write(525,*) "#Gate 1, x, y, Total Charge, Electron Charge, Nuclear Charge", this%counter
       do i=1,this%nxtot
         do j=1,this%nytot
           write(524,40) xg(i), yg(j), this%top(i,j)+sig(i,j,1,5), sig(i,j,1,5), this%top(i,j)
           write(525,40) xg(i), yg(j), this%bottom(i,j)+sig(i,j,this%nztot,6),sig(i,j,this%nztot,6),this%bottom(i,j)
         enddo
           write(524,30)
           write(525,30)
       enddo
       close(524)
       close(525)
      endif
      30 format()
      40 format(5E24.16)

!      WRITE(522,*)this%counter,this%tot_nuc_charge_energy*hartree,temporary*hartree,THIS%ESURF*hartree
!    endif
    SAFE_DEALLOCATE_A(sig)

    POP_SUB(egate)

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

      PUSH_SUB(egate.cap)

      SIG= M_ZERO

      DO I=1,THIS%NXTOT 
        DO J=1,THIS%NYTOT
          DO K=1,THIS%NZTOT
            IF(this%ipio(I,J,K) == 0)THEN ! I,J,K is in Poisson grid
              IF(this%ipio(I-1,J,K) == 1)THEN
                BORDER=2*(THIS%VH_BIG(I,J,K)-this%vbound_lap(I-1,J,K))/DXG(I)
                SIG(I,J,K,1)=SIG(I,J,K,1)-BORDER*DYG(J)*DZ(K)* &
                  this%dielectric(I,J,K)
              end if
              IF(this%ipio(I+1,J,K) == 1)THEN
                BORDER=2*(THIS%VH_BIG(I,J,K)-this%vbound_lap(I+1,J,K))/DXG(I)
                SIG(I,J,K,2)=SIG(I,J,K,2)-BORDER*DYG(J)*DZ(K)* &
                  this%dielectric(I,J,K)
              end if
              IF(this%ipio(I,J-1,K) == 1)THEN
                BORDER=2*(THIS%VH_BIG(I,J,K)-this%vbound_lap(I,J-1,K))/DYG(J)
                SIG(I,J,K,3)=SIG(I,J,K,3)-BORDER*DXG(I)*DZ(K)* &
                  this%dielectric(I,J,K)
              end if
              IF(this%ipio(I,J+1,K) == 1)THEN
                BORDER=2*(THIS%VH_BIG(I,J,K)-this%vbound_lap(I,J+1,K))/DYG(J)
                SIG(I,J,K,4)=SIG(I,J,K,4)-BORDER*DXG(I)*DZ(K)* &
                  this%dielectric(I,J,K)
              end if
              IF(this%ipio(I,J,K-1) == 1)THEN
                BORDER=2*(THIS%VH_BIG(I,J,K)-this%vbound_lap(I,J,K-1))/DZ(K)
                SIG(I,J,K,5)=SIG(I,J,K,5)-BORDER*DXG(I)*DYG(J)* &
                  this%dielectric(I,J,K)
              end if
              IF(this%ipio(I,J,K+1) == 1)THEN
                BORDER=2*(THIS%VH_BIG(I,J,K)-this%vbound_lap(I,J,K+1))/DZ(K)
                SIG(I,J,K,6)=SIG(I,J,K,6)-BORDER*DXG(I)*DYG(J)* &
                  this%dielectric(I,J,K)
              end if
            end if
          end do
        end do
      end do

      POP_SUB(egate.cap)
    end subroutine cap

  end subroutine egate

  FLOAT function poisson_sete_energy(this) result(energy)
    type(poisson_sete_t), intent(in)    :: this

    PUSH_SUB(poisson_sete_energy)

    if (this%free==1) then
     energy = this%esurf
    else
     energy = 0.0
    endif

    POP_SUB(poisson_sete_energy)
  end function poisson_sete_energy

end module poisson_sete_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

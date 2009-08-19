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
    poisson_sete_init,   &
    poisson_sete_solve,  &
    poisson_sete_end,    &
    poisson_sete_energy

  public ::              &
    rho_nuc,             &
    count_atoms

  FLOAT, allocatable :: rho(:, :, :), vh(:, :, :)
  FLOAT, allocatable :: aw(:), awp(:), q2(:), x2(:), qs(:)
  integer, allocatable :: iad(:), jad(:), iy(:), jy(:), ky(:)
  FLOAT :: tol=0.01, hartree=CNST(2.0*13.60569193), bohr=CNST(0.52917720859) ! bohr radius in nm
  integer :: isete_on, nelt, ntot, lenw, leniw, idev, ngates, &
    isym = 0, itol = 2, itmax = 201, itermin = 5, iter, ierr, iunit = 0, &
    nxbot, nybot, nzbot, nxl, nyl, nxtot, nytot, nztot, md
  FLOAT, allocatable :: vbound(:,:,:), dielectric(:,:,:), rhotest(:,:,:), vh_big(:,:,:)
  FLOAT, allocatable :: xg(:), yg(:), zg(:), dxg(:), dyg(:), dz(:), rwork(:), dxl(:), dyl(:), vt(:), vtv(:), adiag(:)
  FLOAT :: xwidth, ywidth, zwidth, dielectric0, pconst
  integer, allocatable :: ipio(:,:,:)
  integer, allocatable :: iwork(:), idiag(:)
  integer :: nx2, ny2, nz2
  FLOAT, allocatable :: sig(:,:,:,:)
  FLOAT :: esurf, charge_top, charge_bot, charge_tot, border, cs1, cs2, cs3, charge_surf
  FLOAT, allocatable :: rho_nuc(:)
  integer :: noatoms, count_atoms

contains

  subroutine poisson_sete_init(nx, ny, nz, xl, yl, zl,number_atoms)
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    FLOAT,   intent(in) :: xl
    FLOAT,   intent(in) :: yl
    FLOAT,   intent(in) :: zl
    integer, intent(in) :: number_atoms

    character(len=40) :: cdum, fil1
    FLOAT :: bohrnm, angsnm
    integer :: i, j
    FLOAT :: zcen, xcen, ycen

    pconst = CNST(4.0)*M_PI ! need 4 pi for Hartrees I think (??)
    bohrnm = bohr/CNST(10.0)
    angsnm = CNST(10.0)/bohr
    noatoms = number_atoms
    count_atoms = 0

    open(56, file = 'poisq',status = 'unknown')
    fil1 = 'test1.pois'
    open(unit = 57, FILE = FIL1, STATUS = 'UNKNOWN') ! status info file

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

      read(57,*) vtv(1), vtv(2) ! voltage on plates in volts and nm?

      vt = vtv/hartree

      read(57,*) dielectric0    ! dielectric constant of region betw plates
      dielectric0 = dielectric0/pconst

      ! mesh related
      read(57,*)xwidth,ywidth ! bounding box around octopus box

      xwidth = xwidth*angsnm
      ywidth = ywidth*angsnm

      read(57,*) xcen,ycen
      xcen = xcen*angsnm
      ycen = ycen*angsnm

      read(57,*) nx2, ny2, nz2   ! additional mesh points in bd box
      read(57,*) nxl, nyl       ! padding points (to sim zone edge)

      SAFE_ALLOCATE(dxl(1:nxl))
      SAFE_ALLOCATE(dyl(1:nyl))

      read(57,*) (dxl(i), i = 1, nxl)
      read(57,*) (dyl(j), j = 1, nyl)

    else
      write(6,*)' idev not supported '
      return
    end if

    call xyzgrid(nx, ny, nz, xl, yl, zl, xcen, ycen, zcen)  ! establish x-y-z grids

    call cbsurf() !Assign boundaries 

    nelt = 32 + 20*((nxtot - 2) + (nytot - 2) +(nztot - 2)) + &
      12*((nxtot - 2)*(nytot - 2) + (nxtot-2)*(nztot-2) + (nytot-2)*(nztot-2)) + &
      7*(nxtot - 2)*(nytot - 2)*(nztot - 2)

    lenw = nelt + 9*nxtot*nytot*nztot
    leniw = nelt + 5*nxtot*nytot*nztot + 12

    SAFE_ALLOCATE(q2(1:ntot))
    SAFE_ALLOCATE(aw(1:nelt))
    SAFE_ALLOCATE(qs(1:ntot))
    SAFE_ALLOCATE(adiag(1:ntot))
    SAFE_ALLOCATE(idiag(1:ntot))
    SAFE_ALLOCATE(iad(1:nelt))
    SAFE_ALLOCATE(jad(1:nelt))
    
    call poissonm()

    qs = q2 ! store Dirichlet BC info in QS

    SAFE_DEALLOCATE_A(q2)

  end subroutine poisson_sete_init

  !---------------------------------------

  subroutine poisson_sete_solve(icase, rho, vh, nx, ny, nz, xl, yl, zl, icalc)
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

    FLOAT   :: bohrnm, angsnm, err
    integer :: j1, k1, m, i, j, k, i1, idum
    FLOAT   :: vhmin, vhmax

    pconst = CNST(4.0)*M_PI ! need 4 pi for Hartrees I think (??)
    bohrnm = bohr/CNST(10.0)
    angsnm = CNST(10.0)/bohr

    SAFE_ALLOCATE(Q2(1:NTOT))
    SAFE_ALLOCATE(rhotest(1:nxtot, 1:nytot, 1:nztot))

    write(message(1),'(a,i7)') "Sete: count atoms ", count_atoms

    if (count_atoms <= noatoms) then
      idum = count_atoms

      SAFE_ALLOCATE(x2(1:ntot))

      do i = 1, ntot
        call quickrnd(idum, x2(i))
        if (x2(i) > M_ZERO) then
          x2(i) = - x2(i)
        endif
      enddo

    else if (count_atoms == noatoms + 1) then

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
      q2(i) = qs(i) ! recall stored BC info
    enddo

    do m = 1, ntot
      i = iy(M)
      j = jy(M)
      k = ky(M)
      if(ipio(i, j, k) == 0)then
        if((i > nxl + nxbot .and. i <= nxl + nxbot + nx).and. &
          (j > nyl + nybot .and. j <= nyl + nybot + ny).and. &
          (k > nzbot .and. k <= nzbot + nz)) then
          i1 = i - nxl + nxbot
          j1 = j - nyl + nybot
          k1 = k - nzbot
          if(i1 > 0.and.j1 > 0.and.k1 > 0)then ! this is redundant
            q2(m) = q2(m) + dyg(j)*dxg(i)*rho(i1, j1, k1)
            rhotest(i, j, k) = rho(i1, j1, k1)
          end if
        end if
      else
        q2(m) = M_ZERO
      end if
    end do

    do i = 1,ntot
      q2(i) = q2(i)/adiag(i)!! Laplacian scaled with diagonal elements
    enddo

    SAFE_ALLOCATE(rwork(1:lenw))
    SAFE_ALLOCATE(iwork(1:leniw))

    call dslucs(ntot,q2,x2,nelt,iad,jad,aw,isym,itol, &
      tol,itmax,iter,itermin,err,ierr,iunit,rwork,lenw, &
      iwork,leniw)

    leniw = iwork(9)
    lenw  = iwork(10) ! new work array lengths

    deallocate(rwork)
    deallocate(iwork)

    do m=1,ntot ! store x2(m) i,j,k wise
      i = iy(m)
      j = jy(m)
      k = ky(m)
      vh_big(i,j,k) = x2(m)
    end do

    do m = 1, ntot
      i = iy(m)
      j = jy(m)
      k = ky(M)
      if(ipio(i, j, k) == 0)then
        if((i > nxl + nxbot .and. i <= nxl + nxbot + nx) .and. &
          (j > nyl + nybot .and. j <= nyl + nybot + ny).and. &
          (k > nzbot .and. k <= nzbot + nz)) then
          i1 = i - nxl + nxbot
          j1 = j - nyl + nybot
          k1 = k - nzbot
          vh(i1, j1, k1) = x2(m)
        end if
      end if
    end do

    call egate()

    write(358,*) "#x,z, vh_big(x,11,k)"

    do i = 1, nxtot
      do k = 1, nztot
        write(358,*) i, k, vh_big(i, 11, k)
      enddo
    enddo

    vhmin = minval(vh)
    vhmax = maxval(vh)

    SAFE_DEALLOCATE_A(Q2)
    SAFE_DEALLOCATE_A(rhotest)

    if (count_atoms <= noatoms ) then
      SAFE_DEALLOCATE_A(x2)
    endif

  end subroutine poisson_sete_solve

  ! ---------------------------------------------

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

    FLOAT   :: av(7)
    integer :: joff(7)
    FLOAT   :: d1, d2, dconst, del
    integer :: i, j, k, ia, iav, icol, idebug, iec, n, ii

    idebug = 0

    if(idebug == 1)then
      open(57,file='aw.dat',status='unknown')
    endif

    joff(1) = -md ! used to compute row indices jad
    joff(2) = -nztot
    joff(3) = -1
    joff(4) = 0
    joff(5) = 1
    joff(6) = nztot
    joff(7) = md

    i = 1 
    j = 1
    k = 1

    iec = 1

    q2(1:ntot) = M_ZERO

    do ia = 1, ntot
      av = M_ZERO

      if(ipio(i, j, k) == 0) then ! the point itself must be
        !                                   in the grid
        do n = 1, 6 ! for all neighbors of this ia
          if(n > 3) then
            iav = n + 1
          else
            iav = n
          end if

          ! set grid intervals based on which direction
          ! this neighbor is located in. (+-i,+-j,+-k).
          select case(iav)
          case(1) ! i-1,j,k
            d1 = dz(k)
            d2 = dyg(j)
            select case(ipio(i - 1, j, k))
            case(1) ! dirichlet boundary point
              del = CNST(0.5)*dxg(i)*dz(k)
              dconst = dielectric(i, j, k)
              av(iav) = M_ZERO
              av(4) = av(4) + dconst*d1*d2/del
              q2(ia) = q2(ia) + vbound(i - 1, j, k)*dconst*d1*d2/del
            case(2) ! neumann boundary point
              av(iav) = M_ZERO
              del = M_ZERO
            case default ! interior point
              del = CNST(0.5)*(dxg(i) + dxg(i - 1))*dz(k)
              dconst = CNST(0.5)*(dielectric(i - 1, j, k) + dielectric(i, j, k))
              av(iav) = -dconst*d1*d2/del
              av(4) = av(4) + dconst*d1*d2/del
            end select
          case(2) ! i,j-1,k
            d1 = dz(k)
            d2 = dxg(i)
            select case(ipio(i, j - 1, k))
            case(1) ! dirichlet
              del = CNST(0.5)*dyg(j)*dz(k)
              dconst = dielectric(i, j, k)
              av(iav) = M_ZERO
              av(4) = av(4) + dconst*d1*d2/del
              q2(ia) = q2(ia) + vbound(i, j - 1, k)*dconst*d1*d2/del
            case(2) ! neumann
              av(iav) = M_ZERO
              del = M_ZERO
            case default! interior point
              del = CNST(0.5)*(dyg(j) + dyg(j - 1))*dz(k)
              dconst = CNST(0.5)*(dielectric(i, j, k) + dielectric(i, j - 1, k))
              av(iav) = -dconst*d1*d2/del
              av(4) = av(4) + dconst*d1*d2/del
            end select
          case(3)
            d1 = dxg(i)
            d2 = dyg(j)
            select case(ipio(i, j, k - 1))
            case(1) ! dirichlet
              del = CNST(0.5)*dz(k)*dz(k)
              dconst = dielectric(i, j, k)
              av(iav) = M_ZERO
              av(4) = av(4) + dconst*d1*d2/del
              q2(ia) = q2(ia) + vbound(i, j, k - 1)*dconst*d1*d2/del
            case(2) ! neumann
              av(iav) = M_ZERO; del = M_ZERO
            case default! interior point
              del = CNST(0.5)*(dz(k)+dz(k - 1))*dz(k)
              dconst = CNST(0.5)*(dielectric(i, j, k) + dielectric(i, j, k - 1))
              av(iav) = -dconst*d1*d2/del
              av(4) = av(4) + dconst*d1*d2/del
            end select
          case(5) ! i,j,k+1
            d1 = dxg(i); d2 = dyg(j)
            select case(ipio(i, j, k + 1))
            case(1) ! dirichlet
              del = CNST(0.5)*dz(k)*dz(k)
              dconst = dielectric(i, j, k)
              av(iav) = M_ZERO
              av(4) = av(4)+dconst*d1*d2/del
              q2(ia) = q2(ia)+vbound(i,j,k+1)*dconst*d1*d2/del
            case(2) ! neumann
              av(iav) = M_ZERO; del = M_ZERO
            case default! interior point
              del = CNST(0.5)*(dz(k)+dz(k+1))*dz(k)
              dconst = CNST(0.5)*(dielectric(i,j,k)+dielectric(i,j,k+1))
              av(iav) = -dconst*d1*d2/del
              av(4) = av(4)+dconst*d1*d2/del
            end select
          case(6) ! i,j+1,k
            d1 = dz(k); d2 = dxg(i)
            select case(ipio(i, j + 1, k))
            case(1) ! dirichlet
              del = CNST(0.5)*dyg(j)*dz(k)
              dconst = dielectric(i, j, k)
              av(iav) = M_ZERO
              av(4) = av(4) + dconst*d1*d2/del
              q2(ia) = q2(ia) + vbound(i, j + 1, k)*dconst*d1*d2/del
            case(2) ! neumann
              av(iav) = M_ZERO; del = M_ZERO
            case default
              del = CNST(0.5)*(dyg(j) + dyg(j+1))*dz(k)
              dconst = CNST(0.5)*(dielectric(i, j, k) + dielectric(i, j+1, k))
              av(iav) = -dconst*d1*d2/del
              av(4) = av(4) + dconst*d1*d2/del
            end select
          case(7) ! i+1,j,k
            d1 = dz(k)
            d2 = dyg(j)
            select case(ipio(i+1,j,k))
            case(1) ! dirichlet
              del = CNST(0.5)*dxg(i)*dz(k)
              dconst = dielectric(i, j, k)
              av(iav) = M_ZERO
              av(4) = av(4) + dconst*d1*d2/del
              q2(ia) = q2(ia) + vbound(i + 1, j, k)*dconst*d1*d2/del
            case(2) ! neumann
              av(iav) = M_ZERO
              del = M_ZERO
            case default
              del = CNST(0.5)*(dxg(i)+dxg(i+1))*dz(k)
              dconst = CNST(0.5)*(dielectric(i,j,k)+dielectric(i+1,j,k))
              av(iav) = -dconst*d1*d2/del
              av(4) = av(4)+dconst*d1*d2/del
            end select
          end select
        end do
      else
        write(57,*) "at else", i, j, k
        av(4) = M_ONE/pconst ! for points that are out of the grid.
      end if

      ! dslucs matrix and normalizations
      adiag(ia) = av(4)

      if(adiag(ia) == M_ZERO) then
        write(message(1), '(a, i8)') 'Sete: diagonal element zero ', ia
        call write_fatal(1)
      end if

      if(idebug == 1) then
        write(57, '(1x,i8,8(1x,e10.4))') iec, (av(ii),ii = 1,7), q2(ia)
      endif

      do icol = 1,7  !  this is slap triad format
        if(av(icol) /= M_ZERO)then
          aw(iec) = av(icol)/adiag(ia)
          iad(iec) = ia
          jad(iec) = ia+joff(icol)
          if(icol == 4)then
            idiag(ia) = iec
          end if
          iec = iec+1
        end if
      end do
      !                                     increment i,j,k
      k = k+1
      if(k == nztot + 1)then
        k = 1
        j = j+1
        if(j == nytot + 1)then
          j = 1
          i = i + 1
        end if
      end if


    end do
    if(idebug == 1)then
      close(57)
    endif

    iec = iec-1 ! overcounted non-zero elements by one.
    write(6,*)' nelt, iec ',nelt,iec
    nelt = iec

    ! in case nelt in static differs from that counted here
!!$        allocate(AWP(NELT))
!!$        DO I = 1,NELT
!!$           AWP(I) = AW(I)
!!$        end do
!!$        deallocate(AW)
!!$        allocate(AW(NELT))
!!$        AW = AWP
!!$        deallocate(AWP)

    if(idebug == 1)then
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
    !  a homogeneous grid approximation is made: derivative
    !  at metal-semi border calculated discretely with 1/DXG (1/DYG,1/DZ)
    !  rather than properly averaged spacings.
    allocate(SIG(NXTOT,NYTOT,NZTOT,6))
    SIG=0.0
    DO I=1,NXTOT 
      DO J=1,NYTOT
        DO K=1,NZTOT
          IF(IPIO(I,J,K) == 0)THEN ! I,J,K is in Poisson grid
            IF(IPIO(I-1,J,K) == 1)THEN
              BORDER=2*(VH_BIG(I,J,K)-VBOUND(I-1,J,K))/DXG(I)
              SIG(I,J,K,1)=SIG(I,J,K,1)-BORDER*DYG(J)*DZ(K)* &
                DIELECTRIC(I,J,K)
            end if
            IF(IPIO(I+1,J,K) == 1)THEN
              BORDER=2*(VH_BIG(I,J,K)-VBOUND(I+1,J,K))/DXG(I)
              SIG(I,J,K,2)=SIG(I,J,K,2)-BORDER*DYG(J)*DZ(K)* &
                DIELECTRIC(I,J,K)
            end if
            IF(IPIO(I,J-1,K) == 1)THEN
              BORDER=2*(VH_BIG(I,J,K)-VBOUND(I,J-1,K))/DYG(J)
              SIG(I,J,K,3)=SIG(I,J,K,3)-BORDER*DXG(I)*DZ(K)* &
                DIELECTRIC(I,J,K)
            end if
            IF(IPIO(I,J+1,K) == 1)THEN
              BORDER=2*(VH_BIG(I,J,K)-VBOUND(I,J+1,K))/DYG(J)
              SIG(I,J,K,4)=SIG(I,J,K,4)-BORDER*DXG(I)*DZ(K)* &
                DIELECTRIC(I,J,K)
            end if
            IF(IPIO(I,J,K-1) == 1)THEN
              BORDER=2*(VH_BIG(I,J,K)-VBOUND(I,J,K-1))/DZ(K)
              SIG(I,J,K,5)=SIG(I,J,K,5)-BORDER*DXG(I)*DYG(J)* &
                DIELECTRIC(I,J,K)
            end if
            IF(IPIO(I,J,K+1) == 1)THEN
              BORDER=2*(VH_BIG(I,J,K)-VBOUND(I,J,K+1))/DZ(K)
              SIG(I,J,K,6)=SIG(I,J,K,6)-BORDER*DXG(I)*DYG(J)* &
                DIELECTRIC(I,J,K)
            end if
          end if
        end do
      end do
    end do

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
    SAFE_ALLOCATE(ipio(0:nxtot + 1, 0:nytot + 1, 0:nztot + 1))
    SAFE_ALLOCATE(vh_big(1:nxtot, 1:nytot, 1:nztot))
    vh_big(1:nxtot, 1:nytot, 1:nztot) = M_ZERO

    SAFE_ALLOCATE(vbound(0:nxtot + 1, 0:nytot + 1, 0:nztot + 1))
    SAFE_ALLOCATE(dielectric(0:nxtot + 1, 0:nytot + 1, 0:nztot + 1))

    ipio = 0
    vbound = M_ZERO

    if(idev == 1) then ! parallel plates
      dielectric = dielectric0
      dielectric(:, :, 0)= CNST(6.5)/pconst  !For platinum
      dielectric(:, :, nztot + 1)=CNST(6.5)/pconst !For platinum
      ipio(:, :, 0) = 1
      ipio(:, :, nztot + 1) = 1
      vbound(:, :, 0) = vt(1)
      vbound(:, :, nztot + 1) = vt(2)

      do k=1, nztot
        if(k /= 0 .and. k /= nztot + 1) then !Redundant?
          ! lateral bc`s
          ipio(0, :, k) = 2
          ipio(nxtot + 1, :, k) = 2
          ipio(:, 0, k) = 2
          ipio(:, nytot + 1, k) = 2
        end if
      end do
    end if

  end subroutine cbsurf

  subroutine poisson_sete_end()
    SAFE_DEALLOCATE_A(IPIO)
    SAFE_DEALLOCATE_A(VBOUND)
    SAFE_DEALLOCATE_A(VT)
    SAFE_DEALLOCATE_A(VTV)
    SAFE_DEALLOCATE_A(QS)
    SAFE_DEALLOCATE_A(IAD)
    SAFE_DEALLOCATE_A(JAD)
    SAFE_DEALLOCATE_A(ADIAG)
    SAFE_DEALLOCATE_A(IDIAG);
    SAFE_DEALLOCATE_A(DXL)
    SAFE_DEALLOCATE_A(DYL)
    SAFE_DEALLOCATE_A(dielectric)
    SAFE_DEALLOCATE_A(zg)
    SAFE_DEALLOCATE_A(xg)
    SAFE_DEALLOCATE_A(yg)
    SAFE_DEALLOCATE_A(dz)
    SAFE_DEALLOCATE_A(dyg)
    SAFE_DEALLOCATE_A(dxg)
    SAFE_DEALLOCATE_A(ky)
    SAFE_DEALLOCATE_A(iy)
    SAFE_DEALLOCATE_A(jy)
    SAFE_DEALLOCATE_A(VH_BIG)
    SAFE_DEALLOCATE_A(aw)
  end subroutine poisson_sete_end

  !--------------------------------------------

  subroutine egate

    integer :: i, j, k

    !
    !     to do energy calculations previously performed in setr
    !  compute capacitances of six surfaces

    call cap()

    esurf = M_ZERO
    charge_surf = M_ZERO
    charge_top = M_ZERO

    charge_bot = M_ZERO
    cs1 = M_ZERO
    cs2 = M_ZERO
    cs3 = M_ZERO

    do i = 1, nxtot
      do j = 1, nytot
        do k = 1, nztot ! n.b. [sig]=charge (not charge/area)
          esurf = esurf + vbound(i - 1, j, k)*sig(i, j, k, 1)
          esurf = esurf + vbound(i + 1, j, k)*sig(i, j, k, 2)
          esurf = esurf + vbound(i, j - 1, k)*sig(i, j, k, 3)
          esurf = esurf + vbound(i, j + 1, k)*sig(i, j, k, 4)
          esurf = esurf + vbound(i, j, k - 1)*sig(i, j, k, 5)
          esurf = esurf + vbound(i, j, k + 1)*sig(i, j, k, 6)
          cs1 = cs1 + sig(i, j, k, 1) + sig(i, j, k, 2)
          cs2 = cs2 + sig(i, j, k, 3) + sig(i, j, k, 4)
          cs3 = cs3 + sig(i, j, k, 5) + sig(i, j, k, 6)
        end do
        charge_top = charge_top + sig(i, j, 1, 5) + sig(i, j, 1, 6)
        charge_bot = charge_bot + sig(i, j, nztot, 5) + sig(i, j, nztot, 6)
      end do
    end do

    SAFE_DEALLOCATE_A(SIG)

    esurf = CNST(0.5)*esurf
    charge_surf = cs1 + cs2 + cs3

    WRITE(520,*)'Esurf', ESURF*27.2
    write(521,*)' CS1 ',CS1,' CS2 ',CS2,' CS3 ',CS3
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

    ywidth_big = sum(DYG)

    yg(1) = -CNST(0.5)*ywidth_big + CNST(0.5)*dyg(1)
    do i = 2, nytot - 1
      yg(i) = yg(i - 1) + CNST(0.5)*(dyg(i - 1) + dyg(i))
    end do
    yg(nytot) = CNST(0.5)*ywidth_big - CNST(0.5)*dyg(nytot)

    ntot = nxtot*nytot*nztot
    md = nztot*nytot

    SAFE_ALLOCATE(iy(1:ntot))
    SAFE_ALLOCATE(jy(1:ntot))
    SAFE_ALLOCATE(ky(1:ntot))

    m=0
    do i = 1, nxtot
      do j = 1, nytot
        do k = 1, nztot
          m = m + 1
          iy(m) = i
          jy(m) = j
          ky(m) = k
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


  FLOAT function poisson_sete_energy() result(energy)
    energy = esurf
  end function poisson_sete_energy

end module poisson_sete_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

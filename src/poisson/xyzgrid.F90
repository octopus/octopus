#include <global.h>

      subroutine xyzgrid(nx, ny, nz, xl, yl, zl)
     ! Input: NX, NY, NZ -> No. of octopus grid points in each dimension
     ! Input: XL, YL, ZL -> Size of the octopus box

	Use pois_data_g ! global Poisson data
        use pois_data_l ! local Poisson data

        implicit none

        integer, intent(in) :: nx, ny, nz
        FLOAT,   intent(in) :: xl, yl, zl

        FLOAT   :: dx1b, dx1t, dx_oct, xbot, xtest, xtop, xcen, xwidth_big
        FLOAT   :: dy1b, dy1t, dy_oct, ybot, ytop, ycen, ywidth_big
        FLOAT   :: dz1b, dz1t, dz_oct, zbot, ztest, ztop, zcen
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
        NZTOT=NZ+NZ2
        allocate(zg(nztot))
        allocate(dz(nztot))
	!top and bottom z coordinates.
	!Depend on POISSON_SETE (ZWIDTH) 
	!and Octopus (ZL)parameters.
        ZTOP=0.5*ZWIDTH-ZCEN-0.5*ZL
        ZBOT=0.5*ZWIDTH+ZCEN-0.5*ZL

        NZTOP=INT(NZ2*ZTOP/(ZTOP+ZBOT))
        NZBOT=NZ2-NZTOP
        DZ1B=ZBOT/TOFLOAT(NZBOT)
        DO K=1,NZBOT
           DZ(K)=DZ1B
        ENDDO
        DZ_OCT=ZL/TOFLOAT(NZ)
        DO K=NZBOT+1,NZBOT+NZ
           DZ(K)=DZ_OCT
        ENDDO
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
101     format(10(1x,e10.4))
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

        RETURN
      END SUBROUTINE XYZGRID

#include <global.h>

!module poisson_sete_m
!  use POIS_DATA_G
!  use POIS_DATA_L
!  use mesh_m

  !implicit none

  !private
  !public ::             &
  !   poisson_sete_init, &
  !   pois

!contains
     subroutine poisson_sete_init(nx,ny,nz,xl,yl,zl)


        USE POIS_DATA_G ! global Poisson data
        USE POIS_DATA_L ! local Poisson data
        IMPLICIT FLOAT (A-H,O-Z)
	FLOAT, intent(in) :: XL,YL,ZL
	INTEGER, intent(in) :: NX,NY,NZ
        CHARACTER*40 :: CDUM,FIL1
        FLOAT :: BOHRNM,ANGSNM
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
           READ(57,*)IDEV,IPOISSON_SETE_ON ! for future use - device configuration
           write(*,*)"IPOISSON_SETE is ON", IPOISSON_SETE_ON
           IF(IDEV.EQ.1)THEN
! physical device characteristics
              READ(57,*)ZWIDTH;ZWIDTH=ZWIDTH*ANGSNM  ! distance between plates (nm)
              READ(57,*)ZCEN;ZCEN=ZCEN*ANGSNM
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
      CALL XYZGRID(nx,ny,nz,xl,yl,zl)  ! establish x-y-z grids
      CALL CBSURF !Assign boundaries 
           NELT=32+20*((NXTOT-2)+(NYTOT-2)+(NZTOT-2))+ &
                12*((NXTOT-2)*(NYTOT-2)+(NXTOT-2)*(NZTOT-2)+(NYTOT-2)*(NZTOT-2)) + &
                7*(NXTOT-2)*(NYTOT-2)*(NZTOT-2)
           LENW=NELT+9*NXTOT*NYTOT*NZTOT
           LENIW=NELT+5*NXTOT*NYTOT*NZTOT+12
           allocate(Q2(NTOT)); allocate(AW(NELT))
           allocate(QS(NTOT))
           allocate(ADIAG(NTOT)); allocate(IDIAG(NTOT))
           allocate(IAD(NELT)); allocate(JAD(NELT))
           allocate(X2(NTOT)); X2=0.0
           CALL POISSONM
           QS=Q2 ! store Dirichlet BC info in QS
           deallocate(Q2)
         RETURN
      END subroutine



     SUBROUTINE POIS(ICASE,RHO,VH,NX,NY,NZ,XL,YL,ZL,ICALC)
     ! Input: NX, NY, NZ -> No. of octopus grid points in each dimension
     ! Input: XL, YL, ZL -> Size of the octopus box

        USE POIS_DATA_G ! global Poisson data
        USE POIS_DATA_L ! local Poisson data

        IMPLICIT FLOAT (A-H,O-Z)
	FLOAT, intent(in) :: XL,YL,ZL
	INTEGER, intent(in) :: NX,NY,NZ,ICASE,ICALC
        
        FLOAT :: RHO(NX,NY,NZ),VH(NX,NY,NZ)
        CHARACTER*40 :: CDUM,FIL1
        FLOAT :: BOHRNM,ANGSNM
	INTEGER :: J1,K1, M, I,J,K,IDUM,I1
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

           allocate(RWORK(LENW));allocate(IWORK(LENIW))
           CALL DSLUCS(NTOT,Q2,X2,NELT,IAD,JAD,AW,ISYM,ITOL, &
                TOL,ITMAX,ITER,ITERMIN,ERR,IERR,IUNIT,RWORK,LENW, &
                IWORK,LENIW)
           LENIW=IWORK(9); LENW=IWORK(10) ! new work array lengths
           deallocate(RWORK); deallocate(IWORK)
           
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
           VHMIN=MINVAL(VH); VHMAX=MAXVAL(VH)
           deallocate(Q2)
           deallocate(rhotest);! deallocate(VH_BIG);

         RETURN
      END subroutine

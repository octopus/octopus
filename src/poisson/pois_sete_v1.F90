#include <global.h>

     SUBROUTINE POIS(ICASE,RHO,VH,NX,NY,NZ,XL,YL,ZL,ICALC)
     ! Input: NX, NY, NZ -> No. of octopus grid points in each dimension
     ! Input: XL, YL, ZL -> Size of the octopus box

        USE POIS_DATA_G ! global Poisson data
        USE POIS_DATA_L ! local Poisson data

        IMPLICIT FLOAT (A-H,O-Z)
	FLOAT, intent(in) :: XL,YL,ZL
	INTEGER, intent(in) :: NX,NY,NZ
        FLOAT :: RHO(NX,NY,NZ),VH(NX,NY,NZ)
        CHARACTER*40 :: CDUM,FIL1
        FLOAT :: BOHRNM,ANGSNM
	INTEGER :: J1,K1, M, I,J,K
	FLOAT :: ZCEN, XCEN, YCEN, VHMIN, VHMAX
	ICALC=0
	write(*,*) "Welcome to SETE ICASE ICALC", ICASE, ICALC
        IF(ICASE.EQ.1)THEN
! INITIALIZATION
! read input file, set up poisson equation
           !write(6,*)' opening poisq '
           open(56,file='poisq',status='unknown')
           PI=DACOS(-1.D0)
           PCONST=4*PI ! need 4 pi for Hartrees I think (??)
	   BOHRNM=BOHR/10.0
	   ANGSNM=10.0/BOHR
           FIL1='test1.pois'
           OPEN(unit=57,FILE=FIL1,STATUS='UNKNOWN') ! status info file
! READ FILE
           READ(57,'(A40)')CDUM
           READ(57,'(A40)')CDUM
           READ(57,'(A40)')CDUM
           READ(57,*)IDEV,ISETE_ON ! for future use - device configuration
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
              !WRITE(6,*)' IDEV not supported '
              RETURN
           ENDIF
	   !XL,YL,ZL are given in terms of
	   !Bohr radii (due to internal usage).
	   !Then, we need to go back to A.
	   !Convert Octopus box from A to nm
	   !XL=XL*BOHRNM;YL=YL*BOHRNM;ZL=ZL*BOHRNM
           CALL XYZGRID(NX,NY,NZ,XL,YL,ZL)  ! establish x-y-z grids
!
           CALL CBSURF !Assign boundaries for plate case scenario.
!
! guess for number of non-zero elements in Laplacian
           NELT=32+20*((NXTOT-2)+(NYTOT-2)+(NZTOT-2))+ &
                12*((NXTOT-2)*(NYTOT-2)+(NXTOT-2)*(NZTOT-2)+(NYTOT-2)*(NZTOT-2)) + &
                7*(NXTOT-2)*(NYTOT-2)*(NZTOT-2)
! guess for input to DSLUCS - refined and output in IWORK(9) and IWORK(10)
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

        ELSEIF(ISETE_ON.EQ.1)THEN
           allocate(Q2(NTOT))
           allocate(rhotest(nxtot,nytot,nztot))
           Q2=QS ! recall stored BC info
           rhotest=0
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

           Q2=Q2/ADIAG ! Laplacian scaled with diagonal elements
           IDUM=1
!           DO I=1,NTOT
!              X2(I)=RAN1(IDUM) ! randomize sol'n vector Pois. Eq.
!           ENDDO
!           X2=0.0
!           WRITE(6,*)' ntot ',ntot
!           write(6,*)' nelt ',nelt
!           write(6,*)' lenw ',lenw
!           write(6,*)' leniw ',leniw
!           write(6,*)' size(aw) ',size(aw),size(rwork),size(iwork)
!           write(6,*)' size(q2,x2) ',size(q2),size(x2),size(iad),size(jad)

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
           VHMIN=MINVAL(VH); VHMAX=MAXVAL(VH)
           deallocate(Q2)
!           open(57,file='vhbig.dat',status='unknown')
!           vh_big(:,:,1)=250.0; vh_big(:,:,nztot)=250.0
!           xx1=-0.367;yy1=0.0;zz1=0.0
!           xx2=0.367;yy2=0.0;zz2=0.0
!           do i=1,nxtot
!              do j=1,nytot
!                 do k=1,nztot
!                    if(i.gt.nxl.and.i.le.nxtot-nxl.and.j.gt.nyl &
!		       .and.j.le.nytot-nyl)then
!                       r1=sqrt((xx1-xg(i))**2+(yy1-yg(j))**2+(zz1-zg(k))**2)
!                       r2=sqrt((xx2-xg(i))**2+(yy2-yg(j))**2+(zz2-zg(k))**2)
!                       if(r1.lt.0.15.or.r2.lt.0.15)then
!                          rhotest(i,j,k)=50.0
!                       endif
                       !write(57,104)xg(i),yg(j),zg(k),vh_big(i,j,k), rhotest(i,j,k)
!104                    format(3(1x,f12.6),3(1x,e12.6))
!                    endif
!                 enddo
!              enddo
!           enddo
!           close(57)

           
!           do k=1,nztot
!	      !write(57,206)zg(k),(vh_big(i,1+(nytot/2),k),i=1,nxtot)
!206           format(1x,f10.5,200(1x,e12.6))
!           enddo
!           close(57)
           deallocate(rhotest);! deallocate(VH_BIG);
	   ICALC=1
      ENDIF

         RETURN
      END

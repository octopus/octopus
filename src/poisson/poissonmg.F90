
#include <global.h>

      SUBROUTINE POISSONM

        use pois_data_g ! global Poisson data
        use pois_data_l ! local Poisson data

        implicit none

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
                    CASE DEFAULT
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
                    CASE DEFAULT
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
                    CASE DEFAULT
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
              AV(4)=1.D0 ! for points that are out of the grid.
           ENDIF
!                                 DSLUCS matrix and normalizations
           ADIAG(IA)=AV(4)
           IF(ADIAG(IA).EQ.0.D0)THEN
              WRITE(56,*)' DIAG ELEM ZERO ',IA
              STOP ' GOTTA QUIT '
           ENDIF
           
           if(idebug.eq.1)then
              write(57,101)iec,(av(ii),ii=1,7),q2(ia)
101           format(1x,i8,8(1x,e10.4))
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
     
        RETURN
      END SUBROUTINE POISSONM

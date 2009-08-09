#include <global.h>

      SUBROUTINE CAP
	USE DRIVER_DATA
        USE POIS_DATA_G ! global Poisson data
        USE POIS_DATA_L ! local Poisson data

      implicit none

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

    END SUBROUTINE CAP


      SUBROUTINE EGATE
        USE POIS_DATA_G ! global Poisson data
        USE POIS_DATA_L ! local Poisson data
      IMPLICIT REAL*8 (A-H,O-Z)
!
!     to do energy calculations previously performed in setr
!  compute capacitances of six surfaces
      ESURF=0.0; CHARGE_SURF=0.0; CHARGE_TOP=0.0
      !write(64,*) "#I Vbound(I,5,5), ESURF"
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
	 !write (64,*) I, ESURF
      ENDDO
      deallocate(SIG)
      ESURF=0.5*ESURF; CHARGE_SURF=CS1+CS2+CS3
      !WRITE(6,*)' CS1 ',CS1,' CS2 ',CS2,' CS3 ',CS3
      RETURN
      END

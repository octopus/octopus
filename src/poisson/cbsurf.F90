#include <global.h>
      SUBROUTINE CBSURF

        USE POIS_DATA_G ! global Poisson data
        USE POIS_DATA_L ! local Poisson data

        implicit none

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
        write(*,*) "Allocating data"
           allocate(VH_BIG(NXTOT,NYTOT,NZTOT))
        allocate(IPIO(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1))
        allocate(VBOUND(0:NXTOT+1,0:NYTOT+1,0:NZTOT+1))
        allocate(DIELECTRIC(NXTOT,NYTOT,NZTOT))
        IPIO=0; VBOUND=0.0
        IF(IDEV.EQ.1)THEN ! parallel plates
           DIELECTRIC=DIELECTRIC0
           IPIO(:,:,0)=1
           IPIO(:,:,NZTOT+1)=1
           VBOUND(:,:,0)=VT(1)
           VBOUND(:,:,NZTOT+1)=VT(2)
           DO K=1,NZTOT
              IF(K.NE.0.AND.K.NE.NZTOT+1)THEN
		! lateral BC's
                 IPIO(0,:,K)=2
                 IPIO(NXTOT+1,:,K)=2
                 IPIO(:,0,K)=2
                 IPIO(:,NYTOT+1,K)=2
              ENDIF
           ENDDO
        ENDIF
        
      END SUBROUTINE CBSURF
      SUBROUTINE CBSURF_end

        USE POIS_DATA_G ! global Poisson data
        USE POIS_DATA_L ! local Poisson data

        implicit none

	deallocate(IPIO)
      	deallocate(VBOUND)
	deallocate(VT);deallocate(VTV);deallocate(QS)
	deallocate(IAD);deallocate(JAD);deallocate(ADIAG)
	deallocate(IDIAG);deallocate(X2)
	deallocate(DXL);deallocate(DYL)
!	deallocate(VH_BIG)
	end subroutine CBSURF_end

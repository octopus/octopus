#include <global.h>

!------------------------------------------------------------------------------
      MODULE DRIVER_DATA
        FLOAT, DIMENSION(:,:,:), ALLOCATABLE :: RHO,VH

      END MODULE DRIVER_DATA
!
!------------------------------------------------------------------------------
      MODULE POIS_DATA_G
!------------------------------------------------------------------------------

        FLOAT, DIMENSION(:), ALLOCATABLE :: AW,AWP,Q2,X2,QS
        INTEGER, DIMENSION(:), ALLOCATABLE :: IAD,JAD,IY,JY,KY
        FLOAT :: TOL=0.01,HARTREE=2.0*13.60569193,BOHR=0.52917720859 ! Bohr radius in nm
!                                                Hartree in eV
        INTEGER :: ISETE_ON,NELT,NTOT,LENW,LENIW,IDEV,NGATES, &
             ISYM=0,ITOL=2,ITMAX=201,ITERMIN=5,ITER,IERR,IUNIT=0, &
             NXBOT,NYBOT,NZBOT,NXL,NYL,NXTOT,NYTOT,NZTOT,MD,PI,PCONST

      END MODULE POIS_DATA_G
!
!------------------------------------------------------------------------------
      MODULE POIS_DATA_L
!------------------------------------------------------------------------------

        FLOAT, DIMENSION(:,:,:), ALLOCATABLE :: VBOUND,DIELECTRIC, &
             rhotest,VH_BIG 
        FLOAT, DIMENSION(:), ALLOCATABLE :: XG,YG,ZG, &
             DXG,DYG,DZ,RWORK,DXL,DYL,VT,VTV,ADIAG
        FLOAT :: XWIDTH,YWIDTH,ZWIDTH,DIELECTRIC0
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IPIO
        INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK,IDIAG
        INTEGER :: NX2,NY2,NZ2
        FLOAT, DIMENSION(:,:,:,:), ALLOCATABLE :: SIG
        FLOAT :: ESURF,CHARGE_TOP,CHARGE_BOT,CHARGE_TOT,BORDER,&
                CS1,CS2,CS3,CHARGE_SURF
	FLOAT, DIMENSION(:), ALLOCATABLE :: rho_nuc(:)!, v_nuc(:)


      END MODULE POIS_DATA_L
!

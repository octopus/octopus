!      SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W) 
!      IMPLICIT REAL*8 (A-H,O-Z) 
!      DIMENSION X(*),W(*) 
      SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,CALFUN) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION X(*),W(*) 
      INTERFACE
        SUBROUTINE CALFUN(N,X,F)
          IMPLICIT REAL*8 (A-H,O-Z)
          DIMENSION X(*)
        END SUBROUTINE CALFUN
      END INTERFACE

!                                                                       
!     This subroutine seeks the least value of a function of many variab
!     by a trust region method that forms quadratic models by interpolat
!     There can be some freedom in the interpolation conditions, which i
!     taken up by minimizing the Frobenius norm of the change to the sec
!     derivative of the quadratic model, beginning with a zero matrix. T
!     arguments of the subroutine are as follows.                       
!                                                                       
!     N must be set to the number of variables and must be at least two.
!     NPT is the number of interpolation conditions. Its value must be i
!       interval [N+2,(N+1)(N+2)/2].                                    
!     Initial values of the variables must be set in X(1),X(2),...,X(N).
!       will be changed to the values that give the least calculated F. 
!     RHOBEG and RHOEND must be set to the initial and final values of a
!       region radius, so both must be positive with RHOEND<=RHOBEG. Typ
!       RHOBEG should be about one tenth of the greatest expected change
!       variable, and RHOEND should indicate the accuracy that is requir
!       the final values of the variables.                              
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls 
!       amount of printing. Specifically, there is no output if IPRINT=0
!       there is output only at the return if IPRINT=1. Otherwise, each 
!       value of RHO is printed, with the best vector of variables so fa
!       the corresponding value of the objective function. Further, each
!       value of F with its variables are output if IPRINT=3.           
!     MAXFUN must be set to an upper bound on the number of calls of CAL
!     The array W will be used for working space. Its length must be at 
!     (NPT+13)*(NPT+N)+3*N*(N+3)/2.                                     
!                                                                       
!     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must se
!     the value of the objective function for the variables X(1),X(2),..
!                                                                       
!     Partition the working space array, so that different parts of it c
!     treated separately by the subroutine that performs the main calcul
!                                                                       
      NP=N+1 
      NPTM=NPT-NP 
      IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN 
          PRINT 10 
   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in',       &
     &      ' the required interval')                                   
          GO TO 20 
      END IF 
      NDIM=NPT+N 
      IXB=1 
      IXO=IXB+N 
      IXN=IXO+N 
      IXP=IXN+N 
      IFV=IXP+N*NPT 
      IGQ=IFV+NPT 
      IHQ=IGQ+N 
      IPQ=IHQ+(N*NP)/2 
      IBMAT=IPQ+NPT 
      IZMAT=IBMAT+NDIM*N 
      ID=IZMAT+NPT*NPTM 
      IVL=ID+N 
      IW=IVL+NDIM 
!                                                                       
!     The above settings provide a partition of W for subroutine NEWUOB.
!     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements 
!     W plus the space that is needed by the last array of NEWUOB.      
!                                                                       
!      CALL NEWUOB (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),          &
!     &  W(IXO),W(IXN),W(IXP),W(IFV),W(IGQ),W(IHQ),W(IPQ),W(IBMAT),      &
!     &  W(IZMAT),NDIM,W(ID),W(IVL),W(IW))                               
      CALL NEWUOB (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),          &
     &  W(IXO),W(IXN),W(IXP),W(IFV),W(IGQ),W(IHQ),W(IPQ),W(IBMAT),      &
     &  W(IZMAT),NDIM,W(ID),W(IVL),W(IW),CALFUN)                               
   20 RETURN 
      END                                           

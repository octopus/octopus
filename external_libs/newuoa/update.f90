      SUBROUTINE UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*) 
!                                                                       
!     The arrays BMAT and ZMAT with IDZ are updated, in order to shift t
!     interpolation point that has index KNEW. On entry, VLAG contains t
!     components of the vector Theta*Wcheck+e_b of the updating formula 
!     (6.11), and BETA holds the value of the parameter that has this na
!     The vector W is used for working space.                           
!                                                                       
!     Set some constants.                                               
!                                                                       
      ONE=1.0D0 
      ZERO=0.0D0 
      NPTM=NPT-N-1 
!                                                                       
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.    
!                                                                       
      JL=1 
      DO 20 J=2,NPTM 
      IF (J .EQ. IDZ) THEN 
          JL=IDZ 
      ELSE IF (ZMAT(KNEW,J) .NE. ZERO) THEN 
          TEMP=DSQRT(ZMAT(KNEW,JL)**2+ZMAT(KNEW,J)**2) 
          TEMPA=ZMAT(KNEW,JL)/TEMP 
          TEMPB=ZMAT(KNEW,J)/TEMP 
          DO 10 I=1,NPT 
          TEMP=TEMPA*ZMAT(I,JL)+TEMPB*ZMAT(I,J) 
          ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,JL) 
   10     ZMAT(I,JL)=TEMP 
          ZMAT(KNEW,J)=ZERO 
      END IF 
   20 END DO 
!                                                                       
!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.             
!                                                                       
      TEMPA=ZMAT(KNEW,1) 
      IF (IDZ .GE. 2) TEMPA=-TEMPA 
      IF (JL .GT. 1) TEMPB=ZMAT(KNEW,JL) 
      DO 30 I=1,NPT 
      W(I)=TEMPA*ZMAT(I,1) 
      IF (JL .GT. 1) W(I)=W(I)+TEMPB*ZMAT(I,JL) 
   30 END DO 
      ALPHA=W(KNEW) 
      TAU=VLAG(KNEW) 
      TAUSQ=TAU*TAU 
      DENOM=ALPHA*BETA+TAUSQ 
      VLAG(KNEW)=VLAG(KNEW)-ONE 
!                                                                       
!     Complete the updating of ZMAT when there is only one nonzero eleme
!     in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to
!     then the first column of ZMAT will be exchanged with another one l
!                                                                       
      IFLAG=0 
      IF (JL .EQ. 1) THEN 
          TEMP=DSQRT(DABS(DENOM)) 
          TEMPB=TEMPA/TEMP 
          TEMPA=TAU/TEMP 
          DO 40 I=1,NPT 
   40     ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I) 
          IF (IDZ .EQ. 1 .AND. TEMP .LT. ZERO) IDZ=2 
          IF (IDZ .GE. 2 .AND. TEMP .GE. ZERO) IFLAG=1 
      ELSE 
!                                                                       
!     Complete the updating of ZMAT in the alternative case.            
!                                                                       
          JA=1 
          IF (BETA .GE. ZERO) JA=JL 
          JB=JL+1-JA 
          TEMP=ZMAT(KNEW,JB)/DENOM 
          TEMPA=TEMP*BETA 
          TEMPB=TEMP*TAU 
          TEMP=ZMAT(KNEW,JA) 
          SCALA=ONE/DSQRT(DABS(BETA)*TEMP*TEMP+TAUSQ) 
          SCALB=SCALA*DSQRT(DABS(DENOM)) 
          DO 50 I=1,NPT 
          ZMAT(I,JA)=SCALA*(TAU*ZMAT(I,JA)-TEMP*VLAG(I)) 
   50     ZMAT(I,JB)=SCALB*(ZMAT(I,JB)-TEMPA*W(I)-TEMPB*VLAG(I)) 
          IF (DENOM .LE. ZERO) THEN 
              IF (BETA .LT. ZERO) IDZ=IDZ+1 
              IF (BETA .GE. ZERO) IFLAG=1 
          END IF 
      END IF 
!                                                                       
!     IDZ is reduced in the following case, and usually the first column
!     of ZMAT is exchanged with a later one.                            
!                                                                       
      IF (IFLAG .EQ. 1) THEN 
          IDZ=IDZ-1 
          DO 60 I=1,NPT 
          TEMP=ZMAT(I,1) 
          ZMAT(I,1)=ZMAT(I,IDZ) 
   60     ZMAT(I,IDZ)=TEMP 
      END IF 
!                                                                       
!     Finally, update the matrix BMAT.                                  
!                                                                       
      DO 70 J=1,N 
      JP=NPT+J 
      W(JP)=BMAT(KNEW,J) 
      TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM 
      TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM 
      DO 70 I=1,JP 
      BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I) 
      IF (I .GT. NPT) BMAT(JP,I-NPT)=BMAT(I,J) 
   70 CONTINUE 
      RETURN 
      END                                           

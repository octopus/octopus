!> @file
!!   Routines to calculate the Faddeeva-function
!! @author
!!      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
!!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!!      VOL. 16, NO. 1, PP. 47.


!>  @brief GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
!!  THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
!!  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I
!!  MEANS SQRT(-1).
!!
!!  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
!!  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
!!  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
!!  OF THE FUNCTION.
!!  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
!!
!!
!!  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
!!     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
!!                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
!!                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
!!                FLOATING-POINT ARITHMETIC
!!     RMAXEXP  = LN(RMAX) - LN(2)
!!     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
!!                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
!!  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL
!!  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
!!
!!
!!  PARAMETER LIST
!!     XI     = REAL      PART OF Z
!!     YI     = IMAGINARY PART OF Z
!!     U      = REAL      PART OF W(Z)
!!     V      = IMAGINARY PART OF W(Z)
!!     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
!!              OCCUR OR NOT; TYPE LOGICAL;
!!              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
!!              MEANING :
!!              FLAG=.FALSE. : NO ERROR CONDITION
!!              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
!!                             BECOMES INACTIVE
!!  XI, YI      ARE THE INPUT-PARAMETERS
!!  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
!!
!!  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
!!
!!  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
!!  PUT TO 0 UPON UNDERFLOW;
!!
!!  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
!!  THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE.
SUBROUTINE WOFZ (XI, YI, U, V, FLAG)

  implicit none
  !IMPLICIT DOUBLE PRECISION (A-H, O-Z)

  real(kind=8), intent(in)      :: xi
  real(kind=8), intent(in)      :: yi
  real(kind=8), intent(out)     :: u
  real(kind=8), intent(out)     :: v
  logical, intent(out)          :: flag
  LOGICAL A, B
  real(kind=8) :: factor, rmaxreal, rmaxexp, rmaxgoni, xabs, yabs, x, y, qrho, xabsq
  real(kind=8) :: xquad, yquad, h, h2, kapn, j, xsum, ysum, xaux, u1, qlambda, rx, v1
  real(kind=8) :: daux, ry, u2, sx, v2, sy, tx, ty, c, w1
  integer :: i, n, nu, np1
  PARAMETER (FACTOR   = 1.12837916709551257388D0,&
       RMAXREAL = 0.5D+154,RMAXEXP  = 708.503061461606D0,&
       RMAXGONI = 3.53711887601422D+15)

  FLAG = .FALSE.

  XABS = DABS(XI)
  YABS = DABS(YI)
  X    = XABS/6.3
  Y    = YABS/4.4


!     THE FOLLOWING IF-STATEMENT PROTECTS
!     QRHO = (X**2 + Y**2) AGAINST OVERFLOW

  IF ((XABS.GT.RMAXREAL).OR.(YABS.GT.RMAXREAL)) then
     FLAG = .TRUE.
     RETURN
  end IF


  QRHO = X**2 + Y**2

  XABSQ = XABS**2
  XQUAD = XABSQ - YABS**2
  YQUAD = 2*XABS*YABS

  A     = QRHO.LT.0.085264D0

  IF (A) THEN
     !IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
     !USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
     !N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
     !ACCURACY
     
     QRHO  = (1-0.85*Y)*DSQRT(QRHO)
     N     = IDNINT(6 + 72*QRHO)
     J     = 2*N+1
     XSUM  = 1.0/J
     YSUM  = 0.0D0
     DO I=N, 1, -1
        J    = J - 2
        XAUX = (XSUM*XQUAD - YSUM*YQUAD)/I
        YSUM = (XSUM*YQUAD + YSUM*XQUAD)/I
        XSUM = XAUX + 1.0/J
     end DO
     U1   = -FACTOR*(XSUM*YABS + YSUM*XABS) + 1.0
     V1   =  FACTOR*(XSUM*XABS - YSUM*YABS)
     DAUX =  DEXP(-XQUAD)
     U2   =  DAUX*DCOS(YQUAD)
     V2   = -DAUX*DSIN(YQUAD)

     U    = U1*U2 - V1*V2
     V    = U1*V2 + V1*U2

  ELSE

     !IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
     !CONTINUED FRACTION
     !NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
     !ACCURACY
     !
     !IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
     !BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
     !IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
     !KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
     !TO OBTAIN THE REQUIRED ACCURACY
     !NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
     !TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY

     IF (QRHO.GT.1.0) THEN
        H    = 0.0D0
        KAPN = 0
        QRHO = DSQRT(QRHO)
        NU   = IDINT(3 + (1442/(26*QRHO+77)))
     ELSE
        QRHO = (1-Y)*DSQRT(1-QRHO)
        H    = 1.88*QRHO
        H2   = 2*H
        KAPN = IDNINT(7  + 34*QRHO)
        NU   = IDNINT(16 + 26*QRHO)
     ENDIF

     B = (H.GT.0.0)

     IF (B) QLAMBDA = H2**KAPN

     RX = 0.0
     RY = 0.0
     SX = 0.0
     SY = 0.0

     DO N=NU, 0, -1
        NP1 = N + 1
        TX  = YABS + H + NP1*RX
        TY  = XABS - NP1*RY
        C   = 0.5/(TX**2 + TY**2)
        RX  = C*TX
        RY  = C*TY
        IF ((B).AND.(N.LE.KAPN)) THEN
           TX = QLAMBDA + SX
           SX = RX*TX - RY*SY
           SY = RY*TX + RX*SY
           QLAMBDA = QLAMBDA/H2
        ENDIF
     end DO

     IF (H.EQ.0.0) THEN
        U = FACTOR*RX
        V = FACTOR*RY
     ELSE
        U = FACTOR*SX
        V = FACTOR*SY
     END IF

     IF (YABS.EQ.0.0) U = DEXP(-XABS**2)

  END IF


  !
  !  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
  !

  IF (YI.LT.0.0) THEN

     IF (A) THEN
        U2    = 2*U2
        V2    = 2*V2
     ELSE
        XQUAD =  -XQUAD

        !         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
        !         AGAINST OVERFLOW

        IF ((YQUAD.GT.RMAXGONI).OR. &
             (XQUAD.GT.RMAXEXP)) then
           FLAG = .TRUE.
           !print *,'BBB',xi,yi,YQUAD,XQUAD,RMAXGONI,RMAXEXP
           RETURN
        end IF

        W1 =  2*DEXP(XQUAD)
        U2  =  W1*DCOS(YQUAD)
        V2  = -W1*DSIN(YQUAD)
     END IF

     U = U2 - U
     V = V2 - V
     IF (XI.GT.0.0) V = -V
  ELSE
     IF (XI.LT.0.0) V = -V
  END IF

END SUBROUTINE WOFZ


!> The same function as above multiplied by exp(-yi**2), to protect against overflows
subroutine wofz_mod(alpha,m,q,jm,u,v,flag)
   use Poisson_Solver, only: dp
   implicit none
   real(dp), intent(in) :: alpha
   integer, intent(in) :: q,m,jm
   logical, intent(out) :: flag
   real(dp), intent(out) :: u,v
   !local variables
   real(dp), parameter :: factor=1.12837916709551257388_dp,rmaxreal=0.5e+154_dp
   real(dp), parameter :: rmaxexp=708.503061461606_dp !n(c) ,rmaxgoni=3.53711887601422e+15_dp
   real(dp), parameter :: pi=3.1415926535897932384_dp
   logical :: a,b
   integer :: j,n,i,kapn,nu,np1
   real(dp) :: xabs,yabs,x,y,qrho,xquad,yquad,xsum,ysum,xaux,u1,v1,u2,v2,daux,h,qlambda,h2
   real(dp) :: rx,ry,sx,sy,tx,ty,c,w1,xabsq,xi,yi,fac,yquadmod
   
   flag = .false.
   
   fac=pi/real(2*m,dp)/alpha
   
   xi=fac*real(q,dp)
   yi=alpha*real(jm,dp)
   xabs = dabs(xi)
   yabs = dabs(yi)
   x    = xabs/6.3_dp
   y    = yabs/4.4_dp
   
   
    ! the following if-statement protects
    ! qrho = (x**2 + y**2) against overflow
   if ((xabs > rmaxreal).or.(yabs > rmaxreal)) then
      flag = .true.
      return
   end if
   
   qrho = x**2 + y**2
   
   xabsq = xabs**2
   xquad = xabsq - yabs**2
   yquad = 2*xabs*yabs
   !yquad to be passed to trigonometric functions
   yquadmod = real(modulo(abs(q*jm),2*m),dp)/real(m,dp)*pi
   
   a = qrho < 0.085264_dp
   
   if (a) then
      !if (qrho < 0.085264d0) then the faddeeva-function is evaluated
      !using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
      !n is the minimum number of terms needed to obtain the required
      !accuracy
      
      qrho  = (1._dp-0.85_dp*y)*sqrt(qrho)
      n     = nint(6._dp + 72._dp*qrho)
      j     = 2*n+1
      xsum  = 1.0_dp/real(j,dp)
      ysum  = 0.0_dp
      do i=n, 1, -1
         j    = j - 2
         xaux = (xsum*xquad - ysum*yquad)/real(i,dp)
         ysum = (xsum*yquad + ysum*xquad)/real(i,dp)
         xsum = xaux + 1.0_dp/real(j,dp)
      end do
      u1   = -factor*(xsum*yabs + ysum*xabs) + 1.0_dp
      v1   =  factor*(xsum*xabs - ysum*yabs)
      !MODIFICATION: the exponent is corrected, yabs disappears
   
      daux = fac**2
      daux = exp(-daux)
      daux = daux**(q**2)
   
      !daux =  exp(-xquad-yabs**2)
      u2   =  daux*dcos(yquadmod)
      v2   = -daux*dsin(yquadmod)
   
      u    = u1*u2 - v1*v2
      v    = u1*v2 + v1*u2
   
   else
   
      !if (qrho.gt.1.0) then w(z) is evaluated using the laplace
      !continued fraction
      !nu is the minimum number of terms needed to obtain the required
      !accuracy
      !
      !if ((qrho.gt.0.085264).and.(qrho.lt.1.0)) then w(z) is evaluated
      !by a truncated taylor expansion, where the laplace continued fraction
      !is used to calculate the derivatives of w(z)
      !kapn is the minimum number of terms in the taylor expansion needed
      !to obtain the required accuracy
      !nu is the minimum number of terms of the continued fraction needed
      !to calculate the derivatives with the required accuracy
   
      if (qrho > 1.0) then
         h    = 0.0_dp
         kapn = 0
         qrho = sqrt(qrho)
         nu   = int(3._dp + (1442._dp/(26._dp*qrho+77._dp)))
      else
         qrho = (1-y)*sqrt(1-qrho)
         h    = 1.88*qrho
         h2   = 2._dp*h
         kapn = nint(7._dp  + 34._dp*qrho)
         nu   = nint(16._dp + 26._dp*qrho)
      endif
   
      b = (h > 0.0_dp)
   
   
   
      if (b) qlambda = h2**kapn
   
      rx = 0.0_dp
      ry = 0.0_dp
      sx = 0.0_dp
      sy = 0.0_dp
   
      do n=nu, 0, -1
         np1 = n + 1
         tx  = yabs + h + real(np1,dp)*rx
         ty  = xabs - real(np1,dp)*ry
         c   = 0.5_dp/(tx**2 + ty**2)
         rx  = c*tx
         ry  = c*ty
         if ((b).and.(n <= kapn)) then
            tx = qlambda + sx
            sx = rx*tx - ry*sy
            sy = ry*tx + rx*sy
            qlambda = qlambda/h2
         endif
         !print *,'nu,yabs,xabs,tx,ty,rx,ry',n,nu,yabs,xabs,tx,ty,rx,ry
      end do
   
      if (h == 0.0_dp) then
         u = factor*rx
         v = factor*ry
      else
         u = factor*sx
         v = factor*sy
      end if
   
      !MODIFICATION: here the exponent is added
      daux=exp(-yabs**2)
      u=daux*u
      v=daux*v
   
      !still do not now what happens to v in that case
      if (yabs == 0.0_dp) u = dexp(-xabs**2)
   
   end if
   
   
   !  evaluation of w(z) in the other quadrants
   if (jm < 0) then
   
      if (a) then
         !no modification is needed
         u2    = 2._dp*u2
         v2    = 2._dp*v2
      else
         xquad =  -xquad
   
         !the following if-statement protects 2*exp(-z**2)
         !against overflow, taking into account the modification
         if (xquad-yabs**2 > rmaxexp) then
            flag = .true.
            return
         end if
         !daux = fac**2
         !daux = exp(-daux)
         !w1 = 2.0_dp*daux**(q**2)
         !avoid floating point exceptions
         daux = real(q,kind=8)*fac
         daux = daux*daux
         daux = dexp(-daux)
         w1 = 2.0_dp*daux    
         !w1 =  2*dexp(xquad-yabs**2)
         u2  =  w1*dcos(yquadmod)
         v2  = -w1*dsin(yquadmod)
         !print *,'w1,u2,v2',yquad,yquadmod,xquad,yabs**2,xquad-yabs**2,w1,u2,v2
      end if
   
      u = u2 - u
      v = v2 - v
      !print *,u,v
      if (q > 0) v = -v
   else
      if (q < 0) v = -v
   end if
END SUBROUTINE wofz_mod



!! The code in this file is taken from the expokit distribution,
!! The subroutines have been modified in order to make them Fortran 90.


!!/*                              NOTICE
!!
!!Permission to use, copy, modify, and distribute EXPOKIT and its
!!supporting documentation for non-commercial purposes, is hereby
!!granted without fee, provided that this permission message and
!!copyright notice appear in all copies. Approval must be sought for
!!commercial purposes.
!!
!!Neither the Institution (University of Queensland) nor the Authors
!!make any representations about the suitability of this software for
!!any purpose.  This software is provided ``as is'' without express or
!!implied warranty.
!!
!!The work resulting from EXPOKIT has been published in ACM-Transactions 
!!on Mathematical Software, 24(1):130-156, 1998.
!!
!!Certain elements of the current software may include inadequacies
!!that may be corrected at any time, as they are discovered. The Web 
!!always contains the latest updates.
!!
!!                             Roger B. Sidje
!!           Department of Mathematics, University of Queensland 
!!     Brisbane, QLD-4072, Australia, (c) 1996-1999 All Rights Reserved */

!----------------------------------------------------------------------|
!----------------------------------------------------------------------|
      subroutine ZGPADM(m, t, H, exph, iflag)

      ! I have changed the interface from the original expokit distribution.
      implicit none
      FLOAT, intent(in) :: t
      integer :: m, iflag
      CMPLX, intent(in)  :: H(m,m)
      CMPLX, intent(out) :: expH(m, m)

      integer, parameter :: ideg = 6
      integer :: ns, ldh, lwsp, iexph
      integer, allocatable :: ipiv(:)
      CMPLX, allocatable :: wsp(:)

!/*-----Purpose----------------------------------------------------------|
!!$
!!$     Computes exp(t*H), the matrix exponential of a general complex 
!!$     matrix in full, using the irreducible rational Pade approximation
!!$     to the exponential exp(z) = r(z) = (+/-)( I + 2*(q(z)/p(z)) ),
!!$     combined with scaling-and-squaring.
!!$
!!$-----Arguments--------------------------------------------------------|
!!$
!!$     ideg      : (input) the degre of the diagonal Pade to be used.
!!$                 a value of 6 is generally satisfactory.
!!$
!!$     m         : (input) order of H.
!!$
!!$     H(ldh,m)  : (input) argument matrix.
!!$
!!$     t         : (input) time-scale (can be < 0).
!!$                  
!!$     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
!!$
!!$     ipiv(m)   : (workspace)
!!$
!!$>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
!!$                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
!!$                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!!$                 NOTE: if the routine was called with wsp(iptr), 
!!$                       then exp(tH) will start at wsp(iptr+iexph-1).
!!$
!!$     ns        : (output) number of scaling-squaring used.
!!$
!!$     iflag     : (output) exit flag.
!!$                       0 - no problem
!!$                      <0 - problem
!!$
!!$----------------------------------------------------------------------|
!!$     Roger B. Sidje (rbs@maths.uq.edu.au)
!!$     EXPOKIT: Software Package for Computing Matrix Exponentials.
!!$     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
!!$----------------------------------------------------------------------|*/

      integer i,j,k,icoef,mm,ih2,iodd,iused,ifree,iq,ip,iput,iget
      double precision hnorm
      CMPLX :: cp, cq, scale, scale2, ZERO, ONE

      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      intrinsic ABS, cmplx, DBLE, INT, LOG, MAX

      ! This I have added to simplify the interface...
      ldh = m
      lwsp = 4*m*m+ideg+1
      allocate(ipiv(m), wsp(lwsp))

!---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
      if ( iflag.ne.0 ) stop 'bad sizes (in input of ZGPADM)'
!
!---  initialise pointers ...
!
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
!
!---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; 
!     and set scale = t/2^ns ...

      do i = 1,m
         wsp(i) = ZERO
      enddo
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
         enddo
      enddo
      hnorm = 0.0d0
      do i = 1,m
         hnorm = MAX( hnorm,DBLE(wsp(i)) )
      enddo
      hnorm = ABS( t*hnorm )
      if ( hnorm.eq.0.0d0 ) stop 'Error - null H in input of ZGPADM.'
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
      scale =  cmplx( t/DBLE(2**ns),0.0d0 )
      scale2 = scale*scale

!---  compute Pade coefficients ...

      i = ideg+1
      j = 2*ideg+1
      wsp(icoef) = ONE
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
      enddo

!---  H2 = scale2*H*H ...
!
!!$      call ZGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,ZERO,wsp(ih2),m )
      call ZGEMM( 'n','n',m,m,m,scale2,H(1, 1),ldh,H(1, 1),ldh,ZERO,wsp(ih2),m )

!---  initialise p (numerator) and q (denominator) ...

      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = ZERO
            wsp(iq + (j-1)*m + i-1) = ZERO
         enddo
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      enddo

!---  Apply Horner rule ...

      iodd = 1
      k = ideg - 1
 100  continue
      iused = iodd*iq + (1-iodd)*ip
!!$      call ZGEMM( 'n','n',m,m,m, ONE,wsp(iused:),m, &
!!$                   wsp(ih2:),m, ZERO,wsp(ifree:),m )
      call blas_gemm( 'n','n',m,m,m, ONE,wsp(iused),m, &
                   wsp(ih2),m, ZERO,wsp(ifree),m )
      do j = 1,m
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
      enddo
      ip = (1-iodd)*ifree + iodd*ip
      iq = iodd*ifree + (1-iodd)*iq
      ifree = iused
      iodd = 1-iodd
      k = k-1
      if ( k.gt.0 )  goto 100

!---  Obtain (+/-)(I + 2*(p\q)) ...

      if ( iodd.ne.0 ) then
!!$         call ZGEMM( 'n','n',m,m,m, scale,wsp(iq),m, &
!!$                      H,ldh, ZERO,wsp(ifree),m )
         call blas_gemm( 'n','n',m,m,m, scale,wsp(iq),m, &
                      H(1, 1),ldh, ZERO,wsp(ifree),m )
         iq = ifree
      else
!!$         call ZGEMM( 'n','n',m,m,m, scale,wsp(ip),m, &
!!$                      H,ldh, ZERO,wsp(ifree),m )
         call blas_gemm( 'n','n',m,m,m, scale,wsp(ip),m, &
                      H(1, 1),ldh, ZERO,wsp(ifree),m )
         ip = ifree
      endif
      call ZAXPY( mm, -ONE,wsp(ip),1, wsp(iq),1 )
      call ZGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
      if ( iflag.ne.0 ) stop 'Problem in ZGESV (within ZGPADM)'
      call ZDSCAL( mm, 2.0d0, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + ONE
      enddo
      iput = ip
      if ( ns.eq.0 .and. iodd.ne.0 ) then
         call ZDSCAL( mm, -1.0d0, wsp(ip), 1 )
         goto 200
      endif

!--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...

      iodd = 1
      do k = 1,ns
         iget = iodd*ip + (1-iodd)*iq
         iput = (1-iodd)*ip + iodd*iq
!!$         call ZGEMM( 'n','n',m,m,m, ONE,wsp(iget),m, wsp(iget),m, &
!!$                      ZERO,wsp(iput),m )
         call blas_gemm( 'n','n',m,m,m, ONE,wsp(iget),m, wsp(iget),m, &
                      ZERO,wsp(iput),m )
         iodd = 1-iodd
      enddo
 200  continue
      iexph = iput

      k = 0
      do i = 1, m
         do j = 1, m
            exph(j, i) = wsp(iexph + k)
            k = k + 1
         enddo
      enddo

      deallocate(ipiv, wsp)
      end subroutine ZGPADM

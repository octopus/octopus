!! $Id$
!!

#if defined(SINGLE_PRECISION)
#  define DLAPACK(x) s ## x
#  define ZLAPACK(x) c ## x
#else
#  define DLAPACK(x) d ## x
#  define ZLAPACK(x) z ## x
#endif

! computes all the eigenvalues and the eigenvectors of a real
! generalized symmetric-definite eigenproblem, of the form  A*x=(lambda)*B*x
! A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
! Here A and B are assumed to be symmetric and B is also positive definite.
subroutine dgeneigensolve(n, a, b, e)
  integer, intent(in)    :: n
  FLOAT,   intent(inout) :: a(n,n)
  FLOAT,   intent(in)    :: b(n,n)
  FLOAT,   intent(out)   :: e(n)

  interface
    subroutine DLAPACK(sygv) (itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: itype, n, lda, ldb, lwork
      FLOAT,        intent(inout) :: a, b    ! a(lda,n), b(ldb,n)
      FLOAT,        intent(out)   :: w, work ! w(n), work(lwork)
      integer,      intent(out)   :: info
    end subroutine DLAPACK(sygv)
  end interface

  integer :: info, lwork
  FLOAT, allocatable :: bp(:,:), work(:)

  lwork = 5*n
  allocate(bp(n, n), work(lwork))
  bp = b
  call DLAPACK(sygv) (1, 'V', 'U', n, a(1, 1), n, bp(1, 1), n, e(1), work(1), lwork, info)
  deallocate(bp, work)

  if(info.ne.0) then
    write(message(1),'(a,i5)') 'In dgeneigensolve, LAPACK dsygv returned error message ', info
    call write_fatal(1)
  end if

end subroutine dgeneigensolve

! computes all the eigenvalues and the eigenvectors of a complex
! generalized Hermitian-definite eigenproblem, of the form  A*x=(lambda)*B*x,
! A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
! Here A and B are assumed to be Hermitian and B is also positive definite.
subroutine zgeneigensolve(n, a, b, e)
  integer, intent(in)    :: n
  CMPLX,   intent(inout) :: a(n,n)
  CMPLX,   intent(in)    :: b(n,n)
  FLOAT,   intent(out)   :: e(n)

  interface
    subroutine ZLAPACK(hegv) (itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, itype, lda, ldb, lwork
      CMPLX,        intent(inout) :: a, b     ! a(lda,n), b(ldb,n)
      FLOAT,        intent(out)   :: w, rwork ! w(n), rwork(max(1,3*n-2))
      CMPLX,        intent(out)   :: work     ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine ZLAPACK(hegv)
  end interface

  integer :: info, lwork
  FLOAT, allocatable :: rwork(:)
  CMPLX, allocatable :: bp(:,:), work(:)

  lwork = 5*n
  allocate(bp(n, n), work(lwork), rwork(max(1, 3*n-2)))
  bp = b
  call ZLAPACK(hegv) (1, 'V', 'U', n, a(1, 1), n, bp(1, 1), n, e(1), work(1), lwork, rwork(1), info)
  deallocate(bp, work)

  if(info.ne.0) then
    write(message(1),'(a,i5)') 'In dgeneigensolve, LAPACK dhegv returned error message ', info
    call write_fatal(1)
  end if

end subroutine zgeneigensolve

! computes all eigenvalues and eigenvectors of a real symmetric square matrix A.
subroutine deigensolve(n, a, b, e)
  integer, intent(in)  :: n
  FLOAT,   intent(in)  :: a(n,n)
  FLOAT,   intent(out) :: b(n,n)
  FLOAT,   intent(out) :: e(n)

  interface
    subroutine DLAPACK(syev) (jobz, uplo, n, a, lda, w, work, lwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, lda, lwork
      FLOAT,        intent(inout) :: a       ! a(lda,n)
      FLOAT,        intent(out)   :: w, work ! w(n), work(lwork)
      integer,      intent(out)   :: info
    end subroutine DLAPACK(syev)
  end interface

  integer :: info, lwork
  FLOAT, allocatable :: work(:)

  lwork = 6*n
  allocate(work(lwork))
  b = a
  call DLAPACK(syev) ('V', 'U', n, b(1,1), n, e(1), work(1), lwork, info)
  if(info.ne.0) then
    write(message(1),'(a,i5)') 'In deigensolve, LAPACK dsyev returned error message ', info
    call write_fatal(1)
  end if
  deallocate(work)

end subroutine deigensolve

! computes all eigenvalues and eigenvectors of a complex Hermitian square matrix A.
subroutine zeigensolve(n, a, b, e)
  integer, intent(in)  :: n
  CMPLX,   intent(in)  :: a(n,n)
  CMPLX,   intent(out) :: b(n,n)
  FLOAT,   intent(out) :: e(n)

  interface
    subroutine ZLAPACK(heev) (jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, lda, lwork
      CMPLX,        intent(inout) :: a        ! a(lda,n)
      FLOAT,        intent(out)   :: w, rwork ! w(n), rwork(max(1,3*n-2))
      CMPLX,        intent(out)   :: work     ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine ZLAPACK(heev)
  end interface

  integer :: info, lwork
  CMPLX, allocatable :: work(:)
  FLOAT, allocatable :: rwork(:)

  lwork = 6*n
  allocate(work(lwork), rwork(max(1,3*n-2)))
  b = a
  call ZLAPACK(heev) ('V','U', n, b(1,1), n, e(1), work(1), lwork, rwork(1), info)
  if(info.ne.0) then
    write(message(1),'(a,i5)') 'In zeigensolve, LAPACK zheev returned error message ', info
    call write_fatal(1)
  end if
  deallocate(work, rwork)

end subroutine zeigensolve

! return the determinant of a general square matrix a of dimensions (n,n)
FLOAT function ddet(a, n)
  integer, intent(in)    :: n
  FLOAT,   intent(inout) :: a(n,n)

  interface
    subroutine DLAPACK(getrf) (m, n, a, lda, ipiv, info)
      integer, intent(in)    :: m, n, lda
      FLOAT,   intent(inout) :: a    ! a(lda,n)
      integer, intent(out)   :: ipiv ! ipiv(min(m,n))
      integer, intent(out)   :: info
    end subroutine DLAPACK(getrf)
  end interface

  integer :: i, info, ipiv(n)

  call DLAPACK(getrf) (n, n, a(1,1), n, ipiv(1), info)
  if(info.ne.0) then
    write(message(1),'(a,i4)') 'In ddet, LAPACK dgetrf returned unsuccesful code ',info
    ddet = M_ZERO
    call write_warning(1)
    return
  endif

  ddet = M_ONE
  do i = 1, n
     if(ipiv(i).ne.i) then
       ddet = - ddet*a(i, i)
     else
       ddet = ddet*a(i, i)
     endif
  end do

end function ddet

! return the determinant of a general square matrix a of dimensions (n,n)
CMPLX function zdet(a, n)
  integer, intent(in)    :: n
  CMPLX,   intent(inout) :: a(n,n)

  interface
    subroutine ZLAPACK(getrf) (m, n, a, lda, ipiv, info)
      integer, intent(in)    :: m, n, lda
      CMPLX,   intent(inout) :: a    ! a(lda,n)
      integer, intent(out)   :: ipiv ! ipiv(min(m,n))
      integer, intent(out)   :: info
    end subroutine ZLAPACK(getrf)
  end interface

  integer :: info, i, ipiv(n)

  call ZLAPACK(getrf) (n, n, a(1,1), n, ipiv(1), info)
  if(info < 0) then
    write(message(1),'(a,i4)') 'In zdet, LAPACK zgetrf returned unsuccesfull code ',info
    zdet = M_ZERO
    call write_warning(1)
    return
  endif

  zdet = M_z1
  do i = 1, n
     if(ipiv(i).ne.i) then
       zdet = - zdet*a(i, i)
     else
       zdet = zdet*a(i, i)
     endif
  end do

end function zdet

! Invert a real symmetric square matrix a
subroutine dinvert(n, a)
  integer, intent(in)    :: n
  FLOAT,   intent(inout) :: a(n,n)

  interface
    subroutine DLAPACK(sytrf) (uplo, n, a, lda, ipiv, work, lwork, info)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda, lwork
      FLOAT,        intent(inout) :: a      ! a(lda,n)
      integer,      intent(out)   :: ipiv   ! ipiv(min(m,n))
      FLOAT,        intent(out)   :: work   ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine DLAPACK(sytrf)

    subroutine DLAPACK(sytri) (uplo, n, a, lda, ipiv, work, info)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda
      FLOAT,        intent(inout) :: a      ! a(lda,n)
      integer,      intent(out)   :: ipiv   ! ipiv(min(m,n))
      FLOAT,        intent(out)   :: work   ! work(n)
      integer,      intent(out)   :: info
    end subroutine DLAPACK(sytri)
  end interface

  integer :: info, i, j
  integer, allocatable :: ipiv(:)
  FLOAT, allocatable :: work(:)

  allocate(work(n), ipiv(n))

  call DLAPACK(sytrf) ('u', n, a(1, 1), n, ipiv(1), work(1), n, info)
  if(info /= 0) then
    write(message(1), '(a, i3)') 'In dinvert, LAPACK dsytrf returned info = ', info
    call write_fatal(1)
  end if

  call DLAPACK(sytri) ('u', n, a(1, 1), n, ipiv(1), work(1), info)
  if(info /= 0) then
    write(message(1), '(a, i3)') 'In dinvert, LAPACK dsytri returned info = ', info
    call write_fatal(1)
  end if

  deallocate(work, ipiv)

  ! complete the matrix
  do i = 1, n
    do j = i + 1, n
      a(j, i) = a(i, j)
    end do
  end do

end subroutine dinvert

! compute the solution to a real system of linear equations A*X = B,
!  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
subroutine dlinsyssolve(n, nhrs, a, b, x)
  integer, intent(in)    :: n, nhrs
  FLOAT,   intent(inout) :: a(n, n), b(n, nhrs)
  FLOAT,   intent(out)   :: x(n, nhrs)

  interface
    subroutine DLAPACK(gesvx) (fact, trans, n, nhrs, a, lda, af, ldaf, ipiv, equed, r, &
                               c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info)
      character(1), intent(in)    :: fact, trans
      integer,      intent(in)    :: n, nhrs, lda, ldaf, ldb, ldx
      FLOAT,        intent(inout) :: a, af, r, c, b      ! a(lda,n), af(ldaf,n), r(n), c(n), b(ldb,nhrs)
      integer,      intent(inout) :: ipiv                ! ipiv(n)
      FLOAT,        intent(out)   :: x, ferr, berr, work ! x(ldx,nhrs), ferr(nhrs), berr(nhrs), work(4*n)
      FLOAT,        intent(out)   :: rcond
      character(1), intent(inout) :: equed
      integer,      intent(out)   :: iwork               ! iwork(n)
      integer,      intent(out)   :: info
    end subroutine DLAPACK(gesvx)
  end interface

  integer :: info
  integer, allocatable :: ipiv(:), iwork(:)
  FLOAT :: rcond
  FLOAT, allocatable :: ferr(:), berr(:), work(:), r(:), c(:), af(:,:)
  character(1) :: equed

  allocate(ipiv(n), iwork(n), ferr(nhrs), berr(nhrs), work(4*n), r(n), c(n), af(n, n))
  call DLAPACK(gesvx) ("N", "N", n, nhrs, a(1, 1), n, af(1, 1), n, ipiv(1), equed, r(1), c(1), b(1, 1), n, x(1, 1), n, &
                       rcond, ferr(1), berr(1), work(1), iwork(1), info)
  if(info /= 0) then
    write(message(1), '(a, i3)') 'In dlinsyssolve, LAPACK dgesvx returned info = ', info
    call write_fatal(1)
  end if
  deallocate(ipiv, iwork, ferr, berr, work, r, c, af)


end subroutine dlinsyssolve

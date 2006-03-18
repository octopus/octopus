! $Id$
!

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
  ALLOCATE(bp(n, n), n*n)
  ALLOCATE(work(lwork), lwork)
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
  ALLOCATE(bp(n, n),    n*n)
  ALLOCATE(work(lwork), lwork)
  ALLOCATE(rwork(max(1, 3*n-2)), max(1, 3*n-2))
  bp = b
  call ZLAPACK(hegv) (1, 'V', 'U', n, a(1, 1), n, bp(1, 1), n, e(1), work(1), lwork, rwork(1), info)
  deallocate(bp, work, rwork)

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
  ALLOCATE(work(lwork), lwork)
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
  ALLOCATE(work(lwork), lwork)
  ALLOCATE(rwork(max(1, 3*n-2)), max(1, 3*n-2))
  b = a
  call ZLAPACK(heev) ('V','U', n, b(1,1), n, e(1), work(1), lwork, rwork(1), info)
  if(info.ne.0) then
    write(message(1),'(a,i5)') 'In zeigensolve, LAPACK zheev returned error message ', info
    call write_fatal(1)
  end if
  deallocate(work, rwork)

end subroutine zeigensolve


! Invert a real symmetric square matrix a
FLOAT function ddeterminant(n, a, invert) result(d)
  integer, intent(in)           :: n
  FLOAT,   intent(inout)        :: a(n,n)
  logical, intent(in), optional :: invert

  interface
    subroutine DLAPACK(getrf) (m, n, a, lda, ipiv, info)
      integer,      intent(in)    :: m, n, lda
      FLOAT,        intent(inout) :: a          ! a(lda, n)
      integer,      intent(out)   :: ipiv       ! ipiv(min(m,n)
      integer,      intent(out)   :: info
    end subroutine DLAPACK(getrf)

    subroutine DLAPACK(getri) (n, a, lda, ipiv, work, lwork, info )
      integer,      intent(in)    :: n, lda, lwork
      FLOAT,        intent(inout) :: a       ! a(lda, n)
      integer,      intent(in)    :: ipiv    ! ipiv(n)
      FLOAT,        intent(inout) :: work    ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine DLAPACK(getri)
  end interface

  integer :: info, i
  integer, allocatable :: ipiv(:)
  FLOAT, allocatable :: work(:)
  logical :: invert_

  ALLOCATE(work(n), n)
  ALLOCATE(ipiv(n), n)

  call DLAPACK(getrf)(n, n, a(1, 1), n, ipiv(1), info)
  if(info /= 0) then
    write(message(1), '(a, i3)') 'In dinvert, LAPACK dgetrf returned info = ', info
    call write_fatal(1)
  end if

  d = M_ONE
  do i = 1, n
    if(ipiv(i).ne.i) then
      d = - d*a(i, i)
    else
      d = d*a(i, i)
    end if
  end do

  invert_ = .true.; if(present(invert)) invert_ = invert
  if(invert_) then
    call DLAPACK(getri)(n, a(1, 1), n, ipiv(1), work(1), n, info)
    if(info /= 0) then
      write(message(1), '(a, i3)') 'In dinvert, LAPACK dgetri returned info = ', info
      call write_fatal(1)
    end if
  end if

  deallocate(work, ipiv)
end function ddeterminant


! Invert a real symmetric square matrix a
CMPLX function zdeterminant(n, a, invert) result(d)
  integer, intent(in)           :: n
  CMPLX,   intent(inout)        :: a(n,n)
  logical, intent(in), optional :: invert

  interface
    subroutine ZLAPACK(getrf) (m, n, a, lda, ipiv, info)
      integer,      intent(in)    :: m, n, lda
      CMPLX,        intent(inout) :: a          ! a(lda, n)
      integer,      intent(out)   :: ipiv       ! ipiv(min(m,n)
      integer,      intent(out)   :: info
    end subroutine ZLAPACK(getrf)

    subroutine ZLAPACK(getri) (n, a, lda, ipiv, work, lwork, info )
      integer,      intent(in)    :: n, lda, lwork
      CMPLX,        intent(inout) :: a       ! a(lda, n)
      integer,      intent(in)    :: ipiv    ! ipiv(n)
      CMPLX,        intent(inout) :: work    ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine ZLAPACK(getri)
  end interface

  integer :: info, i
  integer, allocatable :: ipiv(:)
  CMPLX, allocatable :: work(:)
  logical :: invert_

  ALLOCATE(work(n), n)
  ALLOCATE(ipiv(n), n)

  call ZLAPACK(getrf)(n, n, a(1, 1), n, ipiv(1), info)
  if(info /= 0) then
    write(message(1), '(a, i3)') 'In dinvert, LAPACK dgetrf returned info = ', info
    call write_fatal(1)
  end if

  d = M_ONE
  do i = 1, n
    if(ipiv(i).ne.i) then
      d = - d*a(i, i)
    else
      d = d*a(i, i)
    end if
  end do

  invert_ = .true.; if(present(invert)) invert_ = invert
  if(invert_) then
    call ZLAPACK(getri)(n, a(1, 1), n, ipiv(1), work(1), n, info)
    if(info /= 0) then
      write(message(1), '(a, i3)') 'In dinvert, LAPACK dgetri returned info = ', info
      call write_fatal(1)
    end if
  end if

  deallocate(work, ipiv)
end function zdeterminant


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

  ALLOCATE(ipiv(n),    n)
  ALLOCATE(iwork(n),   n)
  ALLOCATE(ferr(nhrs), nhrs)
  ALLOCATE(berr(nhrs), nhrs)
  ALLOCATE(work(4*n),  4*n)
  ALLOCATE(r(n),       n)
  ALLOCATE(c(n),       n)
  ALLOCATE(af(n, n),   n*n)

  call DLAPACK(gesvx) ("N", "N", n, nhrs, a(1, 1), n, af(1, 1), n, ipiv(1), equed, r(1), c(1), b(1, 1), n, x(1, 1), n, &
    rcond, ferr(1), berr(1), work(1), iwork(1), info)

  if(info /= 0) then
    write(message(1), '(a, i3)') 'In dlinsyssolve, LAPACK dgesvx returned info = ', info
    call write_fatal(1)
  end if

  deallocate(ipiv, iwork, ferr, berr, work, r, c, af)

end subroutine dlinsyssolve

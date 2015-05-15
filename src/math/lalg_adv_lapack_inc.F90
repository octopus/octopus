!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

! ---------------------------------------------------------
!> Compute the Cholesky decomposition of real symmetric or complex Hermitian positive definite
!! matrix a, dim(a) = n x n. On return a = u^T u with u upper triangular matrix.
subroutine X(cholesky)(n, a, bof, err_code)
  integer,           intent(in)    :: n
  R_TYPE,            intent(inout) :: a(:,:)   !< (n,n)
  logical, optional, intent(inout) :: bof      !< Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer :: info

  call profiling_in(cholesky_prof, "CHOLESKY")
  PUSH_SUB(X(cholesky))

  ASSERT(n > 0)

  call lapack_potrf('U', n, a(1, 1), lead_dim(a), info)
  if(info /= 0) then
    if(optional_default(bof, .true.)) then
      write(message(1), '(5a,i5)') 'In ', TOSTRING(X(cholesky)), ' LAPACK ', TOSTRING(X(potrf)), ' returned error message ', info
      call messages_fatal(1)
    else
      if(present(bof)) then
        bof = .true.
      end if
    end if
  else
    if(present(bof)) then
      bof = .false.
    end if
  end if
  if(present(err_code)) then
    err_code = info
  end if

  call profiling_out(cholesky_prof)
  POP_SUB(X(cholesky))
end subroutine X(cholesky)


! ---------------------------------------------------------
!> Computes all the eigenvalues and the eigenvectors of a real symmetric or complex Hermitian
!! generalized definite eigenproblem, of the form \f$ Ax=\lambda Bx \f$. B is also positive definite.
subroutine X(geneigensolve)(n, a, b, e, bof, err_code)
  integer,           intent(in)    :: n
  R_TYPE,            intent(inout) :: a(:,:)   !< (n,n)
  R_TYPE,            intent(inout) :: b(:,:)   !< (n,n)
  FLOAT,             intent(out)   :: e(:)     !< (n)
  logical, optional, intent(inout) :: bof      !< Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer :: info, lwork, ii, jj
  FLOAT, allocatable :: rwork(:)
  R_TYPE, allocatable :: work(:), diag(:)

  call profiling_in(eigensolver_prof, "DENSE_EIGENSOLVER")
  PUSH_SUB(X(geneigensolve))

  ASSERT(n > 0)

  SAFE_ALLOCATE(diag(1:n))

  ! store the diagonal of b  
  forall(ii = 1:n) diag(ii) = b(ii, ii)

  lwork = 5*n
  SAFE_ALLOCATE(work(1:lwork))
#ifdef R_TCOMPLEX
  SAFE_ALLOCATE(rwork(1:max(1, 3*n-2)))
  call lapack_hegv(1, 'V', 'U', n, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), e(1), work(1), lwork, rwork(1), info)
  SAFE_DEALLOCATE_A(rwork)
#else
  call lapack_sygv(1, 'V', 'U', n, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), e(1), work(1), lwork, info)
#endif
  SAFE_DEALLOCATE_A(work)

  ! b was destroyed, so we rebuild it
  do ii = 1, n
    do jj = 1, ii - 1
      b(jj, ii) = b(ii, jj)
    end do
    b(ii, ii) = diag(ii)
  end do

  SAFE_DEALLOCATE_A(diag)

  if(info /= 0) then
    if(optional_default(bof, .true.)) then
      write(message(1),'(3a)') 'In ', TOSTRING(X(geneigensolve)), ' LAPACK '
#ifdef R_TCOMPLEX
      write(message(1),'(3a,i5)') trim(message(1)), TOSTRING(X(hegv)), ' returned error message ', info
#else
      write(message(1),'(3a,i5)') trim(message(1)), TOSTRING(X(sygv)), ' returned error message ', info
#endif
      call messages_fatal(1)
    else
      if(present(bof)) then
        bof = .true.
      end if
    end if
  else
    if(present(bof)) then
      bof = .false.
    end if
  end if
  if(present(err_code)) then
    err_code = info
  end if
  
  call profiling_out(eigensolver_prof)
  POP_SUB(X(geneigensolve))
end subroutine X(geneigensolve)


! ---------------------------------------------------------
!> Computes all the eigenvalues and the right (left) eigenvectors of a real or complex
!! (non-Hermitian) eigenproblem, of the form A*x=(lambda)*x
subroutine X(eigensolve_nonh)(n, a, e, err_code, side, sort_eigenvectors)
  integer,                intent(in)    :: n
  R_TYPE,                 intent(inout) :: a(:, :)   !< (n,n)
  CMPLX,                  intent(out)   :: e(:)      !< (n)
  integer,      optional, intent(out)   :: err_code
  character(1), optional, intent(in)    :: side      !< which eigenvectors ('L' or 'R')
  logical,      optional, intent(in)    :: sort_eigenvectors !< only applies to complex version, sorts by real part

  integer              :: info, lwork, ii
  FLOAT, allocatable   :: rwork(:), re(:)
  R_TYPE, allocatable  :: work(:), vl(:, :), vr(:, :), a_copy(:, :)
  CMPLX, allocatable   :: e_copy(:)
  character(1)         :: side_
  integer, allocatable :: ind(:)

  PUSH_SUB(X(eigensolve_nonh))

  ASSERT(n > 0)

  if (present(side)) then
    side_ = side
  else
    side_ = 'R'
  end if

  lwork = -1
  ! A bug in the query mode of zgeev demands that the working array has to be larger than 1
  ! problem here is that a value is written somewhere into the array whether it is
  ! allocated or not. I noticed that it happens (hopefully) always at an index which
  ! is below the matrix dimension.
  SAFE_ALLOCATE(work(1:n))
  SAFE_ALLOCATE(vl(1:1, 1:1))
  SAFE_ALLOCATE(vr(1:n, 1:n)) ! even in query mode, the size of vr is checked, so we allocate it
  SAFE_ALLOCATE(rwork(1:1))
  call lapack_geev('N', 'V', n, a, lead_dim(a), e, vl, lead_dim(vl), vr, lead_dim(vr), &
    work, lwork, rwork, info)

  lwork = int(work(1))
  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(vl)
  SAFE_DEALLOCATE_A(vr)
  SAFE_DEALLOCATE_A(rwork)

  SAFE_ALLOCATE(work(1:lwork))
  SAFE_ALLOCATE(rwork(1:max(1, 2*n)))
  if(side_ == 'L'.or.side_ == 'l') then
    SAFE_ALLOCATE(vl(1:n, 1:n))
    SAFE_ALLOCATE(vr(1:1, 1:1))
    call lapack_geev('V', 'N', n, a, lead_dim(a), e, vl, lead_dim(vl), vr, lead_dim(vr), &
      work, lwork, rwork, info) ! check info status
    a(1:n, 1:n) = vl(1:n, 1:n)
  else
    SAFE_ALLOCATE(vl(1:1, 1:1))
    SAFE_ALLOCATE(vr(1:n, 1:n))
    call lapack_geev('N', 'V', n, a, lead_dim(a), e, vl, lead_dim(vl), vr, lead_dim(vr), &
      work, lwork, rwork, info)
    a(1:n, 1:n) = vr(1:n, 1:n)
  end if
  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(rwork)
  SAFE_DEALLOCATE_A(vr)
  SAFE_DEALLOCATE_A(vl)

  if(info /= 0) then
    write(message(1),'(5a,i5)') 'In ', TOSTRING(X(eigensolve_nonh)), &
      ' LAPACK ', TOSTRING(X(geev)), ' returned error message ', info
    call messages_fatal(1)
  end if
  if(present(err_code)) then
    err_code = info
  end if

  if(optional_default(sort_eigenvectors, .false.)) then
    SAFE_ALLOCATE(re(1:n))
    SAFE_ALLOCATE(ind(1:n))
    SAFE_ALLOCATE(e_copy(1:n))
    SAFE_ALLOCATE(a_copy(1:n, 1:n))
    re = real(e, REAL_PRECISION)
    e_copy = e
    a_copy = a
    call sort(re, ind)
    forall(ii = 1:n)
      e(ii) = e_copy(ind(ii))
      a(1:n, ii) = a_copy(1:n, ind(ii))
    end forall
    SAFE_DEALLOCATE_A(e_copy)
    SAFE_DEALLOCATE_A(a_copy)
    SAFE_DEALLOCATE_A(re)
    SAFE_DEALLOCATE_A(ind)
  end if

  POP_SUB(X(eigensolve_nonh))
end subroutine X(eigensolve_nonh)


#ifdef R_TREAL
! ---------------------------------------------------------
!> Computes the k lowest eigenvalues and the eigenvectors of a real symmetric or complex Hermitian
!! generalized definite eigenproblem, of the form  A*x=(lambda)*B*x. B is also positive definite.
subroutine dlowest_geneigensolve(k, n, a, b, e, v, bof, err_code)
  integer,           intent(in)    :: k, n
  FLOAT,             intent(inout) :: a(:,:)   !< (n, n)
  FLOAT,             intent(inout) :: b(:,:)   !< (n, n)
  FLOAT,             intent(out)   :: e(:)     !< (n)
  FLOAT,             intent(out)   :: v(:,:)   !< (n, n)
  logical, optional, intent(inout) :: bof      !< Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer            :: m, iwork(5*n), ifail(n), info, lwork
  FLOAT              :: abstol
  FLOAT, allocatable :: work(:)
  
  PUSH_SUB(dlowest_geneigensolve)

  ASSERT(n > 0)

  abstol = 2*sfmin()

  ! Work size query.
  SAFE_ALLOCATE(work(1:1))
  call X(sygvx)(1, 'V', 'I', 'U', n, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), lead_dim(v), work(1), -1, iwork(1), ifail(1), info) ! check info
  lwork = int(work(1))
  SAFE_DEALLOCATE_A(work)

  SAFE_ALLOCATE(work(1:lwork))

  call X(sygvx)(1, 'V', 'I', 'U', n, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), lead_dim(v), work(1), lwork, iwork(1), ifail(1), info)

  SAFE_DEALLOCATE_A(work)

  if(info /= 0) then
    if(optional_default(bof, .true.)) then
      write(message(1),'(3a,i5)') 'In dlowest_geneigensolve, LAPACK ', TOSTRING(X(sygvx)), ' returned error message ', info
      call messages_fatal(1)
    else
      if(present(bof)) then
        bof = .true.
      end if
    end if
  else
    if(present(bof)) then
      bof = .false.
    end if
  end if
  if(present(err_code)) then
    err_code = info
  end if

  POP_SUB(dlowest_geneigensolve)
end subroutine dlowest_geneigensolve

#else

! ---------------------------------------------------------
!> Computes the k lowest eigenvalues and the eigenvectors of a complex
!! generalized Hermitian-definite eigenproblem, of the form  A*x=(lambda)*B*x.
!! Here A and B are assumed to be Hermitian and B is also positive definite.
subroutine zlowest_geneigensolve(k, n, a, b, e, v, bof, err_code)
  integer,           intent(in)    :: k, n
  CMPLX,             intent(inout) :: a(:,:)   !< (n,n)
  CMPLX,             intent(inout) :: b(:,:)   !< (n,n)
  FLOAT,             intent(out)   :: e(:)     !< (n)
  CMPLX,             intent(out)   :: v(:,:)   !< (n)
  logical, optional, intent(inout) :: bof      !< Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer            :: m, iwork(5*n), ifail(n), info, lwork
  FLOAT              :: abstol
  FLOAT              :: rwork(7*n)
  CMPLX, allocatable :: work(:)

  PUSH_SUB(zlowest_geneigensolve)

  ASSERT(n > 0)

  abstol = 2*sfmin()

  ! Work size query.
  SAFE_ALLOCATE(work(1:1))
  call X(hegvx)(1, 'V', 'I', 'U', n, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), lead_dim(v), work(1), -1, rwork(1), iwork(1), ifail(1), info)
  lwork = int(real(work(1)))
  SAFE_DEALLOCATE_A(work)

  SAFE_ALLOCATE(work(1:lwork))
  call X(hegvx)(1, 'V', 'I', 'U', n, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), lead_dim(v), work(1), lwork, rwork(1), iwork(1), ifail(1), info)
  SAFE_DEALLOCATE_A(work)

  if(info /= 0) then
    if(optional_default(bof, .true.)) then
      write(message(1),'(3a,i5)') 'In zlowest_geneigensolve, LAPACK ', TOSTRING(X(hegvx)), ' returned error message ', info
      call messages_fatal(1)
    else
      if(present(bof)) then
        bof = .true.
      end if
    end if
  else
    if(present(bof)) then
      bof = .false.
    end if
  end if
  if(present(err_code)) then
    err_code = info
  end if

  POP_SUB(zlowest_geneigensolve)
end subroutine zlowest_geneigensolve
#endif

#ifdef R_TREAL
! ---------------------------------------------------------
!> Computes all eigenvalues and eigenvectors of a real symmetric square matrix A.
subroutine deigensolve(n, a, e, bof, err_code)
  integer, intent(in)              :: n
  FLOAT,   intent(inout)           :: a(:,:)   !< (n,n)
  FLOAT,   intent(out)             :: e(:)     !< (n)
  logical, optional, intent(inout) :: bof      !< Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer            :: info, lwork
  FLOAT, allocatable :: work(:)

  PUSH_SUB(deigensolve)
  call profiling_in(eigensolver_prof, "DENSE_EIGENSOLVER")

  ASSERT(n > 0)

  lwork = 6*n
  SAFE_ALLOCATE(work(1:lwork))
  call lapack_syev('V', 'U', n, a(1, 1), lead_dim(a), e(1), work(1), lwork, info)
  SAFE_DEALLOCATE_A(work)

  if(info /= 0) then
    if(optional_default(bof, .true.)) then
      write(message(1),'(3a,i5)') 'In deigensolve, LAPACK ', TOSTRING(X(syev)), ' returned error message ', info
      call messages_fatal(1)
    else
      if(present(bof)) then
        bof = .true.
      end if
    end if
  else
    if(present(bof)) then
      bof = .false.
    end if
  end if
  if(present(err_code)) then
    err_code = info
  end if

  call profiling_out(eigensolver_prof)
  POP_SUB(deigensolve)
end subroutine deigensolve

#else

! ---------------------------------------------------------
!> Computes all eigenvalues and eigenvectors of a complex Hermitian square matrix A.
subroutine zeigensolve(n, a, e, bof, err_code)
  integer,           intent(in)    :: n
  CMPLX,             intent(inout) :: a(:,:)   !< (n,n)
  FLOAT,             intent(out)   :: e(:)     !< (n)
  logical, optional, intent(inout) :: bof      !< Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer            :: info, lwork
  CMPLX, allocatable :: work(:)
  FLOAT, allocatable :: rwork(:)

  PUSH_SUB(zeigensolve)
  call profiling_in(eigensolver_prof, "DENSE_EIGENSOLVER")

  ASSERT(n > 0)

  lwork = 6*n
  SAFE_ALLOCATE(work(1:lwork))
  SAFE_ALLOCATE(rwork(1:max(1, 3*n-2)))
  call lapack_heev('V','U', n, a(1, 1), lead_dim(a), e(1), work(1), lwork, rwork(1), info)
  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(rwork)

  if(info /= 0) then
    if(optional_default(bof, .true.)) then
      write(message(1),'(3a,i5)') 'In zeigensolve, LAPACK ', TOSTRING(X(heev)), ' returned error message ', info
      call messages_fatal(1)
    else
      if(present(bof)) then
        bof = .true.
      end if
    end if
  else
    if(present(bof)) then
      bof = .false.
    end if
  end if
  if(present(err_code)) then
    err_code = info
  end if

  call profiling_out(eigensolver_prof)
  POP_SUB(zeigensolve)
end subroutine zeigensolve
#endif

#ifdef R_TREAL
! ---------------------------------------------------------
!> Computes the k lowest eigenvalues and the eigenvectors of a real
!! standard symmetric-definite eigenproblem, of the form  A*x=(lambda)*x.
!! Here A is assumed to be symmetric.
subroutine dlowest_eigensolve(k, n, a, e, v)
  integer, intent(in)  :: k, n
  FLOAT,   intent(in)  :: a(:,:) !< (n, n)
  FLOAT,   intent(out) :: e(:)   !< (n)
  FLOAT,   intent(out) :: v(:,:) !< (n, k)

  integer            :: m, iwork(5*n), ifail(n), info, lwork
  FLOAT              :: abstol
  FLOAT, allocatable :: work(:)
  
  PUSH_SUB(dlowest_eigensolve)

  ASSERT(n > 0)

  abstol = 2*sfmin()

  ! Work size query.
  SAFE_ALLOCATE(work(1:1))
  call X(syevx)('V', 'I', 'U', n, a(1, 1), n, M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), n, work(1), -1, iwork(1), ifail(1), info) ! check info
  lwork = int(work(1))
  SAFE_DEALLOCATE_A(work)

  SAFE_ALLOCATE(work(1:lwork))
  call X(syevx)('V', 'I', 'U', n, a(1, 1), n, M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), n, work(1), lwork, iwork(1), ifail(1), info)
  SAFE_DEALLOCATE_A(work)

  if(info /= 0) then
    write(message(1),'(3a,i5)') &
      'In dlowest_eigensolve, LAPACK ', TOSTRING(X(syevx)), ' returned error message ', info
    call messages_fatal(1)
  end if

  POP_SUB(dlowest_eigensolve)
end subroutine dlowest_eigensolve

#else

! ---------------------------------------------------------
!> Computes the k lowest eigenvalues and the eigenvectors of a complex
!! standard Hermitian-definite eigenproblem, of the form  A*x=(lambda)*x.
!! Here A is assumed to be Hermitian.
subroutine zlowest_eigensolve(k, n, a, e, v)
  integer, intent(in)  :: k, n
  CMPLX,   intent(in)  :: a(:,:) !< (n, n)
  FLOAT,   intent(out) :: e(:)   !< (n)
  CMPLX,   intent(out) :: v(:,:) !< (n, k)

  integer            :: m, iwork(5*n), ifail(n), info, lwork
  FLOAT              :: abstol
  CMPLX, allocatable :: work(:)
   
  PUSH_SUB(zlowest_eigensolve)

  ASSERT(n > 0)

  abstol = 2*sfmin()

  ! Work size query.
  SAFE_ALLOCATE(work(1:1))
  call X(heevx)('V', 'I', 'U', n, a(1, 1), n, M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), n, work(1), -1, iwork(1), ifail(1), info)
  lwork = int(work(1))
  SAFE_DEALLOCATE_A(work)

  SAFE_ALLOCATE(work(1:lwork))
  call X(heevx)('V', 'I', 'U', n, a(1, 1), n, M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), n, work(1), lwork, iwork(1), ifail(1), info)
  SAFE_DEALLOCATE_A(work)

  if(info /= 0) then
    write(message(1),'(3a,i5)') &
      'In zlowest_eigensolve, LAPACK ', TOSTRING(X(heevx)), ' returned error message ', info
    call messages_fatal(1)
  end if

  POP_SUB(zlowest_eigensolve)
end subroutine zlowest_eigensolve
#endif

! ---------------------------------------------------------
!> Invert a real symmetric or complex Hermitian square matrix a
R_TYPE function X(determinant)(n, a, invert) result(d)
  integer, intent(in)           :: n
  R_TYPE,   intent(inout)       :: a(:,:) !< (n,n)
  logical, intent(in), optional :: invert

  interface
    subroutine X(getrf) (m, n, a, lda, ipiv, info)
      implicit none
      integer,      intent(in)    :: m, n, lda
      R_TYPE,        intent(inout) :: a         !< a(lda, n)
      integer,      intent(out)   :: ipiv       !< ipiv(min(m,n)
      integer,      intent(out)   :: info
    end subroutine X(getrf)

    subroutine X(getri) (n, a, lda, ipiv, work, lwork, info )
      implicit none
      integer,      intent(in)    :: n, lda, lwork
      R_TYPE,       intent(inout) :: a       !< a(lda, n)
      integer,      intent(in)    :: ipiv    !< ipiv(n)
      R_TYPE,       intent(inout) :: work    !< work(lwork)
      integer,      intent(out)   :: info
    end subroutine X(getri)
  end interface

  integer :: info, i
  integer, allocatable :: ipiv(:)
  R_TYPE, allocatable :: work(:)

  ! No PUSH_SUB, called too often

  ASSERT(n > 0)

  SAFE_ALLOCATE(work(1:n))
  SAFE_ALLOCATE(ipiv(1:n))

  call X(getrf)(n, n, a(1, 1), n, ipiv(1), info)
  if(info < 0) then
    write(message(1), '(5a, i5)') 'In ', TOSTRING(X(determinant)), ', LAPACK ', TOSTRING(X(getrf)), ' returned info = ', info
    call messages_fatal(1)
  end if

  d = M_ONE
  do i = 1, n
    if(ipiv(i) /= i) then
      d = - d*a(i, i)
    else
      d = d*a(i, i)
    end if
  end do

  if(optional_default(invert, .true.)) then
    call X(getri)(n, a(1, 1), n, ipiv(1), work(1), n, info)
    if(info /= 0) then
      write(message(1), '(5a, i5)') 'In ', TOSTRING(X(determinant)), ', LAPACK ', TOSTRING(X(getri)), ' returned info = ', info
      call messages_fatal(1)
    end if
  end if

  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(ipiv)

end function X(determinant)


! ---------------------------------------------------------
!> Invert a real/complex symmetric square matrix a
subroutine X(sym_inverter)(uplo, n, a)
  character(1), intent(in)      :: uplo
  integer, intent(in)           :: n
  R_TYPE,  intent(inout)        :: a(:,:) !< (n,n)

  interface
    subroutine X(sytrf) (uplo, n, a, lda, ipiv, work, lwork, info)
      implicit none
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda, lwork
      R_TYPE,       intent(inout) :: a
      integer,      intent(out)   :: ipiv
      R_TYPE,       intent(inout) :: work
      integer,      intent(out)   :: info
    end subroutine X(sytrf)

    subroutine X(sytri) (uplo, n, a, lda, ipiv, work, info)
      implicit none
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda
      R_TYPE,       intent(inout) :: a
      integer,      intent(in)    :: ipiv
      R_TYPE,       intent(inout) :: work
      integer,      intent(out)   :: info
    end subroutine X(sytri)
  end interface

  integer :: info
  integer, allocatable :: ipiv(:)
  R_TYPE, allocatable :: work(:)

  PUSH_SUB(X(sym_inverter))

  ASSERT(n > 0)

  SAFE_ALLOCATE(work(1:n))
  SAFE_ALLOCATE(ipiv(1:n))

  call X(sytrf)(uplo, n, a(1, 1), lead_dim(a), ipiv(1), work(1), n, info)
  if(info < 0) then
    write(message(1), '(5a, i5)') 'In ', TOSTRING(X(sym_inverter)), ', LAPACK ', TOSTRING(X(sytrf)), ' returned info = ', info
    call messages_fatal(1)
  end if

  call X(sytri)(uplo, n, a(1, 1), lead_dim(a), ipiv(1), work(1), info)
  if(info /= 0) then
    write(message(1), '(5a, i5)') 'In ', TOSTRING(X(sym_inverter)), ', LAPACK ', TOSTRING(X(sytri)), ' returned info = ', info
    call messages_fatal(1)
  end if

  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(ipiv)
  POP_SUB(X(sym_inverter))
end subroutine X(sym_inverter)

#ifdef R_TREAL
! ---------------------------------------------------------
!> compute the solution to a real system of linear equations A*X = B,
!!  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
subroutine dlinsyssolve(n, nrhs, a, b, x)
  integer, intent(in)    :: n, nrhs
  FLOAT,   intent(inout) :: a(:,:) !< (n, n)
  FLOAT,   intent(inout) :: b(:,:) !< (n, nrhs)
  FLOAT,   intent(out)   :: x(:,:) !< (n, nrhs)

  integer :: info
  integer, allocatable :: ipiv(:), iwork(:)
  FLOAT :: rcond
  FLOAT, allocatable :: ferr(:), berr(:), work(:), r(:), c(:), af(:,:)
  character(1) :: equed

  ! no PUSH_SUB, called too often

  ASSERT(n > 0)

  SAFE_ALLOCATE(ipiv(1:n))
  SAFE_ALLOCATE(iwork(1:n))
  SAFE_ALLOCATE(ferr(1:nrhs))
  SAFE_ALLOCATE(berr(1:nrhs))
  SAFE_ALLOCATE(work(1:4*n))
  SAFE_ALLOCATE(r(1:n))
  SAFE_ALLOCATE(c(1:n))
  SAFE_ALLOCATE(af(1:n, 1:n))

  call X(gesvx) ("N", "N", n, nrhs, a(1, 1), n, af(1, 1), n, ipiv(1), equed, r(1), c(1), b(1, 1), n, x(1, 1), n, &
    rcond, ferr(1), berr(1), work(1), iwork(1), info)

  if(info /= 0) then
    write(message(1), '(3a, i5)') 'In dlinsyssolve, LAPACK ', TOSTRING(X(gesvx)), ' returned info = ', info
    if(info == n+1) then
      message(2) = '(reciprocal of the condition number is less than machine precision)'
      call messages_warning(2)
    else
      call messages_fatal(1)
    end if
  end if

  SAFE_DEALLOCATE_A(ipiv)
  SAFE_DEALLOCATE_A(iwork)
  SAFE_DEALLOCATE_A(ferr)
  SAFE_DEALLOCATE_A(berr)
  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(r)
  SAFE_DEALLOCATE_A(c)
  SAFE_DEALLOCATE_A(af)

end subroutine dlinsyssolve

#else

! ---------------------------------------------------------
!> compute the solution to a complex system of linear equations A*X = B,
!!  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
subroutine zlinsyssolve(n, nrhs, a, b, x)
  integer, intent(in)    :: n, nrhs
  CMPLX,   intent(inout) :: a(:,:) !< (n, n)
  CMPLX,   intent(inout) :: b(:,:) !< (n, nrhs)
  CMPLX,   intent(out)   :: x(:,:) !< (n, nrhs)

  integer              :: info
  integer, allocatable :: ipiv(:)
  FLOAT,   allocatable :: rwork(:), ferr(:), berr(:), r(:), c(:)
  FLOAT                :: rcond
  CMPLX, allocatable   :: work(:), af(:,:)
  character(1)         :: equed

  PUSH_SUB(zlinsyssolve)

  ASSERT(n > 0)

  SAFE_ALLOCATE(ipiv(1:n))
  SAFE_ALLOCATE(rwork(1:2*n))
  SAFE_ALLOCATE(ferr(1:nrhs))
  SAFE_ALLOCATE(berr(1:nrhs))
  SAFE_ALLOCATE(work(1:2*n))
  SAFE_ALLOCATE(r(1:n))
  SAFE_ALLOCATE(c(1:n))
  SAFE_ALLOCATE(af(1:n, 1:n))

  equed = 'N'

  call X(gesvx) ("N", "N", n, nrhs, a(1, 1), n, af(1, 1), n, ipiv(1), equed, r(1), c(1), b(1, 1), n, x(1, 1), n, &
    rcond, ferr(1), berr(1), work(1), rwork(1), info)
  if(info /= 0) then
    write(message(1), '(3a, i5)') 'In zlinsyssolve, LAPACK ', TOSTRING(X(gesvx)), ' returned info = ', info
    call messages_fatal(1)
  end if

  SAFE_DEALLOCATE_A(ipiv)
  SAFE_DEALLOCATE_A(rwork)
  SAFE_DEALLOCATE_A(ferr)
  SAFE_DEALLOCATE_A(berr)
  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(r)
  SAFE_DEALLOCATE_A(c)
  SAFE_DEALLOCATE_A(af)

  POP_SUB(zlinsyssolve)
end subroutine zlinsyssolve
#endif

#ifdef R_TCOMPLEX
! ---------------------------------------------------------
!> computes the singular value decomposition of a complex NxN matrix a(:,:)
subroutine zsingular_value_decomp(n, a, u, vt, sg_values)
  integer, intent(in)    :: n
  CMPLX,   intent(inout) :: a(:,:)          !< (n,n)
  CMPLX,   intent(out)   :: u(:,:), vt(:,:) !< (n,n)  
  FLOAT,   intent(out)   :: sg_values(:)    !< (n)

  interface
    subroutine X(gesvd) ( jobu, jobvt, m, n, a, lda, s, u, ldu, &
      vt, ldvt, work, lwork, rwork, info )
      implicit none
      character(1), intent(in)    :: jobu, jobvt
      integer,      intent(in)    :: m, n
      CMPLX,        intent(inout) :: a, u, vt ! a(lda,n), b(ldu,n), b(ldvt,n)
      CMPLX,        intent(out)   :: work     ! work(lwork)
      integer,      intent(in)    :: lda, ldu, ldvt, lwork
      integer,      intent(out)   :: info
      FLOAT,        intent(out)   :: s        ! s(min(m,n))
      FLOAT,        intent(inout) :: rwork    ! rwork(5*min(m,n))
    end subroutine X(gesvd)
  end interface

  integer :: m, info, lwork
  CMPLX, allocatable :: work(:)
  FLOAT, allocatable :: rwork(:)

  PUSH_SUB(zsingular_value_decomp)

  ASSERT(n > 0)

  ! for now we treat only square matrices
  m = n 

  ! double minimum lwork size to increase performance (see manpage)
  lwork = 2*( 2*min(m, n) + max(m, n) )

  SAFE_ALLOCATE(work(1:lwork))
  SAFE_ALLOCATE(rwork(1:5*min(m, n)))

  call X(gesvd)( &
    'A', 'A', m, n, a(1, 1), m, sg_values(1), u(1, 1), m, vt(1, 1), n, work(1), lwork, rwork(1), info )

  if(info /= 0) then
    write(message(1), '(3a, i5)') 'In zsingular_value_decomp, LAPACK ', TOSTRING(X(gesvd)), ' returned info = ', info
    call messages_fatal(1)
  end if

  SAFE_DEALLOCATE_A(rwork)
  SAFE_DEALLOCATE_A(work)
  POP_SUB(zsingular_value_decomp)
end subroutine zsingular_value_decomp


! ---------------------------------------------------------
!> computes inverse of a complex NxN matrix a(:,:) using the SVD decomposition 
subroutine zsvd_inverse(n, a, threshold)
  integer, intent(in)           :: n
  CMPLX,   intent(inout)        :: a(:,:)    !< (n,n); a will be replaced by its inverse
  FLOAT,   intent(in), optional :: threshold

  CMPLX, allocatable :: u(:,:), vt(:,:)
  FLOAT, allocatable :: sg_values(:)
  CMPLX   :: tmp
  FLOAT   :: sg_inverse, threshold_
  integer :: j, k, l

  ASSERT(n > 0)

  SAFE_ALLOCATE( u(1:n, 1:n))
  SAFE_ALLOCATE(vt(1:n, 1:n))
  SAFE_ALLOCATE(sg_values(1:n))

  PUSH_SUB(zsvd_inverse)

  call zsingular_value_decomp(n, a, u, vt, sg_values)

  threshold_ = CNST(1e-10)
  if(present(threshold)) threshold_ = threshold

  ! build inverse
  do j = 1, n
    do k = 1, n
      tmp = M_ZERO
      do l = 1, n
        if (sg_values(l) < threshold_) then
          write(message(1), '(a)') 'In zsvd_inverse: singular value below threshold.'
          call messages_warning(1)
          sg_inverse = M_ZERO
        else
          sg_inverse = M_ONE/sg_values(l)
        end if
        tmp = tmp + conjg(vt(l, k))*sg_inverse*conjg(u(j, l))
      end do
      a(j, k) = tmp
    end do
  end do

  SAFE_DEALLOCATE_A(sg_values)
  SAFE_DEALLOCATE_A(vt)
  SAFE_DEALLOCATE_A(u)
  POP_SUB(zsvd_inverse)
end subroutine zsvd_inverse
#endif

! ---------------------------------------------------------
!> Calculate the inverse of a real/complex upper triangular matrix (in
!! unpacked storage). (lower triangular would be a trivial variant of this)
subroutine X(invert_upper_triangular)(n, a)
  integer, intent(in)    :: n
  R_TYPE,  intent(inout) :: a(:,:) !< (n,n)

  integer :: info

  interface
    subroutine X(trtri)(uplo, diag, n, a, lda, info)
      implicit none
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: diag
      integer,      intent(in)    :: n
      R_TYPE,       intent(inout) :: a
      integer,      intent(in)    :: lda
      integer,      intent(out)   :: info
    end subroutine X(trtri)
  end interface
  
  PUSH_SUB(X(invert_upper_triangular))

  ASSERT(n > 0)

  call X(trtri)('U', 'N', n, a(1, 1), n, info)

  if(info /= 0) then
    write(message(1), '(5a,i5)') &
      'In ', TOSTRING(Xinvert_upper_triangular), ' LAPACK ', TOSTRING(X(trtri)), ' returned error message ', info
    call messages_fatal(1)
  end if

  POP_SUB(X(invert_upper_triangular))
end subroutine X(invert_upper_triangular)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

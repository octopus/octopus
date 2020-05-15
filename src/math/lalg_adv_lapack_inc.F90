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
      write(message(1), '(5a,i5)') 'In ', TOSTRING(X(cholesky)), ', LAPACK ', TOSTRING(X(potrf)), ' returned error message ', info
! http://www.netlib.org/lapack/explore-3.1.1-html/dpotrf.f.html and zpotrf.f.html
!      *  INFO    (output) INTEGER
!      *          = 0:  successful exit
!      *          < 0:  if INFO = -i, the i-th argument had an illegal value
!      *          > 0:  if INFO = i, the leading minor of order i is not
!      *                positive definite, and the factorization could not be
!      *                completed.
      if(info < 0) then
        write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
      else
        write(message(2), '(a,i5,a)') 'The leading minor of order ', info, ' is not positive definite.'
      end if
      call messages_fatal(2)
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
subroutine X(geneigensolve)(n, a, b, e, preserve_mat, bof, err_code)
  integer,           intent(in)    :: n
  R_TYPE,            intent(inout) :: a(:,:)   !< (n,n)
  R_TYPE,            intent(inout) :: b(:,:)   !< (n,n)
  FLOAT,             intent(out)   :: e(:)     !< (n)
  logical,           intent(in)    :: preserve_mat !< If true, the matrix a and b on exit are the same
                                                   !< as on input
  logical, optional, intent(inout) :: bof      !< Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer :: info, lwork, ii, jj
#ifdef R_TCOMPLEX
  FLOAT, allocatable :: rwork(:)
#endif
  R_TYPE, allocatable :: work(:), diag(:)

  call profiling_in(eigensolver_prof, "DENSE_EIGENSOLVER")
  PUSH_SUB(X(geneigensolve))

  ASSERT(n > 0)

  if(preserve_mat) then
    SAFE_ALLOCATE(diag(1:n))
    ! store the diagonal of b  
    do ii = 1, n
      diag(ii) = b(ii, ii)
    end do
  end if

  lwork = 5*n ! get this from workspace query
  SAFE_ALLOCATE(work(1:lwork))
#ifdef R_TCOMPLEX
  SAFE_ALLOCATE(rwork(1:max(1, 3*n-2)))
  call lapack_hegv(1, 'V', 'U', n, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), e(1), work(1), lwork, rwork(1), info)
  SAFE_DEALLOCATE_A(rwork)
#else
  call lapack_sygv(1, 'V', 'U', n, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), e(1), work(1), lwork, info)
#endif
  SAFE_DEALLOCATE_A(work)

  if(preserve_mat) then
    ! b was destroyed, so we rebuild it
    do ii = 1, n
      do jj = 1, ii - 1
        b(jj, ii) = R_CONJ(b(ii, jj))
      end do
      b(ii, ii) = diag(ii)
    end do

    SAFE_DEALLOCATE_A(diag)
  end if

  if(info /= 0) then
    if(optional_default(bof, .true.)) then
      write(message(1),'(3a)') 'In ', TOSTRING(X(geneigensolve)), ', LAPACK '
#ifdef R_TCOMPLEX
      write(message(1),'(3a,i5)') trim(message(1)), TOSTRING(X(hegv)), ' returned error message ', info
#else
      write(message(1),'(3a,i5)') trim(message(1)), TOSTRING(X(sygv)), ' returned error message ', info
#endif
!*  INFO    (output) INTEGER
!*          = 0:  successful exit
!*          < 0:  if INFO = -i, the i-th argument had an illegal value
!*          > 0:  ZPOTRF or ZHEEV returned an error code:
!*             <= N:  if INFO = i, ZHEEV failed to converge;
!*                    i off-diagonal elements of an intermediate
!*                    tridiagonal form did not converge to zero;
!*             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
!*                    minor of order i of B is not positive definite.
!*                    The factorization of B could not be completed and
!*                    no eigenvalues or eigenvectors were computed.      
      if(info < 0) then
        write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
      else if(info <= n) then
        write(message(2), '(i5,a)') info, ' off-diagonal elements of an intermediate tridiagonal did not converge to zero.'
      else
        write(message(2), '(a,i5,a)') 'The leading minor of order ', info - n, ' of B is not positive definite.'
      end if
      call messages_fatal(2)
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
  ! Initializing info, if not it can cause that the geev query mode fails.
  ! Besides, if info is not initialized valgrind complains about it. 
  info = 0
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
  if(info /= 0) then
    write(message(1),'(5a,i5)') 'In ', TOSTRING(X(eigensolve_nonh)), &
      ', LAPACK ', TOSTRING(X(geev)), ' workspace query returned error message ', info
    call messages_fatal(1)
  end if

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
      work, lwork, rwork, info)
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
      ', LAPACK ', TOSTRING(X(geev)), ' returned error message ', info
!*  INFO    (output) INTEGER
!*          = 0:  successful exit
!*          < 0:  if INFO = -i, the i-th argument had an illegal value.
!*          > 0:  if INFO = i, the QR algorithm failed to compute all the
!*                eigenvalues, and no eigenvectors have been computed;
!*                elements i+1:N of WR and WI contain eigenvalues which
!*                have converged.
    if(info < 0) then
      write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
    else
      write(message(2), '(a,i5,a,i5,a)') 'Only eigenvalues ', info + 1, ' to ', n, ' could be computed.'
    end if
    call messages_fatal(2)
  end if
  if(present(err_code)) then
    err_code = info
  end if

  if(optional_default(sort_eigenvectors, .false.)) then
    SAFE_ALLOCATE(re(1:n))
    SAFE_ALLOCATE(ind(1:n))
    SAFE_ALLOCATE(e_copy(1:n))
    SAFE_ALLOCATE(a_copy(1:n, 1:n))
    re = TOFLOAT(e)
    e_copy = e
    a_copy = a
    call sort(re, ind)
    do ii = 1, n
      e(ii) = e_copy(ind(ii))
      a(1:n, ii) = a_copy(1:n, ind(ii))
    end do
    SAFE_DEALLOCATE_A(e_copy)
    SAFE_DEALLOCATE_A(a_copy)
    SAFE_DEALLOCATE_A(re)
    SAFE_DEALLOCATE_A(ind)
  end if

  POP_SUB(X(eigensolve_nonh))
end subroutine X(eigensolve_nonh)


! ---------------------------------------------------------
!> Computes the k lowest eigenvalues and the eigenvectors of a real symmetric or complex Hermitian
!! generalized definite eigenproblem, of the form  A*x=(lambda)*B*x. B is also positive definite.
subroutine X(lowest_geneigensolve)(k, n, a, b, e, v, preserve_mat, bof, err_code)
  integer,           intent(in)    :: k, n
  R_TYPE,            intent(inout) :: a(:,:)   !< (n, n)
  R_TYPE,            intent(inout) :: b(:,:)   !< (n, n)
  FLOAT,             intent(out)   :: e(:)     !< (n)
  R_TYPE,            intent(out)   :: v(:,:)   !< (n, n)
  logical,           intent(in)    :: preserve_mat !< If true, the matrix a and b on exit are the same
                                                   !< as on input
  logical, optional, intent(inout) :: bof      !< Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer            :: m, iwork(5*n), ifail(n), info, lwork, ii, jj ! allocate me
  FLOAT              :: abstol
  R_TYPE, allocatable :: work(:), diaga(:), diagb(:)
#ifndef R_TREAL
  FLOAT              :: rwork(7*n)
#endif
  PUSH_SUB(X(lowest_geneigensolve))

  ASSERT(n > 0)

  abstol = 2*sfmin()

  if(preserve_mat) then
    SAFE_ALLOCATE(diaga(1:n))
    SAFE_ALLOCATE(diagb(1:n))
  
    ! store the diagonal of a and b  
    do ii = 1, n
      diaga(ii) = a(ii, ii)
      diagb(ii) = b(ii, ii)
    end do
  end if


  ! Work size query.
  SAFE_ALLOCATE(work(1:1))

#ifdef R_TREAL
  call dsygvx(1, 'V', 'I', 'U', n, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), lead_dim(v), work(1), -1, iwork(1), ifail(1), info)
  if(info /= 0) then
    write(message(1),'(3a,i5)') 'In dlowest_geneigensolve, LAPACK ', &
      TOSTRING(dsygvx), ' workspace query returned error message ', info
    call messages_fatal(1)
  end if  
  lwork = int(work(1))
  SAFE_DEALLOCATE_A(work)

  SAFE_ALLOCATE(work(1:lwork))

  call dsygvx(1, 'V', 'I', 'U', n, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), lead_dim(v), work(1), lwork, iwork(1), ifail(1), info)

#else
  call zhegvx(1, 'V', 'I', 'U', n, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), M_ZERO, M_ZERO, &
      1, k, abstol, m, e(1), v(1, 1), lead_dim(v), work(1), -1, rwork(1), iwork(1), ifail(1), info)
    if(info /= 0) then
      write(message(1),'(3a,i5)') 'In zlowest_geneigensolve, LAPACK ', &
        TOSTRING(zhegvx), ' workspace query returned error message ', info
      call messages_fatal(1)
    end if
    lwork = int(real(work(1)))
  SAFE_DEALLOCATE_A(work)
  
  SAFE_ALLOCATE(work(1:lwork))
  call zhegvx(1, 'V', 'I', 'U', n, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), M_ZERO, M_ZERO, &
      1, k, abstol, m, e(1), v(1, 1), lead_dim(v), work(1), lwork, rwork(1), iwork(1), ifail(1), info)
#endif

  if(preserve_mat) then
    ! b was destroyed, so we rebuild it
    do ii = 1, n
      do jj = 1, ii - 1
        a(jj, ii) = R_CONJ(a(ii, jj))
        b(jj, ii) = R_CONJ(b(ii, jj))
      end do
      a(ii, ii) = diaga(ii)
      b(ii, ii) = diagb(ii)
    end do
  
    SAFE_DEALLOCATE_A(diaga)
    SAFE_DEALLOCATE_A(diagb)
  end if


  SAFE_DEALLOCATE_A(work)

  if(info /= 0) then
    if(optional_default(bof, .true.)) then
#ifdef R_TREAL
      write(message(1),'(3a,i5)') 'In dlowest_geneigensolve, LAPACK ', &
        TOSTRING(dsygvx), ' returned error message ', info
#else
      write(message(1),'(3a,i5)') 'In zlowest_geneigensolve, LAPACK ', &
        TOSTRING(zhegvx), ' returned error message ', info
#endif
!        INFO    (output) INTEGER
!*          = 0:  successful exit
!*          < 0:  if INFO = -i, the i-th argument had an illegal value
!*          > 0:  DPOTRF or DSYEVX returned an error code:
!*             <= N:  if INFO = i, DSYEVX failed to converge;
!*                    i eigenvectors failed to converge.  Their indices
!*                    are stored in array IFAIL.
!*             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
!*                    minor of order i of B is not positive definite.
!*                    The factorization of B could not be completed and
!*                    no eigenvalues or eigenvectors were computed.
      if(info < 0) then
        write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
      else if(info <= n) then
        write(message(2), *) info, ' eigenvectors failed to converge: ', ifail(1:info)
      else
        write(message(2), '(a,i5,a)') 'The leading minor of order ', info - n, ' of B is not positive definite.'
      end if
      call messages_fatal(2)
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

  POP_SUB(X(lowest_geneigensolve))
end subroutine X(lowest_geneigensolve)

! ---------------------------------------------------------
!> Computes all eigenvalues and eigenvectors of a real symmetric  or hermitian square matrix A.
subroutine X(eigensolve)(n, a, e, bof, err_code)
  integer, intent(in)              :: n
  R_TYPE,  intent(inout)           :: a(:,:)   !< (n,n)
  FLOAT,   intent(out)             :: e(:)     !< (n)
  logical, optional, intent(inout) :: bof      !< Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer             :: info, lwork
  R_TYPE, allocatable :: work(:)
#ifndef R_TREAL
  FLOAT, allocatable :: rwork(:)
#endif

  PUSH_SUB(X(eigensolve))
  call profiling_in(eigensolver_prof, "DENSE_EIGENSOLVER")

  ASSERT(n > 0)

  lwork = 6*n ! query?
  SAFE_ALLOCATE(work(1:lwork))
#ifdef R_TREAL
  call lapack_syev('V', 'U', n, a(1, 1), lead_dim(a), e(1), work(1), lwork, info)
#else
  SAFE_ALLOCATE(rwork(1:max(1, 3*n-2)))
  call lapack_heev('V','U', n, a(1, 1), lead_dim(a), e(1), work(1), lwork, rwork(1), info)
  SAFE_DEALLOCATE_A(rwork)
#endif
  SAFE_DEALLOCATE_A(work)

  if(info /= 0) then
    if(optional_default(bof, .true.)) then
#ifdef R_TREAL
      write(message(1),'(3a,i5)') 'In eigensolve, LAPACK ', TOSTRING(dsyev), ' returned error message ', info
#else
      write(message(1),'(3a,i5)') 'In eigensolve, LAPACK ', TOSTRING(zheev), ' returned error message   ', info
#endif
!*  INFO    (output) INTEGER
!*          = 0:  successful exit
!*          < 0:  if INFO = -i, the i-th argument had an illegal value
!*          > 0:  if INFO = i, the algorithm failed to converge; i
!*                off-diagonal elements of an intermediate tridiagonal
!*                form did not converge to zero.
      if(info < 0) then
        write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
      else
        write(message(2), '(i5,a)') info, ' off-diagonal elements of an intermediate tridiagonal did not converge to zero.'
      end if
      call messages_fatal(2)
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
  POP_SUB(X(eigensolve))
end subroutine X(eigensolve)


! ---------------------------------------------------------
!> Computes the k lowest eigenvalues and the eigenvectors of a 
!! standard symmetric-definite eigenproblem, of the form  A*x=(lambda)*x.
!! Here A is assumed to be symmetric.
subroutine X(lowest_eigensolve)(k, n, a, e, v, preserve_mat)
  integer, intent(in)    :: k, n
  R_TYPE,  intent(inout) :: a(:,:) !< (n, n)
  FLOAT,   intent(out)   :: e(:)   !< (n)
  R_TYPE,  intent(out)   :: v(:,:) !< (n, k)
  logical, intent(in)    :: preserve_mat !< If true, the matrix a and b on exit are the same
                                         !< as on input

  integer             :: m, iwork(5*n), ifail(n), info, lwork, ii, jj
  FLOAT               :: abstol
  R_TYPE, allocatable :: work(:), diaga(:)
  
  PUSH_SUB(X(lowest_eigensolve))

  ASSERT(n > 0)

  abstol = 2*sfmin()

  if(preserve_mat) then
    SAFE_ALLOCATE(diaga(1:n))

    ! store the diagonal of a and b  
    do ii = 1, n
      diaga(ii) = a(ii, ii)
    end do
  end if


  ! Work size query.
  SAFE_ALLOCATE(work(1:1))
#ifdef R_TREAL
  call dsyevx('V', 'I', 'U', n, a(1, 1), n, M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), n, work(1), -1, iwork(1), ifail(1), info)
  if(info /= 0) then
    write(message(1),'(3a,i5)') 'In dlowest_eigensolve, LAPACK ', &
      TOSTRING(dsyevx), ' workspace query returned error message ', info
    call messages_fatal(1)
  end if
  lwork = int(work(1))
  SAFE_DEALLOCATE_A(work)

  SAFE_ALLOCATE(work(1:lwork))
  call dsyevx('V', 'I', 'U', n, a(1, 1), n, M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), n, work(1), lwork, iwork(1), ifail(1), info)
#else
  call zheevx('V', 'I', 'U', n, a(1, 1), n, M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), n, work(1), -1, iwork(1), ifail(1), info)
  if(info /= 0) then
    write(message(1),'(3a,i5)') 'In zlowest_eigensolve, LAPACK ', &
      TOSTRING(zheevx), ' workspace query returned error message ', info
    call messages_fatal(1)
  end if
  lwork = int(work(1))
  SAFE_DEALLOCATE_A(work)
  
  SAFE_ALLOCATE(work(1:lwork))
  call zheevx('V', 'I', 'U', n, a(1, 1), n, M_ZERO, M_ZERO, &
      1, k, abstol, m, e(1), v(1, 1), n, work(1), lwork, iwork(1), ifail(1), info)
 
#endif

  if(preserve_mat) then
    ! b was destroyed, so we rebuild it
    do ii = 1, n
      do jj = 1, ii - 1
        a(jj, ii) = R_CONJ(a(ii, jj))
      end do
      a(ii, ii) = diaga(ii)
    end do
  
    SAFE_DEALLOCATE_A(diaga)
  end if


  SAFE_DEALLOCATE_A(work)

  if(info /= 0) then
#ifdef R_TREAL
    write(message(1),'(3a,i5)') &
      'In dlowest_eigensolve, LAPACK ', TOSTRING(dsyevx), ' returned error message ', info
#else
    write(message(1),'(3a,i5)') &
      'In zlowest_eigensolve, LAPACK ', TOSTRING(zheevx), ' returned error message ', info
#endif
!    http://www.netlib.org/lapack/explore-3.1.1-html/dsyevx.f.html
!*  INFO    (output) INTEGER
!*          = 0:  successful exit
!*          < 0:  if INFO = -i, the i-th argument had an illegal value
!*          > 0:  if INFO = i, then i eigenvectors failed to converge.
!*                Their indices are stored in array IFAIL.
    if(info < 0) then
      write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
    else
      write(message(2), *) info, ' eigenvectors failed to converge: ', ifail(1:info)
    end if
    call messages_fatal(2)
  end if

  POP_SUB(X(lowest_eigensolve))
end subroutine X(lowest_eigensolve)

! ---------------------------------------------------------
!> Invert a real symmetric or complex Hermitian square matrix a
R_TYPE function X(determinant)(n, a, preserve_mat) result(d)
  integer, intent(in)           :: n
  R_TYPE, target, intent(inout) :: a(:,:) !< (n,n)
  logical, intent(in)           :: preserve_mat

  integer :: info, i
  integer, allocatable :: ipiv(:)
  R_TYPE, pointer :: tmp_a(:,:)

  ! No PUSH_SUB, called too often

  ASSERT(n > 0)

  SAFE_ALLOCATE(ipiv(1:n))

  if(preserve_mat) then
    SAFE_ALLOCATE(tmp_a(1:n, 1:n))
    tmp_a(1:n, 1:n) = a(1:n, 1:n)
  else
    tmp_a => a
  end if

  call lapack_getrf(n, n, tmp_a(1, 1), lead_dim(tmp_a), ipiv(1), info)
  if(info < 0) then
    write(message(1), '(5a, i5)') 'In ', TOSTRING(X(determinant)), ', LAPACK ', TOSTRING(X(getrf)), ' returned info = ', info
    call messages_fatal(1)
  end if

  d = M_ONE
  do i = 1, n
    if(ipiv(i) /= i) then
      d = - d*tmp_a(i, i)
    else
      d = d*tmp_a(i, i)
    end if
  end do

  SAFE_DEALLOCATE_A(ipiv)
  if(preserve_mat) then
    SAFE_DEALLOCATE_P(tmp_a)
  end if

end function X(determinant)

! ---------------------------------------------------------
!> Invert a real symmetric or complex Hermitian square matrix a
subroutine X(inverter)(n, a, det)
  integer,           intent(in)     :: n
  R_TYPE,            intent(inout)  :: a(:,:) !< (n,n)
  R_TYPE,  optional, intent(out)    :: det

  integer :: info, i
  integer, allocatable :: ipiv(:)
  R_TYPE, allocatable :: work(:)

  ! No PUSH_SUB, called too often

  ASSERT(n > 0)

  SAFE_ALLOCATE(work(1:n)) ! query?
  SAFE_ALLOCATE(ipiv(1:n))

  call lapack_getrf(n, n, a(1, 1), lead_dim(a), ipiv(1), info)
  if(info < 0) then
    write(message(1), '(5a, i5)') 'In ', TOSTRING(X(determinant)), ', LAPACK ', TOSTRING(X(getrf)), ' returned info = ', info
    call messages_fatal(1)
  end if

  if(present(det)) then
    det = M_ONE
    do i = 1, n
      if(ipiv(i) /= i) then
        det = - det*a(i, i)
      else
        det = det*a(i, i)
      end if
    end do
  end if


  call lapack_getri(n, a(1, 1), lead_dim(a), ipiv(1), work(1), n, info)
  if(info /= 0) then
    write(message(1), '(5a, i5)') 'In ', TOSTRING(X(determinant)), ', LAPACK ', TOSTRING(X(getri)), ' returned info = ', info
!    http://www.netlib.org/lapack/explore-3.1.1-html/zgetri.f.html
!*  INFO    (output) INTEGER
!*          = 0:  successful exit
!*          < 0:  if INFO = -i, the i-th argument had an illegal value
!*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!*                singular and its inverse could not be computed.
  if(info < 0) then
    write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
  else
    write(message(2), '(a,i5,a)') 'Diagonal element ', info, ' of U is 0; matrix is singular.'
  end if
  call messages_fatal(2)
  end if

  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(ipiv)

end subroutine X(inverter)


! ---------------------------------------------------------
!> Invert a real/complex symmetric square matrix a
subroutine X(sym_inverter)(uplo, n, a)
  character(1), intent(in)      :: uplo
  integer, intent(in)           :: n
  R_TYPE,  intent(inout)        :: a(:,:) !< (n,n)

  integer :: info
  integer, allocatable :: ipiv(:)
  R_TYPE, allocatable :: work(:)

  PUSH_SUB(X(sym_inverter))

  ASSERT(n > 0)

  SAFE_ALLOCATE(work(1:n)) ! query?
  SAFE_ALLOCATE(ipiv(1:n))

  call lapack_sytrf(uplo, n, a(1, 1), lead_dim(a), ipiv(1), work(1), n, info)
  if(info < 0) then
    write(message(1), '(5a, i5)') 'In ', TOSTRING(X(sym_inverter)), ', LAPACK ', TOSTRING(X(sytrf)), ' returned info = ', info
    call messages_fatal(1)
  end if

  call lapack_sytri(uplo, n, a(1, 1), lead_dim(a), ipiv(1), work(1), info)
  if(info /= 0) then
    write(message(1), '(5a, i5)') 'In ', TOSTRING(X(sym_inverter)), ', LAPACK ', TOSTRING(X(sytri)), ' returned info = ', info
!    http://www.netlib.org/lapack/explore-3.1.1-html/dsytri.f.html
!*  INFO    (output) INTEGER
!*          = 0:  successful exit
!*          < 0:  if INFO = -i, the i-th argument had an illegal value
!*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!*                singular and its inverse could not be computed.
    if(info < 0) then
      write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
    else
      write(message(2), '(a,i5,a)') 'Diagonal element ', info, ' of D is 0; matrix is singular.'
    end if
    call messages_fatal(2)
  end if

  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(ipiv)
  POP_SUB(X(sym_inverter))
end subroutine X(sym_inverter)

! MJV 9 nov 2016: why is this stuff explicitly set in cpp instead of using the
! macros X()??? For the moment I have replicated this strategy in svd below.
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
  SAFE_ALLOCATE(iwork(1:n)) ! query?
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
!*  INFO    (output) INTEGER
!*          = 0:  successful exit
!*          < 0:  if INFO = -i, the i-th argument had an illegal value
!*          > 0:  if INFO = i, and i is
!*                <= N:  U(i,i) is exactly zero.  The factorization has
!*                       been completed, but the factor U is exactly
!*                       singular, so the solution and error bounds
!*                       could not be computed. RCOND = 0 is returned.
!*                = N+1: U is nonsingular, but RCOND is less than machine
!*                       precision, meaning that the matrix is singular
!*                       to working precision.  Nevertheless, the
!*                       solution and error bounds are computed because
!*                       there are a number of situations where the
!*                       computed solution can be more accurate than the
!*                       value of RCOND would suggest.
    if(info < 0) then
      write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
      call messages_fatal(2)      
    else if(info == n+1) then
      message(2) = '(reciprocal of the condition number is less than machine precision)'
      call messages_warning(2)
    else
      write(message(2), '(a,i5,a)') 'Diagonal element ', info, ' of U is 0; matrix is singular.'
      call messages_fatal(2)
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

#elif R_TCOMPLEX

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

  SAFE_ALLOCATE(ipiv(1:n)) ! query?
  SAFE_ALLOCATE(rwork(1:2*n))
  SAFE_ALLOCATE(ferr(1:nrhs))
  SAFE_ALLOCATE(berr(1:nrhs))
  SAFE_ALLOCATE(work(1:2*n))
  SAFE_ALLOCATE(r(1:n))
  SAFE_ALLOCATE(c(1:n))
  SAFE_ALLOCATE(af(1:n, 1:n))

  equed = 'N'

  call X(gesvx) ("N", "N", n, nrhs, a(1, 1), lead_dim(a), af(1, 1), lead_dim(af), &
                  ipiv(1), equed, r(1), c(1), b(1, 1), lead_dim(b), x(1, 1), lead_dim(x), &
                  rcond, ferr(1), berr(1), work(1), rwork(1), info)
  if(info /= 0) then
    write(message(1), '(3a, i5)') 'In zlinsyssolve, LAPACK ', TOSTRING(X(gesvx)), ' returned info = ', info
!    http://www.netlib.org/lapack/explore-3.1.1-html/zgesvx.f.html
!*  INFO    (output) INTEGER
!*          = 0:  successful exit
!*          < 0:  if INFO = -i, the i-th argument had an illegal value
!*          > 0:  if INFO = i, and i is
!*                <= N:  U(i,i) is exactly zero.  The factorization has
!*                       been completed, but the factor U is exactly
!*                       singular, so the solution and error bounds
!*                       could not be computed. RCOND = 0 is returned.
!*                = N+1: U is nonsingular, but RCOND is less than machine
!*                       precision, meaning that the matrix is singular
!*                       to working precision.  Nevertheless, the
!*                       solution and error bounds are computed because
!*                       there are a number of situations where the
!*                       computed solution can be more accurate than the
!*                       value of RCOND would suggest.
    if(info < 0) then
      write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
      call messages_fatal(2)      
    else if(info == n+1) then
      message(2) = '(reciprocal of the condition number is less than machine precision)'
      call messages_warning(2)
    else
      write(message(2), '(a,i5,a)') 'Diagonal element ', info, ' of U is 0; matrix is singular.'
      call messages_fatal(2)
    end if
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

#ifdef R_TREAL
! ---------------------------------------------------------
!> computes the singular value decomposition of a real NxN matrix a(:,:)
subroutine dsingular_value_decomp(m, n, a, u, vt, sg_values)
  integer, intent(in)    :: m, n
  FLOAT,   intent(inout) :: a(:,:)          !< (m,n)
  FLOAT,   intent(out)   :: u(:,:), vt(:,:) !< (m,m) (n,n)  
  FLOAT,   intent(out)   :: sg_values(:)    !< (min(m,n))

  interface
    subroutine X(gesvd) ( jobu, jobvt, m, n, a, lda, s, u, ldu, &
      vt, ldvt, work, lwork, info )
      implicit none
      character(1), intent(in)    :: jobu, jobvt
      integer,      intent(in)    :: m, n
      FLOAT,        intent(inout) :: a, u, vt ! a(lda,n), u(ldu,m), vt(ldvt,n)
      FLOAT,        intent(out)   :: work     ! work(lwork)
      integer,      intent(in)    :: lda, ldu, ldvt, lwork
      integer,      intent(out)   :: info
      FLOAT,        intent(out)   :: s        ! s(min(m,n))
    end subroutine X(gesvd)
  end interface

  integer :: info, lwork
  FLOAT, allocatable :: work(:)

  PUSH_SUB(dsingular_value_decomp)

  ASSERT(n > 0)
  ASSERT(m > 0)


  ! double minimum lwork size to increase performance (see manpage)
  lwork = 2*( 2*min(m, n) + max(m, n) )

  SAFE_ALLOCATE(work(1:lwork)) ! query?

  call X(gesvd)( 'A', 'A', m, n, a(1, 1), lead_dim(a), sg_values(1), u(1, 1), lead_dim(u), vt(1, 1), &
                 lead_dim(vt), work(1), lwork, info )

  if(info /= 0) then
    write(message(1), '(3a, i7)') 'In dsingular_value_decomp, LAPACK ', TOSTRING(X(gesvd)), ' returned info = ', info
!    http://www.netlib.org/lapack/explore-3.1.1-html/dgesvd.f.html
!*  INFO    (output) INTEGER
!*          = 0:  successful exit.
!*          < 0:  if INFO = -i, the i-th argument had an illegal value.
!*          > 0:  if ZBDSQR did not converge, INFO specifies how many
!*                superdiagonals of an intermediate bidiagonal form B
!*                did not converge to zero. See the description of WORK
!*                above for details.
    if(info < 0) then
      write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
    else
      write(message(2), '(i5,a)') info, ' superdiagonal elements of an intermediate bidiagonal did not converge to zero.'
    end if
    call messages_fatal(2)
  end if

  SAFE_DEALLOCATE_A(work)
  POP_SUB(dsingular_value_decomp)
end subroutine dsingular_value_decomp


! ---------------------------------------------------------
!> computes inverse of a real NxN matrix a(:,:) using the SVD decomposition 
subroutine dsvd_inverse(m, n, a, threshold)
  integer, intent(in)           :: m, n
  FLOAT,   intent(inout)        :: a(:,:)    !< (m,n); a will be replaced by its inverse
  FLOAT,   intent(in), optional :: threshold

  FLOAT, allocatable :: u(:,:), vt(:,:)
  FLOAT, allocatable :: sg_values(:)
  FLOAT   :: tmp
  FLOAT   :: sg_inverse, threshold_
  integer :: j, k, l, minmn

  ASSERT(n > 0)
  ASSERT(m > 0)
  minmn = min(m,n)

  SAFE_ALLOCATE( u(1:m, 1:m))
  SAFE_ALLOCATE(vt(1:n, 1:n))
  SAFE_ALLOCATE(sg_values(1:minmn))

  PUSH_SUB(dsvd_inverse)

  call dsingular_value_decomp(m, n, a, u, vt, sg_values)

  threshold_ = CNST(1e-10)
  if(present(threshold)) threshold_ = threshold

  ! build inverse
  do j = 1, m
    do k = 1, n
      tmp = M_ZERO
      do l = 1, minmn
        if (sg_values(l) < threshold_) then
          !write(message(1), '(a)') 'In dsvd_inverse: singular value below threshold.'
          !call messages_warning(1)
          sg_inverse = M_ZERO
        else
          sg_inverse = M_ONE/sg_values(l)
        end if
        tmp = tmp + vt(l, k)*sg_inverse*u(j, l)
      end do
      a(j, k) = tmp
    end do
  end do

  SAFE_DEALLOCATE_A(sg_values)
  SAFE_DEALLOCATE_A(vt)
  SAFE_DEALLOCATE_A(u)
  POP_SUB(dsvd_inverse)
end subroutine dsvd_inverse

#elif R_TCOMPLEX
! ---------------------------------------------------------
!> computes the singular value decomposition of a complex MxN matrix a(:,:)
subroutine zsingular_value_decomp(m, n, a, u, vt, sg_values)
  integer, intent(in)    :: m, n
  CMPLX,   intent(inout) :: a(:,:)          !< (m,n)
  CMPLX,   intent(out)   :: u(:,:), vt(:,:) !< (n,n) and (m,m)  
  FLOAT,   intent(out)   :: sg_values(:)    !< (n)

  interface
    subroutine X(gesvd) ( jobu, jobvt, m, n, a, lda, s, u, ldu, &
      vt, ldvt, work, lwork, rwork, info )
      implicit none
      character(1), intent(in)    :: jobu, jobvt
      integer,      intent(in)    :: m, n
      CMPLX,        intent(inout) :: a, u, vt ! a(lda,n), u(ldu,m), vt(ldvt,n)
      CMPLX,        intent(out)   :: work     ! work(lwork)
      integer,      intent(in)    :: lda, ldu, ldvt, lwork
      integer,      intent(out)   :: info
      FLOAT,        intent(out)   :: s        ! s(min(m,n))
      FLOAT,        intent(inout) :: rwork    ! rwork(5*min(m,n))
    end subroutine X(gesvd)
  end interface

  integer :: info, lwork
  CMPLX, allocatable :: work(:)
  FLOAT, allocatable :: rwork(:)

  PUSH_SUB(zsingular_value_decomp)

  ASSERT(n > 0)
  ASSERT(m > 0)

  ! double minimum lwork size to increase performance (see manpage)
  lwork = 2*( 2*min(m, n) + max(m, n) )

  SAFE_ALLOCATE(work(1:lwork)) ! query?
  SAFE_ALLOCATE(rwork(1:5*min(m, n)))

  call X(gesvd)( 'A', 'A', m, n, a(1, 1), lead_dim(a), sg_values(1), u(1, 1), lead_dim(u), &
                 vt(1, 1), lead_dim(vt), work(1), lwork, rwork(1), info )

  if(info /= 0) then
    write(message(1), '(3a, i5)') 'In zsingular_value_decomp, LAPACK ', TOSTRING(X(gesvd)), ' returned info = ', info
!    http://www.netlib.org/lapack/explore-3.1.1-html/zgesvd.f.html
!*  INFO    (output) INTEGER
!*          = 0:  successful exit.
!*          < 0:  if INFO = -i, the i-th argument had an illegal value.
!*          > 0:  if ZBDSQR did not converge, INFO specifies how many
!*                superdiagonals of an intermediate bidiagonal form B
!*                did not converge to zero. See the description of RWORK
!*                above for details.
    if(info < 0) then
      write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
    else
      write(message(2), '(i5,a)') info, ' superdiagonal elements of an intermediate bidiagonal did not converge to zero.'
    end if
    call messages_fatal(2)
  end if

  SAFE_DEALLOCATE_A(rwork)
  SAFE_DEALLOCATE_A(work)
  POP_SUB(zsingular_value_decomp)
end subroutine zsingular_value_decomp


! ---------------------------------------------------------
!> computes inverse of a complex MxN matrix a(:,:) using the SVD decomposition 
subroutine zsvd_inverse(m, n, a, threshold)
  integer, intent(in)           :: m, n
  CMPLX,   intent(inout)        :: a(:,:)    !< (m,n); a will be replaced by its inverse transposed
  FLOAT,   intent(in), optional :: threshold

  CMPLX, allocatable :: u(:,:), vt(:,:)
  FLOAT, allocatable :: sg_values(:)
  CMPLX   :: tmp
  FLOAT   :: sg_inverse, threshold_
  integer :: j, k, l, minmn

  ASSERT(n > 0)
  ASSERT(m > 0)
  minmn = min(m,n)

  SAFE_ALLOCATE( u(1:m, 1:m))
  SAFE_ALLOCATE(vt(1:n, 1:n))
  SAFE_ALLOCATE(sg_values(1:minmn))

  PUSH_SUB(zsvd_inverse)

  call zsingular_value_decomp(m, n, a, u, vt, sg_values)

  threshold_ = CNST(1e-10)
  if(present(threshold)) threshold_ = threshold

  ! build inverse
  do j = 1, m
    do k = 1, n
      tmp = M_ZERO
      do l = 1, minmn
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

  call X(trtri)('U', 'N', n, a(1, 1), lead_dim(a), info)

  if(info /= 0) then
    write(message(1), '(5a,i5)') &
      'In ', TOSTRING(Xinvert_upper_triangular), ', LAPACK ', TOSTRING(X(trtri)), ' returned error message ', info
!http://www.netlib.org/lapack/explore-3.1.1-html/dtrtri.f.html
!*  INFO    (output) INTEGER
!*          = 0: successful exit
!*          < 0: if INFO = -i, the i-th argument had an illegal value
!*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
!*               matrix is singular and its inverse can not be computed.
    if(info < 0) then
      write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
    else
      write(message(2), '(a,i5,a)') 'Diagonal element ', info, ' is 0; matrix is singular.'
    end if
    call messages_fatal(2)
  end if

  POP_SUB(X(invert_upper_triangular))
end subroutine X(invert_upper_triangular)



subroutine X(least_squares_vec)(nn, aa, bb, xx)
  integer, intent(in)    :: nn
  R_TYPE,  intent(inout) :: aa(:, :)
  R_TYPE,  intent(in)    :: bb(:)
  R_TYPE,  intent(out)   :: xx(:)

  R_TYPE :: dlwork
  R_TYPE, allocatable :: work(:)
  integer :: rank, info
  FLOAT, allocatable :: ss(:)
#ifndef R_TREAL
  FLOAT, allocatable :: rwork(:)
#endif 
  
  PUSH_SUB(X(least_squares_vec))


  xx(1:nn) = bb(1:nn)

  SAFE_ALLOCATE(ss(1:nn))

! MJV 2016 11 09 : TODO: this is callable with complex, but does nothing!!!
#ifdef R_TREAL
  
  call lapack_gelss(nn, nn, 1, aa(1, 1), lead_dim(aa), xx(1), nn, ss(1), CNST(-1.0), rank, dlwork, -1, info)

  SAFE_ALLOCATE(work(1:int(dlwork)))

  call lapack_gelss(nn, nn, 1, aa(1, 1), lead_dim(aa), xx(1), nn, ss(1), CNST(-1.0), rank, work(1), int(dlwork), info)
#else
  SAFE_ALLOCATE(rwork(1:5*nn))
  call lapack_gelss(nn, nn, 1, aa(1, 1), lead_dim(aa), xx(1), nn, ss(1), CNST(-1.0), rank, dlwork, -1, rwork(1), info)

  SAFE_ALLOCATE(work(1:int(dlwork)))

  call lapack_gelss(nn, nn, 1, aa(1, 1), lead_dim(aa), xx(1), nn, ss(1), CNST(-1.0), rank, work(1), int(dlwork), rwork(1), info)
  SAFE_DEALLOCATE_A(rwork)
#endif

  if(info /= 0) then
    write(message(1), '(5a,i5)') &
      'In ', TOSTRING(X(lalg_least_squares_vec)), ', LAPACK ', TOSTRING(X(gelss)), ' returned error mess  age ', info
  !https://www.netlib.org/lapack/lapack-3.1.1/html/zgelss.f.html
  !*  INFO    (output) INTEGER
  !*          = 0:  successful exit
  !*          < 0:  if INFO = -i, the i-th argument had an illegal value.
  !*          > 0:  the algorithm for computing the SVD failed to converge;
  !*                if INFO = i, i off-diagonal elements of an intermediate
  !*                bidiagonal form did not converge to zero.
    if(info < 0) then
      write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
    else
      write(message(2), '(a,i5,a)') 'Off-diagonal element ', info, ' of an intermediate bidiagonal form did not converge to zero.'
    end if
    call messages_fatal(2)
  end if

  POP_SUB(X(least_squares_vec))

end subroutine X(least_squares_vec)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#if defined(SINGLE_PRECISION)
#  define DLAPACK(x) s ## x
#  define ZLAPACK(x) c ## x
#else
#  define DLAPACK(x) d ## x
#  define ZLAPACK(x) z ## x
#endif


! ---------------------------------------------------------
! Auxiliary functions.
FLOAT function sfmin()
  interface
    FLOAT function DLAPACK(lamch)(cmach)
      character(1), intent(in) :: cmach
    end function DLAPACK(lamch)
  end interface
  
  sfmin = DLAPACK(lamch)('S')
end function sfmin

  
! ---------------------------------------------------------
!> Compute the Cholesky decomposition of real symmetric positive definite
!! matrix a, dim(a) = n x n. On return a = u^T u with u upper triangular
!! matrix.
subroutine dcholesky(n, a, bof, err_code)
  integer,           intent(in)    :: n
  FLOAT,             intent(inout) :: a(:, :)
  logical, optional, intent(inout) :: bof      ! Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer :: info
  logical :: bof_

  call profiling_in(cholesky_prof, "CHOLESKY")
  PUSH_SUB(dcholesky)

  bof_ = .true.
  if(present(bof)) then
    bof_ = bof
  end if

  call lapack_potrf('U', n, a(1, 1), LD(a), info)
  if(info.ne.0) then
    if(bof_) then
      write(message(1), '(3a,i5)') 'In dcholesky, LAPACK ', TOSTRING(DLAPACK(potrf)), ' returned error message ', info
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
  POP_SUB(dcholesky)
end subroutine dcholesky


! ---------------------------------------------------------
!> Compute the Cholesky decomposition of a complex Hermitian positive definite
!! matrix a, dim(a) = n x n. On return a = u^+ u with u upper triangular
!! matrix.
subroutine zcholesky(n, a, bof, err_code)
  integer,           intent(in)    :: n
  CMPLX,             intent(inout) :: a(:, :)
  logical, optional, intent(inout) :: bof     ! Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer :: info
  logical :: bof_
  call profiling_in(cholesky_prof, "CHOLESKY")
  PUSH_SUB(zcholesky)

  bof_ = .true.
  if(present(bof)) then
    bof_ = bof
  end if

  call lapack_potrf('U', n, a(1, 1), LD(a), info)

  if(info.ne.0) then
    if(bof_) then
      write(message(1), '(3a,i5)') 'In zcholesky, LAPACK ', TOSTRING(ZLAPACK(potrf)), ' returned error message ', info
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
  POP_SUB(zcholesky)
end subroutine zcholesky


! ---------------------------------------------------------
!> Computes all the eigenvalues and the eigenvectors of a real
!! generalized symmetric-definite eigenproblem, of the form  A*x=(lambda)*B*x
!! A*Bx=(lambda)*x, or B*A*x=(lambda)*x.
!! Here A and B are assumed to be symmetric and B is also positive definite.
subroutine dgeneigensolve(n, a, b, e, bof, err_code)
  integer,           intent(in)    :: n
  FLOAT,             intent(inout) :: a(n,n)
  FLOAT,             intent(inout) :: b(n,n)
  FLOAT,             intent(out)   :: e(n)
  logical, optional, intent(inout) :: bof      ! Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer :: info, lwork, ii, jj
  logical :: bof_
  FLOAT, allocatable :: work(:), diag(:)

  call profiling_in(eigensolver_prof, "DENSE_EIGENSOLVER")
  PUSH_SUB(dgeneigensolve)

  bof_ = .true.
  if(present(bof)) then
    bof_ = bof
  end if

  SAFE_ALLOCATE(diag(1:n))

  ! store the diagonal of b  
  forall(ii = 1:n) diag(ii) = b(ii, ii)

  lwork = 5*n
  SAFE_ALLOCATE(work(1:lwork))
  call lapack_sygv(1, 'V', 'U', n, a(1, 1), n, b(1, 1), n, e(1), work(1), lwork, info)
  SAFE_DEALLOCATE_A(work)

  ! b was destroyed, so we rebuild it
  do ii = 1, n
    do jj = 1, ii - 1
      b(jj, ii) = b(ii, jj)
    end do
    b(ii, ii) = diag(ii)
  end do

  SAFE_DEALLOCATE_A(diag)

  if(info.ne.0) then
    if(bof_) then
      write(message(1),'(3a,i5)') 'In dgeneigensolve, LAPACK ', TOSTRING(DLAPACK(sygv)), ' returned error message ', info
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
  POP_SUB(dgeneigensolve)
end subroutine dgeneigensolve


! ---------------------------------------------------------
!> Computes all the eigenvalues and the eigenvectors of a complex
!! generalized Hermitian-definite eigenproblem, of the form  A*x=(lambda)*B*x,
!! A*Bx=(lambda)*x, or B*A*x=(lambda)*x.
!! Here A and B are assumed to be Hermitian and B is also positive definite.
subroutine zgeneigensolve(n, a, b, e, bof, err_code)
  integer,           intent(in)    :: n
  CMPLX,             intent(inout) :: a(n,n)
  CMPLX,             intent(inout) :: b(n,n)
  FLOAT,             intent(out)   :: e(n)
  logical, optional, intent(inout) :: bof      ! Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer            :: info, lwork, ii, jj
  logical            :: bof_
  FLOAT, allocatable :: rwork(:)
  CMPLX, allocatable :: work(:), diag(:)

  call profiling_in(eigensolver_prof, "DENSE_EIGENSOLVER")
  PUSH_SUB(zgeneigensolve)

  bof_ = .true.
  if(present(bof)) then
    bof_ = bof
  end if

  SAFE_ALLOCATE(diag(1:n))
  
  ! store the diagonal of b
  forall(ii = 1:n) diag(ii) = b(ii, ii)

  lwork = 5*n
  SAFE_ALLOCATE(work(1:lwork))
  SAFE_ALLOCATE(rwork(1:max(1, 3*n-2)))
  call lapack_hegv(1, 'V', 'U', n, a(1, 1), n, b(1, 1), n, e(1), work(1), lwork, rwork(1), info)
  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(rwork)

  ! b was destroyed, so we rebuild it
  do ii = 1, n
    do jj = 1, ii - 1
      b(jj, ii) = b(ii, jj)
    end do
    b(ii, ii) = diag(ii)
  end do

  SAFE_DEALLOCATE_A(diag)

  if(info.ne.0) then
    if(bof_) then
      write(message(1),'(3a,i5)') 'In zgeneigensolve, LAPACK ', TOSTRING(ZLAPACK(hegv)), ' returned error message ', info
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
  POP_SUB(zgeneigensolve)
end subroutine zgeneigensolve


! ---------------------------------------------------------
!> Computes all the eigenvalues and the right (left) eigenvectors of a complex
!! (non hermitian) eigenproblem, of the form  A*x=(lambda)*x
subroutine zeigensolve_nonh(n, a, e, err_code, side)
  integer,           intent(in)      :: n
  CMPLX,             intent(inout)   :: a(n,n)
  CMPLX,             intent(out)     :: e(n)
  integer, optional, intent(out)     :: err_code
  character(1), optional, intent(in) :: side ! which eigenvectors ('L' or 'R')

  integer            :: info, lwork
  FLOAT, allocatable :: rwork(:)
  CMPLX, allocatable :: work(:), vl(:, :) ,vr(:, :)
  character(1)       :: side_

  PUSH_SUB(zeigensolve_nonh)

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
  SAFE_ALLOCATE(work(n))
  SAFE_ALLOCATE(vl(1, 1))
  SAFE_ALLOCATE(vr(1, 1))
  SAFE_ALLOCATE(rwork(1))
  call lapack_geev('N', 'V', n, a(1, 1), n, e(1), vl(1, 1), 1, vr(1, 1), n, work(1), lwork, rwork(1), info)

  lwork = int(work(1))
  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(vl)
  SAFE_DEALLOCATE_A(vr)
  SAFE_DEALLOCATE_A(rwork)

  SAFE_ALLOCATE(work(1:lwork))
  SAFE_ALLOCATE(rwork(1:max(1, 2*n)))
  if (side_.eq.'L'.or.side_.eq.'l') then
      SAFE_ALLOCATE(vl(1:n, 1:n))
      SAFE_ALLOCATE(vr(1, 1))
      call lapack_geev('V', 'N', n, a(1, 1), n, e(1), vl(1, 1), n, vr(1,1), 1, work(1), lwork, rwork(1), info)
      a(:, :) = vl(:, :)
  else
      SAFE_ALLOCATE(vl(1, 1))
      SAFE_ALLOCATE(vr(1:n, 1:n))
      call lapack_geev('N', 'V', n, a(1, 1), n, e(1), vl(1, 1), 1, vr(1,1), n, work(1), lwork, rwork(1), info)
      a(:, :) = vr(:, :)
  end if
  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(rwork)
  SAFE_DEALLOCATE_A(vr)
  SAFE_DEALLOCATE_A(vl)

  if(info.ne.0) then
    write(message(1),'(3a,i5)') 'In zeigensolve_nonh, LAPACK ', TOSTRING(ZLAPACK(geev)), ' returned error message ', info
    call messages_fatal(1)
  end if
  if(present(err_code)) then
    err_code = info
  end if

  POP_SUB(zeigensolve_nonh)
end subroutine zeigensolve_nonh


! ---------------------------------------------------------
!> Computes all the eigenvalues and the right (left) eigenvectors of a real
!! generalized (non hermitian) eigenproblem, of the form  A*x=(lambda)*x
subroutine deigensolve_nonh(n, a, e, err_code, side)
  integer,           intent(in)      :: n
  FLOAT,             intent(inout)   :: a(n,n)
  FLOAT,             intent(out)     :: e(n)
  integer, optional, intent(out)     :: err_code
  character(1), optional, intent(in) :: side ! which eigenvectors ('L' or 'R')

  integer            :: info, lwork
  FLOAT, allocatable :: rwork(:), work(:), vl(:, :) ,vr(:, :)
  character(1)       :: side_

  PUSH_SUB(deigensolve_nonh)

  if (present(side)) then
    side_ = side
  else
    side_ = 'R'
  end if

  lwork = -1
  ! A bug in the query mode of dgeev demands that the working array has to be larger than 1
  ! problem here is that a value is written somewhere into the array whether it is
  ! allocated or not. I noticed that it happens (hopefully) always at an index which
  ! is below the matrix dimension.
  SAFE_ALLOCATE(work(n))
  SAFE_ALLOCATE(vl(1, 1))
  SAFE_ALLOCATE(vr(1, 1))
  SAFE_ALLOCATE(rwork(1))
  call lapack_geev('N', 'V', n, a(1, 1), n, e(1), vl(1, 1), 1, vr(1, 1), n, work(1), lwork, rwork(1), info)

  lwork = int(work(1))
  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(vl)
  SAFE_DEALLOCATE_A(vr)
  SAFE_DEALLOCATE_A(rwork)

  SAFE_ALLOCATE(work(1:lwork))
  SAFE_ALLOCATE(rwork(1:max(1, 2*n)))
  if (side_.eq.'L'.or.side_.eq.'l') then
      SAFE_ALLOCATE(vl(1:n, 1:n))
      SAFE_ALLOCATE(vr(1, 1))
      call lapack_geev('V', 'N', n, a(1, 1), n, e(1), vl(1, 1), n, vr(1,1), 1, work(1), lwork, rwork(1), info)
      a(:, :) = vl(:, :)
  else
      SAFE_ALLOCATE(vl(1, 1))
      SAFE_ALLOCATE(vr(1:n, 1:n))
      call lapack_geev('N', 'V', n, a(1, 1), n, e(1), vl(1, 1), 1, vr(1,1), n, work(1), lwork, rwork(1), info)
      a(:, :) = vr(:, :)
  end if
  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(rwork)
  SAFE_DEALLOCATE_A(vr)
  SAFE_DEALLOCATE_A(vl)

  if(info.ne.0) then
    write(message(1),'(3a,i5)') 'In deigensolve_nonh, LAPACK ', TOSTRING(DLAPACK(geev)), ' returned error message ', info
    call messages_fatal(1)
  end if
  if(present(err_code)) then
    err_code = info
  end if

  POP_SUB(deigensolve_nonh)
end subroutine deigensolve_nonh

! ---------------------------------------------------------
!> Computes the k lowest eigenvalues and the eigenvectors of a real
!! generalized symmetric-definite eigenproblem, of the form  A*x=(lambda)*B*x.
!! Here A and B are assumed to be symmetric and B is also positive definite.
subroutine dlowest_geneigensolve(k, n, a, b, e, v, bof, err_code)
  integer,           intent(in)    :: k, n
  FLOAT,             intent(in)    :: a(n,n)
  FLOAT,             intent(in)    :: b(n,n)
  FLOAT,             intent(out)   :: e(n)
  FLOAT,             intent(out)   :: v(n, k)
  logical, optional, intent(inout) :: bof      ! Bomb on failure.
  integer, optional, intent(out)   :: err_code

  interface
    subroutine DLAPACK(sygvx)(itype, jobz, range, uplo, n, a, lda, b, ldb, &
      vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
      integer,      intent(in)  :: itype, n, lda, ldb, il, iu, ldz, lwork
      character(1), intent(in)  :: jobz, range, uplo
      integer,      intent(out) :: m, iwork, ifail, info
      FLOAT,        intent(in)  :: vl, vu, abstol
      FLOAT,        intent(in)  :: a, b
      FLOAT,        intent(out) :: w, z, work
    end subroutine DLAPACK(sygvx)
  end interface
  
  integer            :: m, iwork(5*n), ifail(n), info, lwork
  logical            :: bof_
  FLOAT              :: abstol
  FLOAT, allocatable :: work(:)
  
  PUSH_SUB(dlowest_geneigensolve)

  bof_ = .true.
  if(present(bof)) then
    bof_ = bof
  end if

  abstol = 2*sfmin()

  ! Work size query.
  SAFE_ALLOCATE(work(1))
  call DLAPACK(sygvx)(1, 'V', 'I', 'U', n, a(1, 1), LD(a), b(1, 1), LD(b), M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), LD(v), work(1), -1, iwork(1), ifail(1), info)
  lwork = int(work(1))
  SAFE_DEALLOCATE_A(work)

  SAFE_ALLOCATE(work(1:lwork))

  call DLAPACK(sygvx)(1, 'V', 'I', 'U', n, a(1, 1), LD(a), b(1, 1), LD(b), M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), LD(v), work(1), lwork, iwork(1), ifail(1), info)

  SAFE_DEALLOCATE_A(work)

  if(info.ne.0) then
    if(bof_) then
      write(message(1),'(3a,i5)') 'In dlowest_geneigensolve, LAPACK ', TOSTRING(DLAPACK(sygvx)), ' returned error message ', info
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


! ---------------------------------------------------------
!> Computes the k lowest eigenvalues and the eigenvectors of a complex
!! generalized Hermitian-definite eigenproblem, of the form  A*x=(lambda)*B*x.
!! Here A and B are assumed to be Hermitian and B is also positive definite.
subroutine zlowest_geneigensolve(k, n, a, b, e, v, bof, err_code)
  integer,           intent(in)    :: k, n
  CMPLX,             intent(in)    :: a(n,n)
  CMPLX,             intent(in)    :: b(n,n)
  FLOAT,             intent(out)   :: e(n)
  CMPLX,             intent(out)   :: v(n, k)
  logical, optional, intent(inout) :: bof      ! Bomb on failure.
  integer, optional, intent(out)   :: err_code

  interface
    subroutine ZLAPACK(hegvx)(itype, jobz, range, uplo, n, a, lda, b, ldb, &
      vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
      integer,      intent(in)    :: itype, n, lda, ldb, il, iu, ldz, lwork
      character(1), intent(in)    :: jobz, range, uplo
      integer,      intent(out)   :: m, iwork, ifail, info
      FLOAT,        intent(in)    :: vl, vu, abstol
      FLOAT,        intent(out)   :: w, rwork
      CMPLX,        intent(in)    :: a, b
      CMPLX,        intent(out)   :: z, work
    end subroutine ZLAPACK(hegvx)
  end interface

  integer            :: m, iwork(5*n), ifail(n), info, lwork
  logical            :: bof_
  FLOAT              :: abstol
  FLOAT              :: rwork(7*n)
  CMPLX, allocatable :: work(:)

  PUSH_SUB(zlowest_geneigensolve)

  bof_ = .true.
  if(present(bof)) then
    bof_ = bof
  end if

  abstol = 2*sfmin()

  ! Work size query.
  SAFE_ALLOCATE(work(1))
  call ZLAPACK(hegvx)(1, 'V', 'I', 'U', n, a(1, 1), LD(a), b(1, 1), LD(b), M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), LD(v), work(1), -1, rwork(1), iwork(1), ifail(1), info)
  lwork = int(real(work(1)))
  SAFE_DEALLOCATE_A(work)

  SAFE_ALLOCATE(work(1:lwork))
  call ZLAPACK(hegvx)(1, 'V', 'I', 'U', n, a(1, 1), LD(a), b(1, 1), LD(b), M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), LD(v), work(1), lwork, rwork(1), iwork(1), ifail(1), info)
  SAFE_DEALLOCATE_A(work)

  if(info.ne.0) then
    if(bof_) then
      write(message(1),'(3a,i5)') 'In zlowest_geneigensolve, LAPACK ', TOSTRING(ZLAPACK(hegvx)), ' returned error message ', info
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


! ---------------------------------------------------------
!> Computes all eigenvalues and eigenvectors of a real symmetric square matrix A.
subroutine deigensolve(n, a, e, bof, err_code)
  integer, intent(in)              :: n
  FLOAT,   intent(inout)           :: a(n,n)
  FLOAT,   intent(out)             :: e(n)
  logical, optional, intent(inout) :: bof      !< Bomb on failure.
  integer, optional, intent(out)   :: err_code

  logical            :: bof_
  integer            :: info, lwork
  FLOAT, allocatable :: work(:)

  PUSH_SUB(deigensolve)
  call profiling_in(eigensolver_prof, "DENSE_EIGENSOLVER")

  bof_ = .true.
  if(present(bof)) then
    bof_ = bof
  end if

  lwork = 6*n
  SAFE_ALLOCATE(work(1:lwork))
  call lapack_syev('V', 'U', n, a(1, 1), LD(a), e(1), work(1), lwork, info)
  SAFE_DEALLOCATE_A(work)

  if(info.ne.0) then
    if(bof_) then
      write(message(1),'(3a,i5)') 'In deigensolve, LAPACK ', TOSTRING(DLAPACK(syev)), ' returned error message ', info
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

! ---------------------------------------------------------
!> Computes all eigenvalues and eigenvectors of a complex Hermitian square matrix A.
subroutine zeigensolve(n, a, e, bof, err_code)
  integer, intent(in)    :: n
  CMPLX,   intent(inout) :: a(n,n)
  FLOAT,   intent(out)   :: e(n)
  logical, optional, intent(inout) :: bof      ! Bomb on failure.
  integer, optional, intent(out)   :: err_code

  integer            :: info, lwork
  logical            :: bof_
  CMPLX, allocatable :: work(:)
  FLOAT, allocatable :: rwork(:)

  PUSH_SUB(zeigensolve)
  call profiling_in(eigensolver_prof, "DENSE_EIGENSOLVER")

  bof_ = .true.
  if(present(bof)) then
    bof_ = bof
  end if

  lwork = 6*n
  SAFE_ALLOCATE(work(1:lwork))
  SAFE_ALLOCATE(rwork(1:max(1, 3*n-2)))
  call lapack_heev('V','U', n, a(1, 1), LD(a), e(1), work(1), lwork, rwork(1), info)
  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(rwork)

  if(info.ne.0) then
    if(bof_) then
      write(message(1),'(3a,i5)') 'In zeigensolve, LAPACK ', TOSTRING(ZLAPACK(heev)), ' returned error message ', info
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


! ---------------------------------------------------------
!> Computes the k lowest eigenvalues and the eigenvectors of a real
!! standard symmetric-definite eigenproblem, of the form  A*x=(lambda)*x.
!! Here A is assumed to be symmetric.
subroutine dlowest_eigensolve(k, n, a, e, v)
  integer, intent(in)  :: k, n
  FLOAT,   intent(in)  :: a(n, n)
  FLOAT,   intent(out) :: e(n)
  FLOAT,   intent(out) :: v(n, k)

  interface
    subroutine DLAPACK(syevx)(jobz, range, uplo, n, a, lda, &
      vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
      integer,      intent(in)  :: n, lda, il, iu, ldz, lwork
      character(1), intent(in)  :: jobz, range, uplo
      integer,      intent(out) :: m, iwork, ifail, info
      FLOAT,        intent(in)  :: vl, vu, abstol
      FLOAT,        intent(in)  :: a
      FLOAT,        intent(out) :: w, z, work
    end subroutine DLAPACK(syevx)
  end interface
  
  integer            :: m, iwork(5*n), ifail(n), info, lwork
  FLOAT              :: abstol
  FLOAT, allocatable :: work(:)
  
  PUSH_SUB(dlowest_eigensolve)

  abstol = 2*sfmin()

  ! Work size query.
  SAFE_ALLOCATE(work(1))
  call DLAPACK(syevx)('V', 'I', 'U', n, a(1, 1), n, M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), n, work(1), -1, iwork(1), ifail(1), info)
  lwork = int(work(1))
  SAFE_DEALLOCATE_A(work)

  SAFE_ALLOCATE(work(1:lwork))
  call DLAPACK(syevx)('V', 'I', 'U', n, a(1, 1), n, M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), n, work(1), lwork, iwork(1), ifail(1), info)
  SAFE_DEALLOCATE_A(work)

  if(info.ne.0) then
    write(message(1),'(3a,i5)') &
      'In dlowest_eigensolve, LAPACK ', TOSTRING(DLAPACK(syevx)), ' returned error message ', info
    call messages_fatal(1)
  end if

  POP_SUB(dlowest_eigensolve)
end subroutine dlowest_eigensolve


! ---------------------------------------------------------
!> Computes the k lowest eigenvalues and the eigenvectors of a complex
!! standard Hermitian-definite eigenproblem, of the form  A*x=(lambda)*x.
!! Here A is assumed to be Hermitian.
subroutine zlowest_eigensolve(k, n, a, e, v)
  integer, intent(in)  :: k, n
  CMPLX,   intent(in)  :: a(n, n)
  FLOAT,   intent(out) :: e(n)
  CMPLX,   intent(out) :: v(n, k)

  interface
    subroutine ZLAPACK(heevx)(jobz, range, uplo, n, a, lda, &
      vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
      integer,      intent(in)  :: n, lda, il, iu, ldz, lwork
      character(1), intent(in)  :: jobz, range, uplo
      integer,      intent(out) :: m, iwork, ifail, info
      FLOAT,        intent(in)  :: vl, vu, abstol
      FLOAT,        intent(out) :: w
      CMPLX,        intent(in)  :: a
      CMPLX,        intent(out) :: z, work
    end subroutine ZLAPACK(heevx)
  end interface
  
  integer            :: m, iwork(5*n), ifail(n), info, lwork
  FLOAT              :: abstol
  CMPLX, allocatable :: work(:)
   
  PUSH_SUB(zlowest_eigensolve)

  abstol = 2*sfmin()

  ! Work size query.
  SAFE_ALLOCATE(work(1))
  call ZLAPACK(heevx)('V', 'I', 'U', n, a(1, 1), n, M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), n, work(1), -1, iwork(1), ifail(1), info)
  lwork = int(work(1))
  SAFE_DEALLOCATE_A(work)

  SAFE_ALLOCATE(work(1:lwork))
  call ZLAPACK(heevx)('V', 'I', 'U', n, a(1, 1), n, M_ZERO, M_ZERO, &
    1, k, abstol, m, e(1), v(1, 1), n, work(1), lwork, iwork(1), ifail(1), info)
  SAFE_DEALLOCATE_A(work)

  if(info.ne.0) then
    write(message(1),'(3a,i5)') &
      'In zlowest_eigensolve, LAPACK ', TOSTRING(ZLAPACK(heevx)), ' returned error message ', info
    call messages_fatal(1)
  end if

  POP_SUB(zlowest_eigensolve)
end subroutine zlowest_eigensolve


! ---------------------------------------------------------
!> Invert a real symmetric square matrix a
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

  ! No PUSH_SUB, called to often

  SAFE_ALLOCATE(work(1:n))
  SAFE_ALLOCATE(ipiv(1:n))

  call DLAPACK(getrf)(n, n, a(1, 1), n, ipiv(1), info)
  if(info < 0) then
    write(message(1), '(3a, i3)') 'In ddeterminant, LAPACK ', TOSTRING(DLAPACK(getrf)), ' returned info = ', info
    call messages_fatal(1)
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
      write(message(1), '(3a, i3)') 'In ddeterminant, LAPACK ', TOSTRING(DLAPACK(getri)), ' returned info = ', info
      call messages_fatal(1)
    end if
  end if

  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(ipiv)

end function ddeterminant


! ---------------------------------------------------------
!> Invert a complex Hermitian square matrix a
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

  PUSH_SUB(zdeterminant)

  SAFE_ALLOCATE(work(1:n))
  SAFE_ALLOCATE(ipiv(1:n))

  call ZLAPACK(getrf)(n, n, a(1, 1), n, ipiv(1), info)
  if(info < 0) then
    write(message(1), '(3a, i3)') 'In zdeterminant, LAPACK ', TOSTRING(ZLAPACK(getrf)), ' returned info = ', info
    call messages_fatal(1)
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
      write(message(1), '(3a, i3)') 'In zdeterminant, LAPACK ', TOSTRING(ZLAPACK(getri)), ' returned info = ', info
      call messages_fatal(1)
    end if
  end if

  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(ipiv)
  POP_SUB(zdeterminant)
end function zdeterminant

! ---------------------------------------------------------
!> Invert a real symmetric square matrix a
subroutine dsym_inverter(uplo, n, a)
  character(1), intent(in)      :: uplo
  integer, intent(in)           :: n
  FLOAT,   intent(inout)        :: a(n,n)

  interface
    subroutine DLAPACK(sytrf) (uplo, n, a, lda, ipiv, work, lwork, info)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda, lwork
      FLOAT,        intent(inout) :: a
      integer,      intent(out)   :: ipiv
      FLOAT,        intent(inout) :: work
      integer,      intent(out)   :: info
    end subroutine DLAPACK(sytrf)

    subroutine DLAPACK(sytri) (uplo, n, a, lda, ipiv, work, info )
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda
      FLOAT,        intent(inout) :: a
      integer,      intent(in)    :: ipiv
      FLOAT,        intent(inout) :: work
      integer,      intent(out)   :: info
    end subroutine DLAPACK(sytri)
  end interface

  integer :: info
  integer, allocatable :: ipiv(:)
  FLOAT, allocatable :: work(:)

  PUSH_SUB(dsym_inverter)

  SAFE_ALLOCATE(work(1:n))
  SAFE_ALLOCATE(ipiv(1:n))

  call DLAPACK(sytrf)(uplo, n, a(1, 1), LD(a), ipiv(1), work(1), n, info)
  if(info < 0) then
    write(message(1), '(3a, i3)') 'In dsym_inverter, LAPACK ', TOSTRING(DLAPACK(sytrf)), ' returned info = ', info
    call messages_fatal(1)
  end if

  call DLAPACK(sytri)(uplo, n, a(1, 1), LD(a), ipiv(1), work(1), info)
  if(info /= 0) then
    write(message(1), '(3a, i3)') 'In dsym_inverter, LAPACK ', TOSTRING(DLAPACK(sytri)), ' returned info = ', info
    call messages_fatal(1)
  end if

  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(ipiv)
  POP_SUB(dsym_inverter)
end subroutine dsym_inverter

! ---------------------------------------------------------
!> Invert a complex symmetric square matrix a
subroutine zsym_inverter(uplo, n, a)
  character(1), intent(in)      :: uplo
  integer, intent(in)           :: n
  CMPLX,   intent(inout)        :: a(n,n)

  interface
    subroutine ZLAPACK(sytrf) (uplo, n, a, lda, ipiv, work, lwork, info)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda, lwork
      CMPLX,        intent(inout) :: a 
      integer,      intent(out)   :: ipiv
      CMPLX,        intent(inout) :: work
      integer,      intent(out)   :: info
    end subroutine ZLAPACK(sytrf)

    subroutine ZLAPACK(sytri) (uplo, n, a, lda, ipiv, work, info )
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda
      CMPLX,        intent(inout) :: a
      integer,      intent(in)    :: ipiv
      CMPLX,        intent(inout) :: work
      integer,      intent(out)   :: info
    end subroutine ZLAPACK(sytri)
  end interface

  integer :: info
  integer, allocatable :: ipiv(:)
  CMPLX, allocatable :: work(:)

  PUSH_SUB(zsym_inverter)

  SAFE_ALLOCATE(work(1:n))
  SAFE_ALLOCATE(ipiv(1:n))
  call ZLAPACK(sytrf)(uplo, n, a(1, 1), LD(a), ipiv(1), work(1), n, info)
  if(info < 0) then
    write(message(1), '(3a, i3)') 'In zsym_inverter, LAPACK ', TOSTRING(ZLAPACK(sytrf)), ' returned info = ', info
    call messages_fatal(1)
  end if

  call ZLAPACK(sytri)(uplo, n, a(1, 1), LD(a), ipiv(1), work(1), info)
  if(info /= 0) then
    write(message(1), '(3a, i3)') 'In zsym_inverter, LAPACK ', TOSTRING(ZLAPACK(zsytri)), ' returned info = ', info
    call messages_fatal(1)
  end if

  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(ipiv)
  POP_SUB(zsym_inverter)
end subroutine zsym_inverter

! ---------------------------------------------------------
!> compute the solution to a real system of linear equations A*X = B,
!!  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
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

  ! no PUSH_SUB, called too often

  SAFE_ALLOCATE(ipiv(1:n))
  SAFE_ALLOCATE(iwork(1:n))
  SAFE_ALLOCATE(ferr(1:nhrs))
  SAFE_ALLOCATE(berr(1:nhrs))
  SAFE_ALLOCATE(work(1:4*n))
  SAFE_ALLOCATE(r(1:n))
  SAFE_ALLOCATE(c(1:n))
  SAFE_ALLOCATE(af(1:n, 1:n))

  call DLAPACK(gesvx) ("N", "N", n, nhrs, a(1, 1), n, af(1, 1), n, ipiv(1), equed, r(1), c(1), b(1, 1), n, x(1, 1), n, &
    rcond, ferr(1), berr(1), work(1), iwork(1), info)

  if(info /= 0) then
    write(message(1), '(3a, i3)') 'In dlinsyssolve, LAPACK ', TOSTRING(DLAPACK(gesvx)), ' returned info = ', info
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

! ---------------------------------------------------------
!> compute the solution to a complex system of linear equations A*X = B,
!!  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
subroutine zlinsyssolve(n, nhrs, a, b, x)
  integer, intent(in)    :: n, nhrs
  CMPLX,   intent(inout) :: a(n, n), b(n, nhrs)
  CMPLX,   intent(out)   :: x(n, nhrs)

  interface
    subroutine ZLAPACK(gesvx) (fact, trans, n, nhrs, a, lda, af, ldaf, ipiv, equed, r, &
      c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info)
      character(1), intent(in)    :: fact, trans
      integer,      intent(in)    :: n, nhrs, lda, ldaf, ldb, ldx
      CMPLX,        intent(inout) :: a, af, b            ! a(lda, n), af(ldaf, n), b(ldb, nhrs)
      FLOAT,        intent(inout) :: r, c                ! r(n), c(n)
      integer,      intent(inout) :: ipiv                ! ipiv(n)
      FLOAT,        intent(out)   :: ferr, berr          ! ferr(nhrs), berr(nhrs)
      FLOAT,        intent(out)   :: rcond
      CMPLX,        intent(out)   :: x, work             ! x(ldx, nhrs), work(2*n)
      character(1), intent(inout) :: equed
      FLOAT,        intent(out)   :: rwork               ! rwork(2*n)
      integer,      intent(out)   :: info
    end subroutine ZLAPACK(gesvx)
  end interface

  integer              :: info
  integer, allocatable :: ipiv(:)
  FLOAT,   allocatable :: rwork(:), ferr(:), berr(:), r(:), c(:)
  FLOAT                :: rcond
  CMPLX, allocatable   :: work(:), af(:,:)
  character(1)         :: equed

  PUSH_SUB(zlinsyssolve)

  SAFE_ALLOCATE(ipiv(1:n))
  SAFE_ALLOCATE(rwork(1:2*n))
  SAFE_ALLOCATE(ferr(1:nhrs))
  SAFE_ALLOCATE(berr(1:nhrs))
  SAFE_ALLOCATE(work(1:2*n))
  SAFE_ALLOCATE(r(1:n))
  SAFE_ALLOCATE(c(1:n))
  SAFE_ALLOCATE(af(1:n, 1:n))

  equed = 'N'

  call ZLAPACK(gesvx) ("N", "N", n, nhrs, a(1, 1), n, af(1, 1), n, ipiv(1), equed, r(1), c(1), b(1, 1), n, x(1, 1), n, &
    rcond, ferr(1), berr(1), work(1), rwork(1), info)
  if(info /= 0) then
    write(message(1), '(3a, i3)') 'In zlinsyssolve, LAPACK ', TOSTRING(ZLAPACK(gesvx)), ' returned info = ', info
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


! ---------------------------------------------------------
!> computes the singular value decomposition of a complex NxN matrix a(:,:)
subroutine zsingular_value_decomp(n, a, u, vt, sg_values)
  integer, intent(in)    :: n
  CMPLX,   intent(inout) :: a(n, n)  
  CMPLX,   intent(out)   :: u(n, n), vt(n, n)  
  FLOAT,   intent(out)   :: sg_values(n)

  interface
    subroutine ZLAPACK(gesvd) ( jobu, jobvt, m, n, a, lda, s, u, ldu, &
      vt, ldvt, work, lwork, rwork, info )
      character(1), intent(in)    :: jobu, jobvt
      integer,      intent(in)    :: m, n
      CMPLX,        intent(inout) :: a, u, vt ! a(lda,n), b(ldu,n), b(ldvt,n)
      CMPLX,        intent(out)   :: work     ! work(lwork)
      integer,      intent(in)    :: lda, ldu, ldvt, lwork
      integer,      intent(out)   :: info
      FLOAT,        intent(out)   :: s        ! s(min(m,n))
      FLOAT,        intent(inout) :: rwork    ! rwork(5*min(m,n))
    end subroutine ZLAPACK(gesvd)
  end interface

  integer :: m, info, lwork
  CMPLX, allocatable :: work(:)
  FLOAT, allocatable :: rwork(:)

  PUSH_SUB(zsingular_value_decomp)

  ! for now we treat only square matrices
  m = n 

  ! double minimum lwork size to increase performance (see manpage)
  lwork = 2*( 2*min(m, n) + max(m, n) )

  SAFE_ALLOCATE(work(1:lwork))
  SAFE_ALLOCATE(rwork(1:5*min(m, n)))

  call ZLAPACK(gesvd)( &
    'A', 'A', m, n, a(1, 1), m, sg_values(1), u(1, 1), m, vt(1, 1), n, work(1), lwork, rwork(1), info )

  if(info /= 0) then
    write(message(1), '(3a, i3)') 'In zsingular_value_decomp, LAPACK ', TOSTRING(LAPACK(gesvd)), ' returned info = ', info
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
  CMPLX,   intent(inout)        :: a(n, n)    ! a will be replaced by its inverse
  FLOAT,   intent(in), optional :: threshold

  CMPLX, allocatable :: u(:,:), vt(:,:)
  FLOAT, allocatable :: sg_values(:)
  CMPLX   :: tmp
  FLOAT   :: sg_inverse, threshold_
  integer :: j, k, l

  SAFE_ALLOCATE( u(1:n, 1:n))
  SAFE_ALLOCATE(vt(1:n, 1:n))
  SAFE_ALLOCATE(sg_values(1:n))

  PUSH_SUB(zsvd_inverse)

  call zsingular_value_decomp(n, a, u, vt, sg_values)

  threshold_ = CNST(1e-10); 
  if(present(threshold)) threshold_ = threshold

  ! build inverse
  do j = 1, n
    do k = 1, n
      tmp = M_ZERO
      do l = 1, n
        if (sg_values(l).lt.threshold_) then
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


! ---------------------------------------------------------
!> Calculate the inverse of a real upper triangular matrix (in
!! unpacked storage).
subroutine dinvert_upper_triangular(n, a)
  integer, intent(in)    :: n
  FLOAT,   intent(inout) :: a(n, n)

  integer :: info

  interface
    subroutine DLAPACK(trtri)(uplo, diag, n, a, lda, info)
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: diag
      integer,      intent(in)    :: n
      FLOAT,        intent(inout) :: a
      integer,      intent(in)    :: lda
      integer,      intent(out)   :: info
    end subroutine DLAPACK(trtri)
  end interface
  
  PUSH_SUB(dinvert_upper_triangular)

  call DLAPACK(trtri)('U', 'N', n, a(1, 1), n, info)

  if(info.ne.0) then
    write(message(1), '(3a,i5)') &
      'In dinvert_upper_triangular, LAPACK ', TOSTRING(DLAPACK(trtri)), ' returned error message ', info
    call messages_fatal(1)
  end if

  POP_SUB(dinvert_upper_triangular)
end subroutine dinvert_upper_triangular


! ---------------------------------------------------------
!> Calculate the inverse of a complex upper triangular matrix (in
!! unpacked storage).
subroutine zinvert_upper_triangular(n, a)
  integer, intent(in)    :: n
  CMPLX,   intent(inout) :: a(n, n)

  integer :: info

  interface
    subroutine ZLAPACK(trtri)(uplo, diag, n, a, lda, info)
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: diag
      integer,      intent(in)    :: n
      CMPLX,        intent(inout) :: a
      integer,      intent(in)    :: lda
      integer,      intent(out)   :: info
    end subroutine ZLAPACK(trtri)
  end interface
  
  PUSH_SUB(zinvert_upper_triangular)

  call ZLAPACK(trtri)('U', 'N', n, a(1, 1), n, info)

  if(info.ne.0) then
    write(message(1), '(3a,i5)') &
      'In zinvert_upper_triangular, LAPACK ', TOSTRING(ZLAPACK(trtri)), ' returned error message ', info
    call messages_fatal(1)
  end if

  POP_SUB(zinvert_upper_triangular)
end subroutine zinvert_upper_triangular

! ---------------------------------------------------------
!> Calculate the inverse of a real lower triangular matrix (in
!! unpacked storage).
subroutine dinvert_lower_triangular(n, a)
  integer, intent(in)    :: n
  FLOAT,   intent(inout) :: a(n, n)

  integer :: info

  interface
    subroutine DLAPACK(trtri)(uplo, diag, n, a, lda, info)
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: diag
      integer,      intent(in)    :: n
      FLOAT,        intent(inout) :: a
      integer,      intent(in)    :: lda
      integer,      intent(out)   :: info
    end subroutine DLAPACK(trtri)
  end interface
  
  PUSH_SUB(dinvert_lower_triangular)

  call DLAPACK(trtri)('L', 'N', n, a(1, 1), n, info)

  if(info.ne.0) then
    write(message(1), '(3a,i5)') &
      'In dinvert_lower_triangular, LAPACK ', TOSTRING(DLAPACK(trtri)), ' returned error message ', info
    call messages_fatal(1)
  end if

  POP_SUB(dinvert_lower_triangular)
end subroutine dinvert_lower_triangular

! ---------------------------------------------------------
!> Calculate the inverse of a complex lower triangular matrix (in
!! unpacked storage).
subroutine zinvert_lower_triangular(n, a)
  integer, intent(in)    :: n
  CMPLX,   intent(inout) :: a(n, n)

  integer :: info

  interface
    subroutine ZLAPACK(trtri)(uplo, diag, n, a, lda, info)
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: diag
      integer,      intent(in)    :: n
      CMPLX,        intent(inout) :: a
      integer,      intent(in)    :: lda
      integer,      intent(out)   :: info
    end subroutine ZLAPACK(trtri)
  end interface
  
  PUSH_SUB(zinvert_lower_triangular)

  call ZLAPACK(trtri)('L', 'N', n, a(1, 1), n, info)

  if(info.ne.0) then
    write(message(1), '(3a,i5)') &
      'In zinvert_lower_triangular, LAPACK ', TOSTRING(ZLAPACK(trtri)), ' returned error message ', info
    call messages_fatal(1)
  end if

  POP_SUB(zinvert_lower_triangular)
end subroutine zinvert_lower_triangular


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

!#include "global.h"

! This module is intended to contain "only mathematical" functions and        !
!	procedures.                                                           !

!--------------------------------------------------------------------------------!
! To improve performance only the first element of an array should be passed to  !
! the BLAS and LAPACK routines. The interfaces in this module enforce that. Real !
! rank and dimensions of the arrays are indicated in front of the corresponding  !
! scalar arguments.                                                              !
!--------------------------------------------------------------------------------!
module linalg
  use global

  implicit none
  
  interface

    !--- BLAS interfaces ---!
    
    ! scales a vector by a constant
    subroutine dscal(n, da, dx, incx)
      integer, intent(in)    :: n, incx
      real(8), intent(in)    :: da
      real(8), intent(inout) :: dx ! dx(n)
    end subroutine dscal
    
    ! scales a vector by a constant
    subroutine zscal(n, za, zx, incx)
      integer,    intent(in)    :: n, incx
      complex(8), intent(in)    :: za
      complex(8), intent(inout) :: zx ! zx(n)
    end subroutine zscal
    
    ! constant times a vector plus a vector
    subroutine daxpy(n, da, dx, incx, dy, incy)
      integer, intent(in)    :: n, incx, incy
      real(8), intent(in)    :: da, dx ! dx(n)
      real(8), intent(inout) :: dy     ! dy(n)
    end subroutine daxpy
    
    ! constant times a vector plus a vector
    subroutine zaxpy(n, za, zx, incx, zy, incy)
      integer,    intent(in)    :: n, incx, incy
      complex(8), intent(in)    :: za, zx ! zx(n)
      complex(8), intent(inout) :: zy     ! zy(n)
    end subroutine zaxpy
    
    ! forms the dot product of two vectors
    real(8) function ddot(n, dx, incx, dy, incy)
      integer, intent(in) :: n, incx, incy
      real(8), intent(in) :: dx, dy ! dx(n), dy(n)
    end function ddot
    
    complex(8) function zdotc(n, zx, incx, zy, incy)
      integer,    intent(in) :: n, incx, incy
      complex(8), intent(in) :: zx, zy ! zx(n), zy(n)
    end function zdotc
    
    ! returns the euclidean norm of a vector
    real(8) function dnrm2(n, dx, incx)
      integer, intent(in) :: n, incx
      real(8), intent(in) :: dx ! dx(n)
    end function dnrm2
    
    ! returns the euclidean norm of a vector
    real(8) function dznrm2(n, zx, incx)
      integer, intent(in)    :: n, incx
      complex(8), intent(in) :: zx ! zx(n)
    end function dznrm2
    
    ! matrix-matrix multiplication plus matrix
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      real(8),      intent(in)    :: alpha, beta
      real(8),      intent(in)    :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      real(8),      intent(in)    :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      real(8),      intent(inout) :: c ! c(ldc,n) 
    end subroutine dgemm
    
    ! matrix-matrix multiplication plus matrix
    subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      complex(8),   intent(in)    :: alpha, beta
      complex(8),   intent(in)    :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      complex(8),   intent(in)    :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      complex(8),   intent(inout) :: c ! c(ldc,n)
    end subroutine zgemm
    
    ! matrix-vector multiplication plus vector
    subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      complex(8),   intent(in)    :: alpha, beta
      complex(8),   intent(in)    :: a ! a(lda,n)
      complex(8),   intent(in)    :: x ! x(:)
      complex(8),   intent(inout) :: y ! y(:)
    end subroutine zgemv
    
    ! copies a vector, x, to a vector, y
    subroutine dcopy(n, dx, incx, dy, incy)
      integer, intent(in)  :: n, incx, incy
      real(8), intent(in)  :: dx ! dx(n)
      real(8), intent(out) :: dy ! dy(n)
    end subroutine dcopy

    ! copies a vector, x, to a vector, y
    subroutine zcopy(n, zx, incx, zy, incy)
      integer,    intent(in)  :: n, incx, incy
      complex(8), intent(in)  :: zx ! dz(n)
      complex(8), intent(out) :: zy ! zy(n)
    end subroutine zcopy
    
    
    !--- LAPACK interfaces ---!
    
    ! computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    subroutine dgetrf(m, n, a, lda, ipiv, info)
      integer, intent(in)    :: m, n, lda
      real(8), intent(inout) :: a    ! a(lda,n)
      integer, intent(out)   :: ipiv ! ipiv(min(m,n))
      integer, intent(out)   :: info
    end subroutine dgetrf
    
    ! computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    subroutine zgetrf(m, n, a, lda, ipiv, info)
      integer,    intent(in)    :: m, n, lda
      complex(8), intent(inout) :: a    ! a(lda,n)
      integer,    intent(out)   :: ipiv ! ipiv(min(m,n))
      integer,    intent(out)   :: info
    end subroutine zgetrf
    
    ! computes all eigenvalues and, optionally, eigenvectors of a real symmetric matrix A.
    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, lda, lwork
      real(8),      intent(inout) :: a       ! a(lda,n)
      real(8),      intent(out)   :: w, work ! w(n), work(lwork)
      integer,      intent(out)   :: info 
    end subroutine dsyev
    
    ! computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.
    subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, lda, lwork
      complex(8),   intent(inout) :: a        ! a(lda,n)
      real(8),      intent(out)   :: w, rwork ! w(n), rwork(max(1,3*n-2))
      complex(8),   intent(out)   :: work     ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine zheev
    
    ! computes all the eigenvalues, and optionally, the eigenvectors of a real
    ! generalized symmetric-definite eigenproblem, of the form  A*x=(lambda)*B*x,
    ! A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
    ! Here A and B are assumed to be symmetric and B is also positive definite.
    subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: itype, n, lda, ldb, lwork
      real(8),      intent(inout) :: a, b    ! a(lda,n), b(ldb,n)
      real(8),      intent(out)   :: w, work ! w(n), work(lwork)
      integer,      intent(out)   :: info
    end subroutine dsygv
    
    ! computes all the eigenvalues, and optionally, the eigenvectors of a complex 
    ! generalized Hermitian-definite eigenproblem, of the form  A*x=(lambda)*B*x,  
    ! A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
    ! Here A and B are assumed to be Hermitian and B is also positive definite.
    subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, itype, lda, ldb, lwork
      complex(8),   intent(inout) :: a, b     ! a(lda,n), b(ldb,n)
      real(8),      intent(out)   :: w, rwork ! w(n), rwork(max(1,3*n-2))
      complex(8),   intent(out)   :: work     ! work(lwork)
      integer,      intent(out)   :: info 
    end subroutine zhegv
    
  end interface

contains
  ! ddet and zdet return the determinant of a square matrix a of dimensions (n,n)
  ! (ddet for real and zdet for complex numbers)
  real(r8) function ddet(a, n)
    integer, intent(in) :: n
    real(r8) :: a(n, n)
    
    integer :: i, info, ipiv(n)
    
    call dgetrf(n, n, a(1,1), n, ipiv(1), info)
    if(info < 0) then
      write(message(1),'(a,i4)') 'In det, LAPACK dgetrf returned error code ',info
      call write_fatal(1)
    endif
    ddet = M_ONE
    do i = 1, n
      ddet = ddet*a(i, i)
      if(ipiv(i).ne.i) then
        ipiv(ipiv(i)) = ipiv(i)
        ddet = -ddet
      endif
    enddo
    
  end function ddet

  complex(r8) function zdet(a, n)
    integer, intent(in) :: n
    complex(r8) :: a(n, n)
    
    integer :: info, i, ipiv(n)
    
    call zgetrf(n, n, a(1,1), n, ipiv(1), info)
    if(info < 0) then
      write(message(1),'(a,i4)') 'In zdet, LAPACK zgetrf returned error code ',info
      call write_fatal(1)
    endif
    zdet = M_z1
    do i = 1, n
      zdet = zdet*a(i, i)
      if(ipiv(i).ne.i) then
        ipiv(ipiv(i)) = ipiv(i)
        zdet = -zdet
      endif
    enddo
    
  end function zdet
  
!!! These routines diagonalize a matrix using LAPACK
  subroutine diagonalise(dim, a, b, e)
    integer  :: dim
    real(r8) :: a(dim, dim), b(dim, dim)
    real(r8) :: e(dim)
    
    integer :: info, lwork
    real(r8), allocatable :: work(:)
    
    lwork = 6*dim
    allocate(work(lwork))
    b = a
    call dsyev('V', 'U', dim, b(1,1), dim, e(1), work(1), lwork, info)
    if(info.ne.0) then
      write(message(1),'(a,i5)') 'dsyev returned error message ', info
      call write_fatal(1)
    endif
    deallocate(work)
  end subroutine diagonalise
  
  subroutine ziagonalise(dim, a, b, e)
    integer     :: dim
    complex(r8) :: a(dim, dim), b(dim, dim)
    real(r8)    :: e(dim)
    
    integer :: info, lwork
    complex(r8), allocatable :: work(:)
    real(r8), allocatable :: rwork(:)
    
    lwork = 6*dim
    allocate(work(lwork), rwork(max(1,3*dim-2)))
    b = a
    call zheev('V','U', dim, b(1,1), dim, e(1), work(1), lwork, rwork(1), info)
    if(info.ne.0) then
      write(message(1),'(a,i5)') 'dsyev returned error message ', info
      call write_fatal(1)
    endif
    deallocate(work, rwork)
  end subroutine ziagonalise

end module linalg

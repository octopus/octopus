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

#include "global.h"

! -----------------------------------------------------------------------
!> This module contains interfaces for BLAS routines
!! You should not use these routines directly. Please use the lalg_XXXX
! -----------------------------------------------------------------------
module blas_m

  implicit none

  public ! only interfaces in this module

  ! ---------------------------------------------------------------------
  ! BLAS level I
  ! ---------------------------------------------------------------------

  !> ----------------- swap ------------------
  !! Interchanges two vectors.
  interface blas_swap
    subroutine sswap(n, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      real(4),    intent(inout) :: dx, dy !< dx(n), dy(n)
    end subroutine sswap

    subroutine dswap(n, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      real(8),    intent(inout) :: dx, dy !< dx(n), dy(n)
    end subroutine dswap

    subroutine cswap(n, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      complex(4), intent(inout) :: dx, dy !< dx(n), dy(n)
    end subroutine cswap

    subroutine zswap(n, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      complex(8), intent(inout) :: dx, dy !< dx(n), dy(n)
    end subroutine zswap
  end interface blas_swap

  !> ----------------- scal ------------------
  !! Scales a vector by a constant.
  interface blas_scal
    subroutine sscal(n, da, dx, incx)
      integer,    intent(in)    :: n, incx
      real(4),    intent(in)    :: da
      real(4),    intent(inout) :: dx !< dx(n)
    end subroutine sscal

    subroutine dscal(n, da, dx, incx)
      integer,    intent(in)    :: n, incx
      real(8),    intent(in)    :: da
      real(8),    intent(inout) :: dx !< dx(n)
    end subroutine dscal

    subroutine cscal(n, da, dx, incx)
      integer,    intent(in)    :: n, incx
      complex(4), intent(in)    :: da
      complex(4), intent(inout) :: dx !< dx(n)
    end subroutine cscal

    subroutine zscal(n, da, dx, incx)
      integer,    intent(in)    :: n, incx
      complex(8), intent(in)    :: da
      complex(8), intent(inout) :: dx !< dx(n)
    end subroutine zscal

    subroutine dazscal(n, da, dx)
      integer,    intent(in)    :: n
      real(8),    intent(in)    :: da
      complex(8), intent(inout) :: dx !< dx(n)
    end subroutine dazscal

    subroutine sazscal(n, da, dx)
      integer,    intent(in)    :: n
      real(4),    intent(in)    :: da
      complex(4), intent(inout) :: dx !< dx(n)
    end subroutine sazscal
  end interface blas_scal

  !> ----------------- axpy ------------------
  !! Constant times a vector plus a vector.
  interface blas_axpy
    subroutine saxpy (n, da, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      real(4),    intent(in)    :: da, dx !< dx(n)
      real(4),    intent(inout) :: dy     !< dy(n)
    end subroutine saxpy

    subroutine daxpy (n, da, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      real(8),    intent(in)    :: da, dx !< dx(n)
      real(8),    intent(inout) :: dy     !< dy(n)
    end subroutine daxpy

    subroutine caxpy (n, da, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      complex(4), intent(in)    :: da, dx !< dx(n)
      complex(4), intent(inout) :: dy     !< dy(n)
    end subroutine caxpy

    subroutine zaxpy (n, da, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      complex(8), intent(in)    :: da, dx !< dx(n)
      complex(8), intent(inout) :: dy     !< dy(n)
    end subroutine zaxpy

    subroutine dazaxpy (n, da, dx, dy)
      integer,    intent(in)    :: n
      real(8),    intent(in)    :: da
      complex(8), intent(in)    :: dx     !< dx(n)
      complex(8), intent(inout) :: dy     !< dy(n)
    end subroutine dazaxpy

    subroutine sazaxpy (n, da, dx, dy)
      integer,    intent(in)    :: n
      real(4),    intent(in)    :: da
      complex(4), intent(in)    :: dx     !< dx(n)
      complex(4), intent(inout) :: dy     !< dy(n)
    end subroutine sazaxpy
  end interface blas_axpy

  !> ----------------- copy ------------------
  !! Copies a vector, x, to a vector, y.
  interface blas_copy
    subroutine scopy(n, dx, incx, dy, incy)
      integer,    intent(in)  :: n, incx, incy
      real(4),    intent(in)  :: dx !< dx(n)
      real(4),    intent(out) :: dy !< dy(n)
    end subroutine scopy

    subroutine dcopy(n, dx, incx, dy, incy)
      integer,    intent(in)  :: n, incx, incy
      real(8),    intent(in)  :: dx !< dx(n)
      real(8),    intent(out) :: dy !< dy(n)
    end subroutine dcopy

    subroutine ccopy(n, dx, incx, dy, incy)
      integer,    intent(in)  :: n, incx, incy
      complex(4), intent(in)  :: dx !< dx(n)
      complex(4), intent(out) :: dy !< dy(n)
    end subroutine ccopy

    subroutine zcopy(n, dx, incx, dy, incy)
      integer,    intent(in)  :: n, incx, incy
      complex(8), intent(in)  :: dx !< dx(n)
      complex(8), intent(out) :: dy !< dy(n)
    end subroutine zcopy
  end interface blas_copy

  !> ----------------- dot  ------------------
  !! Forms the dot product of two vectors.
  interface blas_dot
    real(4) function sdot(n, dx, incx, dy, incy)
      integer,    intent(in) :: n, incx, incy
      real(4),    intent(in) :: dx, dy !< dx(n), dy(n)
    end function sdot

    real(8) function ddot(n, dx, incx, dy, incy)
      integer,    intent(in) :: n, incx, incy
      real(8),    intent(in) :: dx, dy !< dx(n), dy(n)
    end function ddot

    complex(4) function cdotc(n, dx, incx, dy, incy)
      integer,    intent(in) :: n, incx, incy
      complex(4), intent(in) :: dx, dy !< dx(n), dy(n)
    end function cdotc

    complex(8) function zdotc(n, dx, incx, dy, incy)
      integer,    intent(in) :: n, incx, incy
      complex(8), intent(in) :: dx, dy !< dx(n), dy(n)
    end function zdotc
  end interface blas_dot

  interface blas_dotu
    complex(4) function cdotu(n, dx, incx, dy, incy)
      integer,    intent(in) :: n, incx, incy
      complex(4), intent(in) :: dx, dy !< dx(n), dy(n)
    end function cdotu

    complex(8) function zdotu(n, dx, incx, dy, incy)
      integer,    intent(in) :: n, incx, incy
      complex(8), intent(in) :: dx, dy !< dx(n), dy(n)
    end function zdotu
  end interface blas_dotu

  !> ----------------- nrm2 ------------------
  !! Returns the euclidean norm of a vector via the function
  !! name, so that
  !! \f[
  !! SNRM2 := sqrt( x'*x )
  !! \f]
  interface blas_nrm2
    real(4) function snrm2(n, dx, incx)
      integer,    intent(in) :: n, incx
      real(4),    intent(in) :: dx !< dx(n)
    end function snrm2

    real(8) function dnrm2(n, dx, incx)
      integer,    intent(in) :: n, incx
      real(8),    intent(in) :: dx !< dx(n)
    end function dnrm2

    real(4) function scnrm2(n, dx, incx)
      integer,    intent(in) :: n, incx
      complex(4), intent(in) :: dx !< dx(n)
    end function scnrm2

    real(8) function dznrm2(n, dx, incx)
      integer,    intent(in) :: n, incx
      complex(8), intent(in) :: dx !< dx(n)
    end function dznrm2
  end interface blas_nrm2


  ! ------------------------------------------------------------------
  ! BLAS level II
  ! ------------------------------------------------------------------

  !> ----------------- symv ------------------
  !! performs the matrix-vector  operation
  !!
  !! \f[
  !!     y := \alpha A x + \beta y
  !! \f]
  !!
  !!  where \f$\alpha\f$ and \f$\beta\f$ are scalars, x and y are n
  !!  element vectors and A is an \f$n\times n\f$ symmetric matrix.
  interface blas_symv
    subroutine ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda, incx, incy
      real(4),      intent(in)    :: alpha, beta
      real(4),      intent(in)    :: a !< a(lda,n)
      real(4),      intent(in)    :: x !< x(:)
      real(4),      intent(inout) :: y !< y(:)
    end subroutine ssymv

    subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda, incx, incy
      real(8),      intent(in)    :: alpha, beta
      real(8),      intent(in)    :: a !< a(lda,n)
      real(8),      intent(in)    :: x !< x(:)
      real(8),      intent(inout) :: y !< y(:)
    end subroutine dsymv

    subroutine csymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda, incx, incy
      complex(4),   intent(in)    :: alpha, beta
      complex(4),   intent(in)    :: a !< a(lda,n)
      complex(4),   intent(in)    :: x !< x(:)
      complex(4),   intent(inout) :: y !< y(:)
    end subroutine csymv

    subroutine zsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda, incx, incy
      complex(8),   intent(in)    :: alpha, beta
      complex(8),   intent(in)    :: a !< a(lda,n)
      complex(8),   intent(in)    :: x !< x(:)
      complex(8),   intent(inout) :: y !< y(:)
    end subroutine zsymv
  end interface blas_symv
  
  !> ----------------- gemv ------------------
  !! SGEMV  performs one of the matrix-vector operations
  !!
  !! \f[
  !!     y := \alpha A x + \beta y,   
  !! \f]
  !! or 
  !! \f[
  !!     y := \alpha A^Tx + \beta y
  !! \f]
  !!  
  !!  where \f$\alpha\f$ and \f$\beta\f$ are scalars, x and y are
  !!  vectors and A is an \f$m\times n\f$ matrix.
  interface blas_gemv
    subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      real(4),      intent(in)    :: alpha, beta
      real(4),      intent(in)    :: a !< a(lda,n)
      real(4),      intent(in)    :: x !< x(:)
      real(4),      intent(inout) :: y !< y(:)
    end subroutine sgemv

    subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      real(8),      intent(in)    :: alpha, beta
      real(8),      intent(in)    :: a !< a(lda,n)
      real(8),      intent(in)    :: x !< x(:)
      real(8),      intent(inout) :: y !< y(:)
    end subroutine dgemv

    subroutine cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      complex(4),   intent(in)    :: alpha, beta
      complex(4),   intent(in)    :: a !< a(lda,n)
      complex(4),   intent(in)    :: x !< x(:)
      complex(4),   intent(inout) :: y !< y(:)
    end subroutine cgemv

    subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      complex(8),   intent(in)    :: alpha, beta
      complex(8),   intent(in)    :: a !< a(lda,n)
      complex(8),   intent(in)    :: x !< x(:)
      complex(8),   intent(inout) :: y !< y(:)
    end subroutine zgemv
  end interface blas_gemv


  ! ------------------------------------------------------------------
  ! BLAS level III
  ! ------------------------------------------------------------------

  !> ----------------- gemm ------------------
  !! performs one of the matrix-matrix operations
  !!
  !! \f[
  !!     C := \alpha op( A ) op( B ) + \beta C,
  !! \f]
  !!
  !!  where  op(X) is one of
  !!
  !! \f[
  !!     op( X ) = X   \mbox{ or }   op( X ) = X^T,
  !! \f]
  !!
  !!  alpha and beta are scalars, and A, B and C are matrices, with op(A)
  !!  an m by k matrix,  op(B)  a  k by n matrix and  C an m by n matrix.
  interface blas_gemm
    subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      real(4),      intent(in)    :: alpha, beta
      real(4),      intent(in)    :: a !< a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      real(4),      intent(in)    :: b !< b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      real(4),      intent(inout) :: c !< c(ldc,n)
    end subroutine sgemm

    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      real(8),      intent(in)    :: alpha, beta
      real(8),      intent(in)    :: a !< a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      real(8),      intent(in)    :: b !< b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      real(8),      intent(inout) :: c !< c(ldc,n)
    end subroutine dgemm

    subroutine cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      complex(4),   intent(in)    :: alpha, beta
      complex(4),   intent(in)    :: a !< a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      complex(4),   intent(in)    :: b !< b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      complex(4),   intent(inout) :: c !< c(ldc,n)
    end subroutine cgemm

    subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      complex(8),   intent(in)    :: alpha, beta
      complex(8),   intent(in)    :: a !< a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      complex(8),   intent(in)    :: b !< b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      complex(8),   intent(inout) :: c !< c(ldc,n)
    end subroutine zgemm
  end interface blas_gemm

  !> ----------------- trmm ------------------ 
  !! Performs one of the matrix-matrix operations
  !!
  !! \f[
  !!     B := \alpha op( A )B,  \mbox{ or }  B := \alpha B op( A ),
  !! \f]
  !!
  !!  where \f$\alpha\f$ is a scalar, B is an \f$m\times n\f$matrix, A
  !!  is a unit, or non-unit, upper or lower triangular matrix and op(
  !!  A ) is one of
  !!
  !! \f[
  !!     op( A ) = A   \mbox{ or }   op( A ) = A^T.
  !! \f]
  interface blas_trmm
    subroutine strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      character(1), intent(in)    :: side, uplo, transa, diag
      integer,      intent(in)    :: m, n, lda, ldb
      real(4),      intent(in)    :: a, alpha
      real(4),      intent(inout) :: b
    end subroutine strmm

    subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      character(1), intent(in)    :: side, uplo, transa, diag
      integer,      intent(in)    :: m, n, lda, ldb
      real(8),      intent(in)    :: a, alpha
      real(8),      intent(inout) :: b
    end subroutine dtrmm

    subroutine ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      character(1), intent(in)    :: side, uplo, transa, diag
      integer,      intent(in)    :: m, n, lda, ldb
      complex(4),   intent(in)    :: a, alpha
      complex(4),   intent(inout) :: b
    end subroutine ctrmm

    subroutine ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      character(1), intent(in)    :: side, uplo, transa, diag
      integer,      intent(in)    :: m, n, lda, ldb
      complex(8),   intent(in)    :: a, alpha
      complex(8),   intent(inout) :: b
    end subroutine ztrmm
  end interface blas_trmm

  !> ----------------- symm, hemm ------------------
  !! performs one of the matrix-matrix operations
  !!
  !! \f[
  !!     C := \alpha A B + \beta C,
  !! \f]
  !!
  !!  or
  !!
  !! \f[
  !!     C := \alpha  B A + \beta C
  !! \f]
  !!
  !!  where \f$\alpha\f$ and \f$\beta\f$ are scalars, A is a symmetric
  !!  matrix and B and C are \f$m\times n\f$ matrices.
  interface blas_symm
    subroutine ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: side, uplo
      integer,      intent(in)    :: m, n, lda, ldb, ldc
      real(4),      intent(in)    :: alpha, beta, a, b
      real(4),      intent(inout) :: c
    end subroutine ssymm

    subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: side, uplo
      integer,      intent(in)    :: m, n, lda, ldb, ldc
      real(8),      intent(in)    :: alpha, beta, a, b
      real(8),      intent(inout) :: c
    end subroutine dsymm

    subroutine csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: side, uplo
      integer,      intent(in)    :: m, n, lda, ldb, ldc
      complex(4),   intent(in)    :: alpha, beta, a, b
      complex(4),   intent(inout) :: c
    end subroutine csymm

    subroutine zsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: side, uplo
      integer,      intent(in)    :: m, n, lda, ldb, ldc
      complex(8),   intent(in)    :: alpha, beta, a, b
      complex(8),   intent(inout) :: c
    end subroutine zsymm
  end interface blas_symm

  !> ----------------- syrk, herk ------------------
  !! performs one of the symmetric rank k operations
  !!
  !! \f[
  !!     C := \alpha A A^T + \beta C,
  !!
  !! \f]
  !! or 
  !! \f[
  !!     C := alpha A^T A + \beta*C
  !! \f]
  !!
  !!  where \f$\alpha\f$ and \f$\beta\f$ are scalars, C is an
  !!  \f$n\times n\f$ symmetric matrix and A is an \f$n\times k\f$
  !!  matrix in the first case and a \f$k\times n\f$ matrix in the
  !!  second case.
  interface blas_herk
    subroutine ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      implicit none

      character(1), intent(in)    :: uplo, trans
      integer,      intent(in)    :: n, k, lda, ldc
      real(4),      intent(in)    :: alpha, beta, a
      real(4),      intent(inout) :: c
    end subroutine ssyrk

    subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      implicit none

      character(1), intent(in)    :: uplo, trans
      integer,      intent(in)    :: n, k, lda, ldc
      real(8),      intent(in)    :: alpha, beta, a
      real(8),      intent(inout) :: c
    end subroutine dsyrk

    subroutine cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      implicit none

      character(1), intent(in)    :: uplo, trans
      integer,      intent(in)    :: n, k, lda, ldc
      real(4),      intent(in)    :: alpha, beta
      complex(4),   intent(in)    :: a
      complex(4),   intent(inout) :: c
    end subroutine cherk

    subroutine zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      implicit none

      character(1), intent(in)    :: uplo, trans
      integer,      intent(in)    :: n, k, lda, ldc
      real(8),      intent(in)    :: alpha, beta
      complex(8),   intent(in)    :: a
      complex(8),   intent(inout) :: c
    end subroutine zherk
  end interface blas_herk

  !> -----------------------trsm-------------------------
  !! Solves one of the matrix equations
  !!
  !! \f[
  !!     op( A )X = \alpha B, 
  !! \f]
  !! or 
  !! \f[
  !! X op( A ) = \alpha B,
  !! \f]
  !!
  !!  where \f$\alpha\f$ is a scalar, X and B are \f$m\times n\f$
  !!  matrices, A is a unit, or non-unit, upper or lower triangular
  !!  matrix and op(A) is one of
  !!
  !! \f[
  !!     op( A ) = A   \mbox{ or }   op( A ) = A^T
  !! \f]
  !!
  !!  The matrix X is overwritten on B.
  interface blas_trsm
    subroutine strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      character(1), intent(in)    :: side
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: transa
      character(1), intent(in)    :: diag
      integer,      intent(in)    :: m
      integer,      intent(in)    :: n
      real(4),      intent(in)    :: alpha
      real(4),      intent(in)    :: a
      integer,      intent(in)    :: lda
      real(4),      intent(inout) :: b
      integer,      intent(in)    :: ldb
    end subroutine strsm

    subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      character(1), intent(in)    :: side
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: transa
      character(1), intent(in)    :: diag
      integer,      intent(in)    :: m
      integer,      intent(in)    :: n
      real(8),      intent(in)    :: alpha
      real(8),      intent(in)    :: a
      integer,      intent(in)    :: lda
      real(8),      intent(inout) :: b
      integer,      intent(in)    :: ldb
    end subroutine dtrsm

    subroutine ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      character(1), intent(in)    :: side
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: transa
      character(1), intent(in)    :: diag
      integer,      intent(in)    :: m
      integer,      intent(in)    :: n
      complex(4),   intent(in)    :: alpha
      complex(4),   intent(in)    :: a
      integer,      intent(in)    :: lda
      complex(4),   intent(inout) :: b
      integer,      intent(in)    :: ldb
    end subroutine ctrsm

    subroutine ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      character(1), intent(in)    :: side
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: transa
      character(1), intent(in)    :: diag
      integer,      intent(in)    :: m
      integer,      intent(in)    :: n
      complex(8),   intent(in)    :: alpha
      complex(8),   intent(in)    :: a
      integer,      intent(in)    :: lda
      complex(8),   intent(inout) :: b
      integer,      intent(in)    :: ldb
    end subroutine ztrsm
  end interface blas_trsm

end module blas_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

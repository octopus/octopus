!! Copyright (C) 2002-2004 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

! -----------------------------------------------------------------------
! This module contains interfaces for BLAS routines
! You should not use these routines directly. Please use the lalg_XXXX
! -----------------------------------------------------------------------
module blas
  implicit none

  ! ---------------------------------------------------------------------
  ! Level 1
  ! ---------------------------------------------------------------------

  ! ----------------- swap ------------------
  interface blas_swap
    subroutine sswap(n, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      real(4),    intent(inout) :: dx, dy ! dx(n), dy(n)
    end subroutine sswap

    subroutine dswap(n, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      real(8),    intent(inout) :: dx, dy ! dx(n), dy(n)
    end subroutine dswap
    
    subroutine cswap(n, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      complex(4), intent(inout) :: dx, dy ! dx(n), dy(n)
    end subroutine cswap

    subroutine zswap(n, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      complex(8), intent(inout) :: dx, dy ! dx(n), dy(n)
    end subroutine zswap
  end interface

  ! ----------------- scal ------------------
  interface blas_scal
    subroutine sscal(n, da, dx, incx)
      integer,    intent(in)    :: n, incx
      real(4),    intent(in)    :: da
      real(4),    intent(inout) :: dx ! dx(n)
    end subroutine sscal
    
    subroutine dscal(n, da, dx, incx)
      integer,    intent(in)    :: n, incx
      real(8),    intent(in)    :: da
      real(8),    intent(inout) :: dx ! dx(n)
    end subroutine dscal

    subroutine cscal(n, da, dx, incx)
      integer,    intent(in)    :: n, incx
      complex(4), intent(in)    :: da
      complex(4), intent(inout) :: dx ! dx(n)
    end subroutine cscal

    subroutine zscal(n, da, dx, incx)
      integer,    intent(in)    :: n, incx
      complex(8), intent(in)    :: da
      complex(8), intent(inout) :: dx ! dx(n)
    end subroutine zscal
  end interface

  ! ----------------- axpy ------------------
  interface blas_axpy
    subroutine saxpy (n, da, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      real(4),    intent(in)    :: da, dx ! dx(n)
      real(4),    intent(inout) :: dy     ! dy(n)
    end subroutine saxpy

    subroutine daxpy (n, da, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      real(8),    intent(in)    :: da, dx ! dx(n)
      real(8),    intent(inout) :: dy     ! dy(n)
    end subroutine daxpy

    subroutine caxpy (n, da, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      complex(4), intent(in)    :: da, dx ! dx(n)
      complex(4), intent(inout) :: dy     ! dy(n)
    end subroutine caxpy

    subroutine zaxpy (n, da, dx, incx, dy, incy)
      integer,    intent(in)    :: n, incx, incy
      complex(8), intent(in)    :: da, dx ! dx(n)
      complex(8), intent(inout) :: dy     ! dy(n)
    end subroutine zaxpy
  end interface

  ! ----------------- copy ------------------
  interface blas_copy
    subroutine scopy(n, dx, incx, dy, incy)
      integer,    intent(in)  :: n, incx, incy
      real(4),    intent(in)  :: dx ! dx(n)
      real(4),    intent(out) :: dy ! dy(n)
    end subroutine scopy

    subroutine dcopy(n, dx, incx, dy, incy)
      integer,    intent(in)  :: n, incx, incy
      real(8),    intent(in)  :: dx ! dx(n)
      real(8),    intent(out) :: dy ! dy(n)
    end subroutine dcopy

    subroutine ccopy(n, dx, incx, dy, incy)
      integer,    intent(in)  :: n, incx, incy
      complex(4), intent(in)  :: dx ! dx(n)
      complex(4), intent(out) :: dy ! dy(n)
    end subroutine ccopy

    subroutine zcopy(n, dx, incx, dy, incy)
      integer,    intent(in)  :: n, incx, incy
      complex(8), intent(in)  :: dx ! dx(n)
      complex(8), intent(out) :: dy ! dy(n)
    end subroutine zcopy
  end interface

  ! ----------------- dot  ------------------
  interface blas_dot
    real(4) function sdot(n, dx, incx, dy, incy)
      integer,    intent(in) :: n, incx, incy
      real(4),    intent(in) :: dx, dy ! dx(n), dy(n)
    end function sdot
    
    real(8) function ddot(n, dx, incx, dy, incy)
      integer,    intent(in) :: n, incx, incy
      real(8),    intent(in) :: dx, dy ! dx(n), dy(n)
    end function ddot

    complex(4) function cdotc(n, dx, incx, dy, incy)
      integer,    intent(in) :: n, incx, incy
      complex(4), intent(in) :: dx, dy ! dx(n), dy(n)
    end function cdotc

    complex(8) function zdotc(n, dx, incx, dy, incy)
      integer,    intent(in) :: n, incx, incy
      complex(8), intent(in) :: dx, dy ! dx(n), dy(n)
    end function zdotc
  end interface

  ! ----------------- nrm2 ------------------
  interface blas_nrm2
    real(4) function snrm2(n, dx, incx)
      integer,    intent(in) :: n, incx
      real(4),    intent(in) :: dx ! dx(n)
    end function snrm2
    
    real(8) function dnrm2(n, dx, incx)
      integer,    intent(in) :: n, incx
      real(8),    intent(in) :: dx ! dx(n)
    end function dnrm2

    real(4) function scnrm2(n, dx, incx)
      integer,    intent(in) :: n, incx
      complex(4), intent(in) :: dx ! dx(n)
    end function scnrm2

    real(8) function dznrm2(n, dx, incx)
      integer,    intent(in) :: n, incx
      complex(8), intent(in) :: dx ! dx(n)
    end function dznrm2
  end interface

  ! ------------------------------------------------------------------
  ! BLAS level II
  ! ------------------------------------------------------------------

  ! ----------------- gemm ------------------
  interface blas_gemm
    subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      real(4),      intent(in)    :: alpha, beta
      real(4),      intent(in)    :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      real(4),      intent(in)    :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      real(4),      intent(inout) :: c ! c(ldc,n) 
    end subroutine sgemm

    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      real(8),      intent(in)    :: alpha, beta
      real(8),      intent(in)    :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      real(8),      intent(in)    :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      real(8),      intent(inout) :: c ! c(ldc,n) 
    end subroutine dgemm

    subroutine cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      complex(4),   intent(in)    :: alpha, beta
      complex(4),   intent(in)    :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      complex(4),   intent(in)    :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      complex(4),   intent(inout) :: c ! c(ldc,n) 
    end subroutine cgemm

    subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      complex(8),   intent(in)    :: alpha, beta
      complex(8),   intent(in)    :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      complex(8),   intent(in)    :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      complex(8),   intent(inout) :: c ! c(ldc,n) 
    end subroutine zgemm
  end interface

  ! ----------------- gemv ------------------
  interface blas_gemv
    subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      real(4),      intent(in)    :: alpha, beta
      real(4),      intent(in)    :: a ! a(lda,n)
      real(4),      intent(in)    :: x ! x(:)
      real(4),      intent(inout) :: y ! y(:)
    end subroutine sgemv

    subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      real(8),      intent(in)    :: alpha, beta
      real(8),      intent(in)    :: a ! a(lda,n)
      real(8),      intent(in)    :: x ! x(:)
      real(8),      intent(inout) :: y ! y(:)
    end subroutine dgemv

    subroutine cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      complex(4),   intent(in)    :: alpha, beta
      complex(4),   intent(in)    :: a ! a(lda,n)
      complex(4),   intent(in)    :: x ! x(:)
      complex(4),   intent(inout) :: y ! y(:)
    end subroutine cgemv

    subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      complex(8),   intent(in)    :: alpha, beta
      complex(8),   intent(in)    :: a ! a(lda,n)
      complex(8),   intent(in)    :: x ! x(:)
      complex(8),   intent(inout) :: y ! y(:)
    end subroutine zgemv
  end interface

end module blas

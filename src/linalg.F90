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

#include "global.h"

#if defined(SINGLE_PRECISION)
#  define DBLAS(x) s ## x
#  define ZBLAS(x) c ## x
#  define ZNRM2    dcnrm2
#else
#  define DBLAS(x) d ## x
#  define ZBLAS(x) z ## x
#  define ZNRM2    dznrm2
#endif

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
  
  !--- BLAS interfaces ---!
  interface la_scal
    ! scales a vector by a constant
    subroutine DBLAS(scal) (n, da, dx, incx)
      integer, intent(in)    :: n, incx
      FLOAT,   intent(in)    :: da
      FLOAT,   intent(inout) :: dx ! dx(n)
    end subroutine DBLAS(scal)
    
    subroutine ZBLAS(scal) (n, za, zx, incx)
      integer, intent(in)    :: n, incx
      CMPLX,   intent(in)    :: za
      CMPLX,   intent(inout) :: zx ! zx(n)
    end subroutine ZBLAS(scal)
  end interface
    
  ! constant times a vector plus a vector
  interface la_axpy
    subroutine DBLAS(axpy) (n, da, dx, incx, dy, incy)
      integer, intent(in)    :: n, incx, incy
      FLOAT,   intent(in)    :: da, dx ! dx(n)
      FLOAT,   intent(inout) :: dy     ! dy(n)
    end subroutine DBLAS(axpy)
    
    subroutine ZBLAS(axpy) (n, za, zx, incx, zy, incy)
      integer, intent(in)    :: n, incx, incy
      CMPLX,   intent(in)    :: za, zx ! zx(n)
      CMPLX,   intent(inout) :: zy     ! zy(n)
    end subroutine ZBLAS(axpy)
  end interface

  ! forms the dot product of two vectors
  interface la_dot
    FLOAT function DBLAS(dot) (n, dx, incx, dy, incy)
      integer, intent(in) :: n, incx, incy
      FLOAT,   intent(in) :: dx, dy ! dx(n), dy(n)
    end function DBLAS(dot)
    
    CMPLX function ZBLAS(dotc) (n, zx, incx, zy, incy)
      integer, intent(in) :: n, incx, incy
      CMPLX,   intent(in) :: zx, zy ! zx(n), zy(n)
    end function ZBLAS(dotc)
  end interface

  ! returns the euclidean norm of a vector
  interface la_nrm2
    FLOAT function DBLAS(nrm2) (n, dx, incx)
      integer, intent(in) :: n, incx
      FLOAT,   intent(in) :: dx ! dx(n)
    end function DBLAS(nrm2)
    
    FLOAT function ZNRM2 (n, zx, incx)
      integer, intent(in)    :: n, incx
      CMPLX,   intent(in) :: zx ! zx(n)
    end function ZNRM2
  end interface
    
  ! matrix-matrix multiplication plus matrix
  interface la_gemm
    subroutine DBLAS(gemm) (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      FLOAT,        intent(in)    :: alpha, beta
      FLOAT,        intent(in)    :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      FLOAT,        intent(in)    :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      FLOAT,        intent(inout) :: c ! c(ldc,n) 
    end subroutine DBLAS(gemm)
    
    subroutine ZBLAS(gemm) (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      CMPLX,        intent(in)    :: alpha, beta
      CMPLX,        intent(in)    :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      CMPLX,        intent(in)    :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      CMPLX,        intent(inout) :: c ! c(ldc,n)
    end subroutine ZBLAS(gemm)
  end interface
    
  ! matrix-vector multiplication plus vector
  interface la_gemv
    subroutine DBLAS(gemv) (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      FLOAT,        intent(in)    :: alpha, beta
      FLOAT,        intent(in)    :: a ! a(lda,n)
      FLOAT,        intent(in)    :: x ! x(:)
      FLOAT,        intent(inout) :: y ! y(:)
    end subroutine DBLAS(gemv)

    subroutine ZBLAS(gemv) (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      CMPLX,        intent(in)    :: alpha, beta
      CMPLX,        intent(in)    :: a ! a(lda,n)
      CMPLX,        intent(in)    :: x ! x(:)
      CMPLX,        intent(inout) :: y ! y(:)
    end subroutine ZBLAS(gemv)
  end interface
    
  ! copies a vector, x, to a vector, y
  interface la_copy
    subroutine DBLAS(copy) (n, dx, incx, dy, incy)
      integer, intent(in)  :: n, incx, incy
      FLOAT,   intent(in)  :: dx ! dx(n)
      FLOAT,   intent(out) :: dy ! dy(n)
    end subroutine DBLAS(copy)

    subroutine ZBLAS(copy) (n, zx, incx, zy, incy)
      integer,    intent(in)  :: n, incx, incy
      CMPLX,      intent(in)  :: zx ! dz(n)
      CMPLX,      intent(out) :: zy ! zy(n)
    end subroutine ZBLAS(copy)
  end interface

contains

#ifdef HAVE_LAPACK
#include "linalg_lapack.F90"
#else
#include "linal_gsl.F90"
#endif

end module linalg

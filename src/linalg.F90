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
      FLOAT, intent(in)    :: da
      FLOAT, intent(inout) :: dx ! dx(n)
    end subroutine dscal
    
    ! scales a vector by a constant
    subroutine zscal(n, za, zx, incx)
      integer,    intent(in)    :: n, incx
      CMPLX, intent(in)    :: za
      CMPLX, intent(inout) :: zx ! zx(n)
    end subroutine zscal
    
    ! constant times a vector plus a vector
    subroutine daxpy(n, da, dx, incx, dy, incy)
      integer, intent(in)    :: n, incx, incy
      FLOAT, intent(in)    :: da, dx ! dx(n)
      FLOAT, intent(inout) :: dy     ! dy(n)
    end subroutine daxpy
    
    ! constant times a vector plus a vector
    subroutine zaxpy(n, za, zx, incx, zy, incy)
      integer,    intent(in)    :: n, incx, incy
      CMPLX, intent(in)    :: za, zx ! zx(n)
      CMPLX, intent(inout) :: zy     ! zy(n)
    end subroutine zaxpy
    
    ! forms the dot product of two vectors
    FLOAT function ddot(n, dx, incx, dy, incy)
      integer, intent(in) :: n, incx, incy
      FLOAT, intent(in) :: dx, dy ! dx(n), dy(n)
    end function ddot
    
    CMPLX function zdotc(n, zx, incx, zy, incy)
      integer,    intent(in) :: n, incx, incy
      CMPLX, intent(in) :: zx, zy ! zx(n), zy(n)
    end function zdotc
    
    ! returns the euclidean norm of a vector
    FLOAT function dnrm2(n, dx, incx)
      integer, intent(in) :: n, incx
      FLOAT, intent(in) :: dx ! dx(n)
    end function dnrm2
    
    ! returns the euclidean norm of a vector
    FLOAT function dznrm2(n, zx, incx)
      integer, intent(in)    :: n, incx
      CMPLX, intent(in) :: zx ! zx(n)
    end function dznrm2
    
    ! matrix-matrix multiplication plus matrix
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      FLOAT,      intent(in)    :: alpha, beta
      FLOAT,      intent(in)    :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      FLOAT,      intent(in)    :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      FLOAT,      intent(inout) :: c ! c(ldc,n) 
    end subroutine dgemm
    
    ! matrix-matrix multiplication plus matrix
    subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      CMPLX,   intent(in)    :: alpha, beta
      CMPLX,   intent(in)    :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      CMPLX,   intent(in)    :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      CMPLX,   intent(inout) :: c ! c(ldc,n)
    end subroutine zgemm
    
    ! matrix-vector multiplication plus vector
    subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      CMPLX,   intent(in)    :: alpha, beta
      CMPLX,   intent(in)    :: a ! a(lda,n)
      CMPLX,   intent(in)    :: x ! x(:)
      CMPLX,   intent(inout) :: y ! y(:)
    end subroutine zgemv
    
    ! copies a vector, x, to a vector, y
    subroutine dcopy(n, dx, incx, dy, incy)
      integer, intent(in)  :: n, incx, incy
      FLOAT, intent(in)  :: dx ! dx(n)
      FLOAT, intent(out) :: dy ! dy(n)
    end subroutine dcopy

    ! copies a vector, x, to a vector, y
    subroutine zcopy(n, zx, incx, zy, incy)
      integer,    intent(in)  :: n, incx, incy
      CMPLX, intent(in)  :: zx ! dz(n)
      CMPLX, intent(out) :: zy ! zy(n)
    end subroutine zcopy
    
  end interface

contains

#ifdef HAVE_LAPACK
#include "linalg_lapack.F90"
#else
#include "linal_gsl.F90"
#endif

end module linalg

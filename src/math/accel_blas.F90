!! Copyright (C) 2016 X. Andrade
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

#include <global.h>

module accel_blas_oct_m
#ifdef HAVE_OPENCL
  use clblas
#endif
  use accel_oct_m
  use global_oct_m
  use iso_c_binding
  use messages_oct_m
  use profiling_oct_m
  use types_oct_m

  implicit none

  private

  public ::                        &
    daccel_dot,                    &
    zaccel_dot,                    &
    daccel_nrm2,                   &
    zaccel_nrm2,                   &
    daccel_herk,                   &
    zaccel_herk,                   &
    daccel_trsm,                   &
    zaccel_trsm,                   &
    daccel_gemm,                   &
    zaccel_gemm

  integer, parameter, public ::                      &
    CUBLAS_DIAG_NON_UNIT = 0,                        &
    CUBLAS_DIAG_UNIT     = 1

  integer, parameter, public ::                      &
    CUBLAS_OP_N = 0,                                 &
    CUBLAS_OP_T = 1,                                 &  
    CUBLAS_OP_C = 2

  integer, parameter, public ::                      &
    CUBLAS_FILL_MODE_LOWER = 0,                      &
    CUBLAS_FILL_MODE_UPPER = 1

  integer, parameter, public ::                      &
    CUBLAS_SIDE_LEFT  = 0,                           &
    CUBLAS_SIDE_RIGHT = 1

#ifdef HAVE_OPENCL
  integer, parameter, public ::                      &
    ACCEL_BLAS_LEFT  = clblasLeft,                   &
    ACCEL_BLAS_RIGHT = clblasRight

  integer, parameter, public ::                      &
    ACCEL_BLAS_LOWER = clblasLower,                  &
    ACCEL_BLAS_UPPER = clblasUpper
  
  integer, parameter, public ::                      &
    ACCEL_BLAS_N = clblasNoTrans,                    &
    ACCEL_BLAS_T = clblasTrans,                      &
    ACCEL_BLAS_C = clblasConjTrans
  
  integer, parameter, public ::                      &
    ACCEL_BLAS_DIAG_NON_UNIT = clblasNonUnit,        &
    ACCEL_BLAS_DIAG_UNIT     = clblasUnit
#else 
  integer, parameter, public ::                      &
    ACCEL_BLAS_LEFT  = CUBLAS_SIDE_LEFT,             &
    ACCEL_BLAS_RIGHT = CUBLAS_SIDE_RIGHT

  integer, parameter, public ::                      &
    ACCEL_BLAS_LOWER = CUBLAS_FILL_MODE_LOWER,       &
    ACCEL_BLAS_UPPER = CUBLAS_FILL_MODE_UPPER
  
  integer, parameter, public ::                      &
    ACCEL_BLAS_N = CUBLAS_OP_N,                      &
    ACCEL_BLAS_T = CUBLAS_OP_T,                      &
    ACCEL_BLAS_C = CUBLAS_OP_C

  integer, parameter, public ::                      &
    ACCEL_BLAS_DIAG_NON_UNIT = CUBLAS_DIAG_NON_UNIT, &
    ACCEL_BLAS_DIAG_UNIT     = CUBLAS_DIAG_UNIT
#endif
  
  ! DOT
  interface
    subroutine cuda_blas_ddot(handle, n, x, offx, incx, y, offy, incy, res, offres)
      use iso_c_binding
      
      implicit none
      
      type(c_ptr),  intent(in)    :: handle
      integer(8),   intent(in)    :: n
      type(c_ptr),  intent(in)    :: x
      integer(8),   intent(in)    :: offx
      integer(8),   intent(in)    :: incx
      type(c_ptr),  intent(in)    :: y
      integer(8),   intent(in)    :: offy
      integer(8),   intent(in)    :: incy
      type(c_ptr),  intent(inout) :: res
      integer(8),   intent(in)    :: offres
    end subroutine cuda_blas_ddot

    subroutine cuda_blas_zdotc(handle, n, x, offx, incx, y, offy, incy, res, offres)
      use iso_c_binding
            
      implicit none
      
      type(c_ptr),  intent(in)    :: handle
      integer(8),   intent(in)    :: n
      type(c_ptr),  intent(in)    :: x
      integer(8),   intent(in)    :: offx
      integer(8),   intent(in)    :: incx
      type(c_ptr),  intent(in)    :: y
      integer(8),   intent(in)    :: offy
      integer(8),   intent(in)    :: incy
      type(c_ptr),  intent(inout) :: res
      integer(8),   intent(in)    :: offres
    end subroutine cuda_blas_zdotc
  end interface

  ! NRM2
  interface
    subroutine cuda_blas_dnrm2(handle, n, x, offx, incx, res, offres)
      use iso_c_binding
      
      implicit none
      
      type(c_ptr),  intent(in)    :: handle
      integer(8),   intent(in)    :: n
      type(c_ptr),  intent(in)    :: x
      integer(8),   intent(in)    :: offx
      integer(8),   intent(in)    :: incx
      type(c_ptr),  intent(inout) :: res
      integer(8),   intent(in)    :: offres
    end subroutine cuda_blas_dnrm2

    subroutine cuda_blas_znrm2(handle, n, x, offx, incx, res, offres)
      use iso_c_binding
            
      implicit none
      
      type(c_ptr),  intent(in)    :: handle
      integer(8),   intent(in)    :: n
      type(c_ptr),  intent(in)    :: x
      integer(8),   intent(in)    :: offx
      integer(8),   intent(in)    :: incx
      type(c_ptr),  intent(inout) :: res
      integer(8),   intent(in)    :: offres
    end subroutine cuda_blas_znrm2
  end interface

  ! GEMM
  interface
    subroutine cuda_blas_dgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
      use iso_c_binding
      
      implicit none
      
      type(c_ptr),  intent(in)    :: handle
      integer,      intent(in)    :: transa
      integer,      intent(in)    :: transb
      integer(8),   intent(in)    :: m
      integer(8),   intent(in)    :: n
      integer(8),   intent(in)    :: k
      type(c_ptr),  intent(in)    :: alpha
      type(c_ptr),  intent(in)    :: A
      integer(8),   intent(in)    :: lda
      type(c_ptr),  intent(in)    :: B
      integer(8),   intent(in)    :: ldb
      type(c_ptr),  intent(in)    :: beta
      type(c_ptr),  intent(inout) :: C
      integer(8),   intent(in)    :: ldc       
    end subroutine cuda_blas_dgemm

    subroutine cuda_blas_zgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
      use iso_c_binding
      
      implicit none

      type(c_ptr),  intent(in)    :: handle
      integer,      intent(in)    :: transa
      integer,      intent(in)    :: transb
      integer(8),   intent(in)    :: m
      integer(8),   intent(in)    :: n
      integer(8),   intent(in)    :: k
      type(c_ptr),  intent(in)    :: alpha
      type(c_ptr),  intent(in)    :: A
      integer(8),   intent(in)    :: lda
      type(c_ptr),  intent(in)    :: B
      integer(8),   intent(in)    :: ldb
      type(c_ptr),  intent(in)    :: beta
      type(c_ptr),  intent(inout) :: C
      integer(8),   intent(in)    :: ldc       
    end subroutine cuda_blas_zgemm
  end interface

  ! SYRK/HERK
  interface
    subroutine cuda_blas_dsyrk(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
      use iso_c_binding
      
      implicit none
      
      type(c_ptr),  intent(in)    :: handle
      integer,      intent(in)    :: uplo
      integer,      intent(in)    :: trans
      integer(8),   intent(in)    :: n
      integer(8),   intent(in)    :: k
      type(c_ptr),  intent(in)    :: alpha
      type(c_ptr),  intent(in)    :: A
      integer(8),   intent(in)    :: lda
      type(c_ptr),  intent(in)    :: beta
      type(c_ptr),  intent(inout) :: C
      integer(8),   intent(in)    :: ldc       
    end subroutine cuda_blas_dsyrk

    subroutine cuda_blas_zherk(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
      use iso_c_binding
      
      implicit none
      
      type(c_ptr),  intent(in)    :: handle
      integer,      intent(in)    :: uplo
      integer,      intent(in)    :: trans
      integer(8),   intent(in)    :: n
      integer(8),   intent(in)    :: k
      type(c_ptr),  intent(in)    :: alpha
      type(c_ptr),  intent(in)    :: A
      integer(8),   intent(in)    :: lda
      type(c_ptr),  intent(in)    :: beta
      type(c_ptr),  intent(inout) :: C
      integer(8),   intent(in)    :: ldc       
    end subroutine cuda_blas_zherk
  end interface

  ! TRSM
  interface
    subroutine cuda_blas_dtrsm(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb)
      use iso_c_binding
      
      implicit none
      
      type(c_ptr),  intent(in)    :: handle
      integer,      intent(in)    :: side
      integer,      intent(in)    :: uplo
      integer,      intent(in)    :: trans
      integer,      intent(in)    :: diag
      integer(8),   intent(in)    :: m
      integer(8),   intent(in)    :: n
      type(c_ptr),  intent(in)    :: alpha
      type(c_ptr),  intent(in)    :: A
      integer(8),   intent(in)    :: lda
      type(c_ptr),  intent(inout) :: B
      integer(8),   intent(in)    :: ldb       
    end subroutine cuda_blas_dtrsm
    
    subroutine cuda_blas_ztrsm(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb)
      use iso_c_binding
      
      implicit none
      
      type(c_ptr),  intent(in)    :: handle
      integer,      intent(in)    :: side
      integer,      intent(in)    :: uplo
      integer,      intent(in)    :: trans
      integer,      intent(in)    :: diag
      integer(8),   intent(in)    :: m
      integer(8),   intent(in)    :: n
      type(c_ptr),  intent(in)    :: alpha
      type(c_ptr),  intent(in)    :: A
      integer(8),   intent(in)    :: lda
      type(c_ptr),  intent(inout) :: B
      integer(8),   intent(in)    :: ldb       
    end subroutine cuda_blas_ztrsm
  end interface
  
contains
  
#include "undef.F90"
#include "complex.F90"
#include "accel_blas_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "accel_blas_inc.F90"

end module accel_blas_oct_m

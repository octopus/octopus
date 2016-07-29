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
    ACCEL_BLAS_DIAG_NON_UNIT = clblasUnit,           &
    ACCEL_BLAS_DIAG_UNIT     = clblasNonUnit
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
  
  public ::                        &
    cuda_blas_ddot,                &
    cuda_blas_zdotc,               &
    cuda_blas_dgemm,               &
    cuda_blas_zgemm,               &
    daccel_syrk,                   &
    zaccel_herk,                   &
    daccel_trsm
  
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

  ! GEMM
  interface
    subroutine cuda_blas_dgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
      use iso_c_binding
      
      implicit none
      
      type(c_ptr),  intent(in)    :: handle
      integer,      intent(in)    :: transa
      integer,      intent(in)    :: transb
      integer,      intent(in)    :: m
      integer,      intent(in)    :: n
      integer,      intent(in)    :: k
      real(8),      intent(in)    :: alpha
      type(c_ptr),  intent(in)    :: A
      integer,      intent(in)    :: lda
      type(c_ptr),  intent(in)    :: B
      integer,      intent(in)    :: ldb
      real(8),      intent(in)    :: beta
      type(c_ptr),  intent(inout) :: C
      integer,      intent(in)    :: ldc       
    end subroutine cuda_blas_dgemm

    subroutine cuda_blas_zgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
      use iso_c_binding
      
      implicit none

      type(c_ptr),  intent(in)    :: handle
      integer,      intent(in)    :: transa
      integer,      intent(in)    :: transb
      integer,      intent(in)    :: m
      integer,      intent(in)    :: n
      integer,      intent(in)    :: k
      complex(8),   intent(in)    :: alpha
      type(c_ptr),  intent(in)    :: A
      integer,      intent(in)    :: lda
      type(c_ptr),  intent(in)    :: B
      integer,      intent(in)    :: ldb
      complex(8),   intent(in)    :: beta
      type(c_ptr),  intent(inout) :: C
      integer,      intent(in)    :: ldc       
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

  end interface
  
contains

  subroutine daccel_syrk(uplo, trans, n, k, alpha, a, offa, lda, beta, c, offc, ldc)
    integer,           intent(in)    :: uplo
    integer,           intent(in)    :: trans
    integer(8),        intent(in)    :: n
    integer(8),        intent(in)    :: k
    real(8),           intent(in)    :: alpha
    type(accel_mem_t), intent(in)    :: a
    integer(8),        intent(in)    :: offa
    integer(8),        intent(in)    :: lda
    real(8),           intent(in)    :: beta
    type(accel_mem_t), intent(inout) :: c
    integer(8),        intent(in)    :: offc
    integer(8),        intent(in)    :: ldc   

    type(accel_mem_t) :: alpha_buffer, beta_buffer
    
    PUSH_SUB(daccel_syrk)

    ASSERT(offa == 0)
    ASSERT(offc == 0)

    call accel_create_buffer(alpha_buffer, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, 1)
    call accel_create_buffer(beta_buffer, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, 1)

    call accel_write_buffer(alpha_buffer, alpha)
    call accel_write_buffer(beta_buffer, beta)

#ifdef HAVE_CUDA    
    call cuda_blas_dsyrk(handle = accel%cublas_handle, uplo = uplo, trans = trans, &
      n = n, k = k, &
      alpha = alpha_buffer%cuda_ptr, &
      A = a%cuda_ptr, lda = lda, &
      beta = beta_buffer%cuda_ptr, &
      C = c%cuda_ptr, ldc = ldc)
#endif

    call accel_finish()
    
    POP_SUB(daccel_syrk)
  end subroutine daccel_syrk

  ! -----------------------------------------------------------------------------------

  subroutine zaccel_herk(uplo, trans, n, k, alpha, a, offa, lda, beta, c, offc, ldc)
    integer,           intent(in)    :: uplo
    integer,           intent(in)    :: trans
    integer(8),        intent(in)    :: n
    integer(8),        intent(in)    :: k
    real(8),           intent(in)    :: alpha
    type(accel_mem_t), intent(in)    :: a
    integer(8),        intent(in)    :: offa
    integer(8),        intent(in)    :: lda
    real(8),           intent(in)    :: beta
    type(accel_mem_t), intent(inout) :: c
    integer(8),        intent(in)    :: offc
    integer(8),        intent(in)    :: ldc   

    type(accel_mem_t) :: alpha_buffer, beta_buffer
    
    PUSH_SUB(zaccel_herk)

    ASSERT(offa == 0)
    ASSERT(offc == 0)

    call accel_create_buffer(alpha_buffer, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, 1)
    call accel_create_buffer(beta_buffer, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, 1)

    call accel_write_buffer(alpha_buffer, alpha)
    call accel_write_buffer(beta_buffer, beta)

#ifdef HAVE_CUDA    
    call cuda_blas_zherk(handle = accel%cublas_handle, uplo = uplo, trans = trans, &
      n = n, k = k, &
      alpha = alpha_buffer%cuda_ptr, &
      A = a%cuda_ptr, lda = lda, &
      beta = beta_buffer%cuda_ptr, &
      C = c%cuda_ptr, ldc = ldc)
#endif

    call accel_finish()
    
    POP_SUB(zaccel_herk)
  end subroutine zaccel_herk

  ! -----------------------------------------------------------------------------------
  
  subroutine daccel_trsm(side, uplo, trans, diag, m, n, alpha, a, offa, lda, b, offb, ldb)
    integer,           intent(in)    :: side
    integer,           intent(in)    :: uplo
    integer,           intent(in)    :: trans
    integer,           intent(in)    :: diag
    integer(8),        intent(in)    :: m
    integer(8),        intent(in)    :: n
    real(8),           intent(in)    :: alpha
    type(accel_mem_t), intent(inout) :: A
    integer(8),        intent(in)    :: offa
    integer(8),        intent(in)    :: lda
    type(accel_mem_t), intent(inout) :: B
    integer(8),        intent(in)    :: offb
    integer(8),        intent(in)    :: ldb
    
    type(accel_mem_t) :: alpha_buffer
    integer :: ierr
    
    PUSH_SUB(daccel_trsm)

    ASSERT(offa == 0)
    ASSERT(offb == 0)

#ifdef HAVE_CUDA    
    call accel_create_buffer(alpha_buffer, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, 1)
    call accel_write_buffer(alpha_buffer, alpha)

    call cuda_blas_dtrsm(handle = accel%cublas_handle, side = side, uplo = uplo, trans = trans, diag = diag, &
      m = m, n = n, alpha = alpha_buffer%cuda_ptr, A = a%cuda_ptr, lda = lda, B = b%cuda_ptr, ldb = ldb)
#endif

#ifdef HAVE_OPENCL
    call clblasDtrsmEx(order = clblasColumnMajor, side = side, uplo = uplo, transA = trans, diag = diag, &
      M = m, N = n, alpha = alpha, A = a%mem, offA = offa, lda = lda, &
      B = b%mem, offB = offb, ldb = ldb, &
      CommandQueue = accel%command_queue, status = ierr)
    if(ierr /= clblasSuccess) call clblas_print_error(ierr, 'clblasDtrsmEx')
#endif
    
    call accel_finish()

#ifdef HAVE_CUDA
    call accel_release_buffer(alpha_buffer)
#endif
    
    POP_SUB(daccel_trsm)
  end subroutine daccel_trsm
  
end module accel_blas_oct_m

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
  use accel_oct_m
  use global_oct_m
  use iso_c_binding
  use messages_oct_m
  use profiling_oct_m
  use types_oct_m

  implicit none

  private

  integer, parameter, public ::  &
    ACCEL_BLAS_L = 0,            &
    ACCEL_BLAS_U = 1
  
  integer, parameter, public ::  &
    ACCEL_BLAS_N = 0,            &
    ACCEL_BLAS_T = 1,            &
    ACCEL_BLAS_C = 2

  integer, parameter, public ::  &
    CUBLAS_OP_N = 0,             &
    CUBLAS_OP_T = 1,             &  
    CUBLAS_OP_C = 2

  integer, parameter, public ::  &
    CUBLAS_FILL_MODE_LOWER = 0,  &
    CUBLAS_FILL_MODE_UPPER = 1

  public ::                        &
    cuda_blas_ddot,                &
    cuda_blas_zdotc,               &
    cuda_blas_dgemm,               &
    cuda_blas_zgemm,               &
    daccel_syrk
  
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

    FLOAT, allocatable :: aa(:, :), cc(:, :)
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
  
end module accel_blas_oct_m

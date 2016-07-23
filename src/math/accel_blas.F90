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

module accel_blas_oct_m
  use iso_c_binding

  implicit none

  private

  integer, parameter, public ::  &
    CUBLAS_OP_N = 0,             &
    CUBLAS_OP_T = 1,             &  
    CUBLAS_OP_C = 2

  public ::                        &
    cuda_blas_ddot,                &
    cuda_blas_zdotc,               &
    cuda_blas_dgemm,               &
    cuda_blas_zgemm

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

end module accel_blas_oct_m

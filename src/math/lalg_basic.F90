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

module lalg_basic_m
  use blas_m
  use global_m
  use messages_m
  use multicomm_m
  use profiling_m
  use utils_m

  implicit none

  private
  public ::      &
    lalg_swap,   &
    lalg_scal,   &
    lalg_axpy,   &
    lalg_copy,   &
    lalg_nrm2,   &
    lalg_symv,   &
    lalg_gemv,   &
    lalg_gemm,   &
    lalg_gemmt,  &
    lalg_trmm,   &
    lalg_symm
  ! ------------------------------------------------------------------
  ! BLAS level I
  ! ------------------------------------------------------------------
  
  !> swap two vectors
  interface lalg_swap
    module procedure swap_1_1
    module procedure swap_2_1
    module procedure swap_3_1
    module procedure swap_4_1
    module procedure swap_1_2
    module procedure swap_2_2
    module procedure swap_3_2
    module procedure swap_4_2
    module procedure swap_1_3
    module procedure swap_2_3
    module procedure swap_3_3
    module procedure swap_4_3
    module procedure swap_1_4
    module procedure swap_2_4
    module procedure swap_3_4
    module procedure swap_4_4
  end interface lalg_swap

  !> scales a vector by a constant
  interface lalg_scal
    module procedure scal_1_1
    module procedure scal_2_1
    module procedure scal_3_1
    module procedure scal_4_1
    module procedure scal_1_2
    module procedure scal_2_2
    module procedure scal_3_2
    module procedure scal_4_2
    module procedure scal_1_3
    module procedure scal_2_3
    module procedure scal_3_3
    module procedure scal_4_3
    module procedure scal_1_4
    module procedure scal_2_4
    module procedure scal_3_4
    module procedure scal_4_4
    module procedure scal_5_3
    module procedure scal_5_4
  end interface lalg_scal

  !> constant times a vector plus a vector
  interface lalg_axpy
    module procedure axpy_1_1
    module procedure axpy_2_1
    module procedure axpy_3_1
    module procedure axpy_4_1
    module procedure axpy_1_2
    module procedure axpy_2_2
    module procedure axpy_3_2
    module procedure axpy_4_2
    module procedure axpy_1_3
    module procedure axpy_2_3
    module procedure axpy_3_3
    module procedure axpy_4_3
    module procedure axpy_1_4
    module procedure axpy_2_4
    module procedure axpy_3_4
    module procedure axpy_4_4
    module procedure axpy_5_3
    module procedure axpy_5_4
    module procedure axpy_6_3
    module procedure axpy_6_4
  end interface lalg_axpy
  
  !> Copies a vector x, to a vector y
  interface lalg_copy
    module procedure copy_1_1
    module procedure copy_2_1
    module procedure copy_3_1
    module procedure copy_4_1
    module procedure copy_1_2
    module procedure copy_2_2
    module procedure copy_3_2
    module procedure copy_4_2
    module procedure copy_1_3
    module procedure copy_2_3
    module procedure copy_3_3
    module procedure copy_4_3
    module procedure copy_1_4
    module procedure copy_2_4
    module procedure copy_3_4
    module procedure copy_4_4
  end interface lalg_copy

  !> Returns the euclidean norm of a vector
  interface lalg_nrm2
    module procedure nrm2_1
    module procedure nrm2_2
    module procedure nrm2_3
    module procedure nrm2_4
  end interface lalg_nrm2
  
  ! ------------------------------------------------------------------
  ! BLAS level II
  ! ------------------------------------------------------------------
    
  !> Matrix-vector multiplication plus vector.
  interface lalg_symv
    module procedure symv_1_1
    module procedure symv_1_2
    module procedure symv_1_3
    module procedure symv_1_4
    module procedure symv_2_1
    module procedure symv_2_2
    module procedure symv_2_3
    module procedure symv_2_4
  end interface lalg_symv

  interface lalg_gemv
    module procedure gemv_1_1
    module procedure gemv_1_2
    module procedure gemv_1_3
    module procedure gemv_1_4
    module procedure gemv_2_1
    module procedure gemv_2_2
    module procedure gemv_2_3
    module procedure gemv_2_4
  end interface lalg_gemv

  ! ------------------------------------------------------------------
  ! BLAS level III
  ! ------------------------------------------------------------------
  
  !> Matrix-matrix multiplication plus matrix.
  interface lalg_gemm
    module procedure gemm_1_1
    module procedure gemm_1_2
    module procedure gemm_1_3
    module procedure gemm_1_4
    module procedure gemm_2_1
    module procedure gemm_2_2
    module procedure gemm_2_3
    module procedure gemm_2_4
  end interface lalg_gemm

  !> The same as above but with (Hermitian) transpose of A.
  interface lalg_gemmt
    module procedure gemmt_1_1
    module procedure gemmt_1_2
    module procedure gemmt_1_3
    module procedure gemmt_1_4
    module procedure gemmt_2_1
    module procedure gemmt_2_2
    module procedure gemmt_2_3
    module procedure gemmt_2_4
  end interface lalg_gemmt
  
  !> The following matrix multiplications all expect upper triangular matrices for a.
  !! For real matrices, \f$A = A^T\f$, for complex matrices \f$A = A^H\f$.
  interface lalg_symm
    module procedure symm_1_1
    module procedure symm_1_2
    module procedure symm_1_3
    module procedure symm_1_4
    module procedure symm_2_1
    module procedure symm_2_2
    module procedure symm_2_3
    module procedure symm_2_4
  end interface lalg_symm

  !> Matrix-matrix multiplication.
  interface lalg_trmm
    module procedure trmm_1_1
    module procedure trmm_1_2
    module procedure trmm_1_3
    module procedure trmm_1_4
  end interface lalg_trmm

  type(profile_t), save :: axpy_profile, copy_profile, gemv_profile, symv_profile

contains

#  define TYPE 1
#  include "lalg_basic_blas_inc.F90"
#  undef TYPE

#  define TYPE 2
#  include "lalg_basic_blas_inc.F90"
#  undef TYPE

#  define TYPE 3
#  include "lalg_basic_blas_inc.F90"
#  undef TYPE

#  define TYPE 4
#  include "lalg_basic_blas_inc.F90"
#  undef TYPE

end module lalg_basic_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

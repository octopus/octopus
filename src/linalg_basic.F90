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
!!
!! $Id$

#include "global.h"

module lib_basic_alg_m
  use blas_m

  implicit none

  private
  public ::    &
    lalg_swap, &
    lalg_scal, &
    lalg_axpy, &
    lalg_copy, &
    lalg_dot,  &
    lalg_nrm2, &
    lalg_gemv, &
    lalg_gemm


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
  end interface

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
  end interface

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
  end interface

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
  end interface

  interface lalg_dot
    module procedure dot_1
    module procedure dot_2
    module procedure dot_3
    module procedure dot_4
  end interface

  interface lalg_nrm2
    module procedure nrm2_1
    module procedure nrm2_2
    module procedure nrm2_3
    module procedure nrm2_4
  end interface

  interface lalg_gemm
    module procedure gemm_1_1
    module procedure gemm_1_2
    module procedure gemm_1_3
    module procedure gemm_1_4
    module procedure gemm_2_1
    module procedure gemm_2_2
    module procedure gemm_2_3
    module procedure gemm_2_4
  end interface

  interface lalg_gemv
    module procedure gemv_1_1
    module procedure gemv_1_2
    module procedure gemv_1_3
    module procedure gemv_1_4
    module procedure gemv_2_1
    module procedure gemv_2_2
    module procedure gemv_2_3
    module procedure gemv_2_4
  end interface

contains

#if defined(HAVE_BLAS)

#  define TYPE 1
#  include "linalg_basic_blas.F90"
#  undef TYPE

#  define TYPE 2
#  include "linalg_basic_blas.F90"
#  undef TYPE

#  define TYPE 3
#  include "linalg_basic_blas.F90"
#  undef TYPE

#  define TYPE 4
#  include "linalg_basic_blas.F90"
#  undef TYPE

#else

#  define TYPE 1
#  include "linalg_basic_int.F90"
#  undef TYPE

#  define TYPE 2
#  include "linalg_basic_int.F90"
#  undef TYPE

#  define TYPE 3
#  include "linalg_basic_int.F90"
#  undef TYPE

#  define TYPE 4
#  include "linalg_basic_int.F90"
#  undef TYPE

#endif

end module lib_basic_alg_m

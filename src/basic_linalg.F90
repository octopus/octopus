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

module lib_basic_alg
  use global

  implicit none

  interface lalg_swap
    module procedure d_swap, z_swap
  end interface

  interface lalg_scal
    module procedure d_scal, z_scal
  end interface
    
  interface lalg_axpy
     module procedure d_axpy, z_axpy
  end interface

  interface lalg_dot
     module procedure d_dot, z_dot
  end interface

  interface lalg_nrm2
     module procedure d_nrm2, dz_nrm2
  end interface

  interface lalg_copy
     module procedure d_copy, z_copy
  end interface

  interface lalg_gemm
     module procedure d_gemm, z_gemm
  end interface

  interface lalg_gemv
     module procedure d_gemv, z_gemv
  end interface

  private
  public :: lalg_swap, lalg_scal, lalg_axpy, lalg_dot, lalg_nrm2, lalg_copy, lalg_gemm, lalg_gemv

contains

#ifdef HAVE_BLAS
#include "basic_linalg_blas.F90"
#else
#include "basic_linalg_int.F90"
#endif

end module lib_basic_alg

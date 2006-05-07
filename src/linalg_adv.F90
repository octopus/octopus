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

module lib_adv_alg_m
  use global_m
  use messages_m

  implicit none

  private
  public ::                     &
    lalg_geneigensolve,         &
    lalg_eigensolve,            &
    lalg_determinant,           &
    lalg_inverter,              &
    lalg_linsyssolve,           &
    lalg_singular_value_decomp, &
    lalg_svd_inverse

  interface lalg_geneigensolve
    module procedure dgeneigensolve, zgeneigensolve
  end interface

  interface lalg_eigensolve
    module procedure deigensolve, zeigensolve
  end interface

  ! Note that lalg_determinant and lalg_inverter are just wrappers
  ! over the same routine.
  interface lalg_determinant
    module procedure ddeterminant, zdeterminant
  end interface

  interface lalg_inverter
    module procedure ddeterminant, zdeterminant
  end interface

  interface lalg_linsyssolve
    module procedure dlinsyssolve
  end interface

  interface lalg_singular_value_decomp
    module procedure zsingular_value_decomp
  end interface

  interface lalg_svd_inverse
    module procedure zsvd_inverse
  end interface

contains

#ifdef HAVE_LAPACK
#include "linalg_adv_lapack.F90"
#endif

end module lib_adv_alg_m

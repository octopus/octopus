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

module lalg_adv_m
  use global_m
  use messages_m

  implicit none

  private
  public ::                       &
    lalg_cholesky,                &
    lalg_geneigensolve,           &
    lalg_eigensolve,              &
    lalg_eigensolve_nonh,         &
    lalg_determinant,             &
    lalg_inverter,                &
    lalg_sym_inverter,            &
    lalg_linsyssolve,             &
    lalg_singular_value_decomp,   &
    lalg_svd_inverse,             &
    lalg_invert_upper_triangular, &
    lalg_invert_lower_triangular, &
    lalg_lowest_geneigensolve,    &
    lalg_lowest_eigensolve

  interface lalg_cholesky
    module procedure dcholesky, zcholesky
  end interface

  interface lalg_geneigensolve
    module procedure dgeneigensolve, zgeneigensolve
  end interface

  interface lalg_eigensolve_nonh
    module procedure zeigensolve_nonh, deigensolve_nonh
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

  interface lalg_sym_inverter
    module procedure dsym_inverter, zsym_inverter
  end interface

  interface lalg_linsyssolve
    module procedure dlinsyssolve, zlinsyssolve
  end interface

  interface lalg_singular_value_decomp
    module procedure zsingular_value_decomp
  end interface

  interface lalg_svd_inverse
    module procedure zsvd_inverse
  end interface

  interface lalg_invert_upper_triangular
    module procedure dinvert_upper_triangular, zinvert_upper_triangular
  end interface
  
  interface lalg_invert_lower_triangular
    module procedure dinvert_lower_triangular, zinvert_lower_triangular
  end interface
  
  interface lalg_lowest_geneigensolve
    module procedure dlowest_geneigensolve, zlowest_geneigensolve
  end interface

  interface lalg_lowest_eigensolve
    module procedure dlowest_eigensolve, zlowest_eigensolve
  end interface
contains

#ifdef HAVE_LAPACK
#include "lalg_adv_lapack.F90"
#endif

end module lalg_adv_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

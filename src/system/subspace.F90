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
!! $Id: eigen.F90 4287 2008-06-15 22:20:10Z xavier $

#include "global.h"

module subspace_m
  use batch_m
  use blas_m
  use blacs_proc_grid_m
  use datasets_m
  use derivatives_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lalg_adv_m
  use lalg_basic_m
  use math_m
  use mesh_m
  use mesh_function_m
  use mesh_batch_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use parser_m
  use pblas_m
  use preconditioners_m
  use profiling_m
  use scalapack_m
  use states_m
  use states_block_m
  use states_calc_m
  use types_m
  use varinfo_m

  implicit none

  private
  public ::             &
    dsubspace_diag,     &
    zsubspace_diag,     &
    dsubspace_test,     &
    zsubspace_test

contains

#include "undef.F90"
#include "real.F90"
#include "subspace_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "subspace_inc.F90"

end module subspace_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

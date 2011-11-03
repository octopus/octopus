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

module mesh_batch_m
  use batch_m
  use blas_m
  use c_pointer_m
  use cl_m
  use cl_kernel_m
  use comm_m
  use global_m
  use hardware_m
  use index_m
  use lalg_basic_m
  use loct_math_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use opencl_m
  use par_vec_m
  use profiling_m
  use types_m

  implicit none

  private
  public ::                         &
    dmesh_batch_dotp_matrix,        &
    zmesh_batch_dotp_matrix,        &
    dmesh_batch_dotp_vector,        &
    zmesh_batch_dotp_vector,        &
    dmesh_batch_dotp_self,          &
    zmesh_batch_dotp_self,          &
    dmesh_batch_rotate,             &
    zmesh_batch_rotate,             &
    dmesh_batch_exchange_points,    &
    zmesh_batch_exchange_points

contains

#include "undef.F90"
#include "real.F90"
#include "mesh_batch_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mesh_batch_inc.F90"

end module mesh_batch_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

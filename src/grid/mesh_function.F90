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

module mesh_function_m
  use batch_m
  use blas_m
  use comm_m
  use cube_function_m
  use global_m
  use hardware_m
  use hypercube_m
  use index_m
  use lalg_basic_m
  use loct_math_m
  use math_m
  use mesh_m
  use messages_m
  use modelmb_1part_m
  use mpi_m
  use par_vec_m
  use profiling_m
  use qshepmod_m
  use splines_m

  implicit none

  private
  public ::                &
    dmf_integrate,         &
    zmf_integrate,         &
    dmf_add,               &
    zmf_add,               &
    dmf_copy,              &
    zmf_copy,              &
    dmf_dotp,              &
    zmf_dotp,              &
    dmf_nrm2,              &
    zmf_nrm2,              &
    dmf_moment,            &
    zmf_moment,            &
    dmf_random,            &
    zmf_random,            &
    dmf_partial_integrate, &
    zmf_partial_integrate, &
    dmf_interpolate,       &
    zmf_interpolate,       &
    dmf_interpolate_points,&
    zmf_interpolate_points,&
    mf_surface_integral,   &
    mf_line_integral,      &
    dmf_put_radial_spline, &
    zmf_put_radial_spline, &
    dmf_dotp_aux,          &
    zmf_dotp_aux,          &
    dmf_multipoles,        &
    zmf_multipoles,        &
    mesh_init_mesh_aux

  interface mf_surface_integral
    module procedure dmf_surface_integral_scalar, dmf_surface_integral_vector, &
                     zmf_surface_integral_scalar, zmf_surface_integral_vector
  end interface mf_surface_integral

  interface mf_line_integral
    module procedure dmf_line_integral_scalar, dmf_line_integral_vector, &
                     zmf_line_integral_scalar, zmf_line_integral_vector
  end interface mf_line_integral

  interface dmf_dotp
    module procedure dmf_dotp_1, dmf_dotp_2
  end interface dmf_dotp

  interface zmf_dotp
    module procedure zmf_dotp_1, zmf_dotp_2
  end interface zmf_dotp

  interface dmf_nrm2
    module procedure dmf_nrm2_1, dmf_nrm2_2
  end interface dmf_nrm2

  interface zmf_nrm2
    module procedure zmf_nrm2_1, zmf_nrm2_2
  end interface zmf_nrm2

  type(mesh_t), pointer :: mesh_aux

  type(profile_t), save ::            &
       C_PROFILING_MF_INTEGRATE,      &
       C_PROFILING_MF_DOTP,           &
       C_PROFILING_MF_REDUCE,         &
       C_PROFILING_MF_NRM2

contains

  subroutine mesh_init_mesh_aux(mesh)
    type(mesh_t), target, intent(in) :: mesh

    PUSH_SUB(mesh_init_mesh_aux)
    mesh_aux => mesh

    POP_SUB(mesh_init_mesh_aux)
  end subroutine mesh_init_mesh_aux

#include "undef.F90"
#include "real.F90"
#include "mesh_function_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mesh_function_inc.F90"

end module mesh_function_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

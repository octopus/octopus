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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module mesh_function_oct_m
  use batch_oct_m
  use blas_oct_m
  use comm_oct_m
  use cube_function_oct_m
  use global_oct_m
  use hardware_oct_m
  use index_oct_m
  use lalg_basic_oct_m
  use loct_math_oct_m
  use math_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use par_vec_oct_m
  use profiling_oct_m
  use quickrnd_oct_m
  use splines_oct_m

  implicit none

  private
  public ::                &
    dmf_integrate,         &
    zmf_integrate,         &
    smf_integrate,         &
    cmf_integrate,         &
    dmf_add,               &
    zmf_add,               &
    smf_add,               &
    cmf_add,               &    
    dmf_dotp,              &
    zmf_dotp,              &
    smf_dotp,              &
    cmf_dotp,              &
    dmf_nrm2,              &
    zmf_nrm2,              &
    smf_nrm2,              &
    cmf_nrm2,              &    
    dmf_moment,            &
    zmf_moment,            &
    smf_moment,            &
    cmf_moment,            &    
    dmf_random,            &
    zmf_random,            &
    dmf_dotp_aux,          &
    zmf_dotp_aux,          &
    dmf_multipoles,        &
    zmf_multipoles,        &
    dmf_local_multipoles,  &
    zmf_local_multipoles,  &
    mesh_init_mesh_aux,    &
    dmf_dotu_aux,          &
    zmf_dotu_aux,          &
    dmf_nrm2_aux,          &
    zmf_nrm2_aux,          &
    dmf_normalize,         &
    zmf_normalize

  ! These variables are to be used by the "distdot" function, that is outside the module
  ! but inside this file.
  ! FIXME: This is very ugly, at least these values should be set by a function.
  public :: mesh_aux
  logical, public :: sp_parallel
  integer, public :: sp_np, sp_dim, sp_st1, sp_st2, sp_kp1, sp_kp2, sp_comm
  integer, public :: sp_distdot_mode
  
  interface dmf_dotp
    module procedure dmf_dotp_1, dmf_dotp_2
  end interface dmf_dotp

  interface zmf_dotp
    module procedure zmf_dotp_1, zmf_dotp_2
  end interface zmf_dotp

  interface smf_dotp
    module procedure smf_dotp_1, smf_dotp_2
  end interface smf_dotp

  interface cmf_dotp
    module procedure cmf_dotp_1, cmf_dotp_2
  end interface cmf_dotp

  interface dmf_nrm2
    module procedure dmf_nrm2_1, dmf_nrm2_2
  end interface dmf_nrm2

  interface zmf_nrm2
    module procedure zmf_nrm2_1, zmf_nrm2_2
  end interface zmf_nrm2

  interface smf_nrm2
    module procedure smf_nrm2_1, smf_nrm2_2
  end interface smf_nrm2

  interface cmf_nrm2
    module procedure cmf_nrm2_1, cmf_nrm2_2
  end interface cmf_nrm2

  type(mesh_t), pointer :: mesh_aux => null()

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

#include "undef.F90"
#include "real_single.F90"
#include "mesh_function_inc.F90"

#include "undef.F90"
#include "complex_single.F90"
#include "mesh_function_inc.F90"

end module mesh_function_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

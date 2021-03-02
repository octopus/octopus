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

module mesh_batch_oct_m
  use accel_oct_m
  use accel_blas_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use blas_oct_m
  use iso_c_binding
  use comm_oct_m
  use global_oct_m
  use hardware_oct_m
  use lalg_basic_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use multicomm_oct_m
#if defined(HAVE_OPENMP)
  use omp_lib
#endif
  use par_vec_oct_m
  use partition_oct_m
  use profiling_oct_m
  use types_oct_m

  implicit none

  private
  public ::                         &
    dmesh_batch_dotp_matrix,        &
    zmesh_batch_dotp_matrix,        &
    dmesh_batch_dotp_vector,        &
    zmesh_batch_dotp_vector,        &
    dmesh_batch_dotp_self,          &
    zmesh_batch_dotp_self,          &
    dmesh_batch_exchange_points,    &
    zmesh_batch_exchange_points,    &
    mesh_batch_nrm2,                &
    dmesh_batch_orthogonalization,  &
    zmesh_batch_orthogonalization,  &
    dmesh_batch_mf_dotp,            &
    zmesh_batch_mf_dotp,            &
    dmesh_batch_codensity,          &
    zmesh_batch_codensity

contains

! -----------------------------------------------------

  subroutine mesh_batch_nrm2(mesh, aa, nrm2, reduce)
    type(mesh_t),            intent(in)    :: mesh
    class(batch_t),          intent(in)    :: aa
    FLOAT,                   intent(out)   :: nrm2(:)
    logical,       optional, intent(in)    :: reduce
    
    PUSH_SUB(mesh_batch_nrm2)
    
    if(aa%type() == TYPE_FLOAT) then
      call dpriv_mesh_batch_nrm2(mesh, aa, nrm2)
    else
      call zpriv_mesh_batch_nrm2(mesh, aa, nrm2)
    end if
    
    if(mesh%parallel_in_domains .and. optional_default(reduce, .true.)) then
      nrm2(1:aa%nst) = nrm2(1:aa%nst)**2
      call mesh%allreduce(nrm2, dim = aa%nst)
      nrm2(1:aa%nst) = sqrt(nrm2(1:aa%nst))
    end if

    POP_SUB(mesh_batch_nrm2)
  end subroutine mesh_batch_nrm2

#include "undef.F90"
#include "real.F90"
#include "mesh_batch_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mesh_batch_inc.F90"

end module mesh_batch_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

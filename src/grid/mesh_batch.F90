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
#ifdef HAVE_OPENCL
  use cl
#ifdef HAVE_CLAMDBLAS
  use clamdblas
#endif
#endif
  use octcl_kernel_m
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
    zmesh_batch_exchange_points,    &
    mesh_batch_nrm2

contains

! -----------------------------------------------------

subroutine mesh_batch_nrm2(mesh, aa, nrm2, reduce)
  type(mesh_t),            intent(in)    :: mesh
  type(batch_t),           intent(in)    :: aa
  FLOAT,                   intent(out)   :: nrm2(:)
  logical,       optional, intent(in)    :: reduce

  PUSH_SUB(mesh_batch_nrm2)

  if(batch_type(aa) == TYPE_FLOAT) then
    call dmesh_batch_nrm2(mesh, aa, nrm2)
  else
    call zmesh_batch_nrm2(mesh, aa, nrm2)
  end if
  
  if(mesh%parallel_in_domains .and. optional_default(reduce, .true.)) then
    nrm2(1:aa%nst) = nrm2(1:aa%nst)**2
    call comm_allreduce(mesh%mpi_grp%comm, nrm2, aa%nst)
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

end module mesh_batch_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

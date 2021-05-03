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

module states_elec_calc_oct_m
  use accel_oct_m
  use accel_blas_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use blas_oct_m
  use blacs_oct_m
  use iso_c_binding
  use comm_oct_m
  use derivatives_oct_m
  use fourier_space_oct_m
  use global_oct_m
  use grid_oct_m
  use hardware_oct_m
  use kpoints_oct_m
  use lalg_adv_oct_m
  use lalg_basic_oct_m
  use messages_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use namespace_oct_m
  use par_vec_oct_m
  use pblas_oct_m
  use physics_op_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use scalapack_oct_m
  use simul_box_oct_m
  use singularity_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_parallel_oct_m
  use types_oct_m
  use wfs_elec_oct_m

  implicit none

  private

  public ::                         &
    states_elec_orthogonalize,           &
    states_elec_rotate,                  &
    dstates_elec_calc_orth_test,         &
    zstates_elec_calc_orth_test,         &
    dstates_elec_orthogonalization,      &
    zstates_elec_orthogonalization,      &
    dstates_elec_orthogonalize_single,   &
    zstates_elec_orthogonalize_single,   &
    dstates_elec_orthogonalize_single_batch,   &
    zstates_elec_orthogonalize_single_batch,   &
    dstates_elec_orthogonalization_full, &
    zstates_elec_orthogonalization_full, &
    dstates_elec_residue,                &
    zstates_elec_residue,                &
    states_elec_calc_momentum,           &
    dstates_elec_angular_momentum,       &
    zstates_elec_angular_momentum,       &
    dstates_elec_matrix,                 &
    zstates_elec_matrix,                 &
    dstates_elec_calc_overlap,           &
    zstates_elec_calc_overlap,           &
    dstates_elec_calc_projections,       &
    zstates_elec_calc_projections,       &
    dstates_elec_me_one_body,            &
    zstates_elec_me_one_body,            &
    dstates_elec_me_two_body,            &
    zstates_elec_me_two_body,            &
    dstates_elec_rrqr_decomposition,     &
    zstates_elec_rrqr_decomposition

  interface states_elec_rotate
    module procedure dstates_elec_rotate, zstates_elec_rotate
  end interface states_elec_rotate
  
contains

  ! ---------------------------------------------------------

  subroutine states_elec_orthogonalize(st, namespace, mesh)
    type(states_elec_t),  intent(inout) :: st
    type(namespace_t),    intent(in)    :: namespace
    type(mesh_t),         intent(in)    :: mesh

    integer :: ik

    PUSH_SUB(states_elec_orthogonalize)

    do ik = st%d%kpt%start, st%d%kpt%end
      if (states_are_real(st)) then
        call dstates_elec_orthogonalization_full(st, namespace, mesh, ik)
      else
        call zstates_elec_orthogonalization_full(st, namespace, mesh, ik)
      end if
    end do

    POP_SUB(states_elec_orthogonalize)
  end subroutine states_elec_orthogonalize

  ! -----------------------------------------------------------------------------

  subroutine states_elec_calc_momentum(st, space, der, kpoints, momentum)
    type(states_elec_t), intent(in)  :: st
    type(space_t),       intent(in)  :: space
    type(derivatives_t), intent(in)  :: der
    type(kpoints_t),     intent(in)  :: kpoints
    FLOAT,               intent(out) :: momentum(:,:,:)

    if (states_are_real(st)) then
      call dstates_elec_calc_momentum(st, space, der, kpoints, momentum)
    else
      call zstates_elec_calc_momentum(st, space, der, kpoints, momentum)
    end if
  end subroutine states_elec_calc_momentum

#include "undef.F90"
#include "real.F90"
#include "states_elec_calc_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_elec_calc_inc.F90"
#include "undef.F90"

end module states_elec_calc_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

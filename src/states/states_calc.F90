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

module states_calc_oct_m
  use accel_oct_m
  use accel_blas_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use blas_oct_m
  use blacs_oct_m
  use iso_c_binding
  use comm_oct_m
  use derivatives_oct_m
  use geometry_oct_m
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
  use pblas_oct_m
  use physics_op_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use scalapack_oct_m
  use simul_box_oct_m
  use states_oct_m
  use states_dim_oct_m
  use states_parallel_oct_m
  use types_oct_m

  implicit none

  private

  public ::                         &
    states_orthogonalize,           &
    states_rotate,                  &
    dstates_calc_orth_test,         &
    zstates_calc_orth_test,         &
    dstates_orthogonalization,      &
    zstates_orthogonalization,      &
    dstates_orthogonalize_single,   &
    zstates_orthogonalize_single,   &
    dstates_orthogonalization_full, &
    zstates_orthogonalization_full, &
    dstates_residue,                &
    zstates_residue,                &
    states_calc_momentum,           &
    dstates_angular_momentum,       &
    zstates_angular_momentum,       &
    dstates_matrix,                 &
    zstates_matrix,                 &
    dstates_calc_overlap,           &
    zstates_calc_overlap,           &
    dstates_calc_projections,       &
    zstates_calc_projections,       &
    dstates_me_one_body,            &
    zstates_me_one_body,            &
    dstates_me_two_body,            &
    zstates_me_two_body

  interface states_rotate
    module procedure dstates_rotate, zstates_rotate
  end interface states_rotate
  
contains

  ! ---------------------------------------------------------

  subroutine states_orthogonalize(st, mesh)
    type(states_t),    intent(inout) :: st
    type(mesh_t),      intent(in)    :: mesh

    integer :: ik

    PUSH_SUB(states_orthogonalize)

    do ik = st%d%kpt%start, st%d%kpt%end
      if (states_are_real(st)) then
        call dstates_orthogonalization_full(st, mesh, ik)
      else
        call zstates_orthogonalization_full(st, mesh, ik)
      end if
    end do

    POP_SUB(states_orthogonalize)
  end subroutine states_orthogonalize

  ! -----------------------------------------------------------------------------

  subroutine states_calc_momentum(st, der, momentum)
    type(states_t),      intent(in)  :: st
    type(derivatives_t), intent(in)  :: der
    FLOAT,               intent(out) :: momentum(:,:,:)

    if (states_are_real(st)) then
      call dstates_calc_momentum(st, der, momentum)
    else
      call zstates_calc_momentum(st, der, momentum)
    end if
  end subroutine states_calc_momentum

#include "undef.F90"
#include "real.F90"
#include "states_calc_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_calc_inc.F90"
#include "undef.F90"

end module states_calc_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

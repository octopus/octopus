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

module states_calc_m
  use batch_m
  use blas_m
  use blacs_m
  use blacs_proc_grid_m
  use c_pointer_m
  use calc_mode_m
  use cl
  use octcl_kernel_m
  use comm_m
  use datasets_m
  use derivatives_m
  use geometry_m
  use global_m
  use grid_m
  use hardware_m
  use io_m
  use kpoints_m
  use lalg_adv_m
  use lalg_basic_m
  use lapack_m
  use loct_m
  use math_m
  use messages_m
  use mesh_m
  use mesh_batch_m
  use mesh_function_m
  use mpi_m
  use mpi_lib_m
  use multicomm_m
  use opencl_m
  use parser_m
  use pblas_m
  use physics_op_m
  use profiling_m
  use scalapack_m
  use simul_box_m
  use smear_m
  use states_m
  use states_block_m
  use states_dim_m
  use unit_m
  use unit_system_m
  use utils_m
  use types_m
  use varinfo_m

  implicit none

  private

  public ::                         &
    states_orthogonalize,           &
    states_degeneracy_matrix,       &
    states_rotate,                  &
    dstates_rotate_in_place,        &
    zstates_rotate_in_place,        &
    dstates_calc_orth_test,         &
    zstates_calc_orth_test,         &
    dstates_orthogonalization,      &
    zstates_orthogonalization,      &
    dstates_orthogonalization_full, &
    zstates_orthogonalization_full, &
    dstates_normalize_orbital,      &
    zstates_normalize_orbital,      &
    dstates_residue,                &
    zstates_residue,                &
    states_calc_momentum,           &
    dstates_angular_momentum,       &
    zstates_angular_momentum,       &
    dstates_matrix,                 &
    zstates_matrix,                 &
    dstates_calc_overlap,           &
    zstates_calc_overlap

contains

  ! ---------------------------------------------------------
  ! This routine transforms the orbitals of state "st", according
  ! to the transformation matrix "uu".
  !
  ! Each row of u contains the coefficients of the new orbitals
  ! in terms of the old ones.
  ! ---------------------------------------------------------
  subroutine states_rotate(mesh, st, stin, uu)
    type(mesh_t),      intent(in)    :: mesh
    type(states_t),    intent(inout) :: st
    type(states_t),    intent(in)    :: stin
    CMPLX,             intent(in)    :: uu(:, :)

    integer :: ik

    PUSH_SUB(states_rotate)

    if(states_are_real(st)) then
      do ik = st%d%kpt%start, st%d%kpt%end
        call lalg_gemm(mesh%np_part*st%d%dim, st%nst, stin%nst, M_ONE, stin%dpsi(:, :, 1:stin%nst, ik), &
          transpose(real(uu(:, :), REAL_PRECISION)), M_ZERO, st%dpsi(:, :, :, ik))
      end do
    else
      do ik = st%d%kpt%start, st%d%kpt%end
        call lalg_gemm(mesh%np_part*st%d%dim, st%nst, stin%nst, M_z1, stin%zpsi(:, :, 1:stin%nst, ik), &
          transpose(uu(:, :)), M_z0, st%zpsi(:, :, :, ik))
      end do
    end if

    POP_SUB(states_rotate)
  end subroutine states_rotate

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

  ! -------------------------------------------------------
  subroutine states_degeneracy_matrix(sb, st)
    type(simul_box_t), intent(in) :: sb
    type(states_t),    intent(in) :: st

    integer :: idir, is, js, inst, inik, iunit
    integer, allocatable :: eindex(:,:), sindex(:)
    integer, allocatable :: degeneracy_matrix(:, :)
    FLOAT,   allocatable :: eigenval_sorted(:)
    FLOAT :: degen_thres, evis, evjs, kpoint(1:MAX_DIM)

    PUSH_SUB(states_degeneracy_matrix)

    SAFE_ALLOCATE(eigenval_sorted(1:st%nst*st%d%nik))
    SAFE_ALLOCATE(         sindex(1:st%nst*st%d%nik))
    SAFE_ALLOCATE(      eindex(1:2, 1:st%nst*st%d%nik))
    SAFE_ALLOCATE(degeneracy_matrix(1:st%nst*st%d%nik, 1:st%nst*st%d%nik))

    ! convert double index "inst, inik" to single index "is"
    ! and keep mapping array
    is = 1
    do inst = 1, st%nst
      do inik = 1, st%d%nik
        eigenval_sorted(is) = st%eigenval(inst, inik)        
        eindex(1, is) = inst
        eindex(2, is) = inik
        is = is + 1
      end do
    end do

    ! sort eigenvalues
    call sort(eigenval_sorted, sindex)

    !%Variable DegeneracyThreshold
    !%Type float
    !%Default 1e-5
    !%Section States
    !%Description
    !% A state j with energy E_j will be considered degenerate with a state
    !% with energy E_i, if  E_i - threshold < E_j < E_i + threshold.
    !%End
    call parse_float(datasets_check('DegeneracyThreshold'), &
       units_from_atomic(units_inp%energy, CNST(1e-5)), degen_thres)
    degen_thres = units_to_atomic(units_inp%energy, degen_thres)

    ! setup degeneracy matrix. the matrix summarizes the degeneracy relations 
    ! among the states
    degeneracy_matrix = 0

    do is = 1, st%nst*st%d%nik
      do js = 1, st%nst*st%d%nik

        ! a state is always degenerate to itself
        if ( is.eq.js ) cycle

        evis = st%eigenval(eindex(1, sindex(is)), eindex(2, sindex(is)))
        evjs = st%eigenval(eindex(1, sindex(js)), eindex(2, sindex(js)))

        ! is evjs in the "evis plus minus threshold" bracket?
        if( (evjs.gt.evis - degen_thres).and.(evjs.lt.evis + degen_thres) ) then
          ! mark forward scattering states with +1 and backward scattering
          ! states with -1
          !WARNING: IS THIS REALLY NECESSARY? - have to calculate momentum
          degeneracy_matrix(is, js) = M_ONE
          !degeneracy_matrix(is, js) = &
          !  sign(M_ONE, st%momentum(1, eindex(1, sindex(js)), eindex(2, sindex(js))))
        end if

      end do
    end do

    if(mpi_grp_is_root(mpi_world)) then

      ! write matrix to "restart/gs" directory
      iunit = io_open(trim(tmpdir)//GS_DIR//'degeneracy_matrix', action='write', is_tmp = .true.)

      write(iunit, '(a)', advance='no') '# index  '
      do idir = 1, sb%dim
        write(iunit, '(2a)', advance='no') 'k', index2axis(idir)
      enddo
      write(iunit, '(a)') ' eigenvalue  degeneracy matrix'

      do is = 1, st%nst*st%d%nik
        write(iunit, '(i6,4e24.16,32767i3)') is, &
          kpoints_get_point(sb%kpoints, states_dim_get_kpoint_index(st%d, eindex(2, sindex(is)))), &
          eigenval_sorted(is), (degeneracy_matrix(is, js), js = 1, st%nst*st%d%nik)
        write(iunit, '(i6)', advance='no') is
        kpoint(1:sb%dim) = kpoints_get_point(sb%kpoints, states_dim_get_kpoint_index(st%d, eindex(2, sindex(is))))
        do idir = 1, sb%dim
          write(iunit, '(e24.16)', advance='no') kpoint(idir)
        enddo
        write(iunit, '(e24.16)', advance='no') eigenval_sorted(is)
        do js = 1, st%nst * st%d%nik
          write(iunit, '(i3)') degeneracy_matrix(is, js)
        enddo
      end do

      call io_close(iunit)

      ! write index vectors to "restart/gs" directory
      iunit = io_open(trim(tmpdir)//GS_DIR//'index_vectors', action='write', is_tmp = .true.)    

      write(iunit, '(a)') '# index  sindex  eindex1 eindex2'

      do is = 1, st%nst*st%d%nik
        write(iunit,'(4i6)') is, sindex(is), eindex(1, sindex(is)), eindex(2, sindex(is))
      end do

      call io_close(iunit)
    end if

    SAFE_DEALLOCATE_A(eigenval_sorted)
    SAFE_DEALLOCATE_A(sindex)
    SAFE_DEALLOCATE_A(eindex)
    SAFE_DEALLOCATE_A(degeneracy_matrix)

    POP_SUB(states_degeneracy_matrix)
  end subroutine states_degeneracy_matrix

  ! -----------------------------------------------------------------------------

  subroutine states_calc_momentum(st, der, momentum)
    type(states_t),      intent(inout) :: st
    type(derivatives_t), intent(inout) :: der
    FLOAT,               intent(out)   :: momentum(:,:,:)

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

end module states_calc_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

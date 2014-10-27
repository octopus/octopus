!! Copyright (C) 2008 X. Andrade
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
!! $Id$

#include "global.h"

module gauge_field_force_m
  use batch_m
  use batch_ops_m
  use datasets_m
  use derivatives_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use lalg_basic_m
  use logrid_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use projector_m
  use ps_m
  use restart_m
  use simul_box_m
  use species_m
  use splines_m
  use states_m
  use states_dim_m
  use submesh_m
  use symmetries_m
  use symmetrizer_m
  use symm_op_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private

  public ::                               &
    gauge_field_get_force

contains

  ! ---------------------------------------------------------
  subroutine gauge_field_get_force(gr, hm, geo, pj, phases, st, force)
    type(grid_t),         intent(inout) :: gr
    type(hamiltonian_t),  intent(in)    :: hm
    type(geometry_t),     intent(in)    :: geo
    type(projector_t),    intent(in)    :: pj(:)
    CMPLX,                intent(in)    :: phases(:, :)
    type(states_t),       intent(inout) :: st
    type(gauge_force_t),  intent(out)   :: force

    integer :: ik, ist, idir, idim, iatom, ip, ib, minst, maxst, ii
    CMPLX, allocatable :: gpsi(:, :, :), psi(:, :), hpsi(:, :), rhpsi(:, :), rpsi(:, :), hrpsi(:, :)
    FLOAT, allocatable :: microcurrent(:, :), symmcurrent(:, :)
    type(profile_t), save :: prof
    type(symmetrizer_t) :: symmetrizer
    type(batch_t) :: hpsib, rhpsib, rpsib, hrpsib
#ifdef HAVE_MPI
    FLOAT :: force_tmp(1:MAX_DIM)
#endif
    logical, parameter :: hamiltonian_current = .false.

    call profiling_in(prof, "GAUGE_FIELD_FORCE")
    PUSH_SUB(gauge_field_get_force)

    SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gpsi(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:st%d%dim))
    SAFE_ALLOCATE(hpsi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(rhpsi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(rpsi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(hrpsi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(microcurrent(1:gr%mesh%np_part, 1:gr%sb%dim))

    microcurrent = M_ZERO

    if(hamiltonian_current) then
      
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end

          call batch_pack(st%group%psib(ib, ik), copy = .true.)

          call batch_copy(st%group%psib(ib, ik), hpsib, reference = .false.)
          call batch_copy(st%group%psib(ib, ik), rhpsib, reference = .false.)
          call batch_copy(st%group%psib(ib, ik), rpsib, reference = .false.)
          call batch_copy(st%group%psib(ib, ik), hrpsib, reference = .false.)

          call zderivatives_batch_set_bc(gr%der, st%group%psib(ib, ik))
          call zhamiltonian_apply_batch(hm, gr%der, st%group%psib(ib, ik), hpsib, ik, set_bc = .false.)

          do idir = 1, gr%sb%dim

            call batch_mul(gr%mesh%np, gr%mesh%x(:, idir), hpsib, rhpsib)
            call batch_mul(gr%mesh%np_part, gr%mesh%x(:, idir), st%group%psib(ib, ik), rpsib)
          
            call zhamiltonian_apply_batch(hm, gr%der, rpsib, hrpsib, ik, set_bc = .false.)

            minst = states_block_min(st, ib)
            maxst = states_block_max(st, ib)
            
            do ist = st%st_start, st%st_end
            
              do idim = 1, st%d%dim
                ii = batch_ist_idim_to_linear(st%group%psib(ib, ik), (/ist, idim/))
                call batch_get_state(st%group%psib(ib, ik), ii, gr%mesh%np, psi(:, idim))
                call batch_get_state(hrpsib, ii, gr%mesh%np, hrpsi(:, idim))
                call batch_get_state(rhpsib, ii, gr%mesh%np, rhpsi(:, idim))
              end do
              
              do idim = 1, st%d%dim
                !$omp parallel do
                do ip = 1, gr%mesh%np
                  microcurrent(ip, idir) = microcurrent(ip, idir) - &
                    CNST(4.0)*M_PI*P_c/gr%sb%rcell_volume*st%d%kweights(ik)*st%occ(ist, ik)*&
                    aimag(conjg(psi(ip, idim))*hrpsi(ip, idim) - conjg(psi(ip, idim))*rhpsi(ip, idim))
                end do
                !$omp end parallel do
              end do
            end do
            
          end do

          call batch_unpack(st%group%psib(ib, ik), copy = .false.)

          call batch_end(hpsib)
          call batch_end(rhpsib)
          call batch_end(rpsib)
          call batch_end(hrpsib)

        end do
      end do
    
    else

      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          
          call states_get_state(st, gr%mesh, ist, ik, psi)
          
          do idim = 1, st%d%dim
            call zderivatives_set_bc(gr%der, psi(:, idim))
          end do

          do idim = 1, st%d%dim

            ! Apply the phase that contains both the k-point and vector-potential terms.
            !$omp parallel do
            do ip = 1, gr%mesh%np_part
              psi(ip, idim) = phases(ip, ik - st%d%kpt%start + 1)*psi(ip, idim)
            end do

            call zderivatives_grad(gr%der, psi(:, idim), gpsi(:, :, idim), set_bc = .false.)
            
          end do
          
          do idir = 1, gr%sb%dim
            do iatom = 1, geo%natoms
              if(species_is_ps(geo%atom(iatom)%spec)) then
                call zprojector_commute_r(pj(iatom), gr, st%d%dim, idir, ik, psi, gpsi(:, idir, :))
              end if
            end do
          end do
          
          do idir = 1, gr%sb%dim
            
            do idim = 1, st%d%dim
              !$omp parallel do
              do ip = 1, gr%mesh%np
                microcurrent(ip, idir) = microcurrent(ip, idir) + &
                  M_FOUR*M_PI*P_c/gr%sb%rcell_volume*st%d%kweights(ik)*st%occ(ist, ik)*&
                  aimag(conjg(psi(ip, idim))*gpsi(ip, idir, idim))
              end do
              !$omp end parallel do
            end do
          end do

        end do
      end do
      
    end if
    
    if(st%symmetrize_density) then
      SAFE_ALLOCATE(symmcurrent(1:gr%mesh%np, 1:gr%sb%dim))
      call symmetrizer_init(symmetrizer, gr%mesh)
      call dsymmetrizer_apply(symmetrizer, field_vector = microcurrent, symmfield_vector = symmcurrent)
      microcurrent(1:gr%mesh%np, 1:gr%sb%dim) = symmcurrent(1:gr%mesh%np, 1:gr%sb%dim)
      call symmetrizer_end(symmetrizer)
      SAFE_DEALLOCATE_A(symmcurrent)
    end if

    do idir = 1, gr%sb%dim
      force%vecpot(idir) = dmf_integrate(gr%mesh, microcurrent(:, idir))
    end do

#ifdef HAVE_MPI
    if(st%parallel_in_states) then
      call MPI_Allreduce(force%vecpot, force_tmp, MAX_DIM, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
      force%vecpot = force_tmp
    end if
    if(st%d%kpt%parallel) then
      call MPI_Allreduce(force%vecpot, force_tmp, MAX_DIM, MPI_FLOAT, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
      force%vecpot = force_tmp
    end if
#endif

    ! The line below should not be added: since the vector potential
    ! is applied as a phase to the states, this term appears
    ! automatically. I keep it with this comment to alert possible
    ! readers of the code who might think that the term is missing.
    !
    !    force%vecpot(1:MAX_DIM) = force%vecpot(1:MAX_DIM) - this%wp2*this%vecpot(1:MAX_DIM)

    SAFE_DEALLOCATE_A(gpsi)

    call profiling_out(prof)
    POP_SUB(gauge_field_get_force)
  end subroutine gauge_field_get_force

end module gauge_field_force_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

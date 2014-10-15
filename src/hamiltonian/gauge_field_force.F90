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
  use datasets_m
  use derivatives_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
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

    integer :: ik, ist, idir, idim, iatom, ip
    CMPLX, allocatable :: gpsi(:, :, :), epsi(:, :)
    FLOAT, allocatable :: microcurrent(:, :), symmcurrent(:, :)
    type(profile_t), save :: prof
    type(symmetrizer_t) :: symmetrizer
#ifdef HAVE_MPI
    FLOAT :: force_tmp(1:MAX_DIM)
#endif

    call profiling_in(prof, "GAUGE_FIELD_FORCE")
    PUSH_SUB(gauge_field_get_force)

    SAFE_ALLOCATE(epsi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gpsi(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:st%d%dim))
    SAFE_ALLOCATE(microcurrent(1:gr%mesh%np_part, 1:gr%sb%dim))

    microcurrent = M_ZERO
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end

        call states_get_state(st, gr%mesh, ist, ik, epsi)

        do idim = 1, st%d%dim
          call zderivatives_set_bc(gr%der, epsi(:, idim))

          ! Apply the phase that contains both the k-point and vector-potential terms.
          forall(ip = 1:gr%mesh%np_part)
            epsi(ip, idim) = phases(ip, ik - st%d%kpt%start + 1)*epsi(ip, idim)
          end forall

          call zderivatives_grad(gr%der, epsi(:, idim), gpsi(:, :, idim), set_bc = .false.)

        end do

        do idir = 1, gr%sb%dim
          do iatom = 1, geo%natoms
            if(species_is_ps(geo%atom(iatom)%spec)) then
              call zprojector_commute_r(pj(iatom), gr, st%d%dim, idir, ik, epsi, gpsi(:, idir, :))
            end if
          end do
        end do
        
        do idir = 1, gr%sb%dim
          do idim = 1, st%d%dim
            microcurrent(1:gr%mesh%np, idir) = microcurrent(1:gr%mesh%np, idir) + &
              M_FOUR*M_PI*P_c/gr%sb%rcell_volume*st%d%kweights(ik)*st%occ(ist, ik)*&
              aimag(conjg(epsi(1:gr%mesh%np, idim))*gpsi(1:gr%mesh%np, idir, idim))
          end do
        end do
        
      end do
    end do

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

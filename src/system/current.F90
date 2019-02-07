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

#include "global.h"

module current_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use boundaries_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use hamiltonian_base_oct_m
  use io_oct_m
  use io_function_oct_m
  use lalg_basic_oct_m
  use logrid_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use projector_oct_m
  use ps_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use species_oct_m
  use splines_oct_m
  use states_oct_m
  use states_dim_oct_m
  use submesh_oct_m
  use symmetries_oct_m
  use symmetrizer_oct_m
  use symm_op_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use xc_oct_m  

  implicit none

  private

  type current_t
    integer :: method
  end type current_t
    

  public ::                               &
    current_t,                            &
    current_init,                         &
    current_end,                          &
    current_calculate

contains

  subroutine current_init(this)
    type(current_t), intent(out)   :: this

    PUSH_SUB(current_init)
    
    POP_SUB(current_init)
  end subroutine current_init

  ! ---------------------------------------------------------

  subroutine current_end(this)
    type(current_t), intent(inout) :: this

    PUSH_SUB(current_end)

    POP_SUB(current_end)
  end subroutine current_end

  ! ---------------------------------------------------------
  subroutine current_calculate(this, der, hm, geo, st, current, current_kpt)
    type(current_t),      intent(in)    :: this
    type(derivatives_t),  intent(inout) :: der
    type(hamiltonian_t),  intent(in)    :: hm
    type(geometry_t),     intent(in)    :: geo
    type(states_t),       intent(inout) :: st
    FLOAT,                intent(out)   :: current(:, :, :) !< current(1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:st%d%nspin)
    FLOAT, pointer,       intent(inout) :: current_kpt(:, :, :) !< current(1:der%mesh%np_part, 1:der%mesh%sb%dim, kpt%start:kpt%end)

    integer :: ik, ist, idir, idim, ip, ib, ii, ispin
    CMPLX, allocatable :: gpsi(:, :, :), psi(:, :), hpsi(:, :), rhpsi(:, :), rpsi(:, :), hrpsi(:, :)
    FLOAT, allocatable :: symmcurrent(:, :)
    type(profile_t), save :: prof
    type(symmetrizer_t) :: symmetrizer
    type(batch_t) :: hpsib, rhpsib, rpsib, hrpsib, epsib
    type(batch_t), allocatable :: commpsib(:)
    logical, parameter :: hamiltonian_current = .false.
    FLOAT :: ww
    CMPLX :: c_tmp

    call profiling_in(prof, "CURRENT")
    PUSH_SUB(current_calculate)

    ASSERT(all(ubound(current) == (/der%mesh%np_part, der%mesh%sb%dim, st%d%nspin/)))
    ASSERT(all(ubound(current_kpt) == (/der%mesh%np_part, der%mesh%sb%dim, st%d%kpt%end/)))
    ASSERT(all(lbound(current_kpt) == (/1, 1, st%d%kpt%start/)))

    SAFE_ALLOCATE(psi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gpsi(1:der%mesh%np, 1:der%mesh%sb%dim, 1:st%d%dim))
    SAFE_ALLOCATE(hpsi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(rhpsi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(rpsi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(hrpsi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(commpsib(1:der%mesh%sb%dim))

    current = M_ZERO
    current_kpt = M_ZERO

    do ik = st%d%kpt%start, st%d%kpt%end
      ispin = states_dim_get_spin_index(st%d, ik)
      do ib = st%group%block_start, st%group%block_end

        call batch_pack(st%group%psib(ib, ik), copy = .true.)
        call batch_copy(st%group%psib(ib, ik), epsib, fill_zeros = .false.)
        call boundaries_set(der%boundaries, st%group%psib(ib, ik))

        if(associated(hm%hm_base%phase)) then
          call zhamiltonian_base_phase(hm%hm_base, der, der%mesh%np_part, ik, &
            conjugate = .false., psib = epsib, src = st%group%psib(ib, ik))
        else
          call batch_copy_data(der%mesh%np_part, st%group%psib(ib, ik), epsib)
        end if

        do idir = 1, der%mesh%sb%dim
          call batch_copy(st%group%psib(ib, ik), commpsib(idir))
          call zderivatives_batch_perform(der%grad(idir), der, epsib, commpsib(idir), set_bc = .false.)
        end do

        call zhamiltonian_base_nlocal_position_commutator(hm%hm_base, der%mesh, st%d, ik, epsib, commpsib)

        do idir = 1, der%mesh%sb%dim

          if(associated(hm%hm_base%phase)) then
            call zhamiltonian_base_phase(hm%hm_base, der, der%mesh%np_part, ik, conjugate = .true., psib = commpsib(idir))
          end if

          do ist = states_block_min(st, ib), states_block_max(st, ib)

            do idim = 1, st%d%dim
              ii = batch_inv_index(st%group%psib(ib, ik), (/ist, idim/))
              call batch_get_state(st%group%psib(ib, ik), ii, der%mesh%np, psi(:, idim))
              call batch_get_state(commpsib(idir), ii, der%mesh%np, hrpsi(:, idim))
            end do

            ww = st%d%kweights(ik)*st%occ(ist, ik) 
            if(st%d%ispin /= SPINORS) then
              !$omp parallel do
              do ip = 1, der%mesh%np
                current_kpt(ip, idir, ik) = &
                  current_kpt(ip, idir, ik) + ww*aimag(conjg(psi(ip, 1))*hrpsi(ip, 1))
              end do
              !$omp end parallel do
            else
              !$omp parallel do private(c_tmp)
              do ip = 1, der%mesh%np
                current(ip, idir, 1) = current(ip, idir, 1) + &
                  ww*aimag(conjg(psi(ip, 1))*hrpsi(ip, 1))
                current(ip, idir, 2) = current(ip, idir, 2) + &
                  ww*aimag(conjg(psi(ip, 2))*hrpsi(ip, 2))
                c_tmp = conjg(psi(ip, 1))*hrpsi(ip, 2) - psi(ip, 2)*conjg(hrpsi(ip, 1))
                current(ip, idir, 3) = current(ip, idir, 3) + ww* real(c_tmp)
                current(ip, idir, 4) = current(ip, idir, 4) + ww*aimag(c_tmp)
              end do
              !$omp end parallel do
            end if

          end do

          call batch_end(commpsib(idir))

        end do

        call batch_end(epsib)
        call batch_unpack(st%group%psib(ib, ik), copy = .false.)

      end do
    end do

    !We sum the current over k-points
    do ik = st%d%kpt%start, st%d%kpt%end
      ispin = states_dim_get_spin_index(st%d, ik)
      do idir = 1, der%mesh%sb%dim
        !$omp parallel do
        do ip = 1, der%mesh%np
          current(ip, idir, ispin) = current(ip, idir, ispin) + current_kpt(ip, idir, ik)
        end do
        !$omp end parallel do
      end do
    end do

    if(st%parallel_in_states .or. st%d%kpt%parallel) then
      ! TODO: this could take dim = (/der%mesh%np, der%mesh%sb%dim, st%d%nspin/)) to reduce the amount of data copied
      call comm_allreduce(st%st_kpt_mpi_grp%comm, current) 
    end if

    if(st%symmetrize_density) then
      SAFE_ALLOCATE(symmcurrent(1:der%mesh%np, 1:der%mesh%sb%dim))
      call symmetrizer_init(symmetrizer, der%mesh)
      do ispin = 1, st%d%nspin
        call dsymmetrizer_apply(symmetrizer, der%mesh%np, field_vector = current(:, :, ispin), &
          symmfield_vector = symmcurrent, suppress_warning = .true.)
        current(1:der%mesh%np, 1:der%mesh%sb%dim, ispin) = symmcurrent(1:der%mesh%np, 1:der%mesh%sb%dim)
      end do
      call symmetrizer_end(symmetrizer)
      SAFE_DEALLOCATE_A(symmcurrent)
    end if

    SAFE_DEALLOCATE_A(gpsi)
    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(rhpsi)
    SAFE_DEALLOCATE_A(rpsi)
    SAFE_DEALLOCATE_A(hrpsi)
    SAFE_DEALLOCATE_A(commpsib)

    call profiling_out(prof)
    POP_SUB(current_calculate)

  end subroutine current_calculate

end module current_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

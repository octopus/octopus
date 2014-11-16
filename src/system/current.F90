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

module current_m
  use batch_m
  use batch_ops_m
  use comm_m
  use datasets_m
  use derivatives_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use hamiltonian_base_m
  use io_m
  use io_function_m
  use lalg_basic_m
  use logrid_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use parser_m
  use poisson_m
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

  type current_t
    integer :: method
  end type current_t
    

  public ::                               &
    current_t,                            &
    current_init,                         &
    current_end,                          &
    current_calculate

  integer, parameter, public ::           &
    CURRENT_GRADIENT           = 1,       &
    CURRENT_HAMILTONIAN        = 2,       &
    CURRENT_POISSON            = 3,       &
    CURRENT_POISSON_CORRECTION = 4

contains

  subroutine current_init(this)
    type(current_t), intent(out)   :: this

    PUSH_SUB(current_init)

    !%Variable CurrentDensity
    !%Default gradient
    !%Type integer
    !%Section Hamiltonian
    !%Description
    !% This variable selects the method used to
    !% calculate the current density. For the moment this variable is
    !% for development purposes and users should not need to use
    !% it.
    !%Option gradient 1
    !% The calculation of current is done using the gradient operator
    !% with additional corrections for non-local operators.
    !%Option hamiltonian 2
    !% The current density is obtained from the commutator of the
    !% Hamiltonian with the position operator. (Experimental)
    !%Option poisson 3
    !% Obtain the current from solving the Poisson equation from the continuity equation. (Experimental)
    !%Option poisson_correction 4
    !% Obtain the current from the hamiltonian and then add a correction term by solving the Poisson equation. (Experimental)
    !%End

    call parse_integer(datasets_check('CurrentDensity'), CURRENT_GRADIENT, this%method)
    if(.not.varinfo_valid_option('CurrentDensity', this%method)) call input_error('CurrentDensity')
    if(this%method /= CURRENT_GRADIENT) &
      call messages_experimental("CurrentDensity /= gradient")

    POP_SUB(current_init)
  end subroutine current_init

  ! ---------------------------------------------------------

  subroutine current_end(this)
    type(current_t), intent(inout) :: this

    PUSH_SUB(current_end)


    POP_SUB(current_end)
  end subroutine current_end

  ! ---------------------------------------------------------
  subroutine current_calculate(this, gr, hm, geo, st, current)
    type(current_t),      intent(in)    :: this
    type(grid_t),         intent(inout) :: gr
    type(hamiltonian_t),  intent(in)    :: hm
    type(geometry_t),     intent(in)    :: geo
    type(states_t),       intent(inout) :: st
    FLOAT,                intent(out)    :: current(:, :, :) !< current(1:gr%mesh%np_part, 1:gr%sb%dim, 1:st%d%nspin)

    integer :: ik, ist, idir, idim, iatom, ip, ib, minst, maxst, ii, ierr
    CMPLX, allocatable :: gpsi(:, :, :), psi(:, :), hpsi(:, :), rhpsi(:, :), rpsi(:, :), hrpsi(:, :)
    FLOAT, allocatable :: symmcurrent(:, :)
    type(profile_t), save :: prof
    type(symmetrizer_t) :: symmetrizer
    type(batch_t) :: hpsib, rhpsib, rpsib, hrpsib
    logical, parameter :: hamiltonian_current = .false.

    call profiling_in(prof, "CURRENT")
    PUSH_SUB(current_calculate)

    ! spin not implemented or tested
    ASSERT(st%d%nspin == 1)
    ASSERT(all(ubound(current) == (/gr%mesh%np_part, gr%sb%dim, st%d%nspin/)))

    SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gpsi(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:st%d%dim))
    SAFE_ALLOCATE(hpsi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(rhpsi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(rpsi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(hrpsi(1:gr%mesh%np_part, 1:st%d%dim))

    current = M_ZERO

    select case(this%method)
    case(CURRENT_HAMILTONIAN)

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
                  current(ip, idir, 1) = current(ip, idir, 1) &
                    - st%d%kweights(ik)*st%occ(ist, ik)&
                    *aimag(conjg(psi(ip, idim))*hrpsi(ip, idim) - conjg(psi(ip, idim))*rhpsi(ip, idim))
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
    
    case(CURRENT_GRADIENT)

      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          
          call states_get_state(st, gr%mesh, ist, ik, psi)
          
          do idim = 1, st%d%dim
            call zderivatives_set_bc(gr%der, psi(:, idim))
          end do

          if(associated(hm%phase)) then 
            ! Apply the phase that contains both the k-point and vector-potential terms.
            do idim = 1, st%d%dim
              !$omp parallel do
              do ip = 1, gr%mesh%np_part
                psi(ip, idim) = hm%phase(ip, ik)*psi(ip, idim)
              end do
              !$omp end parallel do
            end do
          end if

          do idim = 1, st%d%dim
            call zderivatives_grad(gr%der, psi(:, idim), gpsi(:, :, idim), set_bc = .false.)
          end do
          
          do idir = 1, gr%sb%dim
            do iatom = 1, geo%natoms
              if(species_is_ps(geo%atom(iatom)%spec)) then
                call zprojector_commute_r(hm%ep%proj(iatom), gr, st%d%dim, idir, ik, psi, gpsi(:, idir, :))
              end if
            end do
          end do
          
          do idir = 1, gr%sb%dim
            
            do idim = 1, st%d%dim
              !$omp parallel do
              do ip = 1, gr%mesh%np
                current(ip, idir, 1) = current(ip, idir, 1) + &
                  st%d%kweights(ik)*st%occ(ist, ik)*aimag(conjg(psi(ip, idim))*gpsi(ip, idir, idim))
              end do
              !$omp end parallel do
            end do
          end do

        end do
      end do
      
    case(CURRENT_POISSON)

      call calc_current_poisson()

    case(CURRENT_POISSON_CORRECTION)

      call calc_current_poisson_correction()

    case default

      ASSERT(.false.)

    end select

    if(st%parallel_in_states .or. st%d%kpt%parallel) then
      ! TODO: this could take dim = (/gr%mesh%np, gr%sb%dim, st%d%nspin/)) to reduce the amount of data copied
      call comm_allreduce(st%st_kpt_mpi_grp%comm, current) 
    end if
    
    if(st%symmetrize_density) then
      SAFE_ALLOCATE(symmcurrent(1:gr%mesh%np, 1:gr%sb%dim))
      call symmetrizer_init(symmetrizer, gr%mesh)
      call dsymmetrizer_apply(symmetrizer, field_vector = current(:, :, 1), symmfield_vector = symmcurrent)
      current(1:gr%mesh%np, 1:gr%sb%dim, 1) = symmcurrent(1:gr%mesh%np, 1:gr%sb%dim)
      call symmetrizer_end(symmetrizer)
      SAFE_DEALLOCATE_A(symmcurrent)
    end if

    SAFE_DEALLOCATE_A(gpsi)

    call profiling_out(prof)
    POP_SUB(current_calculate)

    contains
      
      subroutine calc_current_poisson()
        
        FLOAT, allocatable :: charge(:), potential(:)
        CMPLX, allocatable :: hpsi(:, :)

        PUSH_SUB(current_calculate.calc_current_poisson)
        
        SAFE_ALLOCATE(charge(1:gr%mesh%np))
        SAFE_ALLOCATE(potential(1:gr%mesh%np_part))
        SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim))

        ASSERT(st%d%dim == 1)

        do ik = st%d%kpt%start, st%d%kpt%end
          do ib = st%group%block_start, st%group%block_end
            
            call batch_pack(st%group%psib(ib, ik), copy = .true.)
            
            call batch_copy(st%group%psib(ib, ik), hpsib, reference = .false.)
            
            call zhamiltonian_apply_batch(hm, gr%der, st%group%psib(ib, ik), hpsib, ik)
            
            minst = states_block_min(st, ib)
            maxst = states_block_max(st, ib)
            
            do ist = st%st_start, st%st_end
              
              do idim = 1, st%d%dim
                ii = batch_ist_idim_to_linear(st%group%psib(ib, ik), (/ist, idim/))
                call batch_get_state(st%group%psib(ib, ik), ii, gr%mesh%np, psi(:, idim))
                call batch_get_state(hpsib, ii, gr%mesh%np, hpsi(:, idim))
              end do
              
              do idim = 1, st%d%dim
                !$omp parallel do
                do ip = 1, gr%mesh%np
                  charge(ip) = charge(ip) &
                    - st%d%kweights(ik)*st%occ(ist, ik)/(CNST(4.0)*M_PI)&
                    *aimag(psi(ip, idim)*conjg(hpsi(ip, idim)) - conjg(psi(ip, idim))*hpsi(ip, idim))
                end do
                !$omp end parallel do
              end do
            end do

          call batch_unpack(st%group%psib(ib, ik), copy = .false.)
          
          call batch_end(hpsib)
          
        end do
      end do

      call dpoisson_solve(psolver, potential, charge)

      call dio_function_output(C_OUTPUT_HOW_PLANE_X, "./continuity", "potential", gr%mesh, potential, unit_one, ierr)
      call dio_function_output(C_OUTPUT_HOW_PLANE_Z, "./continuity", "potential", gr%mesh, potential, unit_one, ierr)
      call dio_function_output(C_OUTPUT_HOW_AXIS_Z, "./continuity", "potential", gr%mesh, potential, unit_one, ierr)

      call dderivatives_grad(gr%der, potential, current(:, :, 1))

      POP_SUB(current_calculate.calc_current_poisson)
    end subroutine calc_current_poisson


    subroutine calc_current_poisson_correction()

      FLOAT, allocatable :: charge(:), potential(:), current2(:, :)
      type(batch_t) :: vpsib

      PUSH_SUB(current_calculate.calc_current_poisson_correction)
      
      SAFE_ALLOCATE(charge(1:gr%mesh%np))
      SAFE_ALLOCATE(potential(1:gr%mesh%np_part))
      SAFE_ALLOCATE(current2(1:gr%mesh%np_part, 1:gr%sb%dim))

      ASSERT(st%d%dim == 1)

      charge = CNST(0.0)
      current = CNST(0.0)

      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end

          call batch_pack(st%group%psib(ib, ik), copy = .true.)

          call batch_copy(st%group%psib(ib, ik), hpsib, reference = .false.)
          call batch_copy(st%group%psib(ib, ik), rhpsib, reference = .false.)
          call batch_copy(st%group%psib(ib, ik), rpsib, reference = .false.)
          call batch_copy(st%group%psib(ib, ik), hrpsib, reference = .false.)
          call batch_copy(st%group%psib(ib, ik), vpsib, reference = .false.)

          call zderivatives_batch_set_bc(gr%der, st%group%psib(ib, ik))

          call zhamiltonian_apply_batch(hm, gr%der, st%group%psib(ib, ik), hpsib, ik, set_bc = .false., &
            terms = TERM_KINETIC)
          call zhamiltonian_apply_batch(hm, gr%der, st%group%psib(ib, ik), vpsib, ik, set_bc = .false., &
            terms = TERM_NON_LOCAL_POTENTIAL)

          minst = states_block_min(st, ib)
          maxst = states_block_max(st, ib)
            
          do ist = st%st_start, st%st_end
            
            do idim = 1, st%d%dim
              ii = batch_ist_idim_to_linear(st%group%psib(ib, ik), (/ist, idim/))
              call batch_get_state(st%group%psib(ib, ik), ii, gr%mesh%np, psi(:, idim))
              call batch_get_state(vpsib, ii, gr%mesh%np, hpsi(:, idim))
            end do
            
            do idim = 1, st%d%dim
              !$omp parallel do
              do ip = 1, gr%mesh%np
                charge(ip) = charge(ip) &
                  - st%d%kweights(ik)*st%occ(ist, ik)/(CNST(4.0)*M_PI)&
                  *aimag(psi(ip, idim)*conjg(hpsi(ip, idim)) - conjg(psi(ip, idim))*hpsi(ip, idim))
              end do
              !$omp end parallel do
            end do
          end do
          
          do idir = 1, gr%sb%dim

            call batch_mul(gr%mesh%np, gr%mesh%x(:, idir), hpsib, rhpsib)
            call batch_mul(gr%mesh%np_part, gr%mesh%x(:, idir), st%group%psib(ib, ik), rpsib)

            call zhamiltonian_apply_batch(hm, gr%der, rpsib, hrpsib, ik, set_bc = .false.)

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
                  current(ip, idir, 1) = current(ip, idir, 1) &
                    - st%d%kweights(ik)*st%occ(ist, ik)&
                    *aimag(conjg(psi(ip, idim))*hrpsi(ip, idim) - conjg(psi(ip, idim))*rhpsi(ip, idim))
                end do
                !$omp end parallel do
              end do
            end do

          end do
          call batch_unpack(st%group%psib(ib, ik), copy = .false.)

          call batch_end(hpsib)
          call batch_end(vpsib)
          call batch_end(rhpsib)
          call batch_end(rpsib)
          call batch_end(hrpsib)

        end do
      end do

      call dpoisson_solve(psolver, potential, charge)

!      call dio_function_output(C_OUTPUT_HOW_PLANE_X, "./continuity", "potential", gr%mesh, potential, unit_one, ierr)
!      call dio_function_output(C_OUTPUT_HOW_PLANE_Z, "./continuity", "potential", gr%mesh, potential, unit_one, ierr)
!      call dio_function_output(C_OUTPUT_HOW_AXIS_Z, "./continuity", "potential", gr%mesh, potential, unit_one, ierr)

      call dderivatives_grad(gr%der, potential, current2)

      do idir = 1, gr%sb%dim
        do ip = 1, gr%mesh%np
          current(ip, idir, 1) = current(ip, idir, 1) + current2(ip, idir)
        end do
      end do

      POP_SUB(current_calculate.calc_current_poisson_correction)
    end subroutine calc_current_poisson_correction

  end subroutine current_calculate

end module current_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

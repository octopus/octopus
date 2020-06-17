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

module system_oct_m
  use accel_oct_m
  use calc_mode_par_oct_m
  use density_oct_m
  use elf_oct_m
  use energy_calc_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use mesh_oct_m
  use messages_oct_m
  use modelmb_particles_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use space_oct_m
  use simul_box_oct_m
  use sort_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use v_ks_oct_m
  use xc_oct_m
  use xc_oep_oct_m

  implicit none

  private
  public ::               &
    electrons_t,          &
    system_init,          &
    system_h_setup

  type :: electrons_t
    ! Components are public by default
    type(space_t)                :: space
    type(geometry_t)             :: geo
    type(grid_t),        pointer :: gr    !< the mesh
    type(states_elec_t), pointer :: st    !< the states
    type(v_ks_t)                 :: ks    !< the Kohn-Sham potentials
    type(output_t)               :: outp  !< the output
    type(multicomm_t)            :: mc    !< index and domain communicators
    type(namespace_t)            :: namespace
    type(hamiltonian_elec_t)     :: hm
  contains
    final :: system_finalize
  end type electrons_t
  
contains
  
  !----------------------------------------------------------
  function system_init(namespace) result(sys)
    class(electrons_t), pointer    :: sys
    type(namespace_t),  intent(in) :: namespace

    type(profile_t), save :: prof

    PUSH_SUB(system_init)
    call profiling_in(prof,"SYSTEM_INIT")

    SAFE_ALLOCATE(sys)

    SAFE_ALLOCATE(sys%gr)
    SAFE_ALLOCATE(sys%st)

    sys%namespace = namespace

    call messages_obsolete_variable(sys%namespace, 'SystemName')

    call space_init(sys%space, sys%namespace)
    
    call geometry_init(sys%geo, sys%namespace, sys%space)
    call grid_init_stage_0(sys%gr, sys%namespace, sys%geo, sys%space)
    call states_elec_init(sys%st, sys%namespace, sys%gr, sys%geo)
    call sys%st%write_info(sys%namespace)
    call grid_init_stage_1(sys%gr, sys%namespace, sys%geo)
    ! if independent particles in N dimensions are being used, need to initialize them
    !  after masses are set to 1 in grid_init_stage_1 -> derivatives_init
    call modelmb_copy_masses (sys%st%modelmbparticles, sys%gr%der%masses)

    call parallel_init()

    call geometry_partition(sys%geo, sys%mc)
    call kpoints_distribute(sys%st%d, sys%mc)
    call states_elec_distribute_nodes(sys%st, sys%namespace, sys%mc)
    call grid_init_stage_2(sys%gr, sys%namespace, sys%mc, sys%geo)
    if(sys%st%symmetrize_density) call mesh_check_symmetries(sys%gr%mesh, sys%gr%sb)

    call v_ks_nullify(sys%ks)
    call output_init(sys%outp, sys%namespace, sys%gr%sb, sys%st, sys%st%nst, sys%ks)
    call states_elec_densities_init(sys%st, sys%gr, sys%geo)
    call states_elec_exec_init(sys%st, sys%namespace, sys%mc)
    call elf_init(sys%namespace)

    call v_ks_init(sys%ks, sys%namespace, sys%gr, sys%st, sys%geo, sys%mc)

    call hamiltonian_elec_init(sys%hm, sys%namespace, sys%gr, sys%geo, sys%st, sys%ks%theory_level, &
      sys%ks%xc, sys%mc, need_exchange = output_need_exchange(sys%outp) .or. sys%ks%oep%level /= XC_OEP_NONE)
    
    if(poisson_is_multigrid(sys%hm%psolver)) call grid_create_multigrid(sys%gr, sys%namespace, sys%geo, sys%mc)

    if (sys%hm%pcm%run_pcm .and. sys%mc%par_strategy /= P_STRATEGY_SERIAL .and. sys%mc%par_strategy /= P_STRATEGY_STATES) then
      call messages_experimental('Parallel in domain calculations with PCM')
    end if

    call profiling_out(prof)
    POP_SUB(system_init)

  contains

    ! ---------------------------------------------------------
    subroutine parallel_init()
      integer :: index_range(4)

      PUSH_SUB(system_init.parallel_init)

      ! store the ranges for these two indices (serves as initial guess
      ! for parallelization strategy)
      index_range(1) = sys%gr%mesh%np_global  ! Number of points in mesh
      index_range(2) = sys%st%nst             ! Number of states
      index_range(3) = sys%st%d%nik           ! Number of k-points
      index_range(4) = 100000                 ! Some large number

      ! create index and domain communicators
      call multicomm_init(sys%mc, sys%namespace, mpi_world, calc_mode_par_parallel_mask(), calc_mode_par_default_parallel_mask(), &
        mpi_world%size, index_range, (/ 5000, 1, 1, 1 /))

      POP_SUB(system_init.parallel_init)
    end subroutine parallel_init

  end function system_init

  !----------------------------------------------------------
  subroutine system_finalize(sys)
    type(electrons_t), intent(inout) :: sys

    PUSH_SUB(system_finalize)

    call hamiltonian_elec_end(sys%hm)

    call multicomm_end(sys%mc)

    call v_ks_end(sys%ks, sys%gr)
    
    call output_end(sys%outp)
    
    if(associated(sys%st)) then
      call states_elec_end(sys%st)
      SAFE_DEALLOCATE_P(sys%st)
    end if

    call geometry_end(sys%geo)
    call simul_box_end(sys%gr%sb)
    call grid_end(sys%gr)

    call space_end(sys%space)

    SAFE_DEALLOCATE_P(sys%gr)

    POP_SUB(system_finalize)
  end subroutine system_finalize


  !----------------------------------------------------------
  subroutine system_h_setup(sys, calc_eigenval, calc_current)
    type(electrons_t), intent(inout) :: sys
    logical, optional, intent(in)    :: calc_eigenval !< default is true
    logical, optional, intent(in)    :: calc_current !< default is true

    integer, allocatable :: ind(:)
    integer :: ist, ik
    FLOAT, allocatable :: copy_occ(:)
    logical :: calc_eigenval_
    logical :: calc_current_

    PUSH_SUB(system_h_setup)

    calc_eigenval_ = optional_default(calc_eigenval, .true.)
    calc_current_ = optional_default(calc_current, .true.)
    call states_elec_fermi(sys%st, sys%namespace, sys%gr%mesh)
    call density_calc(sys%st, sys%gr, sys%st%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%hm, sys%st, sys%geo, calc_eigenval = calc_eigenval_, calc_current = calc_current_) ! get potentials

    if(sys%st%restart_reorder_occs .and. .not. sys%st%fromScratch) then
      message(1) = "Reordering occupations for restart."
      call messages_info(1)

      SAFE_ALLOCATE(ind(1:sys%st%nst))
      SAFE_ALLOCATE(copy_occ(1:sys%st%nst))

      do ik = 1, sys%st%d%nik
        call sort(sys%st%eigenval(:, ik), ind)
        copy_occ(1:sys%st%nst) = sys%st%occ(1:sys%st%nst, ik)
        do ist = 1, sys%st%nst
          sys%st%occ(ist, ik) = copy_occ(ind(ist))
        end do
      end do

      SAFE_DEALLOCATE_A(ind)
      SAFE_DEALLOCATE_A(copy_occ)
    end if

    if(calc_eigenval_) call states_elec_fermi(sys%st, sys%namespace, sys%gr%mesh) ! occupations
    call energy_calc_total(sys%namespace, sys%hm, sys%gr, sys%st)

    POP_SUB(system_h_setup)
  end subroutine system_h_setup

end module system_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

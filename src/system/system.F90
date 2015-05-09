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
!! $Id$

#include "global.h"

module system_m
  use calc_mode_par_m
  use density_m
  use elf_m
  use energy_calc_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_function_m
  use json_m
  use mesh_m
  use messages_m
  use modelmb_particles_m
  use mpi_m
  use multicomm_m
  use octcl_kernel_m
  use opencl_m
  use output_m
  use parser_m
  use pcm_m
  use poisson_m
  use profiling_m
  use space_m
  use species_m
  use simul_box_m
  use sort_om
  use ssys_config_m
  use ssys_handle_m
  use ssys_model_m
  use ssys_states_m
  use ssys_system_m
  use states_m
  use states_dim_m
  use unit_m
  use unit_system_m
  use v_ks_m

  implicit none

  private
  public ::               &
    system_t,             &
    system_init,          &
    system_end,           &
    system_h_setup

  type system_t
    type(space_t)               :: space
    type(geometry_t)            :: geo
    type(grid_t),       pointer :: gr    !< the mesh
    type(states_t),     pointer :: st    !< the states
    type(ssys_model_t), pointer :: subsys_model !< the Subsystems model pointer.
    type(v_ks_t)                :: ks    !< the Kohn-Sham potentials
    type(output_t)              :: outp  !< the output
    type(multicomm_t)           :: mc    !< index and domain communicators
  end type system_t

contains

  !----------------------------------------------------------
  subroutine system_init(sys, subsys_handle, config)
    type(system_t),                intent(out)   :: sys
    type(ssys_handle_t), optional, intent(inout) :: subsys_handle
    type(json_object_t), optional, intent(inout) :: config

    type(ssys_system_t), pointer :: subsys_system
    type(ssys_states_t), pointer :: subsys_states

    type(profile_t), save :: prof
    PUSH_SUB(system_init)
    call profiling_in(prof,"SYSTEM_INIT")
    
    SAFE_ALLOCATE(sys%gr)
    SAFE_ALLOCATE(sys%st)

    call opencl_init(mpi_world)
    call octcl_kernel_global_init()

    call messages_obsolete_variable('SystemName')

    call space_init(sys%space)
    
    call geometry_init(sys%geo, sys%space)
    call grid_init_stage_0(sys%gr, sys%geo, sys%space)
    call states_init(sys%st, sys%gr, sys%geo)
    call states_write_info(sys%st)
    call grid_init_stage_1(sys%gr, sys%geo)
    ! if independent particles in N dimensions are being used, need to initialize them
    !  after masses are set to 1 in grid_init_stage_1 -> derivatives_init
    call modelmb_copy_masses (sys%st%modelmbparticles, sys%gr%der%masses)

    call parallel_init()

    call geometry_partition(sys%geo, sys%mc)
    call kpoints_distribute(sys%st%d, sys%mc)
    call states_distribute_nodes(sys%st, sys%mc)
    call grid_init_stage_2(sys%gr, sys%mc, sys%geo)
    call output_init(sys%outp, sys%gr%sb, sys%st%nst)

    nullify(sys%subsys_model, subsys_system, subsys_states)
    if(present(subsys_handle).and.present(config))then
      call ssys_config_parse(config, sys%geo, sys%space%dim, sys%st%d%nspin)
      call ssys_handle_init(subsys_handle, config)
      call ssys_handle_start(subsys_handle, sys%gr)
      call ssys_handle_get(subsys_handle, sys%subsys_model)
      ASSERT(associated(sys%subsys_model))
      call ssys_model_get(sys%subsys_model, subsys_system)
      ASSERT(associated(subsys_system))
      call ssys_system_get(subsys_system, subsys_states)
      ASSERT(associated(subsys_states))
      call states_add_substates(sys%st, subsys_states)
      nullify(subsys_system, subsys_states)
    end if

    call states_densities_init(sys%st, sys%gr, sys%geo)
    call states_exec_init(sys%st, sys%mc)
    call elf_init()

    call poisson_init(psolver, sys%gr%der, sys%mc, theta = sys%st%cmplxscl%theta)
    if(poisson_is_multigrid(psolver)) call grid_create_multigrid(sys%gr, sys%geo)

    call v_ks_init(sys%ks, sys%gr, sys%st, sys%geo, sys%mc)

    call profiling_out(prof)
    POP_SUB(system_init)

  contains

    ! ---------------------------------------------------------
    subroutine parallel_init()
      integer :: index_dim, index_range(4)

      PUSH_SUB(system_init.parallel_init)

      ! for the moment we need communicators for domains plus three
      ! index-dimensions: states, k-points, and others
      index_dim = 4

      ! store the ranges for these two indices (serves as initial guess
      ! for parallelization strategy)
      index_range(1) = sys%gr%mesh%np_global  ! Number of points in mesh
      index_range(2) = sys%st%nst             ! Number of states
      index_range(3) = sys%st%d%nik           ! Number of k-points
      index_range(4) = 100000                 ! Some large number

      ! create index and domain communicators
      call multicomm_init(sys%mc, mpi_world, calc_mode_par_parallel_mask(), calc_mode_par_default_parallel_mask(), &
        mpi_world%size, index_dim, index_range, (/ 5000, 1, 1, 1 /))

      POP_SUB(system_init.parallel_init)
    end subroutine parallel_init

  end subroutine system_init


  !----------------------------------------------------------
  subroutine system_end(sys)
    type(system_t), intent(inout) :: sys

    PUSH_SUB(system_end)

    call multicomm_end(sys%mc)

    call poisson_end(psolver)
    call v_ks_end(sys%ks)
    
    call output_end(sys%outp)
    
    nullify(sys%subsys_model)

    if(associated(sys%st)) then
      call states_end(sys%st)
      SAFE_DEALLOCATE_P(sys%st)
    end if

    call geometry_end(sys%geo)
    call simul_box_end(sys%gr%sb)
    call grid_end(sys%gr)

    call space_end(sys%space)

    call octcl_kernel_global_end()
    call opencl_end()

    SAFE_DEALLOCATE_P(sys%gr)

    POP_SUB(system_end)
  end subroutine system_end


  !----------------------------------------------------------
  subroutine system_h_setup(sys, hm, calc_eigenval)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    logical,   optional, intent(in)    :: calc_eigenval !< default is true

    integer, allocatable :: ind(:)
    integer :: ist, ik
    FLOAT, allocatable :: copy_occ(:)
    logical :: calc_eigenval_

    PUSH_SUB(system_h_setup)

    calc_eigenval_ = optional_default(calc_eigenval, .true.)
    call states_fermi(sys%st, sys%gr%mesh)
    if(hm%cmplxscl%space) then
      call density_calc(sys%st, sys%gr, sys%st%zrho%Re, sys%st%zrho%Im)
    else
      call density_calc(sys%st, sys%gr, sys%st%rho)
    end if
    call v_ks_calc(sys%ks, hm, sys%st, sys%geo, calc_eigenval = calc_eigenval_) ! get potentials

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
        enddo
      enddo

      SAFE_DEALLOCATE_A(ind)
      SAFE_DEALLOCATE_A(copy_occ)
    endif

    call states_fermi(sys%st, sys%gr%mesh) ! occupations
    call energy_calc_total(hm, sys%gr, sys%st)

    POP_SUB(system_h_setup)
  end subroutine system_h_setup

end module system_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

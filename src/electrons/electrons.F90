!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2009 X. Andrade
!! Copyright (C) 2020 M. Oliveira
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

module electrons_oct_m
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
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use unit_system_oct_m
  use v_ks_oct_m
  use xc_oct_m
  use xc_oep_oct_m

  implicit none

  private
  public ::               &
    electrons_t

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
    final :: electrons_finalize
  end type electrons_t
  
  interface electrons_t
    procedure electrons_constructor
  end interface electrons_t

contains
  
  !----------------------------------------------------------
  function electrons_constructor(namespace, generate_epot) result(sys)
    class(electrons_t), pointer    :: sys
    type(namespace_t),  intent(in) :: namespace
    logical,  optional, intent(in) :: generate_epot

    type(profile_t), save :: prof
    FLOAT :: mesh_global, mesh_local, wfns

    PUSH_SUB(electrons_constructor)
    call profiling_in(prof,"ELECTRONS_CONSTRUCTOR")

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

    ! Print memory requirements
    call messages_print_stress(stdout, 'Approximate memory requirements', namespace=sys%namespace)

    mesh_global = mesh_global_memory(sys%gr%mesh)
    mesh_local  = mesh_local_memory(sys%gr%mesh)

    call messages_write('Mesh')
    call messages_new_line()
    call messages_write('  global  :')
    call messages_write(mesh_global, units = unit_megabytes, fmt = '(f10.1)')
    call messages_new_line()
    call messages_write('  local   :')
    call messages_write(mesh_local, units = unit_megabytes, fmt = '(f10.1)')
    call messages_new_line()
    call messages_write('  total   :')
    call messages_write(mesh_global + mesh_local, units = unit_megabytes, fmt = '(f10.1)')
    call messages_new_line()
    call messages_info()

    wfns = states_elec_wfns_memory(sys%st, sys%gr%mesh)
    call messages_write('States')
    call messages_new_line()
    call messages_write('  real    :')
    call messages_write(wfns, units = unit_megabytes, fmt = '(f10.1)')
    call messages_write(' (par_kpoints + par_states + par_domains)')
    call messages_new_line()
    call messages_write('  complex :')
    call messages_write(2.0_8*wfns, units = unit_megabytes, fmt = '(f10.1)')
    call messages_write(' (par_kpoints + par_states + par_domains)')
    call messages_new_line()
    call messages_info()

    call messages_print_stress(stdout)

    if (optional_default(generate_epot, .false.)) then
      message(1) = "Info: Generating external potential"
      call messages_info(1)
      call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, sys%gr, sys%geo, sys%st)
      message(1) = "      done."
      call messages_info(1)
    end if

    call profiling_out(prof)
    POP_SUB(electrons_constructor)
  contains

    ! ---------------------------------------------------------
    subroutine parallel_init()
      integer :: index_range(4)

      PUSH_SUB(electrons_constructor.parallel_init)

      ! store the ranges for these two indices (serves as initial guess
      ! for parallelization strategy)
      index_range(1) = sys%gr%mesh%np_global  ! Number of points in mesh
      index_range(2) = sys%st%nst             ! Number of states
      index_range(3) = sys%st%d%nik           ! Number of k-points
      index_range(4) = 100000                 ! Some large number

      ! create index and domain communicators
      call multicomm_init(sys%mc, sys%namespace, mpi_world, calc_mode_par_parallel_mask(), calc_mode_par_default_parallel_mask(), &
        mpi_world%size, index_range, (/ 5000, 1, 1, 1 /))

      POP_SUB(electrons_constructor.parallel_init)
    end subroutine parallel_init

  end function electrons_constructor

  !----------------------------------------------------------
  subroutine electrons_finalize(sys)
    type(electrons_t), intent(inout) :: sys

    PUSH_SUB(electrons_finalize)

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

    POP_SUB(electrons_finalize)
  end subroutine electrons_finalize

end module electrons_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

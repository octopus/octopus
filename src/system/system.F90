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

module system_m
  use calc_mode_m
  use density_m
  use elf_m
  use energy_m
  use geometry_m
  use global_m
  use grid_m
  use output_m
  use hamiltonian_m
  use io_function_m
  use mesh_m
  use messages_m
  use modelmb_particles_m
  use mpi_m
  use multicomm_m
  use opencl_m
  use poisson_m
  use profiling_m
  use space_m
  use simul_box_m
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
    type(space_t)             :: space
    type(geometry_t)          :: geo
    type(grid_t),     pointer :: gr    !< the mesh
    type(states_t),   pointer :: st    !< the states
    type(v_ks_t)              :: ks    !< the Kohn-Sham potentials
    type(output_t)      :: outp  !< the output
    type(multicomm_t)         :: mc    !< index and domain communicators
  end type system_t

contains

  !----------------------------------------------------------
  subroutine system_init(sys)
    type(system_t), intent(out) :: sys

    PUSH_SUB(system_init)

    SAFE_ALLOCATE(sys%gr)
    SAFE_ALLOCATE(sys%st)

#ifdef HAVE_OPENCL
    call opencl_init()
#endif

    call messages_obsolete_variable('SystemName')

    call space_init(sys%space)

    call geometry_init(sys%geo, sys%space)
    call grid_init_stage_0(sys%gr, sys%geo, sys%space)
    call states_init(sys%st, sys%gr, sys%geo)
    call grid_init_stage_1(sys%gr, sys%geo)
    ! if independent particles in N dimensions are being used, need to initialize them
    !  after masses are set to 1 in grid_init_stage_1 -> derivatives_init
    call modelmb_copy_masses (sys%st%modelmbparticles, sys%gr%der%masses)

    call parallel_init()

    call geometry_partition(sys%geo, sys%mc)
    call kpoints_distribute(sys%st%d, sys%mc)
    call states_distribute_nodes(sys%st, sys%mc)
    call grid_init_stage_2(sys%gr, sys%mc, sys%geo)
    call output_init(sys%gr%sb, sys%st%nst, sys%outp)
    call states_densities_init(sys%st, sys%gr, sys%geo, sys%mc)
    call states_lead_densities_init(sys%st, sys%gr)
    call elf_init()

    call poisson_init(psolver, sys%gr%der, sys%geo, sys%mc%master_comm)
    if(poisson_is_multigrid(psolver)) call grid_create_multigrid(sys%gr, sys%geo)

    call v_ks_init(sys%ks, sys%gr, sys%st%d, sys%geo, sys%mc, sys%st%qtot)

    !print the mesh information if it is required
    call print_r()

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
      call multicomm_init(sys%mc, mpi_world, calc_mode_parallel_mask(), calc_mode_default_parallel_mask(), &
        mpi_world%size, index_dim, index_range, (/ 5000, 1, 1, 1 /))

      POP_SUB(system_init.parallel_init)
    end subroutine parallel_init

    subroutine print_r()
      integer :: i, ierr
      character(len=80) :: fname

      PUSH_SUB(system_init.print_r)

      if(iand(sys%outp%what, C_OUTPUT_R).ne.0) then

        do i=1, sys%gr%mesh%sb%dim
          write(fname, '(a,i1)') 'r-', i
          call dio_function_output(sys%outp%how, 'exec/', fname, sys%gr%mesh, sys%gr%mesh%x(:,i), &
            units_out%length, ierr, geo = sys%geo)
        end do
      end if

      POP_SUB(system_init.print_r)
    end subroutine print_r

  end subroutine system_init


  !----------------------------------------------------------
  subroutine system_end(sys)
    type(system_t), intent(inout) :: sys

    PUSH_SUB(system_end)

    call multicomm_end(sys%mc)

    call poisson_end(psolver)
    call v_ks_end(sys%ks, sys%gr, sys%geo)

    if(associated(sys%st)) then
      call states_lead_densities_end(sys%st, sys%gr)
      call states_end(sys%st)
      SAFE_DEALLOCATE_P(sys%st)
    end if

    call geometry_end(sys%geo)
    call simul_box_end(sys%gr%sb)
    call grid_end(sys%gr)

    call space_end(sys%space)

#ifdef HAVE_OPENCL
    call opencl_end()
#endif

    SAFE_DEALLOCATE_P(sys%gr)

    POP_SUB(system_end)
  end subroutine system_end


  !----------------------------------------------------------
  subroutine system_h_setup(sys, hm)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm

    PUSH_SUB(system_h_setup)

    call states_fermi(sys%st, sys%gr%mesh)
    call density_calc(sys%st, sys%gr, sys%st%rho)

    call v_ks_calc(sys%ks, hm, sys%st, sys%geo, calc_eigenval = .true.) ! get potentials
    call states_fermi(sys%st, sys%gr%mesh)                             ! occupations
    call total_energy(hm, sys%gr, sys%st, -1)

    POP_SUB(system_h_setup)
  end subroutine system_h_setup

end module system_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

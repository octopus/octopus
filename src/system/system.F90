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
  use energy_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use mesh_m
  use messages_m
  use modelmb_particles_m
  use mpi_m
  use multicomm_m
  use h_sys_output_m
  use io_function_m
  use profiling_m
  use poisson_m
  use simul_box_m
  use states_m
  use states_dim_m
  use units_m
  use v_ks_m
  use elf_m

  implicit none

  private
  public ::               &
    system_t,             &
    system_init,          &
    system_end,           &
    system_h_setup

  type system_t
    type(geometry_t)          :: geo
    type(grid_t),     pointer :: gr    ! the mesh
    type(states_t),   pointer :: st    ! the states
    type(v_ks_t)              :: ks    ! the Kohn-Sham potentials
    type(h_sys_output_t)      :: outp  ! the output
    type(multicomm_t)         :: mc    ! index and domain communicators
  end type system_t

contains

  !----------------------------------------------------------
  subroutine system_init(sys)
    type(system_t), intent(out) :: sys
    integer :: ii

    call push_sub('system.system_init')

    SAFE_ALLOCATE(sys%gr)
    SAFE_ALLOCATE(sys%st)
    
    call obsolete_variable('SystemName')

    call geometry_init(sys%geo)
    call simul_box_init(sys%gr%sb, sys%geo)
    call states_init(sys%st, sys%gr, sys%geo)
    call grid_init_stage_1(sys%gr, sys%geo)
! if independent particles in N dimensions are being used, need to initialize them
!  after masses are set to 1 in grid_init_stage_1 -> derivatives_init 
    call modelmb_copy_masses (sys%st%modelmbparticles,sys%gr%der%masses)

    call parallel_init()

    call geometry_partition(sys%geo, sys%mc)
    call kpoints_distribute(sys%st%d, sys%mc)
    call grid_init_stage_2(sys%gr, sys%mc, sys%geo)
    call states_densities_init(sys%st, sys%gr, sys%geo, sys%mc)
    call h_sys_output_init(sys%gr%sb, sys%outp)
    call elf_init()
    call poisson_init(sys%gr, sys%geo)
    call v_ks_init(sys%gr, sys%ks, sys%st%d, sys%st%qtot)


    !print the mesh information if it is required
    call print_r()

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine parallel_init()
      integer :: index_dim, index_range(4)

      call push_sub('system.system_init.parallel_init')

      ! for the moment we need communicators for domains plus three
      ! index-dimensions: states, kpoints, and others
      index_dim = 4

      ! store the ranges for these two indices (serves as initial guess
      ! for parallelization strategy)
      index_range(1) = sys%gr%mesh%np_global  ! Number of points in mesh
      index_range(2) = sys%st%nst             ! Number of states
      index_range(3) = sys%st%d%nik           ! Number of kpoints
      index_range(4) = 100000                 ! Some large number

      ! create index and domain communicators
      call multicomm_init(sys%mc, calc_mode_parallel_mask(), calc_mode_default_parallel_mask(), mpi_world%size, index_dim, &
         index_range, (/ 5000, 1, 1, 1 /))

      call pop_sub()
    end subroutine parallel_init

    subroutine print_r()
      FLOAT :: u
      integer :: i, ierr
      character(len=80) :: fname

      call push_sub('system.system_init.print_r')      
      
      if(iand(sys%outp%what, output_r).ne.0) then
        u = M_ONE/units_out%length%factor
        
        do i=1, sys%gr%mesh%sb%dim
          write(fname, '(a,i1)') 'r-', i
          call doutput_function(sys%outp%how, 'exec/', fname, sys%gr%mesh, sys%gr%sb, sys%gr%mesh%x(:,i), u, ierr, geo = sys%geo)
        end do
      end if

      call pop_sub()
    end subroutine print_r
    
  end subroutine system_init


  !----------------------------------------------------------
  subroutine system_end(s)
    type(system_t), intent(inout) :: s

    call push_sub('system.system_end')

#if defined(HAVE_MPI)
    call multicomm_end(s%mc)
#endif

    call poisson_end()
    call v_ks_end(s%ks)

    if(associated(s%st)) then
      call states_end(s%st)
      SAFE_DEALLOCATE_P(s%st); nullify(s%st)
    end if

    call geometry_end(s%geo)
    call simul_box_end(s%gr%sb)
    call grid_end(s%gr)
    SAFE_DEALLOCATE_P(s%gr);  nullify(s%gr)

    call pop_sub()
  end subroutine system_end


  !----------------------------------------------------------
  subroutine system_h_setup(sys, hm)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm

    call push_sub('system.system_h_setup')

    call states_fermi(sys%st, sys%gr%mesh)
    call states_calc_dens(sys%st, sys%gr)

    call v_ks_calc(sys%gr, sys%ks, hm, sys%st, calc_eigenval=.true.) ! get potentials
    call states_fermi(sys%st, sys%gr%mesh)                           ! occupations
    call total_energy(hm, sys%gr, sys%st, -1)

    call pop_sub()
  end subroutine system_h_setup

end module system_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

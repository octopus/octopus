!! Copyright (C) 2002-2006 F. Bonafe, H. Appel, R. Jestaedt
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

module system_mxll_oct_m
  use calc_mode_par_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_mxll_oct_m
  use mesh_oct_m
  use messages_oct_m
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
  use states_mxll_oct_m
  use system_oct_m
  use distributed_oct_m
  
  implicit none

  private
  public ::               &
    system_mxll_t,        &
    system_mxll_init,     &
    system_mxll_end

  integer, parameter, public ::           &
    MULTIGRID_MX_TO_MA_EQUAL   = 1,       &
    MULTIGRID_MX_TO_MA_LARGE   = 2

  type system_mxll_t
    ! Components are public by default
    type(space_t)                :: space
    type(geometry_t)             :: geo
    type(grid_t),        pointer :: gr    !< the mesh
    type(states_mxll_t), pointer :: st    !< the states
    type(output_t)               :: outp  !< the output
    type(multicomm_t)            :: mc    !< index and domain communicators
    type(namespace_t)            :: namespace
    type(hamiltonian_mxll_t)     :: hm
  end type system_mxll_t

contains
  
  !----------------------------------------------------------
  subroutine system_mxll_init(sys, namespace)
    type(system_mxll_t), intent(inout) :: sys
    type(namespace_t), intent(in)  :: namespace

!    type(base_states_t), pointer :: subsys_states

    type(profile_t), save :: prof

    PUSH_SUB(system_mxll_init)
    call profiling_in(prof,"SYSTEM_INIT")

    SAFE_ALLOCATE(sys%gr)
    SAFE_ALLOCATE(sys%st)
    sys%namespace = namespace

!    sys%gr%maxwell_grid           = .true.
!    sys%gr%sb%maxwell_simul_box   = .true.
!    sys%gr%mesh%maxwell_mesh      = .true.
!    sys%st%maxwell_states         = .true.
!    sys%space%maxwell_space       = .true. 
    
    call messages_obsolete_variable(sys%namespace, 'SystemName')

    call space_init(sys%space, sys%namespace)

    ! The geometry needs to be nullified in order to be able to call grid_init_stage_*
    !call geometry_nullify(sys%geo)

    nullify(sys%geo%space, sys%geo%atom, sys%geo%catom, sys%geo%species)
    sys%geo%natoms=0
    sys%geo%ncatoms=0
    sys%geo%nspecies=0
    sys%geo%only_user_def=.false.
    sys%geo%kinetic_energy=M_ZERO
    sys%geo%nlpp=.false.
    sys%geo%nlcc=.false.
    call distributed_nullify(sys%geo%atoms_dist, 0)
    sys%geo%reduced_coordinates=.false.
    sys%geo%periodic_dim=0
    sys%geo%lsize=M_ZERO

    
    call grid_init_stage_0(sys%gr, sys%namespace, sys%geo, sys%space)
    call states_mxll_init(sys%st, sys%namespace, sys%gr, sys%geo)

    call grid_init_stage_1(sys%gr, sys%namespace, sys%geo)
    
    call parallel_mxll_init(sys)

    call grid_init_stage_2(sys%gr, sys%namespace, sys%mc, sys%geo)

    call output_mxll_init(sys%outp, sys%namespace, sys%gr%sb)

    call hamiltonian_mxll_init(sys%hm, sys%namespace, sys%gr, sys%st)
    
    call profiling_out(prof)
    POP_SUB(system_mxll_init)    

  contains

    ! ---------------------------------------------------------
    subroutine parallel_mxll_init(sys)
      type(system_mxll_t), intent(inout) :: sys      

      integer :: index_range(4)

      PUSH_SUB(system_mxll_init.parallel_init)

      ! store the ranges for these two indices (serves as initial guess
      ! for parallelization strategy)
      index_range(1) = sys%gr%mesh%np_global  ! Number of points in mesh
      index_range(2) = sys%st%nst             ! Number of states
      index_range(3) = sys%st%d%nik           ! Number of k-points
      index_range(4) = 100000                 ! Some large number

      !sys%mc%maxwell_mc = .true.
      
      ! create index and domain communicators
      call multicomm_init(sys%mc, sys%namespace, mpi_world, calc_mode_par_parallel_mask(), &
           &calc_mode_par_default_parallel_mask(),mpi_world%size, index_range, (/ 5000, 1, 1, 1 /))

      POP_SUB(system_mxll_init.parallel_init)
    end subroutine parallel_mxll_init

  end subroutine system_mxll_init

  !----------------------------------------------------------
  subroutine system_mxll_end(sys)
    type(system_mxll_t), intent(inout) :: sys

    PUSH_SUB(system_mxll_end)

    call hamiltonian_mxll_end(sys%hm)

    call multicomm_end(sys%mc)

    !call output_mxll_end(sys%outp)
    
    if(associated(sys%st)) then
      call states_mxll_end(sys%st)
      SAFE_DEALLOCATE_P(sys%st)
    end if

    call simul_box_end(sys%gr%sb)
    call grid_end(sys%gr)

    call space_end(sys%space)

    SAFE_DEALLOCATE_P(sys%gr)

    POP_SUB(system_mxll_end)
  end subroutine system_mxll_end


    !----------------------------------------------------------
  subroutine system_check_match_maxwell_and_matter(sys_elec, sys_mxll, multigrid_mode)
    type(system_t),           intent(in)  :: sys_elec
    type(system_mxll_t),      intent(in)  :: sys_mxll
    integer,        optional, intent(out) :: multigrid_mode

    integer :: err

    PUSH_SUB(system_check_match_maxwell_and_matter)

    err = 0
    if (sys_mxll%gr%sb%box_shape /= PARALLELEPIPED)  then
      write(message(1), '(a)') 'The MaxwellBoxShape has to be set as:'
      write(message(2), '(a)') 'MaxwellBoxShape == PARALLELEPIPED'
      call messages_fatal(2)
    end if
    if (sys_mxll%gr%mesh%spacing(1) /= sys_elec%gr%mesh%spacing(1)) then
      if (sys_mxll%gr%mesh%spacing(1) < M_FOUR * sys_elec%gr%sb%lsize(1)) err = err + 1
    end if
    if (sys_mxll%gr%mesh%spacing(2) /= sys_elec%gr%mesh%spacing(2)) then
      if (sys_mxll%gr%mesh%spacing(2) < M_FOUR * sys_elec%gr%sb%lsize(2)) err = err + 1
    end if
    if (sys_mxll%gr%mesh%spacing(3) /= sys_elec%gr%mesh%spacing(3)) then
      if (sys_mxll%gr%mesh%spacing(3) < M_FOUR * sys_elec%gr%sb%lsize(3)) err = err + 1
    end if
    if (err /= 0) then
      write(message(1), '(a)') 'The "MaxwellSpacing" for the Maxwell grid has to be equal to'
      write(message(2), '(a)') 'the "Spacing" of the matter grid for all 3 dimensions.'
      write(message(3), '(a)') 'or'
      write(message(4), '(a)') 'The "MaxwellSpacing" for the Maxwell grid hast to be at least'
      write(message(5), '(a)') 'four times larger than the matter box size "Lsize" for all'
      write(message(6), '(a)') '3 dimensions.'
      call messages_fatal(6)
    end if
    if ( (sys_mxll%gr%mesh%spacing(1) == sys_elec%gr%mesh%spacing(1)) .and. &
         (sys_mxll%gr%mesh%spacing(2) == sys_elec%gr%mesh%spacing(2)) .and. &
         (sys_mxll%gr%mesh%spacing(3) == sys_elec%gr%mesh%spacing(3)) ) then
      multigrid_mode = MULTIGRID_MX_TO_MA_EQUAL 
    else if ( (sys_mxll%gr%mesh%spacing(1) >= M_FOUR * sys_elec%gr%sb%lsize(1)) .and. &
              (sys_mxll%gr%mesh%spacing(2) >= M_FOUR * sys_elec%gr%sb%lsize(2)) .and. &
              (sys_mxll%gr%mesh%spacing(3) >= M_FOUR * sys_elec%gr%sb%lsize(3)) ) then
      multigrid_mode = MULTIGRID_MX_TO_MA_LARGE
    else
      write(message(1), '(a)') 'There is no valid multigird option for current Maxwell and Matter grids.'
      call messages_fatal(1)
    end if

    POP_SUB(system_check_match_maxwell_and_matter)
  end subroutine system_check_match_maxwell_and_matter


end module system_mxll_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

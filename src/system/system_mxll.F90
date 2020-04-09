!! Copyright (C) 2019-2020 F. Bonafe, H. Appel, R. Jestaedt
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
  use clock_oct_m
  use distributed_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_mxll_oct_m
  use interaction_abst_oct_m
  use iso_c_binding
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m
  use quantity_oct_m
  use simul_box_oct_m
  use sort_oct_m
  use space_oct_m
  use system_abst_oct_m
  use states_mxll_oct_m
  use system_oct_m
  
  implicit none

  private
  public ::               &
    system_mxll_t

  integer, parameter, public ::           &
    MULTIGRID_MX_TO_MA_EQUAL   = 1,       &
    MULTIGRID_MX_TO_MA_LARGE   = 2

  type, extends(system_abst_t) :: system_mxll_t
    type(states_mxll_t), pointer :: st    !< the states
    type(hamiltonian_mxll_t)     :: hm
    type(geometry_t)             :: geo
    type(grid_t),        pointer :: gr    !< the mesh
    type(output_t)               :: outp  !< the output
    type(multicomm_t)            :: mc    !< index and domain communicators

    type(c_ptr) :: output_handle
  contains
    procedure :: add_interaction_partner => system_mxll_add_interaction_partner
    procedure :: has_interaction => system_mxll_has_interaction
    procedure :: do_td_operation => system_mxll_do_td
    procedure :: write_td_info => system_mxll_write_td_info
    procedure :: td_write_init => system_mxll_td_write_init
    procedure :: td_write_iter => system_mxll_td_write_iter
    procedure :: td_write_end => system_mxll_td_write_end
    procedure :: is_tolerance_reached => system_mxll_is_tolerance_reached
    procedure :: store_current_status => system_mxll_store_current_status
    procedure :: update_quantity => system_mxll_update_quantity
    procedure :: update_exposed_quantity => system_mxll_update_exposed_quantity
    procedure :: set_pointers_to_interaction => system_mxll_set_pointers_to_interaction
    final :: system_mxll_finalize
  end type system_mxll_t

  interface system_mxll_t
    procedure system_mxll_init
  end interface system_mxll_t

contains

  ! ---------------------------------------------------------
  function system_mxll_init(namespace) result(sys)
    class(system_mxll_t), pointer    :: sys
    type(namespace_t),    intent(in) :: namespace

    type(profile_t), save :: prof

    PUSH_SUB(system_mxll_init)
    call profiling_in(prof,"SYSTEM_INIT")

    SAFE_ALLOCATE(sys)

    sys%namespace = namespace

    SAFE_ALLOCATE(sys%gr)
    SAFE_ALLOCATE(sys%st) ! memleak

    call messages_obsolete_variable(sys%namespace, 'SystemName')

    call space_init(sys%space, sys%namespace)

    ! The geometry needs to be nullified in order to be able to call grid_init_stage_*

    nullify(sys%geo%space, sys%geo%atom, sys%geo%catom, sys%geo%species)
    sys%geo%natoms = 0
    sys%geo%ncatoms = 0
    sys%geo%nspecies = 0
    sys%geo%only_user_def = .false.
    sys%geo%kinetic_energy = M_ZERO
    sys%geo%nlpp = .false.
    sys%geo%nlcc = .false.
    call distributed_nullify(sys%geo%atoms_dist, 0)
    sys%geo%reduced_coordinates = .false.
    sys%geo%periodic_dim = 0
    sys%geo%lsize = M_ZERO
    
    call grid_init_stage_0(sys%gr, sys%namespace, sys%geo, sys%space)  ! memleak in kpoints_init, symmetries_init
    call states_mxll_init(sys%st, sys%namespace, sys%gr, sys%geo)

    call grid_init_stage_1(sys%gr, sys%namespace, sys%geo) ! memleak in kpoints_init, double_grid_init
                                             ! sencil_allocate, derivatives_init, mesh_init_stage_2 
    
    call parallel_mxll_init(sys)

    call grid_init_stage_2(sys%gr, sys%namespace, sys%mc, sys%geo) ! memleak in mem_init_stage_3,
                                                                   ! nl_operator_build, mesh_cube_map

    call output_mxll_init(sys%outp, sys%namespace, sys%gr%sb)

    call hamiltonian_mxll_init(sys%hm, sys%namespace, sys%gr, sys%st)
    
    call profiling_out(prof)

    ! Initialize the propagator
    !call sys%init_propagator()

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

      ! create index and domain communicators
      call multicomm_init(sys%mc, sys%namespace, mpi_world, calc_mode_par_parallel_mask(), &
           &calc_mode_par_default_parallel_mask(),mpi_world%size, index_range, (/ 5000, 1, 1, 1 /))

      POP_SUB(system_mxll_init.parallel_init)
    end subroutine parallel_mxll_init

  end function system_mxll_init

  !----------------------------------------------------------
  subroutine system_mxll_end(sys)
    type(system_mxll_t), intent(inout) :: sys

    PUSH_SUB(system_mxll_end)

    call hamiltonian_mxll_end(sys%hm)

    call multicomm_end(sys%mc)

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


  ! ---------------------------------------------------------
  subroutine system_mxll_add_interaction_partner(this, partner)
    class(system_mxll_t), target, intent(inout) :: this
    class(system_abst_t),         intent(inout) :: partner

    PUSH_SUB(system_mxll_add_interaction_partner)

    POP_SUB(system_mxll_add_interaction_partner)
  end subroutine system_mxll_add_interaction_partner

  ! ---------------------------------------------------------
  logical function system_mxll_has_interaction(this, interaction)
    class(system_mxll_t),      intent(in) :: this
    class(interaction_abst_t), intent(in) :: interaction

    PUSH_SUB(system_mxll_has_interaction)

    POP_SUB(system_mxll_has_interaction)
  end function system_mxll_has_interaction

  ! ---------------------------------------------------------
  subroutine system_mxll_do_td(this, operation)
    class(system_mxll_t), intent(inout) :: this
    integer,              intent(in)    :: operation

    PUSH_SUB(system_mxll_do_td)

    select case(operation)
    case (VERLET_UPDATE_POS)

    case (VERLET_COMPUTE_ACC)

    case (VERLET_COMPUTE_VEL)

    case (BEEMAN_PREDICT_POS)

    case (BEEMAN_PREDICT_VEL)

    case( BEEMAN_CORRECT_POS)

    case (BEEMAN_CORRECT_VEL)

    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

   POP_SUB(system_mxll_do_td)
  end subroutine system_mxll_do_td

  ! ---------------------------------------------------------
  logical function system_mxll_is_tolerance_reached(this, tol) result(converged)
    class(system_mxll_t),   intent(in)    :: this
    FLOAT,                  intent(in)    :: tol

    PUSH_SUB(system_mxll_is_tolerance_reached)

    converged = .false.

    POP_SUB(system_mxll_is_tolerance_reached)
   end function system_mxll_is_tolerance_reached

   ! ---------------------------------------------------------
   subroutine system_mxll_store_current_status(this)
     class(system_mxll_t),   intent(inout)    :: this

     PUSH_SUB(system_mxll_store_current_status)

     POP_SUB(system_mxll_store_current_status)
   end subroutine system_mxll_store_current_status

  ! ---------------------------------------------------------
  subroutine system_mxll_write_td_info(this)
    class(system_mxll_t), intent(in) :: this

    PUSH_SUB(system_mxll_write_td_info)

    POP_SUB(system_mxll_write_td_info)
  end subroutine system_mxll_write_td_info

  ! ---------------------------------------------------------
  subroutine system_mxll_td_write_init(this, dt)
    class(system_mxll_t), intent(inout) :: this
    FLOAT,                intent(in)    :: dt

    PUSH_SUB(system_mxll_td_write_init)


    POP_SUB(system_mxll_td_write_init)
  end subroutine system_mxll_td_write_init

  ! ---------------------------------------------------------
  subroutine system_mxll_td_write_iter(this, iter)
    class(system_mxll_t), intent(inout) :: this
    integer,              intent(in)    :: iter

    integer :: idir
    character(len=50) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(system_mxll_td_write_iter)

    if(iter == 0) then
      ! header
      call write_iter_clear(this%output_handle)
      call write_iter_string(this%output_handle,'################################################################################')
      call write_iter_nl(this%output_handle)
      call write_iter_string(this%output_handle,'# HEADER')
      call write_iter_nl(this%output_handle)

      ! first line: column names
      call write_iter_header_start(this%output_handle)

      do idir = 1, this%space%dim
        write(aux, '(a2,i3,a1)') 'x(', idir, ')'
        call write_iter_header(this%output_handle, aux)
      end do
      do idir = 1, this%space%dim
        write(aux, '(a2,i3,a1)') 'v(', idir, ')'
        call write_iter_header(this%output_handle, aux)
      end do
      call write_iter_nl(this%output_handle)


      call write_iter_string(this%output_handle,'################################################################################')
      call write_iter_nl(this%output_handle)
    end if

    call write_iter_start(this%output_handle)
    call write_iter_nl(this%output_handle)

    POP_SUB(system_mxll_td_write_iter)
  end subroutine system_mxll_td_write_iter

  ! ---------------------------------------------------------
  subroutine system_mxll_td_write_end(this)
    class(system_mxll_t), intent(inout) :: this

    PUSH_SUB(system_mxll_td_write_end)

    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_end(this%output_handle)
    end if

    POP_SUB(system_mxll_td_write_end)
  end subroutine system_mxll_td_write_end

  ! ---------------------------------------------------------
  subroutine system_mxll_update_quantity(this, iq, clock)
    class(system_mxll_t),   intent(inout) :: this
    integer,                   intent(in)    :: iq
    class(clock_t),            intent(in)    :: clock

    PUSH_SUB(system_mxll_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case (MASS)
      !The celestial body has a mass, but it is not necessary to update it, as it does not change with time.
      call this%quantities(iq)%clock%set_time(clock)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(system_mxll_update_quantity)
  end subroutine system_mxll_update_quantity

 ! ---------------------------------------------------------
 logical function system_mxll_update_exposed_quantity(this, iq, clock) result(updated)
    class(system_mxll_t),   intent(inout) :: this
    integer,                   intent(in)    :: iq
    class(clock_t),            intent(in)    :: clock

    PUSH_SUB(system_mxll_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case (MASS)
      !The celestial body has a mass, but it does not require any update, as it does not change with time.
      updated = .true.
      call this%quantities(iq)%clock%set_time(this%clock)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(system_mxll_update_exposed_quantity)
  end function system_mxll_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine system_mxll_set_pointers_to_interaction(this, inter)
    class(system_mxll_t), target,  intent(inout) :: this
    class(interaction_abst_t),        intent(inout) :: inter

    PUSH_SUB(system_mxll_set_pointers_to_interaction)

!    select type(inter)
!    type is(interaction_gravity_t)
!    class default
!      message(1) = "Unsupported interaction."
!      call messages_fatal(1)
!    end select

    POP_SUB(system_mxll_set_pointers_to_interaction)
  end subroutine system_mxll_set_pointers_to_interaction

  ! ---------------------------------------------------------
  subroutine system_mxll_finalize(this)
    type(system_mxll_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter
    class(interaction_abst_t), pointer :: interaction

    PUSH_SUB(system_mxll_finalize)

    !deallocate(this%prop)

    call iter%start(this%interactions)
    do while (iter%has_next())
      interaction => iter%get_next_interaction()
      SAFE_DEALLOCATE_P(interaction)
    end do

    POP_SUB(system_mxll_finalize)
  end subroutine system_mxll_finalize


end module system_mxll_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

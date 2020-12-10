!! Copyright (C) 2019-2020 Franco Bonafe, Heiko Appel, Rene Jestaedt
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
  use algorithm_oct_m
  use calc_mode_par_oct_m
  use clock_oct_m
  use distributed_oct_m
  use external_densities_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_mxll_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use iso_c_binding
  use loct_oct_m
  use lorentz_force_oct_m
  use maxwell_boundary_op_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mesh_interpolation_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use propagator_exp_mid_oct_m
  use propagator_oct_m
  use propagator_mxll_oct_m
  use quantity_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use sort_oct_m
  use space_oct_m
  use system_oct_m
  use states_mxll_oct_m
  use states_mxll_restart_oct_m
  use electrons_oct_m
  use td_write_oct_m
  use unit_oct_m
  use unit_system_oct_m


  implicit none

  private
  public ::               &
    system_mxll_t,        &
    system_mxll_init

  integer, parameter, public ::           &
    MULTIGRID_MX_TO_MA_EQUAL   = 1,       &
    MULTIGRID_MX_TO_MA_LARGE   = 2

  type, extends(system_t) :: system_mxll_t
    type(states_mxll_t), pointer :: st    !< the states
    type(hamiltonian_mxll_t)     :: hm
    type(geometry_t)             :: geo
    type(grid_t),        pointer :: gr    !< the mesh
    type(output_t)               :: outp  !< the output
    type(multicomm_t)            :: mc    !< index and domain communicators

    type(mesh_interpolation_t)   :: mesh_interpolate

    type(propagator_mxll_t)      :: tr_mxll   !< contains the details of the Maxwell time-evolution
    type(td_write_t)             :: write_handler
    type(c_ptr)                  :: output_handle

    CMPLX, allocatable           :: ff_rs_inhom_t1(:,:), ff_rs_inhom_t2(:,:)
    CMPLX, allocatable           :: rs_state_init(:,:)
    FLOAT                        :: bc_bounds(2,MAX_DIM), dt_bounds(2,MAX_DIM)
    FLOAT                        :: etime
    integer                      :: energy_update_iter
    integer                      :: mxll_td_relax_iter
    integer                      :: mxll_ks_relax_iter

  contains
    procedure :: init_interaction => system_mxll_init_interaction
    procedure :: init_parallelization => system_mxll_init_parallelization
    procedure :: initial_conditions => system_mxll_initial_conditions
    procedure :: do_td_operation => system_mxll_do_td
    procedure :: iteration_info => system_mxll_iteration_info
    procedure :: is_tolerance_reached => system_mxll_is_tolerance_reached
    procedure :: update_quantity => system_mxll_update_quantity
    procedure :: update_exposed_quantity => system_mxll_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => system_mxll_copy_quantities_to_interaction
    procedure :: output_start => system_mxll_output_start
    procedure :: output_write => system_mxll_output_write
    procedure :: output_finish => system_mxll_output_finish
    final :: system_mxll_finalize
  end type system_mxll_t

  interface system_mxll_t
    procedure system_mxll_constructor
  end interface system_mxll_t

contains

  ! ---------------------------------------------------------
  function system_mxll_constructor(namespace) result(sys)
    class(system_mxll_t), pointer  :: sys
    type(namespace_t),  intent(in) :: namespace

    PUSH_SUB(system_mxll_constructor)

    SAFE_ALLOCATE(sys)

    call system_mxll_init(sys, namespace)

    POP_SUB(system_mxll_constructor)
  end function system_mxll_constructor


  ! ---------------------------------------------------------
  subroutine system_mxll_init(this, namespace)
    class(system_mxll_t), intent(inout) :: this
    type(namespace_t),    intent(in)    :: namespace

    type(profile_t), save :: prof

    PUSH_SUB(system_mxll_init)

    call profiling_in(prof,"SYSTEM_MXLL_INIT")

    this%namespace = namespace

    SAFE_ALLOCATE(this%gr)
    SAFE_ALLOCATE(this%st)

    call messages_obsolete_variable(this%namespace, 'SystemName')
    call space_init(this%space, this%namespace)

    ! The geometry needs to be nullified in order to be able to call grid_init_stage_*

    nullify(this%geo%space, this%geo%atom, this%geo%catom, this%geo%species)
    this%geo%natoms = 0
    this%geo%ncatoms = 0
    this%geo%nspecies = 0
    this%geo%only_user_def = .false.
    this%geo%kinetic_energy = M_ZERO
    this%geo%nlpp = .false.
    this%geo%nlcc = .false.
    call distributed_nullify(this%geo%atoms_dist, 0)
    this%geo%reduced_coordinates = .false.
    this%geo%periodic_dim = 0
    this%geo%lsize = M_ZERO

    call simul_box_init(this%gr%sb, this%namespace, this%geo, this%space)
    call states_mxll_init(this%st, this%namespace, this%gr, this%geo)
    call grid_init_stage_1(this%gr, this%namespace, this%geo)

    this%quantities(E_FIELD)%required = .true.
    this%quantities(E_FIELD)%protected = .true.
    this%quantities(B_FIELD)%required = .true.
    this%quantities(B_FIELD)%protected = .true.

    call mesh_interpolation_init(this%mesh_interpolate, this%gr%mesh)

    call this%supported_interactions_as_partner%add(LORENTZ_FORCE)

    call profiling_out(prof)

    POP_SUB(system_mxll_init)
  end subroutine system_mxll_init

  ! ---------------------------------------------------------
  subroutine system_mxll_init_interaction(this, interaction)
    class(system_mxll_t), target, intent(inout) :: this
    class(interaction_t),         intent(inout) :: interaction

    PUSH_SUB(system_mxll_init_interaction)

    select type (interaction)
    class default
      ! Currently Maxwell system does not know any type of interaction
      message(1) = "Trying to initialize an unsupported interaction by Maxwell."
      call messages_fatal(1)
    end select

    POP_SUB(system_mxll_init_interaction)
  end subroutine system_mxll_init_interaction

  ! ---------------------------------------------------------
  subroutine system_mxll_init_parallelization(this, grp)
    class(system_mxll_t), intent(inout) :: this
    type(mpi_grp_t),      intent(in)    :: grp

    integer :: index_range(4)

    PUSH_SUB(system_mxll_init_parallelization)

    call mpi_grp_copy(this%grp, grp)

    ! store the ranges for these two indices (serves as initial guess
    ! for parallelization strategy)
    index_range(1) = this%gr%mesh%np_global  ! Number of points in mesh
    index_range(2) = this%st%nst             ! Number of states
    index_range(3) = 1                      ! Number of k-points
    index_range(4) = 100000                 ! Some large number

    ! create index and domain communicators
    call multicomm_init(this%mc, this%namespace, mpi_world, calc_mode_par_parallel_mask(), &
         &calc_mode_par_default_parallel_mask(),mpi_world%size, index_range, (/ 5000, 1, 1, 1 /))

    call grid_init_stage_2(this%gr, this%namespace, this%mc)
    call output_mxll_init(this%outp, this%namespace, this%gr%sb)
    call hamiltonian_mxll_init(this%hm, this%namespace, this%gr, this%st)

    POP_SUB(system_mxll_init_parallelization)
  end subroutine system_mxll_init_parallelization

  ! ---------------------------------------------------------
  subroutine system_mxll_initial_conditions(this, from_scratch)
    class(system_mxll_t), intent(inout) :: this
    logical,              intent(in)    :: from_scratch

    integer :: inter_steps_default
    integer :: relax_iter_default
    FLOAT   :: courant

    type(profile_t), save :: prof

    PUSH_SUB(system_mxll_initial_conditions)

    call profiling_in(prof,"SYSTEM_MXLL_INITIAL_CONDITIONS")

    courant = M_ONE/(P_c * sqrt(M_ONE/this%gr%mesh%spacing(1)**2 + M_ONE/this%gr%mesh%spacing(2)**2 + &
         M_ONE/this%gr%mesh%spacing(3)**2) )

    if(this%prop%dt > M_TWO * courant) then
      write(message(1),'(a,es9.2)') 'Time step seems too large, check this value'
      call messages_warning(1, namespace=this%namespace)
    end if

    !%Variable TDMaxwellTDRelaxationSteps
    !%Type integer
    !%Default 0
    !%Section Time-Dependent::Propagation
    !%Description
    !% Time steps needed to relax the Maxwell states in the presence of a matter system, to avoid
    !% spurious relaxation effects. After these steps, the Maxwell-matter coupling can be switched on.
    !% of the relaxation dynamics.
    !%End
    relax_iter_default = 0
    call parse_variable(this%namespace, 'TDMaxwellTDRelaxationSteps', relax_iter_default, this%mxll_td_relax_iter)

    !%Variable TDMaxwellKSRelaxationSteps
    !%Type integer
    !%Default 0
    !%Section Time-Dependent::Propagation
    !%Description
    !% Time steps in which the coupled Maxwell-matter system  relax the Maxwell states evolves under
    !% free dynamics conditions. After these many steps, the external fields and currents are
    !% switched on. The full requested simulation effectively states after this value.
    !%End
    relax_iter_default = 0
    call parse_variable(this%namespace, 'TDMaxwellKSRelaxationSteps', relax_iter_default, this%mxll_ks_relax_iter)

    ! maxwell delay time
    this%tr_mxll%delay_time = this%mxll_ks_relax_iter * this%prop%dt

    if ( (this%mxll_td_relax_iter /= 0) .and. (this%mxll_ks_relax_iter /= 0) ) then
      if ( .not. ( (this%mxll_td_relax_iter < this%mxll_ks_relax_iter) ) ) then
        call messages_write('TDMaxwellTDRelaxationSteps ')
        call messages_write(' has to be smaller than ')
        call messages_write(' TDMaxwellKSRelaxationSteps. ')
        call messages_write(' TDMaxwellTDRelaxationSteps ')
        call messages_write(' and ')
        call messages_write(' TDMaxwellKSRelaxationSteps ')
        call messages_write(' have to be smaller than TDMaxSteps.')
        call messages_fatal()
      end if
    end if

    !%Variable MaxwellTDIntervalSteps
    !%Type integer
    !%Section Time-Dependent::Propagation
    !%Description
    !% This variable determines how many intervall steps the Maxwell field propagation
    !% does until it reaches the matter time step. In case that MaxwellTDIntervalSteps is
    !% equal to one, the Maxwell time step is equal to the matter one. The default value is 1.
    !%End
    inter_steps_default = 1
    call parse_variable(this%namespace, 'MaxwellTDIntervalSteps', inter_steps_default, this%tr_mxll%inter_steps)

    if (this%tr_mxll%inter_steps < 1) then
      call messages_write('MaxwellTDIntervalSteps hast to be larger than 0 !')
      call messages_fatal()
    end if

    this%st%energy_rate = M_ZERO
    this%st%delta_energy = M_ZERO
    this%st%energy_via_flux_calc = M_ZERO
    this%st%trans_energy_rate = M_ZERO
    this%st%trans_delta_energy = M_ZERO
    this%st%trans_energy_via_flux_calc = M_ZERO
    this%st%plane_waves_energy_rate = M_ZERO
    this%st%plane_waves_delta_energy = M_ZERO
    this%st%plane_waves_energy_via_flux_calc = M_ZERO

    SAFE_ALLOCATE(this%rs_state_init(1:this%gr%mesh%np_part, 1:this%st%dim))
    this%rs_state_init(:,:) = M_z0

    this%energy_update_iter = 1

    call propagator_mxll_init(this%gr, this%namespace, this%st, this%hm, this%tr_mxll)
    call states_mxll_allocate(this%st, this%gr%mesh)
    call external_current_init(this%st, this%namespace, this%gr%mesh)
    this%hm%propagation_apply = .true.

    if (parse_is_defined(this%namespace, 'MaxwellIncidentWaves') .and. (this%tr_mxll%bc_plane_waves)) then
      this%st%rs_state_plane_waves(:,:) = M_z0
    end if

    this%hm%plane_waves_apply = .true.
    this%hm%spatial_constant_apply = .true.
    call bc_mxll_init(this%hm%bc, this%namespace, this%gr, this%st, this%gr%sb, this%geo, this%prop%dt/this%tr_mxll%inter_steps)
    this%bc_bounds(:,1:3) = this%hm%bc%bc_bounds(:,1:3)
    call inner_and_outer_points_mapping(this%gr%mesh, this%st, this%bc_bounds)
    this%dt_bounds(2, 1:3) = this%bc_bounds(1, 1:3)
    this%dt_bounds(1, 1:3) = this%bc_bounds(1, 1:3) - this%gr%der%order * this%gr%mesh%spacing(1:3)
    call surface_grid_points_mapping(this%gr%mesh, this%st, this%dt_bounds)

    if (parse_is_defined(this%namespace, 'UserDefinedInitialMaxwellStates')) then
      call states_mxll_read_user_def(this%gr%mesh, this%st, this%rs_state_init, this%namespace)
      call messages_print_stress(stdout, "Setting initial EM field inside box")
      ! TODO: add consistency check that initial state fulfills Gauss laws
      this%st%rs_state(:,:) = this%st%rs_state + this%rs_state_init
      if (this%tr_mxll%bc_plane_waves) then
        this%st%rs_state_plane_waves(:,:) = this%rs_state_init
      end if
    end if

    ! initialize the spatial constant field according to the conditions set in the
    ! UserDefinedConstantSpatialMaxwellField block
    if (this%tr_mxll%bc_constant) then
      call spatial_constant_calculation(this%tr_mxll%bc_constant, this%st, this%gr, this%hm, M_ZERO, &
           this%prop%dt/this%tr_mxll%inter_steps, this%tr_mxll%delay_time, this%st%rs_state, &
           set_initial_state = .true.)
      this%st%rs_state_const(:) = this%st%rs_state(this%gr%mesh%idx%lxyz_inv(0,0,0),:)
    end if

    if (parse_is_defined(this%namespace, 'UserDefinedInitialMaxwellStates')) then
      SAFE_DEALLOCATE_A(this%rs_state_init)
    end if

    call hamiltonian_mxll_update(this%hm, time = M_ZERO)

    ! calculate Maxwell energy
    call energy_mxll_calc(this%gr, this%st, this%hm, this%hm%energy, this%st%rs_state, &
         this%st%rs_state_plane_waves)

    this%st%rs_state_trans(:,:) = this%st%rs_state

    call get_rs_state_at_point(this%st%selected_points_rs_state(:,:), this%st%rs_state, &
      this%st%selected_points_coordinate(:,:), this%st, this%gr%mesh)

    call profiling_out(prof)

    POP_SUB(system_mxll_initial_conditions)
  end subroutine system_mxll_initial_conditions

  ! ---------------------------------------------------------
  subroutine system_mxll_do_td(this, operation)
    class(system_mxll_t),           intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    CMPLX, allocatable :: current_density_ext(:,:), charge_density_ext(:)
    type(profile_t), save :: prof

    PUSH_SUB(system_mxll_do_td)

    select case (operation%id)
    case (SKIP)
      ! Do nothing

    case (STORE_CURRENT_STATUS)
      ! For the moment we do nothing

    case (EXPMID_START)
      ! This variable is used to compute the elapsed time during the time-step.
      ! This is incorrect when there is more than one system, as the operations for the different systems
      ! are intermingled. Therefore it needs to be changed (maybe have the propagator handle it?)
      this%etime = loct_clock()

      SAFE_ALLOCATE(this%ff_rs_inhom_t1(1:this%gr%mesh%np_part, 1:this%hm%dim))
      SAFE_ALLOCATE(this%ff_rs_inhom_t2(1:this%gr%mesh%np_part, 1:this%hm%dim))

    case (EXPMID_FINISH)

      SAFE_DEALLOCATE_A(this%ff_rs_inhom_t1)
      SAFE_DEALLOCATE_A(this%ff_rs_inhom_t2)

    case (EXPMID_PREDICT_DT_2)  ! predict: psi(t+dt/2) = 0.5*(U_H(dt) psi(t) + psi(t)) or via extrapolation
      ! Empty for the moment

      this%quantities(E_FIELD)%clock = this%quantities(E_FIELD)%clock + CLOCK_TICK
      this%quantities(B_FIELD)%clock = this%quantities(B_FIELD)%clock + CLOCK_TICK
    case (UPDATE_HAMILTONIAN)   ! update: H(t+dt/2) from psi(t+dt/2)
      ! Empty for the moment
    case (EXPMID_PREDICT_DT)    ! predict: psi(t+dt) = U_H(t+dt/2) psi(t)

      call profiling_in(prof, "SYSTEM_MXLL_DO_TD")

      ! Propagation

      !We first compute thre external charge and current densities and we convert them as RS vectors
      SAFE_ALLOCATE(current_density_ext(1:this%gr%mesh%np, 1:this%gr%sb%dim))
      SAFE_ALLOCATE(charge_density_ext(1:this%gr%mesh%np))

      ! calculation of external RS density at time (time-dt)
      if (this%hm%current_density_ext_flag) then
        call get_rs_density_ext(this%st, this%gr%mesh, this%clock%time(), current_density_ext)
      else
        current_density_ext = M_z0
      end if
      !No charge density at the moment
      charge_density_ext = M_z0

      call transform_rs_densities(this%hm, this%gr%mesh, charge_density_ext, &
                     current_density_ext, this%ff_rs_inhom_t1, RS_TRANS_FORWARD)

      ! calculation of external RS density at time (time)
      if (this%hm%current_density_ext_flag) then
        call get_rs_density_ext(this%st, this%gr%mesh, this%clock%time()+this%prop%dt, current_density_ext)
      end if
      !No charge density at the moment

      call transform_rs_densities(this%hm, this%gr%mesh, charge_density_ext, &
                     current_density_ext, this%ff_rs_inhom_t2, RS_TRANS_FORWARD)

      SAFE_DEALLOCATE_A(current_density_ext)
      SAFE_DEALLOCATE_A(charge_density_ext)

      ! Propagation dt with H_maxwell
      call mxll_propagation_step(this%hm, this%namespace, this%gr, this%st, this%tr_mxll,&
          this%st%rs_state, this%ff_rs_inhom_t1, this%ff_rs_inhom_t2, this%clock%time(), this%prop%dt)

      this%st%rs_state_trans(:,:) = this%st%rs_state

      ! calculate Maxwell energy
      call energy_mxll_calc(this%gr, this%st, this%hm, this%hm%energy, this%st%rs_state, &
           this%st%rs_state_plane_waves)

      ! get RS state values for selected points
      call get_rs_state_at_point(this%st%selected_points_rs_state(:,:), this%st%rs_state, this%st%selected_points_coordinate(:,:), &
        this%st, this%gr%mesh)

      this%quantities(E_FIELD)%clock = this%quantities(E_FIELD)%clock + CLOCK_TICK
      this%quantities(B_FIELD)%clock = this%quantities(B_FIELD)%clock + CLOCK_TICK

      this%etime = loct_clock()

      call profiling_out(prof)

    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select

   POP_SUB(system_mxll_do_td)
  end subroutine system_mxll_do_td

  ! ---------------------------------------------------------
  subroutine system_mxll_iteration_info(this)
    class(system_mxll_t), intent(in) :: this

    write(message(1), '(i8,1x,f13.6,2x,f13.6,6x,f13.6)') this%clock%get_tick(), &
      units_from_atomic(units_out%time, this%clock%time()),                     &
      units_from_atomic(units_out%energy, this%hm%energy%energy),               &
      loct_clock() - this%etime
    call messages_info(1)

  end subroutine system_mxll_iteration_info

  ! ---------------------------------------------------------
  logical function system_mxll_is_tolerance_reached(this, tol) result(converged)
    class(system_mxll_t),   intent(in)    :: this
    FLOAT,                  intent(in)    :: tol

    PUSH_SUB(system_mxll_is_tolerance_reached)

    converged = .false.

    POP_SUB(system_mxll_is_tolerance_reached)
   end function system_mxll_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine system_mxll_update_quantity(this, iq)
    class(system_mxll_t),      intent(inout) :: this
    integer,                   intent(in)    :: iq

    PUSH_SUB(system_mxll_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(system_mxll_update_quantity)
  end subroutine system_mxll_update_quantity

 ! ---------------------------------------------------------
 subroutine system_mxll_update_exposed_quantity(partner, iq)
    class(system_mxll_t),      intent(inout) :: partner
    integer,                   intent(in)    :: iq

    PUSH_SUB(system_mxll_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(system_mxll_update_exposed_quantity)
  end subroutine system_mxll_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine system_mxll_copy_quantities_to_interaction(partner, interaction)
    class(system_mxll_t),       intent(inout) :: partner
    class(interaction_t),       intent(inout) :: interaction

    CMPLX :: interpolated_value(3)
    FLOAT :: e_field(3)
    FLOAT :: b_field(3)
    type(profile_t), save :: prof

    PUSH_SUB(system_mxll_copy_quantities_to_interaction)

    call profiling_in(prof, "SYSTEM_MXLL_COPY_QUANTITIES_TO_INTERACTION")

    select type (interaction)
    type is (lorentz_force_t)
      call mesh_interpolation_evaluate(partner%mesh_interpolate, partner%st%rs_state(:,1), &
        interaction%system_pos, interpolated_value(1))
      call mesh_interpolation_evaluate(partner%mesh_interpolate, partner%st%rs_state(:,2), &
        interaction%system_pos, interpolated_value(2))
      call mesh_interpolation_evaluate(partner%mesh_interpolate, partner%st%rs_state(:,3), &
        interaction%system_pos, interpolated_value(3))
      call get_electric_field_vector(interpolated_value, e_field)
      call get_magnetic_field_vector(interpolated_value, 1, b_field)
      interaction%partner_E_field = e_field
      interaction%partner_B_field = b_field
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    call profiling_out(prof)

    POP_SUB(system_mxll_copy_quantities_to_interaction)
  end subroutine system_mxll_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine system_mxll_output_start(this)
    class(system_mxll_t), intent(inout) :: this

    type(profile_t), save :: prof

    PUSH_SUB(system_mxll_output_start)

    call profiling_in(prof, "SYSTEM_MXLL_OUTPUT_START")

    call td_write_mxll_init(this%write_handler, this%namespace, this%gr, this%st, this%hm, 0, this%prop%dt)
    call td_write_mxll_iter(this%write_handler, this%gr, this%st, this%hm, this%prop%dt, 0)
    call td_write_mxll_free_data(this%write_handler, this%namespace, this%gr, this%st, this%hm, this%geo, this%outp, this%clock)

    ! Currently we print this header here, but this needs to be changed.
    write(message(1), '(a10,1x,a10,1x,a20,1x,a18)') 'Iter ', 'Time ',  'Maxwell energy', 'Elapsed Time'
    call messages_info(1)
    call messages_print_stress(stdout)

    call profiling_out(prof)

    POP_SUB(system_mxll_output_start)
  end subroutine system_mxll_output_start

  ! ---------------------------------------------------------
  subroutine system_mxll_output_write(this)
    class(system_mxll_t), intent(inout) :: this

    logical :: stopping
    type(profile_t), save :: prof

    PUSH_SUB(system_mxll_output_write)

    call profiling_in(prof, "SYSTEM_MXLL_OUTPUT_WRITE")

    stopping = clean_stop(this%mc%master_comm)

    call td_write_mxll_iter(this%write_handler, this%gr, this%st, this%hm, this%prop%dt, this%clock%get_tick())

    if ((this%outp%output_interval > 0 .and. mod(this%clock%get_tick(), this%outp%output_interval) == 0) .or. stopping) then
      call td_write_mxll_free_data(this%write_handler, this%namespace, this%gr, this%st, this%hm, this%geo, this%outp, this%clock)
    end if

    call profiling_out(prof)

    POP_SUB(system_mxll_output_write)
  end subroutine system_mxll_output_write

  ! ---------------------------------------------------------
  subroutine system_mxll_output_finish(this)
    class(system_mxll_t), intent(inout) :: this

    type(profile_t), save :: prof

    PUSH_SUB(system_mxll_output_finish)

    call profiling_in(prof, "SYSTEM_MXLL_OUTPUT_FINISH")

    call td_write_mxll_end(this%write_handler)

    call profiling_out(prof)

    POP_SUB(system_mxll_output_finish)
  end subroutine system_mxll_output_finish

  ! ---------------------------------------------------------
  subroutine system_mxll_finalize(this)
    type(system_mxll_t), intent(inout) :: this

    type(profile_t), save :: prof

    PUSH_SUB(system_mxll_finalize)

    call profiling_in(prof, "SYSTEM_MXLL_FINALIZE")

    call system_end(this)

    ! free memory
    SAFE_DEALLOCATE_A(this%rs_state_init)

    call hamiltonian_mxll_end(this%hm)

    call multicomm_end(this%mc)

    if(associated(this%st)) then
      call states_mxll_end(this%st)
      SAFE_DEALLOCATE_P(this%st)
    end if

    call simul_box_end(this%gr%sb)
    call grid_end(this%gr)

    call space_end(this%space)

    SAFE_DEALLOCATE_P(this%gr)

    call profiling_out(prof)

    POP_SUB(system_mxll_finalize)
  end subroutine system_mxll_finalize

end module system_mxll_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

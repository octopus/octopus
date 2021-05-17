!! Copyright (C) 2019 R. Jestaedt, F. Bonafe, H. Appel
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

#include "global.h"

module propagator_mxll_oct_m
  use boundary_op_oct_m
  use batch_ops_oct_m
  use batch_oct_m
  use comm_oct_m
  use cube_function_oct_m
  use cube_oct_m
  use derivatives_oct_m
  use density_oct_m
  use energy_calc_oct_m
  use energy_mxll_oct_m
  use exponential_oct_m
  use external_densities_oct_m
  use fft_oct_m
  use fourier_space_oct_m
  use grid_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_mxll_oct_m
  use index_oct_m
  use io_oct_m
  use io_function_oct_m
  use ions_oct_m
  use lalg_adv_oct_m
  use lalg_basic_oct_m
  use loct_math_oct_m
  use math_oct_m
  use maxwell_boundary_op_oct_m
  use maxwell_function_oct_m
  use mesh_oct_m
  use mesh_cube_parallel_map_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use output_oct_m
  use parser_oct_m
  use par_vec_oct_m
  use profiling_oct_m
  use poisson_oct_m
  use states_elec_oct_m
  use states_mxll_oct_m
  use string_oct_m
  use tdfunction_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                                  &
    propagator_mxll_t,                       &
    propagator_mxll_init,                    &
    mxll_propagation_step,                   &
    transform_rs_densities,                  &
    energy_mxll_calc,                        &
    spatial_constant_calculation,            &
    set_medium_rs_state

  ! The following routines are currently unused, but will be used in the near future.
  ! In order not to generate warnings about them, we declared them as public
  public ::                                  &
    calculate_matter_longitudinal_field,     &
    get_vector_pot_and_transverse_field,     &
    calculate_vector_potential

  type propagator_mxll_t
    integer             :: op_method
    logical             :: bc_add_ab_region  = .false.
    logical             :: bc_zero           = .false.
    logical             :: bc_constant       = .false.
    logical             :: bc_mirror_pec     = .false.
    logical             :: bc_mirror_pmc     = .false.
    logical             :: bc_periodic       = .false.
    logical             :: bc_plane_waves    = .false.
    logical             :: bc_medium         = .false.
    type(exponential_t) :: te
    integer             :: inter_steps
    FLOAT               :: delay_time
    FLOAT               :: scf_threshold
    logical             :: plane_waves_in_box
    integer             :: tr_etrs_approx
  end type propagator_mxll_t

  integer, public, parameter ::   &
     RS_TRANS_FORWARD  = 1,       &
     RS_TRANS_BACKWARD = 2

  integer, parameter ::    & 
    MXWLL_ETRS_FULL  = 0,  &
    MXWLL_ETRS_CONST = 1

contains

  ! ---------------------------------------------------------
  subroutine propagator_mxll_init(gr, namespace, st, hm, tr)
    type(grid_t),                 intent(in)    :: gr
    type(namespace_t),            intent(in)    :: namespace
    type(states_mxll_t),          intent(inout) :: st
    type(hamiltonian_mxll_t),     intent(inout) :: hm
    type(propagator_mxll_t),      intent(inout) :: tr

    integer :: default_propagator, nlines, ncols, icol
    type(block_t) :: blk
    character(len=256) :: string
    logical :: plane_waves_set = .false.
    type(profile_t), save :: prof

    PUSH_SUB(propagator_mxll_init)

    call profiling_in(prof,"PROPAGATOR_MXLL_INIT")

    hm%bc%bc_type(:) = MXLL_BC_ZERO ! default boundary condition is zero

    !%Variable MaxwellBoundaryConditions
    !%Type block
    !%Section Time-Dependent::Propagation
    !%Description
    !% Defines boundary conditions for the electromagnetic field propagation.
    !%
    !% Example:
    !%
    !% <tt>%MaxwellBoundaryConditions
    !% <br>&nbsp;&nbsp;   zero | mirror_pec | consant
    !% <br>%</tt>
    !%
    !%
    !%Option zero 0
    !% Boundaries are set to zero.
    !%Option constant 1
    !% Boundaries are set to a constant.
    !%Option mirror_pec 2
    !% Perfect electric conductor.
    !%Option mirror_pmc 3
    !% Perfect magnetic conductor.
    !%Option plane_waves 4
    !% Boundaries feed in plane waves.
    !%Option periodic 5
    !% Periodic boundary conditions (not yet implemented).
    !%Option medium 6
    !% Boundaries as linear medium (not yet implemented).
    !%End
    if(parse_block(namespace, 'MaxwellBoundaryConditions', blk) == 0) then

      call messages_print_stress(stdout, trim('Maxwell boundary conditions:'))

      ! find out how many lines (i.e. states) the block has
      nlines = parse_block_n(blk)
      if (nlines /= 1) then
        call messages_input_error(namespace, 'MaxwellBoundaryConditions', 'should consist of one line')
      end if
      ncols = parse_block_cols(blk, 0)
      if (ncols /= 3) then
        call messages_input_error(namespace, 'MaxwellBoundaryConditions', 'should consist of three columns')
      end if
      do icol = 1, ncols
        call parse_block_integer(blk, 0, icol-1, hm%bc%bc_type(icol))
      end do
      call parse_block_end(blk)
      call messages_print_stress(stdout, namespace=namespace)
    end if

    do icol = 1, 3
      select case (hm%bc%bc_type(icol))
      case (MXLL_BC_ZERO)
        string = 'Zero'
        hm%bc_zero = .true.
        tr%bc_zero = .true.
      case (MXLL_BC_CONSTANT)
        string = 'Constant'
        tr%bc_constant = .true.
        tr%bc_add_ab_region = .true.
        hm%bc_constant = .true.
        hm%bc_add_ab_region = .true.
      case (MXLL_BC_MIRROR_PEC)
        string = 'PEC Mirror'
        tr%bc_mirror_pec = .true.
        hm%bc_mirror_pec = .true.
      case (MXLL_BC_MIRROR_PMC)
        string = 'PMC Mirror'
        tr%bc_mirror_pmc = .true.
        hm%bc_mirror_pmc = .true.
      case (MXLL_BC_PERIODIC)
        string = 'Periodic'
        tr%bc_periodic = .true.
        hm%bc_periodic = .true.
      case (MXLL_BC_PLANE_WAVES)
        string = 'Plane waves'
        plane_waves_set = .true.
        tr%bc_plane_waves = .true.
        tr%bc_add_ab_region = .true.
        hm%plane_waves = .true.
        hm%bc_plane_waves = .true.
        hm%bc_add_ab_region = .true.
      case (MXLL_BC_MEDIUM)
        string = 'Medium boundary'
      case default
        write(message(1),'(a)') 'Unknown Maxwell boundary condition'
        call messages_fatal(1, namespace=namespace)
      end select
      write(message(1),'(a,I1,a,a)') 'Maxwell boundary condition in direction ', icol, ': ', trim(string)
      call messages_info(1)
      if (plane_waves_set .and. .not. (parse_is_defined(namespace, 'MaxwellIncidentWaves')) ) then
        write(message(1),'(a)') 'Input: Maxwell boundary condition option is set to "plane_waves".'
        write(message(2),'(a)') 'Input: User defined Maxwell plane waves have to be defined!'
        call messages_fatal(2, namespace=namespace)
      end if
    end do

    if (any(hm%bc%bc_type(1:3) == MXLL_BC_CONSTANT)) then
      call td_function_mxll_init(st, namespace, hm)
      SAFE_ALLOCATE(st%rs_state_const(1:st%dim))
      st%rs_state_const = M_z0
    end if

    !%Variable MaxwellTDETRSApprox
    !%Type integer
    !%Default no
    !%Section Time-Dependent::Propagation
    !%Description
    !% Whether to perform  aproximations to the ETRS propagator.
    !%Option no 0
    !% No approximations.
    !%Option const_steps 1
    !% Use constant current density.
    !%End
    call parse_variable(namespace, 'MaxwellTDETRSApprox', MXWLL_ETRS_FULL, tr%tr_etrs_approx)
    call messages_print_var_option(stdout, 'MaxwellTDETRSApprox', tr%tr_etrs_approx)

    !%Variable MaxwellTDOperatorMethod
    !%Type integer
    !%Default op_fd
    !%Section Time-Dependent::Propagation
    !%Description
    !% The Maxwell Operator e.g. the curl operation can be obtained by
    !% two different methods, the finid-difference or the fast fourier
    !% transform.
    !%Option op_fd 1
    !% Maxwell operator calculated by finite differnce method
    !%Option op_fft 2
    !% Maxwell operator calculated by fast fourier transform
    !%End
    default_propagator = OPTION__MAXWELLTDOPERATORMETHOD__OP_FD
    call parse_variable(namespace, 'MaxwellTDOperatorMethod', default_propagator, tr%op_method)
    call messages_print_var_option(stdout, 'MaxwellTDOperatorMethod', tr%op_method)
    hm%op_method = tr%op_method

    !%Variable MaxwellTDSCFThreshold
    !%Type float
    !%Default 1.0e-6
    !%Section Time-Dependent::Propagation
    !%Description
    !% Since the Maxwell-KS propagator is non-linear, each propagation step
    !% should be performed self-consistently.  In practice, for most
    !% purposes this is not necessary, except perhaps in the first
    !% iterations. This variable holds the number of propagation steps
    !% for which the propagation is done self-consistently.
    !%
    !% This variable controls the accuracy threshold for the self consistency.
    !%End
    call parse_variable(namespace, 'MaxwellTDSCFThreshold', CNST(1.0e-6), tr%scf_threshold)

    !%Variable MaxwellPlaneWavesInBox
    !%Type logical
    !%Default no
    !%Section States
    !%Description
    !% Analytic evaluation of the incoming waves inside the box,
    !% not doing any numerical propagation of Maxwells equations.
    !%End
    call parse_variable(namespace, 'MaxwellPlaneWavesInBox', .false., tr%plane_waves_in_box)

    call derivatives_boundary_mask(hm%bc, gr%mesh, hm)

    !tr%te%exp = .true.
    call exponential_init(tr%te, namespace) ! initialize Maxwell propagator

    call profiling_out(prof)

    POP_SUB(propagator_mxll_init)
  end subroutine propagator_mxll_init

  ! ---------------------------------------------------------
  subroutine mxll_propagation_step(hm, namespace, gr, st, tr, rs_state, ff_rs_inhom_t1, ff_rs_inhom_t2, time, dt)
    type(hamiltonian_mxll_t),   intent(inout) :: hm
    type(namespace_t),          intent(in)    :: namespace
    type(grid_t),               intent(inout) :: gr
    type(states_mxll_t),        intent(inout) :: st
    type(propagator_mxll_t),    intent(inout) :: tr
    CMPLX,                      intent(inout) :: rs_state(:,:)
    CMPLX,                      intent(in)    :: ff_rs_inhom_t1(:,:) !> Inhomogeneous term at t
    CMPLX,                      intent(in)    :: ff_rs_inhom_t2(:,:) !> Inhomogeneous term at t+dt
    FLOAT,                      intent(in)    :: time
    FLOAT,                      intent(in)    :: dt

    integer            :: ii, inter_steps, ff_dim, idim, istate
    FLOAT              :: inter_dt, inter_time, delay


    logical            :: pml_check = .false.
    type(profile_t), save :: prof
    type(batch_t) :: ff_rs_stateb, ff_rs_state_pmlb
    type(batch_t) :: ff_rs_inhom_1b, ff_rs_inhom_2b, ff_rs_inhom_meanb

    PUSH_SUB(mxll_propagation_step)

    call profiling_in(prof, 'MXLL_PROPAGATOR_STEP')

    if (hm%ma_mx_coupling_apply) then
      message(1) = "Maxwell-matter coupling not implemented yet"
      call messages_fatal(1, namespace=namespace)
    end if

    if (tr%plane_waves_in_box) then
      call plane_waves_in_box_calculation(hm%bc, time+dt, gr, st, hm, rs_state)
      POP_SUB(mxll_propagation_step)
      return
    end if

    do idim = 1, 3
      if (hm%bc%bc_ab_type(idim) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML) then
        pml_check = .true.
      end if
    end do

    ff_dim = hm%dim

    ! intermediate step variables
    inter_steps   = tr%inter_steps
    inter_dt      = M_ONE / inter_steps * dt

    ! delay time
    delay = tr%delay_time

    call zbatch_init(ff_rs_stateb, 1, 1, hm%dim, gr%mesh%np_part)
    if (st%pack_states) call ff_rs_stateb%do_pack()

    if (pml_check) then
      call ff_rs_stateb%copy_to(ff_rs_state_pmlb)
    end if

    ! first step of Maxwell inhomogeneity propagation with constant current density
    if ((hm%ma_mx_coupling_apply .or. hm%current_density_ext_flag) .and. &
        tr%tr_etrs_approx == MXWLL_ETRS_CONST) then
      call ff_rs_stateb%copy_to(ff_rs_inhom_1b)
      call ff_rs_stateb%copy_to(ff_rs_inhom_2b)
      call ff_rs_stateb%copy_to(ff_rs_inhom_meanb)

      do istate = 1, hm%dim
        call batch_set_state(ff_rs_inhom_meanb, istate, gr%mesh%np, ff_rs_inhom_t1(:, istate))
        call batch_set_state(ff_rs_inhom_2b, istate, gr%mesh%np, ff_rs_inhom_t2(:, istate))
      end do
      call batch_axpy(gr%mesh%np, M_ONE, ff_rs_inhom_2b, ff_rs_inhom_meanb)
      call batch_scal(gr%mesh%np, M_HALF, ff_rs_inhom_meanb)

      ! inhomogeneity propagation
      call ff_rs_inhom_meanb%copy_data_to(gr%mesh%np, ff_rs_inhom_1b)
      call ff_rs_inhom_meanb%copy_data_to(gr%mesh%np, ff_rs_inhom_2b)

      call hamiltonian_mxll_update(hm, time=time)
      hm%cpml_hamiltonian = .false.
      call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, ff_rs_inhom_2b, inter_dt)

      ! add term U(time+dt,time)J(time)
      call batch_axpy(gr%mesh%np, M_ONE, ff_rs_inhom_2b, ff_rs_inhom_1b)
      call ff_rs_inhom_meanb%copy_data_to(gr%mesh%np, ff_rs_inhom_2b)
      call hamiltonian_mxll_update(hm, time=time)
      call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, ff_rs_inhom_2b, inter_dt*M_HALF)
      ! add term U(time+dt/2,time)J(time)
      call batch_axpy(gr%mesh%np, M_ONE, ff_rs_inhom_2b, ff_rs_inhom_1b)
      call ff_rs_inhom_meanb%copy_data_to(gr%mesh%np, ff_rs_inhom_2b)
      call hamiltonian_mxll_update(hm, time=time)
      call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, ff_rs_inhom_2b, -inter_dt*M_HALF)
      ! add term U(time,time+dt/2)J(time)
      call batch_axpy(gr%mesh%np, M_ONE, ff_rs_inhom_2b, ff_rs_inhom_1b)
      call ff_rs_inhom_2b%end()
      call ff_rs_inhom_meanb%end()
    end if

    do ii = 1, inter_steps

      ! intermediate time
      inter_time = time + inter_dt * (ii-1)

      ! transformation of RS state into 3x3 or 4x4 representation
      call transform_rs_state_batch(hm, gr, st, rs_state, ff_rs_stateb, RS_TRANS_FORWARD)

      ! RS state propagation
      call hamiltonian_mxll_update(hm, time=inter_time)
      if (pml_check) then
        call pml_propagation_stage_1_batch(hm, gr, st, tr, ff_rs_stateb, ff_rs_state_pmlb)
      end if

      hm%cpml_hamiltonian = pml_check
      call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, ff_rs_stateb, dt)
      hm%cpml_hamiltonian = .false.

      if (pml_check) then
        call pml_propagation_stage_2_batch(hm, namespace, gr, st, tr, inter_time, inter_dt, delay, ff_rs_state_pmlb, ff_rs_stateb)
      end if

      !Below we add the contribution from the inhomogeneous terms
      if ((hm%ma_mx_coupling_apply) .or. hm%current_density_ext_flag) then
        if (tr%tr_etrs_approx == MXWLL_ETRS_FULL) then
          call ff_rs_stateb%copy_to(ff_rs_inhom_1b)
          call ff_rs_stateb%copy_to(ff_rs_inhom_2b)
          call ff_rs_stateb%copy_to(ff_rs_inhom_meanb)

          ! Interpolation of the external current
          do istate = 1, hm%dim
            call batch_set_state(ff_rs_inhom_meanb, istate, gr%mesh%np, ff_rs_inhom_t2(:, istate))
            call batch_set_state(ff_rs_inhom_1b, istate, gr%mesh%np, ff_rs_inhom_t1(:, istate))
          end do
          ! store t1 - t2 for the interpolation in mean
          call batch_axpy(gr%mesh%np, -M_ONE, ff_rs_inhom_1b, ff_rs_inhom_meanb)
          call ff_rs_inhom_1b%copy_data_to(gr%mesh%np, ff_rs_inhom_2b)
          call batch_axpy(gr%mesh%np, ii / TOFLOAT(inter_steps), ff_rs_inhom_meanb, ff_rs_inhom_2b)
          call batch_axpy(gr%mesh%np, (ii-1) / TOFLOAT(inter_steps), ff_rs_inhom_meanb, ff_rs_inhom_1b)

          hm%cpml_hamiltonian = .false.
          call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, ff_rs_inhom_1b, inter_dt)
          ! add terms U(time+dt,time)J(time) and J(time+dt)
          call batch_axpy(gr%mesh%np, -M_FOURTH * inter_dt, ff_rs_inhom_1b, ff_rs_stateb)
          call batch_axpy(gr%mesh%np, -M_FOURTH * inter_dt, ff_rs_inhom_2b, ff_rs_stateb)

          do istate = 1, hm%dim
            call batch_set_state(ff_rs_inhom_1b, istate, gr%mesh%np, ff_rs_inhom_t1(:, istate))
            call batch_set_state(ff_rs_inhom_2b, istate, gr%mesh%np, ff_rs_inhom_t2(:, istate))
          end do
          call batch_axpy(gr%mesh%np, M_ONE, ff_rs_inhom_2b, ff_rs_inhom_1b)
          call batch_scal(gr%mesh%np, M_HALF, ff_rs_inhom_1b)
          call ff_rs_inhom_1b%copy_data_to(gr%mesh%np, ff_rs_inhom_2b)

          call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, ff_rs_inhom_1b, inter_dt/M_TWO)
          call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, ff_rs_inhom_2b, -inter_dt/M_TWO)

          ! add terms U(time+dt/2,time)J(time) and U(time,time+dt/2)J(time+dt)
          call batch_axpy(gr%mesh%np, -M_FOURTH * inter_dt, ff_rs_inhom_1b, ff_rs_stateb)
          call batch_axpy(gr%mesh%np, -M_FOURTH * inter_dt, ff_rs_inhom_2b, ff_rs_stateb)

          call ff_rs_inhom_1b%end()
          call ff_rs_inhom_2b%end()
          call ff_rs_inhom_meanb%end()
        else if (tr%tr_etrs_approx == MXWLL_ETRS_CONST) then
          call batch_axpy(gr%mesh%np, -M_FOURTH * inter_dt, ff_rs_inhom_1b, ff_rs_stateb)
        end if
      end if

      ! PML convolution function update
      if (pml_check) then
        call cpml_conv_function_update(hm, gr, ff_rs_state_pmlb)
      end if

      ! back tranformation of RS state representation
      call transform_rs_state_batch(hm, gr, st, rs_state, ff_rs_stateb, RS_TRANS_BACKWARD)

      if (tr%bc_constant) then
        ! Propagation dt with H(inter_time+inter_dt) for constant boundaries
        if (st%rs_state_const_external) then
          call spatial_constant_calculation(tr%bc_constant, st, gr, hm, inter_time, inter_dt, delay, rs_state)
        end if
        call constant_boundaries_calculation(tr%bc_constant, hm%bc, hm, st, rs_state)
      end if

      ! Propagation dt with H(inter_time+inter_dt) for PEC mirror boundaries
      call mirror_pec_boundaries_calculation(hm%bc, st, rs_state)

      ! Propagation dt with H(inter_time+inter_dt) for PMC mirror boundaries
      call mirror_pmc_boundaries_calculation(hm%bc, st, rs_state)

      ! Apply mask absorbing boundaries
      call mask_absorbing_boundaries(namespace, gr, hm, st, tr, inter_time, inter_dt, delay, rs_state)

      if (tr%bc_plane_waves) then
        ! Propagation dt with H(inter_time+inter_dt) for plane waves boundaries
        call plane_waves_boundaries_calculation(hm, st, gr%mesh, inter_time+inter_dt, delay ,rs_state)
      end if

    end do

    if (tr%tr_etrs_approx == OPTION__MAXWELLTDETRSAPPROX__CONST_STEPS) then
      call ff_rs_inhom_1b%end()
    end if

    call ff_rs_stateb%end()

    if (pml_check) then
      call ff_rs_state_pmlb%end()
    end if

    call profiling_out(prof)

    POP_SUB(mxll_propagation_step)
  end subroutine mxll_propagation_step

  ! ---------------------------------------------------------
  subroutine set_medium_rs_state(st, gr, hm)
    type(states_mxll_t),      intent(inout) :: st
    type(grid_t),             intent(in)    :: gr
    type(hamiltonian_mxll_t), intent(in)    :: hm

    integer :: ip, ip_in, il, idim
    type(profile_t), save :: prof

    PUSH_SUB(set_medium_rs_state)

    call profiling_in(prof, 'SET_MEDIUM_RS_STATE')

    SAFE_ALLOCATE(st%ep(1:gr%mesh%np_part))
    SAFE_ALLOCATE(st%mu(1:gr%mesh%np_part))
    st%ep = P_ep
    st%mu = P_mu
    if (hm%calc_medium_box) then
      do il = 1, size(hm%medium_boxes)
        do ip_in = 1, hm%medium_boxes(il)%points_number
          ip = hm%medium_boxes(il)%points_map(ip_in)
          st%ep(ip) = hm%medium_boxes(il)%ep(ip_in)
          st%mu(ip) = hm%medium_boxes(il)%mu(ip_in)
        end do
      end do
    end if

    do idim = 1, st%dim
      if (hm%bc%bc_type(idim) == MXLL_BC_MEDIUM) then
        do ip_in = 1, hm%bc%medium(idim)%points_number
          ip = hm%bc%medium(idim)%points_map(ip_in)
          st%ep(ip) = hm%bc%medium(idim)%ep(ip_in)
          st%mu(ip) = hm%bc%medium(idim)%mu(ip_in)
        end do
      end if
    end do

    call profiling_out(prof)

    POP_SUB(set_medium_rs_state)
  end subroutine set_medium_rs_state

  ! ---------------------------------------------------------
  subroutine transform_rs_state_batch(hm, gr, st, rs_state, ff_rs_stateb, sign)
    type(hamiltonian_mxll_t), intent(in)         :: hm
    type(grid_t),             intent(in)         :: gr
    type(states_mxll_t),      intent(in)         :: st
    CMPLX,                    intent(inout)      :: rs_state(:,:)
    type(batch_t),            intent(inout)      :: ff_rs_stateb
    integer,                  intent(in)         :: sign

    CMPLX, allocatable :: rs_state_tmp(:,:)
    type(profile_t), save :: prof
    integer :: ii, np

    PUSH_SUB(transform_rs_state_batch)

    call profiling_in(prof, 'TRANSFORM_RS_STATE')

    ASSERT(sign == RS_TRANS_FORWARD .or. sign == RS_TRANS_BACKWARD)

    np = gr%mesh%np

    if (hm%operator == FARADAY_AMPERE_MEDIUM) then
      if (sign == RS_TRANS_FORWARD) then
        ! 3 to 6
        do ii = 1, 3
          call batch_set_state(ff_rs_stateb, ii, np, rs_state(:, ii))
          call batch_set_state(ff_rs_stateb, ii+3, np, conjg(rs_state(:, ii)))
        end do
      else
        ! 6 to 3
        SAFE_ALLOCATE(rs_state_tmp(gr%mesh%np, st%dim))
        do ii = 1, 3
          call batch_get_state(ff_rs_stateb, ii, np, rs_state(:, ii))
          call batch_get_state(ff_rs_stateb, ii+3, np, rs_state_tmp(:, ii))
          rs_state(1:np, ii) = M_HALF * (rs_state(1:np, ii) + conjg(rs_state_tmp(1:np, ii)))
        end do
        SAFE_DEALLOCATE_A(rs_state_tmp)
      end if

    else if (hm%operator == FARADAY_AMPERE_GAUSS) then
      if (sign == RS_TRANS_FORWARD) then
        ! 3 to 4
        call batch_set_state(ff_rs_stateb, 1, np, -rs_state(1:np, 1)+rs_state(1:np, 2))
        call batch_set_state(ff_rs_stateb, 2, np, rs_state(1:np, 3))
        call batch_set_state(ff_rs_stateb, 3, np, rs_state(1:np, 3))
        call batch_set_state(ff_rs_stateb, 4, np, rs_state(1:np, 1)+rs_state(1:np, 2))
      else
        ! 4 to 3
        SAFE_ALLOCATE(rs_state_tmp(gr%mesh%np, 2))
        call batch_get_state(ff_rs_stateb, 1, np, rs_state_tmp(:, 1))
        call batch_get_state(ff_rs_stateb, 4, np, rs_state_tmp(:, 2))
        rs_state(1:np, 1) = M_HALF * (-rs_state_tmp(1:np, 1) + rs_state_tmp(1:np, 2))
        rs_state(1:np, 2) = M_HALF * (-rs_state_tmp(1:np, 1) - rs_state_tmp(1:np, 2))
        call batch_get_state(ff_rs_stateb, 2, np, rs_state_tmp(:, 1))
        call batch_get_state(ff_rs_stateb, 3, np, rs_state_tmp(:, 2))
        rs_state(1:np, 3) = M_HALF * (rs_state_tmp(1:np, 1) + rs_state_tmp(1:np, 2))
        SAFE_DEALLOCATE_A(rs_state_tmp)
      end if

    else
      if (sign == RS_TRANS_FORWARD) then
        ! 3 to 3
        do ii = 1, 3
          call batch_set_state(ff_rs_stateb, ii, np, rs_state(:, ii))
        end do
      else
        ! 3 to 3
        do ii = 1, 3
          call batch_get_state(ff_rs_stateb, ii, np, rs_state(:, ii))
        end do
      end if
    end if

    call profiling_out(prof)

    POP_SUB(transform_rs_state_batch)

  end subroutine transform_rs_state_batch

  ! ---------------------------------------------------------
  subroutine transform_rs_densities(hm, mesh, rs_charge_density, rs_current_density, ff_density, sign)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    type(mesh_t),             intent(in)    :: mesh
    CMPLX,                    intent(inout) :: rs_charge_density(:)
    CMPLX,                    intent(inout) :: rs_current_density(:,:)
    CMPLX,                    intent(inout) :: ff_density(:,:)
    integer,                  intent(in)    :: sign

    type(profile_t), save :: prof

    ASSERT(sign == RS_TRANS_FORWARD .or. sign == RS_TRANS_BACKWARD)
    ASSERT(size(rs_charge_density) == mesh%np .or. size(rs_charge_density) == mesh%np_part)
    ASSERT(size(rs_current_density, dim=1) == size(rs_charge_density))
    ASSERT(size(rs_current_density, dim=2) == 3)

    PUSH_SUB(transform_rs_densities)

    call profiling_in(prof, 'TRANSFORM_RS_DENSITIES')

    if (hm%operator == FARADAY_AMPERE_MEDIUM) then
      if (sign == RS_TRANS_FORWARD) then
        call transform_rs_densities_to_6x6_rs_densities_forward(mesh, rs_charge_density, &
                               rs_current_density, ff_density)
      else
        call transform_rs_densities_to_6x6_rs_densities_backward(mesh, ff_density, &
                               rs_charge_density, rs_current_density)
      end if
    else if (hm%operator == FARADAY_AMPERE_GAUSS) then
      if (sign == RS_TRANS_FORWARD) then
        call transform_rs_densities_to_4x4_rs_densities_forward(mesh, rs_charge_density,&
            rs_current_density, ff_density)
      else
        call transform_rs_densities_to_4x4_rs_densities_backward(mesh, ff_density, rs_charge_density,&
            rs_current_density)
      end if
    else
      if (sign == RS_TRANS_FORWARD) then
        ff_density(1:mesh%np, 1:3) = rs_current_density(1:mesh%np, 1:3)
      else
        rs_current_density(1:mesh%np, 1:3) = ff_density(1:mesh%np, 1:3)
      end if
    end if

    call profiling_out(prof)

    POP_SUB(transform_rs_densities)

  end subroutine transform_rs_densities

  !----------------------------------------------------------
  subroutine transform_rs_densities_to_6x6_rs_densities_forward(mesh, rs_charge_density, rs_current_density, rs_density_6x6)
    type(mesh_t),             intent(in)    :: mesh
    CMPLX,                    intent(in)    :: rs_charge_density(:)
    CMPLX,                    intent(in)    :: rs_current_density(:,:)
    CMPLX,                    intent(inout) :: rs_density_6x6(:,:)

    integer :: ii

    ASSERT(size(rs_current_density, dim=2) == 3)
    ASSERT(size(rs_density_6x6, dim=2) == 6)

    ! no push_sub, called to frequently
    do ii = 1, 3
      rs_density_6x6(1:mesh%np, ii) = rs_current_density(1:mesh%np, ii)
      rs_density_6x6(1:mesh%np, ii+3) = rs_current_density(1:mesh%np, ii)
    end do

  end subroutine transform_rs_densities_to_6x6_rs_densities_forward

  !----------------------------------------------------------
  subroutine transform_rs_densities_to_6x6_rs_densities_backward(mesh, rs_density_6x6, rs_charge_density, rs_current_density)
    type(mesh_t),             intent(in)    :: mesh
    CMPLX,                    intent(in)    :: rs_density_6x6(:,:)
    CMPLX,                    intent(inout) :: rs_charge_density(:)
    CMPLX,                    intent(inout) :: rs_current_density(:,:)

    integer :: ii

    ASSERT(size(rs_current_density, dim=2) == 3)
    ASSERT(size(rs_density_6x6, dim=2) == 6)

    ! no push_sub, called to frequently
    do ii = 1, 3
      rs_current_density(1:mesh%np, ii) = M_HALF * &
                TOFLOAT(rs_density_6x6(1:mesh%np, ii) + rs_density_6x6(1:mesh%np, ii+3))
    end do

  end subroutine transform_rs_densities_to_6x6_rs_densities_backward

  !----------------------------------------------------------
  subroutine transform_rs_densities_to_4x4_rs_densities_forward(mesh, rs_charge_density, rs_current_density, rs_density_4x4)
    type(mesh_t),             intent(in)    :: mesh
    CMPLX,                    intent(in)    :: rs_charge_density(:)
    CMPLX,                    intent(in)    :: rs_current_density(:,:)
    CMPLX,                    intent(inout) :: rs_density_4x4(:,:)

    ! no push_sub, called to frequently
    rs_density_4x4(1:mesh%np, 1) = M_z1 * (-rs_current_density(1:mesh%np, 1) + rs_current_density(1:mesh%np, 2))
    rs_density_4x4(1:mesh%np, 2) = M_z1 * (rs_current_density(1:mesh%np, 3) - rs_charge_density(1:mesh%np))
    rs_density_4x4(1:mesh%np, 3) = M_z1 * (rs_current_density(1:mesh%np, 3) + rs_charge_density(1:mesh%np))
    rs_density_4x4(1:mesh%np, 4) = M_z1 * (rs_current_density(1:mesh%np, 1) + rs_current_density(1:mesh%np, 2))

  end subroutine transform_rs_densities_to_4x4_rs_densities_forward

  !----------------------------------------------------------
  subroutine transform_rs_densities_to_4x4_rs_densities_backward(mesh, rs_density_4x4, rs_charge_density, rs_current_density)
    type(mesh_t),             intent(in)    :: mesh
    CMPLX,                    intent(in)    :: rs_density_4x4(:,:)
    CMPLX,                    intent(inout) :: rs_charge_density(:)
    CMPLX,                    intent(inout) :: rs_current_density(:,:)

    ! no push_sub, called to frequently
    rs_charge_density(1:mesh%np)     = M_z1 * M_HALF * (-rs_density_4x4(1:mesh%np,2) + rs_density_4x4(1:mesh%np,3))
    rs_current_density(1:mesh%np, 1) = M_z1 * M_HALF * (-rs_density_4x4(1:mesh%np,1) + rs_density_4x4(1:mesh%np,4))
    rs_current_density(1:mesh%np, 2) = M_zI * M_HALF * (-rs_density_4x4(1:mesh%np,1) - rs_density_4x4(1:mesh%np,4))
    rs_current_density(1:mesh%np, 3) = M_z1 * M_HALF * (-rs_density_4x4(1:mesh%np,2) + rs_density_4x4(1:mesh%np,3))

  end subroutine transform_rs_densities_to_4x4_rs_densities_backward

  !----------------------------------------------------------
  subroutine calculate_matter_longitudinal_field(gr_mxll, st_mxll, hm_mxll, gr_elec, st_elec, hm_elec, rs_state_matter, ions)
    type(grid_t),                  intent(in)    :: gr_mxll
    type(states_mxll_t),           intent(in)    :: st_mxll
    type(hamiltonian_mxll_t),      intent(in)    :: hm_mxll
    type(grid_t),                  intent(in)    :: gr_elec
    type(states_elec_t),           intent(in)    :: st_elec
    type(hamiltonian_elec_t),      intent(in)    :: hm_elec
    CMPLX,                         intent(inout) :: rs_state_matter(:,:)
    type(ions_t),        optional, intent(in)    :: ions

    CMPLX, allocatable :: tmp_pot_mx_gr(:,:), tmp_grad_mx_gr(:,:)

    SAFE_ALLOCATE(tmp_pot_mx_gr(1:gr_mxll%mesh%np_part,1))
    SAFE_ALLOCATE(tmp_grad_mx_gr(1:gr_mxll%mesh%np,1:gr_mxll%sb%dim))
    ! this subroutine needs the matter part

    PUSH_SUB(calculate_matter_longitudinal_field)

    tmp_pot_mx_gr(:,:) = M_ZERO
    tmp_grad_mx_gr(:,:) = M_ZERO
    call zderivatives_grad(gr_mxll%der, tmp_pot_mx_gr(:,1), tmp_grad_mx_gr(:,:), set_bc = .false.)
    tmp_grad_mx_gr = - tmp_grad_mx_gr

    rs_state_matter = M_z0
    call build_rs_state(real(tmp_grad_mx_gr(1:gr_mxll%mesh%np,:)), aimag(tmp_grad_mx_gr(1:gr_mxll%mesh%np,:)), st_mxll%rs_sign, &
      rs_state_matter(1:gr_mxll%mesh%np,:), gr_mxll%mesh, st_mxll%ep(1:gr_mxll%mesh%np), st_mxll%mu(1:gr_mxll%mesh%np), &
      gr_mxll%mesh%np)

    SAFE_DEALLOCATE_A(tmp_pot_mx_gr)
    SAFE_DEALLOCATE_A(tmp_grad_mx_gr)

    POP_SUB(calculate_matter_longitudinal_field)
  end subroutine calculate_matter_longitudinal_field

  !----------------------------------------------------------
  subroutine get_vector_pot_and_transverse_field(trans_calc_method, gr_mxll, hm_mxll, st_mxll, tr_mxll, hm, st, &
    poisson_solver, time, field, transverse_field, vector_potential)
    integer,                    intent(in)    :: trans_calc_method
    type(grid_t),               intent(in)    :: gr_mxll
    type(hamiltonian_mxll_t),   intent(in)    :: hm_mxll
    type(states_mxll_t),        intent(in)    :: st_mxll
    type(propagator_mxll_t),    intent(in)    :: tr_mxll
    type(hamiltonian_elec_t),   intent(in)    :: hm
    type(states_elec_t),        intent(in)    :: st
    type(poisson_t),            intent(in)    :: poisson_solver
    FLOAT,                      intent(in)    :: time
    CMPLX,                      intent(inout) :: field(:,:)
    CMPLX,                      intent(inout) :: transverse_field(:,:)
    FLOAT,                      intent(inout) :: vector_potential(:,:)

    integer            :: np

    PUSH_SUB(get_vector_pot_and_transverse_field)

    transverse_field = M_z0
    vector_potential = M_ZERO

    np = gr_mxll%mesh%np

    if (hm_mxll%ma_mx_coupling) then

       ! check what other transverse field methods are needed

       ! trans_calc_method == OPTION__MAXWELLTRANSFIELDCALCULATIONMETHOD__TRANS_FIELD_POISSON

        ! plane waves subtraction
        if (tr_mxll%bc_plane_waves .and. hm_mxll%plane_waves_apply) then
          transverse_field(1:np,:) = field(1:np,:) - st_mxll%rs_state_plane_waves(1:np,:)
        else
          transverse_field(1:np,:) = field(1:np,:)
        end if
        ! apply helmholtz decomposition for transverse field
        call maxwell_helmholtz_decomposition_trans_field(poisson_solver, gr_mxll, hm_mxll, hm, st_mxll, transverse_field)
        ! plane waves addition
        if (tr_mxll%bc_plane_waves .and. hm_mxll%plane_waves_apply) then
          transverse_field(1:np,:) = transverse_field(1:np,:) + st_mxll%rs_state_plane_waves(1:np,:)
        end if

    else

      transverse_field(1:np,:) = field

    end if


    POP_SUB(get_vector_pot_and_transverse_field)

  end subroutine  get_vector_pot_and_transverse_field

  ! ---------------------------------------------------------
  subroutine calculate_vector_potential(poisson_solver, gr, st, field, vector_potential)
    type(poisson_t),            intent(in)    :: poisson_solver
    type(grid_t),               intent(in)    :: gr
    type(states_mxll_t),        intent(in)    :: st
    CMPLX,                      intent(in)    :: field(:,:)
    FLOAT,                      intent(inout) :: vector_potential(:,:)

    integer :: idim
    FLOAT, allocatable :: dtmp(:,:)

    SAFE_ALLOCATE(dtmp(1:gr%mesh%np_part,1:3))

    dtmp = M_ZERO

    call get_magnetic_field_state(field, gr%mesh, st%rs_sign, vector_potential, st%mu, gr%mesh%np_part)
    dtmp = vector_potential
    call dderivatives_curl(gr%der, dtmp, vector_potential, set_bc = .false.)
    do idim=1, st%dim
      call dpoisson_solve(poisson_solver, dtmp(:,idim), vector_potential(:,idim), .true.)
    end do
    vector_potential = M_ONE / (M_FOUR * M_PI) * vector_potential

    SAFE_DEALLOCATE_A(dtmp)

  end subroutine calculate_vector_potential

  ! ---------------------------------------------------------
  subroutine derivatives_boundary_mask(bc, mesh, hm)
    type(bc_mxll_t),  intent(inout)      :: bc
    type(mesh_t),        intent(in)      :: mesh
    type(hamiltonian_mxll_t), intent(in) :: hm

    integer :: ip, ip_in, point_info, idim, dim
    FLOAT   :: bounds(2, mesh%sb%dim), xx(mesh%sb%dim)
    FLOAT   :: ddv(mesh%sb%dim), tmp(mesh%sb%dim), width(mesh%sb%dim)
    FLOAT, allocatable :: mask(:)
    type(profile_t), save :: prof

    PUSH_SUB(derivatives_boundary_mask)

    call profiling_in(prof, 'DERIVATIVES_BOUNDARY_MASK')
    dim = mesh%sb%dim

    if (hm%bc_zero .or. hm%bc_constant .or. hm%bc_plane_waves) then
      bounds(1,1:dim) = ( mesh%idx%nr(2,1:dim) - 2 * mesh%idx%enlarge(1:dim) ) * mesh%spacing(1:dim)
      bounds(2,1:dim) = ( mesh%idx%nr(2,1:dim) -     mesh%idx%enlarge(1:dim) ) * mesh%spacing(1:dim)
    end if

    ip_in=0
    do ip=1, mesh%np
      xx(1:dim) = mesh%x(ip,1:dim)
      if ((abs(xx(1))<=bounds(2,1)) .and. (abs(xx(2))<=bounds(2,2)) .and. (abs(xx(3))<=bounds(2,3))) then
        if ((abs(xx(1))>bounds(1,1)) .or. (abs(xx(2))>bounds(1,2)) .or. (abs(xx(3))>bounds(1,3))) then
          point_info = 1
        else
          point_info = 0
        end if
      else
        point_info = -1
      end if
      if (point_info == 1) then
        ip_in = ip_in + 1
      end if
    end do
    bc%der_bndry_mask_points_number = ip_in
    SAFE_ALLOCATE(bc%der_bndry_mask(1:ip_in))
    SAFE_ALLOCATE(bc%der_bndry_mask_points_map(1:ip_in))

    ip_in=0
    do ip=1, mesh%np
      xx(1:dim) = mesh%x(ip,1:dim)
      if ((abs(xx(1))<=bounds(2,1)) .and. (abs(xx(2))<=bounds(2,2)) .and. (abs(xx(3))<=bounds(2,3))) then
        if ((abs(xx(1))>bounds(1,1)) .or. (abs(xx(2))>bounds(1,2)) .or. (abs(xx(3))>bounds(1,3))) then
          point_info = 1
        else
          point_info = 0
        end if
      else
        point_info = -1
      end if
      if (point_info == 1) then
        ip_in = ip_in + 1
        bc%der_bndry_mask_points_map(ip_in) = ip
      end if
    end do

    SAFE_ALLOCATE(mask(mesh%np))
    mask(:) = M_ONE
    width(:) = bounds(2,:) - bounds(1,:)
    tmp(:)   = M_ZERO

    do ip=1, mesh%np
      tmp = M_ONE
      mask(ip) = M_ONE
      ddv(1:dim) = abs(mesh%x(ip,1:dim)) - bounds(1,1:dim)
      do idim=1, mesh%sb%dim
        if(ddv(idim) >= M_ZERO ) then
          if (ddv(idim)  <=  width(idim)) then
            tmp(idim) = M_ONE - sin(ddv(idim) * M_PI / (M_TWO * (width(idim)) ))**2
          else
            tmp(idim) = M_ONE
          end if
        end if
        mask(ip) = mask(ip) * tmp(idim)
      end do
    end do

    do idim=1, mesh%sb%dim
      do ip_in=1, bc%der_bndry_mask_points_number
        ip = bc%der_bndry_mask_points_map(ip_in)
        bc%der_bndry_mask(ip_in) = mask(ip)
      end do
    end do

    SAFE_DEALLOCATE_A(mask)
    call profiling_out(prof)

    POP_SUB(derivatives_boundary_mask)
  end subroutine derivatives_boundary_mask

  !----------------------------------------------------------
  subroutine energy_density_calc(gr, st, rs_field, energy_dens, e_energy_dens, b_energy_dens, plane_waves_check, &
    rs_field_plane_waves, energy_dens_plane_waves)
    type(grid_t),        intent(in)    :: gr
    type(states_mxll_t), intent(in)    :: st
    CMPLX,               intent(in)    :: rs_field(:,:)
    FLOAT,               intent(inout) :: energy_dens(:)
    FLOAT,               intent(inout) :: e_energy_dens(:)
    FLOAT,               intent(inout) :: b_energy_dens(:)
    logical,   optional, intent(in)    :: plane_waves_check
    CMPLX,     optional, intent(in)    :: rs_field_plane_waves(:,:)
    FLOAT,     optional, intent(inout) :: energy_dens_plane_waves(:)

    integer            :: idim, ip
    type(profile_t), save :: prof

    PUSH_SUB(energy_density_calc)

    call profiling_in(prof, 'ENERGY_DENSITY_CALC')

    e_energy_dens(:) = M_ZERO
    b_energy_dens(:) = M_ZERO
    do ip = 1, gr%mesh%np
      do idim = 1, st%dim
        e_energy_dens(ip) = e_energy_dens(ip) + real(rs_field(ip,idim))**2
        b_energy_dens(ip) = b_energy_dens(ip) + aimag(rs_field(ip,idim))**2
      end do
      energy_dens(ip) = e_energy_dens(ip) + b_energy_dens(ip)
    end do

    if (present(rs_field_plane_waves) .and. present(energy_dens_plane_waves) .and. plane_waves_check) then
      energy_dens_plane_waves(:) = M_ZERO
      do ip = 1, gr%mesh%np
        do idim = 1, st%dim
          energy_dens_plane_waves(ip) = energy_dens_plane_waves(ip) &
                + real(conjg(rs_field_plane_waves(ip,idim)) * rs_field_plane_waves(ip,idim))
        end do
      end do
    end if

    call profiling_out(prof)

    POP_SUB(energy_density_calc)
  end subroutine energy_density_calc

  !----------------------------------------------------------
  subroutine energy_mxll_calc(gr, st, hm, energy_mxll, rs_field, rs_field_plane_waves)
    type(grid_t),             intent(in)  :: gr
    type(states_mxll_t),      intent(in)  :: st
    type(hamiltonian_mxll_t), intent(in)  :: hm
    type(energy_mxll_t),      intent(inout) :: energy_mxll
    CMPLX,                    intent(in)  :: rs_field(:,:)
    CMPLX, optional,          intent(in)  :: rs_field_plane_waves(:,:)

    type(profile_t), save :: prof

    PUSH_SUB(energy_mxll_calc)

    call profiling_in(prof, 'ENERGY_MXLL_CALC')

    call energy_density_calc(gr, st, rs_field, energy_mxll%energy_density, energy_mxll%e_energy_density, &
         energy_mxll%b_energy_density, hm%plane_waves, rs_field_plane_waves, energy_mxll%energy_density_plane_waves)

    energy_mxll%energy    = dmf_integrate(gr%mesh, energy_mxll%energy_density, mask=st%inner_points_mask)
    energy_mxll%e_energy  = dmf_integrate(gr%mesh, energy_mxll%e_energy_density, mask=st%inner_points_mask)
    energy_mxll%b_energy  = dmf_integrate(gr%mesh, energy_mxll%b_energy_density, mask=st%inner_points_mask)
    if (present(rs_field_plane_waves) .and. hm%plane_waves) then
      energy_mxll%energy_plane_waves = dmf_integrate(gr%mesh, energy_mxll%energy_density_plane_waves, mask=st%inner_points_mask)
    else
      energy_mxll%energy_plane_waves = M_ZERO
    end if

    energy_mxll%boundaries = dmf_integrate(gr%mesh, energy_mxll%energy_density, mask=st%boundary_points_mask)

    call profiling_out(prof)

    POP_SUB(energy_mxll_calc)
  end subroutine energy_mxll_calc

  ! ---------------------------------------------------------
  subroutine mask_absorbing_boundaries(namespace, gr, hm, st, tr, time, dt, time_delay, rs_state)
    type(namespace_t),          intent(in)    :: namespace
    type(grid_t),               intent(in)    :: gr
    type(hamiltonian_mxll_t),   intent(inout) :: hm
    type(states_mxll_t),        intent(inout) :: st
    type(propagator_mxll_t),    intent(inout) :: tr
    FLOAT,                      intent(in)    :: time
    FLOAT,                      intent(in)    :: dt
    FLOAT,                      intent(in)    :: time_delay
    CMPLX,                      intent(inout) :: rs_state(:,:)

    integer            :: ip, ip_in, idim
    logical            :: mask_check = .false.
    type(profile_t), save :: prof

    PUSH_SUB(mask_absorbing_boundaries)

    call profiling_in(prof, 'MASK_ABSORBING_BOUNDARIES')

    do idim = 1, 3
      if (hm%bc%bc_ab_type(idim) == OPTION__MAXWELLABSORBINGBOUNDARIES__MASK) then
        mask_check = .true.
      end if
    end do

    if (mask_check) then
      if (tr%bc_plane_waves .and. hm%plane_waves_apply) then
        call plane_waves_propagation(hm, tr, namespace, st, gr, time, dt, time_delay)
        rs_state = rs_state - st%rs_state_plane_waves
        call maxwell_mask(hm, rs_state)
        rs_state = rs_state + st%rs_state_plane_waves
      else if (tr%bc_constant .and. hm%spatial_constant_apply) then
        !call constant_at_absorbing_boundaries_calculation(st, hm%bc)
        call constant_boundaries_calculation(tr%bc_constant, hm%bc, hm, st, rs_state)
        do ip_in=1, hm%bc%constant_points_number
          ip = hm%bc%constant_points_map(ip_in)
          rs_state(ip,:) = rs_state(ip,:) - st%rs_state_const(:)
        end do
        call maxwell_mask(hm, rs_state)
        do ip_in=1, hm%bc%constant_points_number
          ip = hm%bc%constant_points_map(ip_in)
          rs_state(ip,:) = rs_state(ip,:) + st%rs_state_const(:)
        end do
      else
        call maxwell_mask(hm, rs_state)
      end if
    end if

    call profiling_out(prof)

    POP_SUB(mask_absorbing_boundaries)
  end subroutine mask_absorbing_boundaries

  ! ---------------------------------------------------------
  subroutine maxwell_mask(hm, rs_state)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    CMPLX,                    intent(inout) :: rs_state(:,:)

    integer :: ip, ip_in, idim
    type(profile_t), save :: prof

    PUSH_SUB(maxwell_mask)

    call profiling_in(prof, 'MAXWELL_MASK')

    do idim = 1, 3
      if (hm%bc%bc_ab_type(idim) == OPTION__MAXWELLABSORBINGBOUNDARIES__MASK) then
        do ip_in = 1, hm%bc%mask_points_number(idim)
          ip = hm%bc%mask_points_map(ip_in,idim)
          rs_state(ip,:) = rs_state(ip,:) * hm%bc%mask(ip_in,idim)
        end do
      end if
    end do

    call profiling_out(prof)

    POP_SUB(maxwell_mask)
  end subroutine maxwell_mask

  ! ---------------------------------------------------------
  subroutine pml_propagation_stage_1_batch(hm, gr, st, tr, ff_rs_stateb, ff_rs_state_pmlb)
    type(hamiltonian_mxll_t),   intent(inout) :: hm
    type(grid_t),               intent(in)    :: gr
    type(states_mxll_t),        intent(inout) :: st
    type(propagator_mxll_t),    intent(inout) :: tr
    type(batch_t),              intent(in)    :: ff_rs_stateb
    type(batch_t),              intent(inout) :: ff_rs_state_pmlb

    integer            :: ii
    CMPLX, allocatable :: rs_state_constant(:,:)
    type(profile_t), save :: prof

    PUSH_SUB(pml_propagation_stage_1_batch)

    call profiling_in(prof, 'PML_PROP_STAGE_1_BATCH')

    if (tr%bc_plane_waves .and. hm%plane_waves_apply) then
      call transform_rs_state_batch(hm, gr, st, st%rs_state_plane_waves, &
        ff_rs_state_pmlb, RS_TRANS_FORWARD)
      call batch_xpay(gr%mesh%np, ff_rs_stateb, -M_ONE, ff_rs_state_pmlb)
    else if (tr%bc_constant .and. hm%spatial_constant_apply) then
      ! this could be optimized: right now we broadcast the constant value
      !  to the full mesh to be able to use the batch functions easily.
      ! in principle, we would need to do the transform only for one point
      !  and then subtract that value from all points of the state
      SAFE_ALLOCATE(rs_state_constant(1:gr%mesh%np,1:3))
      do ii = 1, 3
        rs_state_constant(1:gr%mesh%np, ii) = st%rs_state_const(ii)
      end do

      call transform_rs_state_batch(hm, gr, st, rs_state_constant, &
        ff_rs_state_pmlb, RS_TRANS_FORWARD)
      call batch_xpay(gr%mesh%np, ff_rs_stateb, -M_ONE, ff_rs_state_pmlb)

      SAFE_DEALLOCATE_A(rs_state_constant)
    else
      ! this copy should not be needed
      call ff_rs_stateb%copy_data_to(gr%mesh%np, ff_rs_state_pmlb)
    end if

    call profiling_out(prof)

    POP_SUB(pml_propagation_stage_1_batch)
  end subroutine pml_propagation_stage_1_batch

  ! ---------------------------------------------------------
  subroutine pml_propagation_stage_2_batch(hm, namespace, gr, st, tr, time, dt, time_delay, ff_rs_state_pmlb, ff_rs_stateb)
    type(hamiltonian_mxll_t),   intent(inout) :: hm
    type(namespace_t),          intent(in)    :: namespace
    type(grid_t),               intent(in)    :: gr
    type(states_mxll_t),        intent(inout) :: st
    type(propagator_mxll_t),    intent(inout) :: tr
    FLOAT,                      intent(in)    :: time
    FLOAT,                      intent(in)    :: dt
    FLOAT,                      intent(in)    :: time_delay
    type(batch_t),              intent(inout) :: ff_rs_state_pmlb
    type(batch_t),              intent(inout) :: ff_rs_stateb

    integer            :: ip, ip_in, ii, ff_dim
    CMPLX, allocatable :: rs_state_constant(:,:), ff_rs_state_constant(:,:)
    type(profile_t), save :: prof
    type(batch_t) :: ff_rs_state_plane_wavesb, ff_rs_constantb

    PUSH_SUB(pml_propagation_stage_2_batch)

    call profiling_in(prof, 'PML_PROP_STAGE_2_BATCH')

    if (tr%bc_plane_waves .and. hm%plane_waves_apply) then
      hm%cpml_hamiltonian = .true.
      call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, ff_rs_state_pmlb, dt)
      hm%cpml_hamiltonian = .false.
      call plane_waves_propagation(hm, tr, namespace, st, gr, time, dt, time_delay)

      call ff_rs_stateb%copy_to(ff_rs_state_plane_wavesb)
      call transform_rs_state_batch(hm, gr, st, st%rs_state_plane_waves, ff_rs_state_plane_wavesb, RS_TRANS_FORWARD)

      select case(ff_rs_stateb%status())
      case(BATCH_NOT_PACKED)
        do ii = 1, ff_rs_stateb%nst_linear
          do ip_in = 1, hm%bc%plane_wave%points_number
            ip = hm%bc%plane_wave%points_map(ip_in)
            ff_rs_stateb%zff_linear(ip, ii) = ff_rs_state_pmlb%zff_linear(ip, ii) + &
              ff_rs_state_plane_wavesb%zff_linear(ip, ii)
          end do
        end do
      case(BATCH_PACKED)
        do ip_in = 1, hm%bc%plane_wave%points_number
          ip = hm%bc%plane_wave%points_map(ip_in)
          do ii = 1, ff_rs_stateb%nst_linear
            ff_rs_stateb%zff_pack(ii, ip) = ff_rs_state_pmlb%zff_pack(ii, ip) + &
              ff_rs_state_plane_wavesb%zff_pack(ii, ip)
          end do
        end do
      end select

    else if (tr%bc_constant .and. hm%spatial_constant_apply) then
      hm%cpml_hamiltonian = .true.
      call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, ff_rs_state_pmlb, dt)
      hm%cpml_hamiltonian = .false.

      call ff_rs_stateb%copy_to(ff_rs_constantb)
      ff_dim = ff_rs_stateb%nst_linear
      SAFE_ALLOCATE(rs_state_constant(1:gr%mesh%np, 1:st%dim))
      ! copy the value to the full mesh to be able to use batches
      ! this is in principle unneeded, but otherwise we could not use batches...
      do ii = 1, st%dim
        rs_state_constant(1:gr%mesh%np, ii) = st%rs_state_const(ii)
      end do

      call transform_rs_state_batch(hm, gr, st, rs_state_constant, ff_rs_constantb, RS_TRANS_FORWARD)
      select case(ff_rs_stateb%status())
      case(BATCH_NOT_PACKED)
        do ii = 1, ff_rs_stateb%nst_linear
          do ip_in = 1, hm%bc%constant_points_number
            ip = hm%bc%constant_points_map(ip_in)
            ff_rs_stateb%zff_linear(ip, ii) = ff_rs_state_pmlb%zff_linear(ip, ii) + &
              ff_rs_constantb%zff_linear(ip, ii)
          end do
        end do
      case(BATCH_PACKED)
        do ip_in = 1, hm%bc%constant_points_number
          ip = hm%bc%constant_points_map(ip_in)
          do ii = 1, ff_rs_stateb%nst_linear
            ff_rs_stateb%zff_pack(ii, ip) = ff_rs_state_pmlb%zff_pack(ii, ip) + &
              ff_rs_constantb%zff_pack(ii, ip)
          end do
        end do
      end select
 
      SAFE_DEALLOCATE_A(rs_state_constant)
      SAFE_DEALLOCATE_A(ff_rs_state_constant)
    end if

    call profiling_out(prof)

    POP_SUB(pml_propagation_stage_2_batch)
  end subroutine pml_propagation_stage_2_batch

  ! ---------------------------------------------------------
  subroutine cpml_conv_function_update(hm, gr, ff_rs_state_pmlb)
    type(hamiltonian_mxll_t), intent(inout) :: hm
    type(grid_t),             intent(in)    :: gr
    type(batch_t),            intent(inout) :: ff_rs_state_pmlb

    integer :: ff_dim, istate
    type(profile_t), save :: prof
    CMPLX, allocatable :: ff_rs_state_pml(:,:)

    PUSH_SUB(cpml_conv_function_update)

    call profiling_in(prof, 'CPML_CONV_FUNCTION_UPDATE')

    if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RS) then
      call cpml_conv_function_update_via_riemann_silberstein(hm, gr, ff_rs_state_pmlb)
    else if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__EM) then
      ff_dim = hm%dim
      SAFE_ALLOCATE(ff_rs_state_pml(1:gr%mesh%np_part, ff_dim))
      do istate = 1, hm%dim
        call batch_get_state(ff_rs_state_pmlb, istate, gr%mesh%np_part, ff_rs_state_pml(:, istate))
      end do
      call cpml_conv_function_update_via_e_b_fields(hm, gr, ff_rs_state_pml)
      SAFE_DEALLOCATE_A(ff_rs_state_pml)
    end if

    call profiling_out(prof)

    POP_SUB(cpml_conv_function_update)
  end subroutine cpml_conv_function_update

  ! ---------------------------------------------------------
  subroutine cpml_conv_function_update_via_riemann_silberstein(hm, gr, ff_rs_state_pmlb)
    type(hamiltonian_mxll_t), intent(inout) :: hm
    type(grid_t),             intent(in)    :: gr
    type(batch_t),            intent(inout) :: ff_rs_state_pmlb

    integer :: ip, ip_in, np_part, rs_sign
    CMPLX :: pml_a, pml_b, pml_g, grad
    FLOAT :: g_real, g_imag
    type(profile_t), save :: prof
    integer :: pml_dir, field_dir, ifield, idir
    integer, parameter :: field_dirs(3, 2) = reshape([2, 3, 1, 3, 1, 2], [3, 2])
    logical :: with_medium
    type(batch_t), allocatable :: gradb(:)

    PUSH_SUB(cpml_conv_function_update_via_riemann_silberstein)

    call profiling_in(prof, 'CPML_CONV_FUN_UPDATE_VIA_RS')

    ASSERT(hm%dim == 3 .or. hm%dim == 6)

    np_part = gr%mesh%np_part
    rs_sign = hm%rs_sign

    SAFE_ALLOCATE(gradb(1:gr%der%dim))
    do idir = 1, gr%der%dim
      call ff_rs_state_pmlb%copy_to(gradb(idir))
    end do
    call zderivatives_batch_grad(gr%der, ff_rs_state_pmlb, gradb)

    with_medium = hm%dim == 6

    do pml_dir = 1, hm%st%dim
      select case(gradb(pml_dir)%status())
      case(BATCH_NOT_PACKED)
        do ip_in=1, hm%bc%pml%points_number
          ip       = hm%bc%pml%points_map(ip_in)
          pml_a = hm%bc%pml%a(ip_in,pml_dir)
          pml_b = hm%bc%pml%b(ip_in,pml_dir)
          do ifield = 1, 2
            field_dir = field_dirs(pml_dir, ifield)
            grad = gradb(pml_dir)%zff_linear(ip, field_dir)
            pml_g = hm%bc%pml%conv_plus(ip_in, pml_dir, field_dir)
            g_real = TOFLOAT(pml_a) * TOFLOAT(grad) + TOFLOAT(pml_b) * TOFLOAT(pml_g)
            g_imag = aimag(pml_a) * aimag(grad) + aimag(pml_b) * aimag(pml_g)
            hm%bc%pml%conv_plus(ip_in, pml_dir, field_dir) = TOCMPLX(g_real, g_imag)
            if (with_medium) then
              grad = gradb(pml_dir)%zff_linear(ip, field_dir+3)
              pml_g = hm%bc%pml%conv_minus(ip_in, pml_dir, field_dir)
              g_real = TOFLOAT(pml_a) * TOFLOAT(grad) + TOFLOAT(pml_b) * TOFLOAT(pml_g)
              g_imag = aimag(pml_a) * aimag(grad) + aimag(pml_b) * aimag(pml_g)
              hm%bc%pml%conv_minus(ip_in, pml_dir, field_dir) = TOCMPLX(g_real, g_imag)
            end if
          end do
        end do
      case(BATCH_PACKED)
        do ip_in=1, hm%bc%pml%points_number
          ip       = hm%bc%pml%points_map(ip_in)
          pml_a = hm%bc%pml%a(ip_in,pml_dir)
          pml_b = hm%bc%pml%b(ip_in,pml_dir)
          do ifield = 1, 2
            field_dir = field_dirs(pml_dir, ifield)
            grad = gradb(pml_dir)%zff_pack(field_dir, ip)
            pml_g = hm%bc%pml%conv_plus(ip_in, pml_dir, field_dir)
            g_real = TOFLOAT(pml_a) * TOFLOAT(grad) + TOFLOAT(pml_b) * TOFLOAT(pml_g)
            g_imag = aimag(pml_a) * aimag(grad) + aimag(pml_b) * aimag(pml_g)
            hm%bc%pml%conv_plus(ip_in, pml_dir, field_dir) = TOCMPLX(g_real, g_imag)
            if (with_medium) then
              grad = gradb(pml_dir)%zff_pack(field_dir+3, ip)
              pml_g = hm%bc%pml%conv_minus(ip_in, pml_dir, field_dir)
              g_real = TOFLOAT(pml_a) * TOFLOAT(grad) + TOFLOAT(pml_b) * TOFLOAT(pml_g)
              g_imag = aimag(pml_a) * aimag(grad) + aimag(pml_b) * aimag(pml_g)
              hm%bc%pml%conv_minus(ip_in, pml_dir, field_dir) = TOCMPLX(g_real, g_imag)
            end if
          end do
        end do
      case(BATCH_DEVICE_PACKED)
        call messages_not_implemented("PML on GPU")
      end select
    end do

    do idir = 1, gr%der%dim
      call gradb(idir)%end()
    end do
    SAFE_DEALLOCATE_A(gradb)

    call profiling_out(prof)

    POP_SUB(cpml_conv_function_update_via_riemann_silberstein)
  end subroutine cpml_conv_function_update_via_riemann_silberstein

  ! ---------------------------------------------------------
  subroutine cpml_conv_function_update_via_e_b_fields(hm, gr, ff_rs_state_pml)
    type(hamiltonian_mxll_t), intent(inout) :: hm
    type(grid_t),        intent(in)         :: gr
    CMPLX,               intent(in)         :: ff_rs_state_pml(:,:)

    integer :: ip, ip_in, np_part
    FLOAT, allocatable :: tmp_e(:,:), tmp_b(:,:), tmp_partial_e(:), tmp_partial_b(:)
    CMPLX              :: pml_a(3), pml_b(3), pml_g(3)
    type(profile_t), save :: prof

    PUSH_SUB(cpml_conv_function_update_via_e_b_fields)

    call profiling_in(prof, 'CPML_CONV_FUNC_UPD_VIA_E_B')

    np_part = gr%mesh%np_part
    SAFE_ALLOCATE(tmp_e(np_part,3))
    SAFE_ALLOCATE(tmp_partial_e(np_part))
    SAFE_ALLOCATE(tmp_b(np_part,3))
    SAFE_ALLOCATE(tmp_partial_b(np_part))

    call get_electric_field_state(ff_rs_state_pml, gr%mesh, tmp_e)
    call get_magnetic_field_state(ff_rs_state_pml, gr%mesh, hm%rs_sign, tmp_b)

    ! calculation g(1,2)
    call dderivatives_partial(gr%der, tmp_e(:,2), tmp_partial_e(:), 1, set_bc = .false.)
    call dderivatives_partial(gr%der, tmp_b(:,2), tmp_partial_b(:), 1, set_bc = .false.)
    tmp_partial_e(:) = sqrt(P_ep/M_TWO) * tmp_partial_e(:)
    tmp_partial_b(:) = sqrt(M_ONE/(M_TWO*P_mu)) * tmp_partial_b(:)
    do ip_in=1, hm%bc%pml%points_number
      ip       = hm%bc%pml%points_map(ip_in)
      pml_a(:) = hm%bc%pml%a(ip_in,:)
      pml_b(:) = hm%bc%pml%b(ip_in,:)
      pml_g(:) = hm%bc%pml%conv_plus(ip_in,1,:)
      pml_g(2) = real(pml_a(1)) * tmp_partial_e(ip) + real(pml_b(1)) * real(pml_g(2)) + &
                 M_zI * ( aimag(pml_a(1)) * tmp_partial_b(ip) + aimag(pml_b(1)) * aimag(pml_g(2)) )
      hm%bc%pml%conv_plus(ip_in,1,:) = pml_g(:)
    end do

    ! calculation g(2,1)
    call dderivatives_partial(gr%der, tmp_e(:,1), tmp_partial_e(:), 2, set_bc = .false.)
    call dderivatives_partial(gr%der, tmp_b(:,1), tmp_partial_b(:), 2, set_bc = .false.)
    tmp_partial_e(:) = sqrt(P_ep/M_TWO) * tmp_partial_e(:)
    tmp_partial_b(:) = sqrt(M_ONE/(M_TWO*P_mu)) * tmp_partial_b(:)
    do ip_in=1, hm%bc%pml%points_number
      ip       = hm%bc%pml%points_map(ip_in)
      pml_a(:) = hm%bc%pml%a(ip_in,:)
      pml_b(:) = hm%bc%pml%b(ip_in,:)
      pml_g(:) = hm%bc%pml%conv_plus(ip_in,2,:)
      pml_g(1) = real(pml_a(2)) * tmp_partial_e(ip) + real(pml_b(2)) * real(pml_g(1)) + &
                 M_zI * ( aimag(pml_a(2)) * tmp_partial_b(ip) + aimag(pml_b(2)) * aimag(pml_g(1)) )
      hm%bc%pml%conv_plus(ip_in,2,:) = pml_g(:)
    end do

    ! calculation g(1,3)
    call dderivatives_partial(gr%der, tmp_e(:,3), tmp_partial_e(:), 1, set_bc = .false.)
    call dderivatives_partial(gr%der, tmp_b(:,3), tmp_partial_b(:), 1, set_bc = .false.)
    tmp_partial_e(:) = sqrt(P_ep/M_TWO) * tmp_partial_e(:)
    tmp_partial_b(:) = sqrt(M_ONE/(M_TWO*P_mu)) * tmp_partial_b(:)
    do ip_in=1, hm%bc%pml%points_number
      ip       = hm%bc%pml%points_map(ip_in)
      pml_a(:) = hm%bc%pml%a(ip_in,:)
      pml_b(:) = hm%bc%pml%b(ip_in,:)
      pml_g(:) = hm%bc%pml%conv_plus(ip_in,1,:)
      pml_g(3) = real(pml_a(1)) * tmp_partial_e(ip) + real(pml_b(1)) * real(pml_g(3)) + &
                 M_zI * ( aimag(pml_a(1)) * tmp_partial_b(ip) + aimag(pml_b(1)) * aimag(pml_g(3)) )
      hm%bc%pml%conv_plus(ip_in,1,:) = pml_g(:)
    end do

    ! calculation g(3,1)
    call dderivatives_partial(gr%der, tmp_e(:,1), tmp_partial_e(:), 3, set_bc = .false.)
    call dderivatives_partial(gr%der, tmp_b(:,1), tmp_partial_b(:), 3, set_bc = .false.)
    tmp_partial_e(:) = sqrt(P_ep/M_TWO) * tmp_partial_e(:)
    tmp_partial_b(:) = sqrt(M_ONE/(M_TWO*P_mu)) * tmp_partial_b(:)
    do ip_in=1, hm%bc%pml%points_number
      ip       = hm%bc%pml%points_map(ip_in)
      pml_a(:) = hm%bc%pml%a(ip_in,:)
      pml_b(:) = hm%bc%pml%b(ip_in,:)
      pml_g(:) = hm%bc%pml%conv_plus(ip_in,3,:)
      pml_g(1) = real(pml_a(3)) * tmp_partial_e(ip) + real(pml_b(3)) * real(pml_g(1)) + &
                 M_zI * ( aimag(pml_a(3)) * tmp_partial_b(ip) + aimag(pml_b(3)) * aimag(pml_g(1)) )
      hm%bc%pml%conv_plus(ip_in,3,:) = pml_g(:)
    end do

    ! calculation g(2,3)
    call dderivatives_partial(gr%der, tmp_e(:,3), tmp_partial_e(:), 2, set_bc = .false.)
    call dderivatives_partial(gr%der, tmp_b(:,3), tmp_partial_b(:), 2, set_bc = .false.)
    tmp_partial_e(:) = sqrt(P_ep/M_TWO) * tmp_partial_e(:)
    tmp_partial_b(:) = sqrt(M_ONE/(M_TWO*P_mu)) * tmp_partial_b(:)
    do ip_in=1, hm%bc%pml%points_number
      ip       = hm%bc%pml%points_map(ip_in)
      pml_a(:) = hm%bc%pml%a(ip_in,:)
      pml_b(:) = hm%bc%pml%b(ip_in,:)
      pml_g(:) = hm%bc%pml%conv_plus(ip_in,2,:)
      pml_g(3) = real(pml_a(2)) * tmp_partial_e(ip) + real(pml_b(2)) * real(pml_g(3)) + &
                 M_zI * ( aimag(pml_a(2)) * tmp_partial_b(ip) + aimag(pml_b(2)) * aimag(pml_g(3)) )
      hm%bc%pml%conv_plus(ip_in,2,:) = pml_g(:)
    end do

    ! calculation g(3,2)
    call dderivatives_partial(gr%der, tmp_e(:,2), tmp_partial_e(:), 3, set_bc = .false.)
    call dderivatives_partial(gr%der, tmp_b(:,2), tmp_partial_b(:), 3, set_bc = .false.)
    tmp_partial_e(:) = sqrt(P_ep/M_TWO) * tmp_partial_e(:)
    tmp_partial_b(:) = sqrt(M_ONE/(M_TWO*P_mu)) * tmp_partial_b(:)
    do ip_in=1, hm%bc%pml%points_number
      ip       = hm%bc%pml%points_map(ip_in)
      pml_a(:) = hm%bc%pml%a(ip_in,:)
      pml_b(:) = hm%bc%pml%b(ip_in,:)
      pml_g(:) = hm%bc%pml%conv_plus(ip_in,3,:)
      pml_g(2) = real(pml_a(3)) * tmp_partial_e(ip) + real(pml_b(3)) * real(pml_g(2)) + &
                 M_zI * ( aimag(pml_a(3)) * tmp_partial_b(ip) + aimag(pml_b(3)) * aimag(pml_g(2)) )
      hm%bc%pml%conv_plus(ip_in,3,:) = pml_g(:)
    end do

    SAFE_DEALLOCATE_A(tmp_e)
    SAFE_DEALLOCATE_A(tmp_partial_e)
    SAFE_DEALLOCATE_A(tmp_b)
    SAFE_DEALLOCATE_A(tmp_partial_b)

    call profiling_out(prof)

    POP_SUB(cpml_conv_function_update_via_e_b_fields)
  end subroutine cpml_conv_function_update_via_e_b_fields

  ! ---------------------------------------------------------
  subroutine td_function_mxll_init(st, namespace, hm)
    type(states_mxll_t),      intent(inout) :: st
    type(namespace_t),        intent(in)    :: namespace
    type(hamiltonian_mxll_t), intent(inout)    :: hm

    type(block_t)        :: blk
    integer              :: il, nlines, idim, ncols, ierr
    FLOAT                :: e_field(st%dim), b_field(st%dim)
    character(len=1024)  :: mxf_expression
    type(profile_t), save :: prof

    PUSH_SUB(td_function_mxll_init)

    call profiling_in(prof, 'TD_FUNCTION_MXLL_INIT')

    !%Variable UserDefinedConstantSpatialMaxwellField
    !%Type block
    !%Section MaxwellStates
    !%Description
    !% Define parameters of spatially constant field.
    !%
    !% Example:
    !%
    !% <tt>%UserDefinedConstantSpatialMaxwellFields
    !% <br>&nbsp;&nbsp;   plane_wave_parser      | E_x | E_y | E_z | B_x | B_y | B_z | "tdf_function"
    !% <br>%</tt>
    !%
    !% This block defines three components of E field, three components of B field, and reference to
    !% the TD function.
    !%
    !%End

    if (parse_block(namespace, 'UserDefinedConstantSpatialMaxwellField', blk) == 0) then
      st%rs_state_const_external = .true.
      nlines = parse_block_n(blk)
      SAFE_ALLOCATE(st%rs_state_const_td_function(nlines))
      SAFE_ALLOCATE(st%rs_state_const_amp(st%dim, nlines))
      ! read all lines
      do il = 1, nlines
        e_field = M_ZERO
        b_field = M_ZERO
        ! Check that number of columns is five or six.
        ncols = parse_block_cols(blk, il - 1)
        if (ncols  /= 7) then
          message(1) = 'Each line in the UserDefinedConstantSpatialMaxwellField block must have'
          message(2) = 'seven columns.'
          call messages_fatal(2)
        end if
        do idim = 1, st%dim
          call parse_block_float( blk, il - 1, idim-1, e_field(idim))
        end do
        do idim = 1, st%dim
          call parse_block_float( blk, il - 1, idim+2, b_field(idim))
        end do
        call parse_block_string( blk, il - 1, 6, mxf_expression)
        call build_rs_vector(e_field, b_field, st%rs_sign, st%rs_state_const_amp(:,il))
        call tdf_read(st%rs_state_const_td_function(il), namespace, trim(mxf_expression), ierr)
      end do
    end if
    call parse_block_end(blk)

    !%Variable PropagateSpatialMaxwellField
    !%Type logical
    !%Default yes
    !%Section MaxwellStates
    !%Description
    !% Allow for numerical propagation of Maxwells equations of spatially constant field.
    !% If set to no, do only analytic evaluation of the field inside the box.
    !%End

    call parse_variable(namespace, 'PropagateSpatialMaxwellField', .true., hm%spatial_constant_propagate)

    call profiling_out(prof)

    POP_SUB(td_function_mxll_init)
  end subroutine td_function_mxll_init

  ! ---------------------------------------------------------
  subroutine spatial_constant_calculation(constant_calc, st, gr, hm, time, dt, delay, rs_state, set_initial_state)
    logical,                  intent(in)    :: constant_calc
    type(states_mxll_t),      intent(inout) :: st
    type(grid_t),             intent(in)    :: gr
    type(hamiltonian_mxll_t), intent(in)    :: hm
    FLOAT,                    intent(in)    :: time
    FLOAT,                    intent(in)    :: dt
    FLOAT,                    intent(in)    :: delay
    CMPLX,                    intent(inout) :: rs_state(:,:)
    logical,        optional, intent(in)    :: set_initial_state

    integer :: ip, ic, icn
    FLOAT   :: tf_old, tf_new
    logical :: set_initial_state_
    type(profile_t), save :: prof

    PUSH_SUB(spatial_constant_calculation)

    call profiling_in(prof, 'SPATIAL_CONSTANT_CALCULATION')

    set_initial_state_ = .false.
    if (present(set_initial_state)) set_initial_state_ = set_initial_state

    if (hm%spatial_constant_apply) then
      if (constant_calc) then
        icn = size(st%rs_state_const_td_function(:))
        st%rs_state_const(:) = M_z0
        do ic = 1, icn
          tf_old = tdf(st%rs_state_const_td_function(ic), time-delay-dt)
          tf_new = tdf(st%rs_state_const_td_function(ic), time-delay)
          do ip = 1, gr%mesh%np
            if (set_initial_state_ .or. (.not. hm%spatial_constant_propagate)) then
              rs_state(ip,:) = st%rs_state_const_amp(:,ic) * tf_new
            else
              rs_state(ip,:) = rs_state(ip,:) + st%rs_state_const_amp(:,ic) * (tf_new - tf_old)
            end if
          end do
          st%rs_state_const(:) = st%rs_state_const(:) + st%rs_state_const_amp(:, ic)
        end do
       st%rs_state_const(:) = st%rs_state_const(:) * tf_new
      end if
    end if

    call profiling_out(prof)

    POP_SUB(spatial_constant_calculation)
  end subroutine spatial_constant_calculation

  ! ---------------------------------------------------------
  subroutine constant_boundaries_calculation(constant_calc, bc, hm, st, rs_state)
    logical,                   intent(in)    :: constant_calc
    type(bc_mxll_t),           intent(inout) :: bc
    type(states_mxll_t),       intent(in)    :: st
    type(hamiltonian_mxll_t),  intent(in)    :: hm
    CMPLX,                     intent(inout) :: rs_state(:,:)

    integer :: ip_in, ip
    type(profile_t), save :: prof

    PUSH_SUB(constant_boundaries_calculation)
    call profiling_in(prof, 'CONSTANT_BOUNDARIES_CALC')

    if (hm%spatial_constant_apply) then
      if (constant_calc) then
        do ip_in = 1, bc%constant_points_number
          ip = bc%constant_points_map(ip_in)
          rs_state(ip,:) = st%rs_state_const(:)
          bc%constant_rs_state(ip_in,:) = st%rs_state_const(:)
        end do
      end if
    end if

    call profiling_out(prof)

    POP_SUB(constant_boundaries_calculation)
  end subroutine constant_boundaries_calculation

  ! ---------------------------------------------------------
  subroutine mirror_pec_boundaries_calculation(bc, st, rs_state)
    type(bc_mxll_t),     intent(in)    :: bc
    type(states_mxll_t), intent(in)    :: st
    CMPLX,               intent(inout) :: rs_state(:,:)

    integer                    :: ip, ip_in, idim
    FLOAT                      :: e_field(st%dim), b_field(st%dim)

    PUSH_SUB(mirror_pec_boundaries_calculation)

    do idim = 1, 3
      if (bc%bc_type(idim) == MXLL_BC_MIRROR_PEC) then
        do ip_in = 1, bc%mirror_points_number(idim)
          ip = bc%mirror_points_map(ip_in, idim)
          e_field(:) = M_ZERO
          call get_magnetic_field_vector(rs_state(ip,:), st%rs_sign, b_field(:), st%mu(ip))
          call build_rs_vector(e_field(:), b_field(:), st%rs_sign, rs_state(ip,:), st%ep(ip), st%mu(ip))
        end do
      end if
    end do

    POP_SUB(mirror_pec_boundaries_calculation)
  end subroutine mirror_pec_boundaries_calculation

  ! ---------------------------------------------------------
  subroutine mirror_pmc_boundaries_calculation(bc, st, rs_state)
    type(bc_mxll_t),     intent(in)    :: bc
    type(states_mxll_t), intent(in)    :: st
    CMPLX,               intent(inout) :: rs_state(:,:)

    integer                    :: ip, ip_in, idim
    FLOAT                      :: e_field(st%dim), b_field(st%dim)

    PUSH_SUB(mirror_pmc_boundaries_calculation)

    do idim = 1, 3
      if (bc%bc_type(idim) == MXLL_BC_MIRROR_PMC) then
        do ip_in = 1, bc%mirror_points_number(idim)
          ip = bc%mirror_points_map(ip_in,idim)
          b_field(:) = M_ZERO
          call get_electric_field_vector(rs_state(ip,:), e_field(:), st%ep(ip))
          call build_rs_vector(e_field(:), b_field(:), st%rs_sign, rs_state(ip,:), st%ep(ip), st%mu(ip))
        end do
      end if
    end do

    POP_SUB(mirror_pmc_boundaries_calculation)
  end subroutine mirror_pmc_boundaries_calculation

  ! ---------------------------------------------------------
  subroutine plane_waves_boundaries_calculation(hm, st, mesh, time, time_delay, rs_state)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    type(states_mxll_t),      intent(in)    :: st
    type(mesh_t),             intent(in)    :: mesh
    FLOAT,                    intent(in )   :: time
    FLOAT,                    intent(in)    :: time_delay
    CMPLX,                    intent(inout) :: rs_state(:,:)

    integer                    :: ip, ip_in, wn
    FLOAT                      :: x_prop(mesh%sb%dim), rr, vv(mesh%sb%dim), k_vector(mesh%sb%dim)
    FLOAT                      :: k_vector_abs, nn
    CMPLX                      :: e0(mesh%sb%dim)
    FLOAT                      :: e_field(mesh%sb%dim), b_field(mesh%sb%dim)
    CMPLX                      :: rs_state_add(mesh%sb%dim)
    type(profile_t), save :: prof

    PUSH_SUB(plane_waves_boundaries_calculation)

    call profiling_in(prof, 'PLANE_WAVES_BOUNDARIES_CALC')

    if (hm%plane_waves_apply) then
      do wn = 1, hm%bc%plane_wave%number
         k_vector(:) = hm%bc%plane_wave%k_vector(1:mesh%sb%dim, wn)
         k_vector_abs = sqrt(sum(k_vector(1:mesh%sb%dim)**2))
         vv(:) = hm%bc%plane_wave%v_vector(1:mesh%sb%dim, wn)
         e0(:) = hm%bc%plane_wave%e_field(1:mesh%sb%dim, wn)
         do ip_in = 1, hm%bc%plane_wave%points_number
          ip = hm%bc%plane_wave%points_map(ip_in)
          if (wn == 1) rs_state(ip,:) = M_Z0
          nn = sqrt(st%ep(ip)/P_ep*st%mu(ip)/P_mu)
          x_prop(1:mesh%sb%dim) = mesh%x(ip,1:mesh%sb%dim) - vv(1:mesh%sb%dim) * (time - time_delay)
          rr = sqrt(sum(x_prop(1:mesh%sb%dim)**2))
          if (hm%bc%plane_wave%modus(wn) == OPTION__MAXWELLINCIDENTWAVES__PLANE_WAVE_MX_FUNCTION) then
            e_field(1:mesh%sb%dim) = TOFLOAT(e0(1:mesh%sb%dim) * mxf(hm%bc%plane_wave%mx_function(wn), x_prop(1:mesh%sb%dim)))
          end if
          b_field(1:3) = dcross_product(k_vector, e_field) / P_c / k_vector_abs
          call build_rs_vector(e_field, b_field, st%rs_sign, rs_state_add, st%ep(ip), st%mu(ip))
          rs_state(ip, :) =  rs_state(ip, :) + rs_state_add(:)
        end do
      end do
      else
        do ip_in = 1, hm%bc%plane_wave%points_number
          ip             = hm%bc%plane_wave%points_map(ip_in)
          rs_state(ip,:) = M_z0
        end do
      end if

      call profiling_out(prof)

    POP_SUB(plane_waves_boundaries_calculation)
  end subroutine plane_waves_boundaries_calculation

  ! ---------------------------------------------------------
  subroutine plane_waves_propagation(hm, tr, namespace, st, gr, time, dt, time_delay)
    type(hamiltonian_mxll_t), intent(inout) :: hm
    type(propagator_mxll_t),  intent(inout) :: tr
    type(namespace_t),        intent(in)    :: namespace
    type(states_mxll_t),      intent(inout) :: st
    type(grid_t),             intent(in)    :: gr
    FLOAT,                    intent(in)    :: time
    FLOAT,                    intent(in)    :: dt
    FLOAT,                    intent(in)    :: time_delay

    type(batch_t) :: ff_rs_stateb
    integer :: ff_dim
    type(profile_t), save :: prof

    PUSH_SUB(plane_waves_propagation)

    call profiling_in(prof, 'PLANE_WAVES_PROPAGATION')

    ff_dim = hm%dim
    call zbatch_init(ff_rs_stateb, 1, 1, hm%dim, gr%mesh%np_part)

    call transform_rs_state_batch(hm, gr, st, st%rs_state_plane_waves, ff_rs_stateb, RS_TRANS_FORWARD)

    ! Time evolution of RS plane waves state without any coupling with H(inter_time)
    call hamiltonian_mxll_update(hm, time=time)
    hm%cpml_hamiltonian = .false.
    call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, ff_rs_stateb, dt)

    call transform_rs_state_batch(hm, gr, st, st%rs_state_plane_waves, ff_rs_stateb, RS_TRANS_BACKWARD)
    call plane_waves_boundaries_calculation(hm, st, gr%mesh, time+dt, time_delay, st%rs_state_plane_waves)

    call profiling_out(prof)
    POP_SUB(plane_waves_propagation)
  end subroutine plane_waves_propagation

  ! ---------------------------------------------------------
  subroutine plane_waves_in_box_calculation(bc, time, gr, st, hm, rs_state)
    type(bc_mxll_t),           intent(in)    :: bc
    FLOAT,                     intent(in)    :: time
    type(grid_t),              intent(in)    :: gr
    type(states_mxll_t),       intent(in)    :: st
    type(hamiltonian_mxll_t),  intent(in)    :: hm
    CMPLX,                     intent(inout) :: rs_state(:,:)

    integer              :: ip, wn, idim, np
    FLOAT                :: x_prop(gr%sb%dim), rr, vv(gr%sb%dim), k_vector(gr%sb%dim), k_vector_abs, nn
    FLOAT                :: e_field(gr%sb%dim), b_field(gr%sb%dim), dummy(gr%sb%dim)
    CMPLX                :: rs_state_add(st%dim), e0(gr%sb%dim)
    type(profile_t), save :: prof

    PUSH_SUB(plane_waves_in_box_calculation)

    call profiling_in(prof, 'PLANE_WAVES_IN_BOX_CALCULATION')

    np = gr%mesh%np_part
    do wn = 1, bc%plane_wave%number
      vv(:) = hm%bc%plane_wave%v_vector(1:gr%sb%dim, wn)
      k_vector(:) = hm%bc%plane_wave%k_vector(1:gr%sb%dim, wn)
      k_vector_abs = sqrt(sum(k_vector(1:gr%sb%dim)**2))
      e0(:) = hm%bc%plane_wave%e_field(1:gr%sb%dim, wn)
      do ip = 1, gr%mesh%np
        if (wn == 1) rs_state(ip,:) = M_Z0
        nn = sqrt(st%ep(ip)/P_ep*st%mu(ip)/P_mu)
        x_prop = gr%mesh%x(ip,:) - vv * time
        rr = sqrt(sum(x_prop(1:gr%sb%dim)**2))
        select case (bc%plane_wave%modus(wn))
        case (OPTION__MAXWELLINCIDENTWAVES__PLANE_WAVE_PARSER)
          do idim = 1, gr%sb%dim
            call parse_expression(e_field(idim), dummy(idim), gr%sb%dim, x_prop, rr, M_ZERO, &
              bc%plane_wave%e_field_string(idim,wn))
            e_field(idim) = units_to_atomic(units_inp%energy/units_inp%length, e_field(idim))
          end do
        case (OPTION__MAXWELLINCIDENTWAVES__PLANE_WAVE_MX_FUNCTION)
          e_field(:) = TOFLOAT(e0(1:gr%sb%dim) * mxf(bc%plane_wave%mx_function(wn), x_prop(1:gr%sb%dim)))
        end select
        b_field(1:3) = M_ONE/(P_c * k_vector_abs) * dcross_product(k_vector, e_field)
        call build_rs_vector(e_field, b_field, st%rs_sign, rs_state_add, st%ep(ip), st%mu(ip))
        rs_state(ip,:) = rs_state(ip,:) + rs_state_add(:)
      end do
    end do

    call profiling_out(prof)

    POP_SUB(plane_waves_in_box_calculation)
  end subroutine plane_waves_in_box_calculation

end module propagator_mxll_oct_m

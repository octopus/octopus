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
!!
!! $Id: propagator.F90 13908 2015-05-05 06:02:30Z xavier $

#include "global.h"

module propagator_mxll_oct_m
  use boundary_op_oct_m
  use batch_ops_oct_m
  use batch_oct_m
  use comm_oct_m
  use cube_function_oct_m
  use cube_oct_m
  use current_oct_m
  use derivatives_oct_m
  use density_oct_m
  use energy_calc_oct_m
  use exponential_oct_m
  use fft_oct_m
  use fourier_space_oct_m
  use grid_oct_m
  use geometry_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_mxll_oct_m
  use index_oct_m
  use io_oct_m
  use io_function_oct_m
  use io_oct_m
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
  use propagator_base_oct_m
  use poisson_oct_m
  use states_elec_oct_m
  use states_mxll_oct_m
  use string_oct_m
  use tdfunction_oct_m
  use varinfo_oct_m
  use xyz_adjust_oct_m

  implicit none

  private
  public ::                                          &
    propagator_mxll_init,                         &
    propagation_mxll_etrs,                        &
    transform_rs_state_forward,              &
    transform_rs_state_backward,             &
    transform_rs_densities_forward,          &
    transform_rs_densities_backward,         &
    rs_state_to_cube_map,                    &
    calculate_matter_longitudinal_field,     &
    get_vector_pot_and_transverse_field,     &
    inner_and_outer_points_mapping,          &
    surface_grid_points_mapping,             &
    energy_density_calc,                     &
    energy_mxll_calc,                             &
    energy_rate_calc,                        &
    poynting_vector_through_box_surfaces,    &
    poynting_vector_through_box_surfaces_plane_waves, &
    fields_through_box_surfaces,             &
    fields_through_box_surfaces_plane_waves, &
    constant_boundaries_calculation,         &
    maxwell_mask,                            &
    td_function_mxll_init,                   &
    plane_waves_in_box_calculation
    ! maxwell_matter_mesh_mapping,                     &
    ! maxwell_grid_points_coupling_points_mapping,     &
    ! get_mx_ma_coupling_points
    
contains

  ! ---------------------------------------------------------
  subroutine propagator_mxll_init(gr, namespace, st, hm, tr)
    type(grid_t),                 intent(in)    :: gr
    type(namespace_t),   intent(in)    :: namespace
    type(states_mxll_t),               intent(inout) :: st
    type(hamiltonian_mxll_t),          intent(inout) :: hm
    type(propagator_mxll_t),   intent(inout) :: tr

    integer :: default_propagator, il, nlines, ncols, icol, idim
    type(block_t) :: blk
    character(len=256) :: string
    logical :: plane_waves_set = .false.

    PUSH_SUB(propagator_mxll_init)
   
    !%Variable MaxwellTDPropagator
    !%Type integer
    !%Default etrs
    !%Section Time-Dependent::Propagation
    !%Description
    !% There are several time-evolution methods for the Maxwell propagation
    !% similar to the methods for the time-evolution of Kohn-Sham orbitals.
    !%Option maxwell_etrs 1
    !% Enforced time-reversal-symmetry propagation (etrs)
    !%End
    default_propagator = OPTION__MAXWELLTDPROPAGATOR__MAXWELL_ETRS
    call parse_variable(namespace, 'MaxwellTDPropagator', default_propagator, tr%tr_method)
    call messages_print_var_option(stdout, 'MaxwellTDPropagator', tr%tr_method)

    !%Variable MaxwellBoundaryConditions
    !%Type block
    !%Section Time-Dependent::Propagation
    !%Description
    !% Follows
    !%
    !% Example:
    !%
    !% <tt>%MaxwellBoundaryConditions
    !% <br>&nbsp;&nbsp;   zero | mirror_pec | consant 
    !% <br>%</tt>
    !%
    !% Description follows
    !%
    !%Option zero 0
    !% follows ...
    !%Option constant 2
    !% follows ...
    !%Option mirror_pec 3
    !% follows ...
    !%Option mirror_pmc 4
    !% follows ...
    !%Option plane_waves 5
    !% follows ...
    !%Option periodic 6
    !% follows ...
    !%Option medium 7
    !% follows ...
    !%Option lossy_layer 8
    !% follows ...
    !%End
    if(parse_block(namespace, 'MaxwellBoundaryConditions', blk) == 0) then

      call messages_print_stress(stdout, trim('Maxwell boundary conditions:'))

      ! find out how many lines (i.e. states) the block has
      nlines = parse_block_n(blk)
      if (nlines /= 1) then
        message(1) = 'MaxwellBoundaryConditions has to consist of one line!'
        call messages_fatal(1)
      end if
      ncols = parse_block_cols(blk, 0)
      if (ncols /= 3) then
        message(1) = 'MaxwellBoundaryConditions has to consist of three columns!'
        call messages_fatal(1)
      end if
      do icol=1, ncols
        call parse_block_integer(blk, 0, icol-1, hm%bc%bc_type(icol))
        if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__ZERO) then
          string = 'Zero'
          tr%bc_zero = .true.
        else if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__CONSTANT) then
          string = 'Constant'
          tr%bc_constant = .true.
          tr%bc_add_ab_region = .true.
          hm%bc_constant = .true.
          hm%bc_add_ab_region = .true.
          SAFE_ALLOCATE(st%rs_state_const(1:st%d%dim))
          st%rs_state_const = M_z0
          call td_function_mxll_init(st, namespace, gr, hm)
        else if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__MIRROR_PEC) then
          string = 'PEC Mirror'
          tr%bc_mirror_pec = .true.
          hm%bc_mirror_pec = .true.
        else if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__MIRROR_PMC) then
          string = 'PMC Mirror'
          tr%bc_mirror_pmc = .true.
          hm%bc_mirror_pmc = .true.
        else if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__PERIODIC) then
          string = 'Periodic'
          tr%bc_periodic = .true.
          hm%bc_periodic = .true.
        else if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__PLANE_WAVES) then
          string = 'Plane waves'
          plane_waves_set = .true.
          tr%bc_plane_waves = .true.
          tr%bc_add_ab_region = .true.
          hm%plane_waves = .true.
          hm%bc_plane_waves = .true.
          hm%bc_add_ab_region = .true.
          SAFE_ALLOCATE(st%rs_state_plane_waves(1:gr%mesh%np_part, 1:st%d%dim))
        else if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__MEDIUM) then
          string = 'Medium boundary'
        end if
        write(message(1),'(a,I1,a,a)') 'Maxwell boundary condition in direction ', icol, ': ', trim(string)
        call messages_info(1)
        if (plane_waves_set .and. .not. (parse_is_defined(namespace, 'UserDefinedMaxwellIncidentWaves')) ) then
          write(message(1),'(a)') 'Input: Maxwell boundary condition option is set to "plane_waves".'
          write(message(2),'(a)') 'Input: User defined Maxwell plane waves have to be defined!'
          call messages_fatal(2)
        end if
      end do

      call messages_print_stress(stdout)

    end if

    !%Variable MaxwellMediumBox
    !%Type block
    !%Section Time-Dependent::Propagation
    !%Description
    !% Follows
    !%
    !% Example:
    !%
    !% <tt>%MaxwellMediumBox
    !% <br>&nbsp;&nbsp;   center_x | center_y | center_z | x_length | y_length | z_length | epsilon_factor | mu_factor | sigma_e | sigma_m | edged/smooth
    !% <br>%</tt>
    !%
    !% Description about MaxwellMediumBox follows
    !%
    !%Option edged 1
    !% Follows
    !%Option smooth 2
    !% Follows
    !%End
    if(parse_block(namespace, 'MaxwellMediumBox', blk) == 0) then

      call messages_print_stress(stdout, trim('Maxwell Medium box:'))
      hm%medium_box = .true.

      ! find out how many lines (i.e. states) the block has
      nlines = parse_block_n(blk)
      SAFE_ALLOCATE(hm%medium_box_center(1:3,1:nlines))
      SAFE_ALLOCATE(hm%medium_box_size(1:3,1:nlines))
      SAFE_ALLOCATE(hm%medium_box_ep_factor(1:nlines))
      SAFE_ALLOCATE(hm%medium_box_mu_factor(1:nlines))
      SAFE_ALLOCATE(hm%medium_box_sigma_e_factor(1:nlines))
      SAFE_ALLOCATE(hm%medium_box_sigma_m_factor(1:nlines))
      SAFE_ALLOCATE(hm%medium_box_shape(1:nlines))
      do il=1, nlines
        ncols = parse_block_cols(blk, il-1)
        if (ncols /= 11) then
          message(1) = 'MaxwellMedium has to consist of eleven columns!'
          call messages_fatal(1)
        end if
        do idim=1,3
          call parse_block_float(blk, il-1, idim-1, hm%medium_box_center(idim,il))
          call parse_block_float(blk, il-1, idim+2, hm%medium_box_size(idim,il))
        end do
        call parse_block_float(blk, il-1, 6, hm%medium_box_ep_factor(il))
        call parse_block_float(blk, il-1, 7, hm%medium_box_mu_factor(il))
        call parse_block_float(blk, il-1, 8, hm%medium_box_sigma_e_factor(il))
        call parse_block_float(blk, il-1, 9, hm%medium_box_sigma_m_factor(il))
        call parse_block_integer(blk, il-1, 10, hm%medium_box_shape(il))
        if (il > 1) then
          write(message(1),'(a)') ""
          write(message(2),'(a,I1)')    'Medium box number:  ', il
          write(message(3),'(a,es9.2,a,es9.2,a,es9.2)') 'Box center:         ', hm%medium_box_center(1,il), ' | ',&
                hm%medium_box_center(2,il), ' | ', hm%medium_box_center(3,il)
          write(message(4),'(a,es9.2,a,es9.2,a,es9.2)') 'Box size  :         ', hm%medium_box_size(1,il), ' | ', &
                hm%medium_box_size(2,il), ' | ', hm%medium_box_size(3,il)
          write(message(5),'(a,es9.2)') 'Box epsilon factor: ', hm%medium_box_ep_factor(il)
          write(message(6),'(a,es9.2)') 'Box mu factor:      ', hm%medium_box_mu_factor(il)
          write(message(7),'(a,es9.2)') 'Box electric sigma: ', hm%medium_box_sigma_e_factor(il)
          write(message(8),'(a,es9.2)') 'Box magnetic sigma: ', hm%medium_box_sigma_m_factor(il)
          if (hm%medium_box_shape(il) == OPTION__MAXWELLMEDIUMBOX__EDGED) then
            write(message(9),'(a,a)')   'Box shape:          ', 'edged'
          else if (hm%medium_box_shape(il) == OPTION__MAXWELLMEDIUMBOX__SMOOTH) then
            write(message(9),'(a,a)')   'Box shape:          ', 'smooth'
          end if
          call messages_info(9)
        else
          write(message(1),'(a,I1)')    'Medium box number:  ', il
          write(message(2),'(a,es9.2,a,es9.2,a,es9.2)') 'Box center:         ', hm%medium_box_center(1,il), ' | ',&
                hm%medium_box_center(2,il), ' | ', hm%medium_box_center(3,il)
          write(message(3),'(a,es9.2,a,es9.2,a,es9.2)') 'Box size  :         ', hm%medium_box_size(1,il), ' | ', &
                hm%medium_box_size(2,il), ' | ', hm%medium_box_size(3,il)
          write(message(4),'(a,es9.2)') 'Box epsilon factor: ', hm%medium_box_ep_factor(il)
          write(message(5),'(a,es9.2)') 'Box mu factor:      ', hm%medium_box_mu_factor(il)
          write(message(6),'(a,es9.2)') 'Box electric sigma: ', hm%medium_box_sigma_e_factor(il)
          write(message(7),'(a,es9.2)') 'Box magnetic sigma: ', hm%medium_box_sigma_m_factor(il)
          if (hm%medium_box_shape(il) == OPTION__MAXWELLMEDIUMBOX__EDGED) then
            write(message(8),'(a,a)')   'Box shape:          ', 'edged'
          else if (hm%medium_box_shape(il) == OPTION__MAXWELLMEDIUMBOX__SMOOTH) then
            write(message(8),'(a,a)')   'Box shape:          ', 'smooth'
          end if
          call messages_info(8)
        end if
      end do

      call generate_medium_boxes(hm, gr, nlines)

      call messages_print_stress(stdout)

    end if

    !%Variable MaxwellTDETRSApprox
    !%Type integer
    !%Default no
    !%Section Time-Dependent::Propagation
    !%Description
    !% follows
    !%Option no 0
    !% follows
    !%Option no_etrs 1
    !% follows
    !%Option const_steps 2
    !% follows
    !%End
    call parse_variable(namespace, 'MaxwellTDETRSApprox', OPTION__MAXWELLTDETRSAPPROX__NO, tr%tr_etrs_approx)
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
    !% not doing any numerical propagation of Maxwell's equations.
    !%End
    call parse_variable(namespace, 'MaxwellPlaneWavesInBox', .false., tr%plane_waves_in_box)
    call set_medium_rs_state(st, gr, hm)

    call derivatives_boundary_mask(hm%bc, gr%mesh, hm)

    !tr%te%exp = .true.
    call exponential_init(tr%te, namespace) ! initialize Maxwell propagator

    POP_SUB(propagator_mxll_init)
  end subroutine propagator_mxll_init


  ! ---------------------------------------------------------
  subroutine propagation_mxll_etrs(hm, namespace, gr, st, tr, rs_state, &
                                      rs_current_density_t1, rs_current_density_t2,             & 
                                      rs_charge_density_t1, rs_charge_density_t2, time, dt,     &
                                      rs_state_pml_predict)
    type(hamiltonian_mxll_t),        intent(inout) :: hm
    type(namespace_t),           intent(in)    :: namespace
    type(grid_t),               intent(inout) :: gr
    type(states_mxll_t),             intent(inout) :: st
    type(propagator_mxll_t), intent(inout) :: tr
    CMPLX,                      intent(inout) :: rs_state(:,:)
    CMPLX,                      intent(inout) :: rs_current_density_t1(:,:)
    CMPLX,                      intent(inout) :: rs_current_density_t2(:,:)
    CMPLX,                      intent(inout) :: rs_charge_density_t1(:)
    CMPLX,                      intent(inout) :: rs_charge_density_t2(:)
    FLOAT,                      intent(in)    :: time
    FLOAT,                      intent(in)    :: dt
    CMPLX,                      intent(in)    :: rs_state_pml_predict(:,:)

    integer            :: ii, iborder, inter_steps, ff_dim, ip, ip_in, idim, istate
    FLOAT              :: inter_dt, inter_time, delay
    CMPLX, allocatable :: ff_rs_state(:,:), ff_rs_inhom_1(:,:), ff_rs_inhom_2(:,:), ff_rs_inhom_mean(:,:)
    CMPLX, allocatable :: ff_rs_state_pml(:,:) !, ff_rs_state_pml_old(:,:), ff_rs_state_pml_predict(:,:), ff_rs_state_pml_slope(:,:)
    logical            :: pml_check = .false.

    PUSH_SUB(propagation_mxll_etrs)


            
    if (tr%plane_waves_in_box) then
      call plane_waves_in_box_calculation(hm%bc, time+dt, gr, st, hm, rs_state)
      POP_SUB(propagation_mxll_etrs)
      return
    end if

    do idim=1, 3
      if (hm%bc%bc_ab_type(idim) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML) then
        pml_check = .true.
      end if
    end do

    if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then
      ff_dim = 6
    else if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_GAUSS) then
      ff_dim = 4
    else 
      ff_dim = 3
    end if

    ! intermediate step variables
    inter_steps   = tr%inter_steps
    inter_dt      = M_ONE / inter_steps * dt

    ! delay time
    delay = tr%delay_time

    SAFE_ALLOCATE(ff_rs_state(1:gr%mesh%np_part,ff_dim))

    if (pml_check) then
      SAFE_ALLOCATE(ff_rs_state_pml(1:gr%mesh%np_part,ff_dim))
    end if

    ! first step of Maxwell inhomogeneity propagation with constant current density
    if ((hm%ma_mx_coupling_apply .or. hm%current_density_ext_flag) .and. &
      (tr%tr_etrs_approx == OPTION__MAXWELLTDETRSAPPROX__CONST_STEPS)) then

      message(1) = "Maxwell-matter coupling or external current not implemented yet"
      call messages_fatal(1)
      
!       !AFE_ALLOCATE(ff_rs_inhom_1(1:gr%mesh%np_part,ff_dim))
!       !AFE_ALLOCATE(ff_rs_inhom_2(1:gr%mesh%np_part,ff_dim))
!       !AFE_ALLOCATE(ff_rs_inhom_mean(1:gr%mesh%np_part,ff_dim))
!       ! inhomogeneity propagation
! 4      call transform_rs_densities_forward(hm, rs_charge_density_t1, rs_current_density_t1, ff_rs_inhom_1)
!       call transform_rs_densities_forward(hm, rs_charge_density_t2, rs_current_density_t2, ff_rs_inhom_2)
!       ff_rs_inhom_mean = (ff_rs_inhom_1 + ff_rs_inhom_2)/M_TWO
!       ! add term J(time)
!       ff_rs_inhom_1 = ff_rs_inhom_mean
!       ff_rs_inhom_2 = ff_rs_inhom_mean
!       call hamiltonian_mxll_update(hm, time=time)
!       call exponential_mxll_apply(hm, gr, st, tr, time, inter_dt, ff_rs_inhom_2)
!       ! add term U(time+dt,time)J(time)
!       ff_rs_inhom_1 = ff_rs_inhom_1 + ff_rs_inhom_2
!       ff_rs_inhom_2 = ff_rs_inhom_mean
!       call hamiltonian_mxll_update(hm, time=time)
!       call exponential_mxll_apply(hm, gr, st, tr, time, inter_dt/M_TWO, ff_rs_inhom_2)
!       ! add term U(time+dt/2,time)J(time)
!       ff_rs_inhom_1 = ff_rs_inhom_1 + ff_rs_inhom_2
!       ff_rs_inhom_2 = ff_rs_inhom_mean
!       call hamiltonian_mxll_update(hm, time=time)
!       call exponential_mxll_apply(hm, gr, st, tr, time, -inter_dt/M_TWO, ff_rs_inhom_2)
!       ! add term U(time,time+dt/2)J(time)
!       ff_rs_inhom_1 = ff_rs_inhom_1 + ff_rs_inhom_2
!       !AFE_DEALLOCATE_A(ff_rs_inhom_2)
!       !AFE_DEALLOCATE_A(ff_rs_inhom_mean)
    end if

    do ii=1, inter_steps

      ! intermediate time
      inter_time = time + inter_dt * (ii-1)

      ! transformation of RS state into 3x3 or 4x4 representation
      ! call transform_rs_state_forward(hm, gr, st, rs_state, ff_rs_state)

      if ((hm%ma_mx_coupling_apply) .or. hm%current_density_ext_flag) then

        message(1) = "Maxwell-matter coupling or external current not implemented yet"
        call messages_fatal(1)

      ! ! Maxwell etrs propagation without any approximations
      !   if (tr%tr_etrs_approx == OPTION__MAXWELLTDETRSAPPROX__NO) then
      !     !AFE_ALLOCATE(ff_rs_inhom_1(1:gr%mesh%np_part,ff_dim))
      !     !AFE_ALLOCATE(ff_rs_inhom_2(1:gr%mesh%np_part,ff_dim))
      !     ! RS state propagation
      !     call hamiltonian_mxll_update(hm, time=inter_time)
      !     if (pml_check) then
      !       call pml_propagation_stage_1(hm, gr, st, tr, ff_rs_state, ff_rs_state_pml)
      !       hm%cpml_hamiltonian = .true.
      !     end if
      !     call exponential_mxll_apply(hm, gr, st, tr, inter_time, inter_dt, ff_rs_state)
      !     if (pml_check) then
      !       hm%cpml_hamiltonian = .false.
      !       call pml_propagation_stage_2(hm, gr, st, tr, inter_time, inter_dt, delay, &
      !                                           ff_rs_state_pml, ff_rs_state)
      !     end if
      !     ! inhomogeneity propagation
      !     call transform_rs_densities_forward(hm, rs_charge_density_t1, rs_current_density_t1, ff_rs_inhom_1)
      !     call transform_rs_densities_forward(hm, rs_charge_density_t2, rs_current_density_t2, ff_rs_inhom_2)
      !     ff_rs_inhom_1 = ff_rs_inhom_1 + (ff_rs_inhom_2 - ff_rs_inhom_1) / real(inter_steps) * inter_dt * (ii-1)
      !     ff_rs_inhom_2 = ff_rs_inhom_1 + (ff_rs_inhom_2 - ff_rs_inhom_1) / real(inter_steps) * inter_dt * (ii)
      !     call exponential_mxll_apply(hm, gr, st, tr, inter_time, inter_dt, ff_rs_inhom_1)
      !     ! add terms U(time+dt,time)J(time) and J(time+dt)
      !     ff_rs_state = ff_rs_state + M_FOURTH * inter_dt * (ff_rs_inhom_1 + ff_rs_inhom_2)
      !     call transform_rs_densities_forward(hm, rs_charge_density_t1, rs_current_density_t1, ff_rs_inhom_1)
      !     call transform_rs_densities_forward(hm, rs_charge_density_t2, rs_current_density_t2, ff_rs_inhom_2)
      !     ff_rs_inhom_1 = M_HALF * (ff_rs_inhom_1 + ff_rs_inhom_2)
      !     ff_rs_inhom_2 = M_HALF * (ff_rs_inhom_1 + ff_rs_inhom_2)
      !     call exponential_mxll_apply(hm, gr, st, tr, inter_time, inter_dt/M_TWO, ff_rs_inhom_1)
      !     call exponential_mxll_apply(hm, gr, st, tr, inter_time, -inter_dt/M_TWO, ff_rs_inhom_2)
      !     ! add terms U(time+dt/2,time)J(time) and U(time,time+dt/2)J(time+dt)
      !     ff_rs_state = ff_rs_state + M_FOURTH * inter_dt * (ff_rs_inhom_1 + ff_rs_inhom_2)
      !     !AFE_DEALLOCATE_A(ff_rs_inhom_1)
      !     !AFE_DEALLOCATE_A(ff_rs_inhom_2)

      !   ! Maxwell no etrs propagation, just straight forward with the time evolution operator
      !   else if (tr%tr_etrs_approx == OPTION__MAXWELLTDETRSAPPROX__NO_ETRS) then
      !     !AFE_ALLOCATE(ff_rs_inhom_1(1:gr%mesh%np_part,ff_dim))
      !     !AFE_ALLOCATE(ff_rs_inhom_2(1:gr%mesh%np_part,ff_dim))
      !     call transform_rs_densities_forward(hm, rs_charge_density_t1, rs_current_density_t1, ff_rs_inhom_1)
      !     call transform_rs_densities_forward(hm, rs_charge_density_t2, rs_current_density_t2, ff_rs_inhom_2)
      !     ! RS state propagation
      !     call hamiltonian_mxll_update(hm, time=inter_time)
      !     if (pml_check) then
      !       call pml_propagation_stage_1(hm, gr, st, tr, ff_rs_state, ff_rs_state_pml)
      !       hm%cpml_hamiltonian = .true.
      !     end if
      !     call exponential_mxll_apply(hm, gr, st, tr, inter_time, inter_dt, ff_rs_state)
      !     if (pml_check) then
      !       hm%cpml_hamiltonian = .false.
      !       call pml_propagation_stage_2(hm, gr, st, tr, inter_time, inter_dt, delay, ff_rs_state_pml, ff_rs_state)
      !     end if
      !     ! inhomogeneity propagation
      !     ff_rs_inhom_1 = ff_rs_inhom_1 + (ff_rs_inhom_2 - ff_rs_inhom_1) / real(inter_steps) * inter_dt * (ii-1)
      !     ff_rs_inhom_2 = ff_rs_inhom_1 + (ff_rs_inhom_2 - ff_rs_inhom_1) / real(inter_steps) * inter_dt * (ii)
      !     call exponential_mxll_apply(hm, gr, st, tr, inter_time, inter_dt, ff_rs_inhom_1)
      !     ! add terms U(time+dt,time)J(time) and J(time+dt)
      !     ff_rs_state = ff_rs_state + M_HALF * inter_dt * (ff_rs_inhom_1 + ff_rs_inhom_2)
      !     !AFE_DEALLOCATE_A(ff_rs_inhom_1)
      !     !AFE_DEALLOCATE_A(ff_rs_inhom_2)

      !   ! Maxwell etrs propagation for small current density changes in matter time to assume them as constant
      !   else if (tr%tr_etrs_approx == OPTION__MAXWELLTDETRSAPPROX__CONST_STEPS) then
      !     ! RS state propagation
      !     call hamiltonian_mxll_update(hm, time=inter_time)
      !     if (pml_check) then
      !       call pml_propagation_stage_1(hm, gr, st, tr, ff_rs_state, ff_rs_state_pml)
      !       hm%cpml_hamiltonian = .true.
      !     end if
      !     call exponential_mxll_apply(hm, gr, st, tr, inter_time, inter_dt, ff_rs_state)
      !     if (pml_check) then
      !       hm%cpml_hamiltonian = .false.
      !       call pml_propagation_stage_2(hm, gr, st, tr, inter_time, inter_dt, delay, &
      !                                           ff_rs_state_pml, ff_rs_state)
      !     end if
      !     ff_rs_state = ff_rs_state + M_FOURTH * inter_dt * ff_rs_inhom_1

      !   end if

      else

        ! RS state propagation
        call hamiltonian_mxll_update(hm, time=inter_time)
        if (pml_check) then
          call pml_propagation_stage_1(hm, gr, st, tr, ff_rs_state, ff_rs_state_pml)
          hm%cpml_hamiltonian = .true.
        end if

        call zbatch_init(st%rsb, 1, 1, st%d%dim, gr%mesh%np_part)

        do istate = 1, st%d%dim
          call batch_set_state(st%rsb, istate, gr%mesh%np_part, rs_state(:, istate))
        end do

        call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, st%rsb, inter_dt)

        do istate = 1, st%d%dim
          call batch_get_state(st%rsb, istate, gr%mesh%np_part, rs_state(:, istate))
        end do

        call st%rsb%end()

        if (pml_check) then
          hm%cpml_hamiltonian = .false.
          call pml_propagation_stage_2(hm, gr, st, tr, inter_time, inter_dt, delay, ff_rs_state_pml, ff_rs_state)
        end if
      end if

      ! PML convolution function update
      if (pml_check) then
        call cpml_conv_function_update(hm, gr, ff_rs_state_pml, ff_dim)
      end if

      ! back tranformation of RS state representation
!      call transform_rs_state_backward(hm, gr, st, ff_rs_state, rs_state)

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
      call mask_absorbing_boundaries(gr, hm, st, tr, inter_time, inter_dt, delay, rs_state)

      if (tr%bc_plane_waves) then
        ! Propagation dt with H(inter_time+inter_dt) for plane waves boundaries
        call plane_waves_boundaries_calculation(hm, tr, st, gr%mesh, inter_time+inter_dt, delay ,rs_state)
      end if

    end do

    if (tr%tr_etrs_approx == OPTION__MAXWELLTDETRSAPPROX__CONST_STEPS) then
      SAFE_DEALLOCATE_A(ff_rs_inhom_1)
    end if

    SAFE_DEALLOCATE_A(ff_rs_state)

    if (pml_check) then
      SAFE_DEALLOCATE_A(ff_rs_state_pml)
    end if

    POP_SUB(propagation_mxll_etrs)
  end subroutine propagation_mxll_etrs


  ! ---------------------------------------------------------
  subroutine exponential_mxll_apply(hm, gr, st, tr, time, dt, ff)
    type(hamiltonian_mxll_t),        intent(in)    :: hm
    type(grid_t),               intent(in)    :: gr
    type(states_mxll_t),             intent(inout) :: st
    type(propagator_mxll_t), intent(inout) :: tr
    FLOAT,                      intent(in)    :: time
    FLOAT,                      intent(in)    :: dt
    CMPLX,                      intent(inout) :: ff(:,:)

    PUSH_SUB(exponential_mxll_apply)

    select case (tr%op_method)
      case(OPTION__MAXWELLTDOPERATORMETHOD__OP_FD)
!        call exponential_apply(tr%te, gr%mesh, gr%der, hm, ff, 1, 1, dt, time)
      case(OPTION__MAXWELLTDOPERATORMETHOD__OP_FFT)
!        call exponential_mxll_fft_apply(hm, gr, st, tr, time, dt, ff)
    end select

    POP_SUB(exponential_mxll_apply)
  end subroutine exponential_mxll_apply


  !----------------------------------------------------------
 ! subroutine exponential_mxll_fft_apply(hm, gr, st, tr, time, dt, ff)
 ! end subroutine exponential_mxll_fft_apply


  !----------------------------------------------------------
!  subroutine exponential_fft_taylor_series(tr, mesh, hm, st, dt, k_vec, cube, dim, ff_dim, fft_func_in, fft_func_out)
!  end subroutine exponential_fft_taylor_series


  ! ---------------------------------------------------------
  subroutine set_medium_rs_state(st, gr, hm)
    type(states_mxll_t),      intent(inout) :: st
    type(grid_t),        intent(in)    :: gr
    type(hamiltonian_mxll_t), intent(in)    :: hm

    integer :: ip, ip_in, il, idim

    PUSH_SUB(set_medium_rs_state)

    SAFE_ALLOCATE(st%ep(1:gr%mesh%np_part))
    SAFE_ALLOCATE(st%mu(1:gr%mesh%np_part))
    st%ep = P_ep
    st%mu = P_mu
    if (hm%medium_box) then
      do il=1, hm%medium_box_number
        do ip_in=1, hm%medium_box_points_number(il)
          ip = hm%medium_box_points_map(ip_in,il)
          st%ep(ip) = hm%medium_box_ep(ip_in,il)
          st%mu(ip) = hm%medium_box_mu(ip_in,il)
        end do
      end do
    end if

    do idim=1, st%d%dim
      if (hm%bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MEDIUM) then
        do ip_in=1, hm%bc%medium_points_number(idim)
          ip = hm%bc%medium_points_map(ip_in,idim)
          st%ep(ip) = hm%bc%medium_ep(ip_in,idim)
          st%mu(ip) = hm%bc%medium_mu(ip_in,idim)
        end do
      end if
    end do

    POP_SUB(set_medium_rs_state)
  end subroutine set_medium_rs_state


  ! ---------------------------------------------------------
  subroutine rs_state_to_cube_map(gr, hm, st)
    type(grid_t),        intent(in)    :: gr
    type(hamiltonian_mxll_t), intent(in)    :: hm
    type(states_mxll_t),      intent(inout) :: st

    integer :: ipx, ipy, ipz, idx, np(st%d%dim), idim

    PUSH_SUB(rs_state_to_cube_map)

    do idim=1, st%d%dim
      ASSERT(hm%cube%rs_n(idim) == hm%cube%fs_n(idim))
      np(idim) = hm%cube%rs_n(idim)
    end do

    SAFE_ALLOCATE(st%rs_state_fft_map(1:np(1),1:np(2),1:np(3)))
    SAFE_ALLOCATE(st%rs_state_fft_map_inv(1:gr%mesh%np_global,1:3))

    idx = 0
    do ipx=1,np(1)
      do ipy=1,np(2)
        do ipz=1,np(3)
          idx = idx+1
          st%rs_state_fft_map(ipx,ipy,ipz) = idx
          st%rs_state_fft_map_inv(idx,1) = ipx
          st%rs_state_fft_map_inv(idx,2) = ipy
          st%rs_state_fft_map_inv(idx,3) = ipz
        end do
      end do
    end do

    POP_SUB(rs_state_to_cube_map)
  end subroutine rs_state_to_cube_map


  ! ---------------------------------------------------------
  subroutine transform_rs_state_forward(hm, gr, st, rs_state, ff_rs_state)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    type(grid_t),        intent(in)    :: gr
    type(states_mxll_t),      intent(in)    :: st
    CMPLX,               intent(in)    :: rs_state(:,:)
    CMPLX,               intent(inout) :: ff_rs_state(:,:)

    CMPLX, allocatable :: rs_state_plus(:,:), rs_state_minus(:,:) 

    if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then
      message(1) = "Maxwell solver in linear media not yet implemented"
      call messages_fatal(1)
      ! !AFE_ALLOCATE(rs_state_plus(1:gr%mesh%np_part,1:st%d%dim))
      ! !AFE_ALLOCATE(rs_state_minus(1:gr%mesh%np_part,1:st%d%dim))
      ! rs_state_plus  = rs_state
      ! rs_state_minus = real(rs_state) - M_zI * aimag(rs_state)
      ! call transform_rs_state_to_6x6_rs_state_forward(rs_state_plus, rs_state_minus, ff_rs_state)
      ! !AFE_DEALLOCATE_A(rs_state_plus)
      ! !AFE_DEALLOCATE_A(rs_state_minus)
    else if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_GAUSS) then
      call transform_rs_state_to_4x4_rs_state_forward(rs_state, ff_rs_state)
    else 
      ff_rs_state(:,:) = rs_state(:,:)
    end if

  end subroutine transform_rs_state_forward


  ! ---------------------------------------------------------
  subroutine transform_rs_state_backward(hm, gr, st, ff_rs_state, rs_state)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    type(grid_t),        intent(in)    :: gr
    type(states_mxll_t),      intent(in)    :: st
    CMPLX,               intent(in)    :: ff_rs_state(:,:)
    CMPLX,               intent(inout) :: rs_state(:,:)

    CMPLX, allocatable :: rs_state_plus(:,:), rs_state_minus(:,:) 

    if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then
      call transform_rs_state_to_6x6_rs_state_backward(ff_rs_state, rs_state)
    else if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then
      call transform_rs_state_to_4x4_rs_state_backward(ff_rs_state, rs_state)
    else 
      rs_state(:,:) = ff_rs_state(:,:)
    end if

  end subroutine transform_rs_state_backward


  ! -------------------------------------------------------
  subroutine transform_rs_densities_forward(hm, rs_charge_density, rs_current_density, ff_density)

    type(hamiltonian_mxll_t), intent(in)    :: hm
    CMPLX,               intent(in)    :: rs_charge_density(:)
    CMPLX,               intent(in)    :: rs_current_density(:,:)
    CMPLX,               intent(inout) :: ff_density(:,:)

    if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then
      call transform_rs_densities_to_6x6_rs_densities_forward(rs_charge_density, rs_current_density, ff_density)
    else if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then
      call transform_rs_densities_to_4x4_rs_densities_forward(rs_charge_density, rs_current_density, ff_density)
    else 
      ff_density(:,:) = rs_current_density(:,:)
    end if

  end subroutine transform_rs_densities_forward


  ! -------------------------------------------------------
  subroutine transform_rs_densities_backward(hm, ff_density, rs_charge_density, rs_current_density)

    type(hamiltonian_mxll_t), intent(in)    :: hm
    CMPLX,               intent(in)    :: ff_density(:,:)
    CMPLX,               intent(inout) :: rs_charge_density(:)
    CMPLX,               intent(inout) :: rs_current_density(:,:)

    if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then
      call transform_rs_densities_to_6x6_rs_densities_backward(ff_density, rs_charge_density, rs_current_density)
    else if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then
      call transform_rs_densities_to_4x4_rs_densities_backward(ff_density, rs_charge_density, rs_current_density)
    else 
      rs_current_density = ff_density(:,:)
    end if

  end subroutine transform_rs_densities_backward


  !----------------------------------------------------------
  subroutine transform_rs_state_to_6x6_rs_state_forward(rs_state_3x3_plus, rs_state_3x3_minus, rs_state_6x6)
    CMPLX, intent(in)    :: rs_state_3x3_plus(:,:)
    CMPLX, intent(in)    :: rs_state_3x3_minus(:,:)
    CMPLX, intent(inout) :: rs_state_6x6(:,:)

    ! no push_sub, called to frequently
    rs_state_6x6(:,1) = rs_state_3x3_plus(:,1)
    rs_state_6x6(:,2) = rs_state_3x3_plus(:,2)
    rs_state_6x6(:,3) = rs_state_3x3_plus(:,3)
    rs_state_6x6(:,4) = rs_state_3x3_minus(:,1)
    rs_state_6x6(:,5) = rs_state_3x3_minus(:,2)
    rs_state_6x6(:,6) = rs_state_3x3_minus(:,3)

  end subroutine transform_rs_state_to_6x6_rs_state_forward


  !----------------------------------------------------------
  subroutine transform_rs_state_to_6x6_rs_state_backward(rs_state_6x6, rs_state)
    CMPLX, intent(in)    :: rs_state_6x6(:,:)
    CMPLX, intent(inout) :: rs_state(:,:)

    ! no push_sub, called to frequently
    rs_state(:,1) = M_HALF * real(rs_state_6x6(:,1)+rs_state_6x6(:,4)) + M_HALF * M_zI * aimag(rs_state_6x6(:,1)-rs_state_6x6(:,4))
    rs_state(:,2) = M_HALF * real(rs_state_6x6(:,2)+rs_state_6x6(:,5)) + M_HALF * M_zI * aimag(rs_state_6x6(:,2)-rs_state_6x6(:,5))
    rs_state(:,3) = M_HALF * real(rs_state_6x6(:,3)+rs_state_6x6(:,6)) + M_HALF * M_zI * aimag(rs_state_6x6(:,3)-rs_state_6x6(:,6))

  end subroutine transform_rs_state_to_6x6_rs_state_backward


  !----------------------------------------------------------
  subroutine transform_rs_densities_to_6x6_rs_densities_forward(rs_charge_density, rs_current_density, rs_density_6x6)
    CMPLX, intent(in)    :: rs_charge_density(:)
    CMPLX, intent(in)    :: rs_current_density(:,:)
    CMPLX, intent(inout) :: rs_density_6x6(:,:)

    ! no push_sub, called to frequently
    rs_density_6x6(:,1) = rs_current_density(:,1)
    rs_density_6x6(:,2) = rs_current_density(:,2)
    rs_density_6x6(:,3) = rs_current_density(:,3)
    rs_density_6x6(:,4) = rs_current_density(:,1)
    rs_density_6x6(:,5) = rs_current_density(:,2)
    rs_density_6x6(:,6) = rs_current_density(:,3)

  end subroutine transform_rs_densities_to_6x6_rs_densities_forward


  !----------------------------------------------------------
  subroutine transform_rs_densities_to_6x6_rs_densities_backward(rs_density_6x6, rs_charge_density, rs_current_density)
    CMPLX, intent(in)    :: rs_density_6x6(:,:)
    CMPLX, intent(inout) :: rs_charge_density(:)
    CMPLX, intent(inout) :: rs_current_density(:,:)

    ! no push_sub, called to frequently
    rs_current_density(:,1) = M_HALF * real( rs_density_6x6(:,1) + rs_density_6x6(:,4) )
    rs_current_density(:,2) = M_HALF * real( rs_density_6x6(:,2) + rs_density_6x6(:,5) )
    rs_current_density(:,3) = M_HALF * real( rs_density_6x6(:,3) + rs_density_6x6(:,6) )

  end subroutine transform_rs_densities_to_6x6_rs_densities_backward


  !----------------------------------------------------------
  subroutine transform_rs_state_to_4x4_rs_state_forward(rs_state_3x3, rs_state_4x4)
    CMPLX, intent(in)    :: rs_state_3x3(:,:)
    CMPLX, intent(inout) :: rs_state_4x4(:,:)

    ! no push_sub, called to frequently
    rs_state_4x4(:,1) = - M_z1 * rs_state_3x3(:,1) + M_zI * rs_state_3x3(:,2)
    rs_state_4x4(:,2) =   M_z1 * rs_state_3x3(:,3)
    rs_state_4x4(:,3) =   M_z1 * rs_state_3x3(:,3)
    rs_state_4x4(:,4) =   M_z1 * rs_state_3x3(:,1) + M_zI * rs_state_3x3(:,2)

  end subroutine transform_rs_state_to_4x4_rs_state_forward


  !----------------------------------------------------------
  subroutine transform_rs_state_to_4x4_rs_state_backward(rs_state_4x4, rs_state_3x3)
    CMPLX, intent(in)    :: rs_state_4x4(:,:)
    CMPLX, intent(inout) :: rs_state_3x3(:,:)

    ! no push_sub, called to frequently
    rs_state_3x3(:,1) = - M_z1 * M_HALF * rs_state_4x4(:,1) + M_z1 * M_HALF * rs_state_4x4(:,4)
    rs_state_3x3(:,2) = - M_zI * M_HALF * rs_state_4x4(:,1) - M_zI * M_HALF * rs_state_4x4(:,4)
    rs_state_3x3(:,3) =   M_z1 * M_HALF * rs_state_4x4(:,2) + M_z1 * M_HALF * rs_state_4x4(:,3)

  end subroutine transform_rs_state_to_4x4_rs_state_backward


  !----------------------------------------------------------
  subroutine transform_rs_densities_to_4x4_rs_densities_forward(rs_charge_density, rs_current_density, rs_density_4x4)
    CMPLX, intent(in)    :: rs_charge_density(:)
    CMPLX, intent(in)    :: rs_current_density(:,:)
    CMPLX, intent(inout) :: rs_density_4x4(:,:)

    ! no push_sub, called to frequently
    rs_density_4x4(:,1) = - M_z1 * rs_current_density(:,1) + M_zI * rs_current_density(:,2)
    rs_density_4x4(:,2) =   M_z1 * rs_current_density(:,3) - M_z1 * rs_charge_density(:)
    rs_density_4x4(:,3) =   M_z1 * rs_current_density(:,3) + M_z1 * rs_charge_density(:)
    rs_density_4x4(:,4) =   M_z1 * rs_current_density(:,1) + M_zI * rs_current_density(:,2)

  end subroutine transform_rs_densities_to_4x4_rs_densities_forward


  !----------------------------------------------------------
  subroutine transform_rs_densities_to_4x4_rs_densities_backward(rs_density_4x4, rs_charge_density, rs_current_density)
    CMPLX, intent(in)    :: rs_density_4x4(:,:)
    CMPLX, intent(inout) :: rs_charge_density(:)
    CMPLX, intent(inout) :: rs_current_density(:,:)

    ! no push_sub, called to frequently
    rs_charge_density(:)    = - M_z1 * M_HALF * rs_density_4x4(:,2) + M_z1 * M_HALF * rs_density_4x4(:,3)
    rs_current_density(:,1) = - M_z1 * M_HALF * rs_density_4x4(:,1) + M_z1 * M_HALF * rs_density_4x4(:,4)
    rs_current_density(:,2) = - M_zI * M_HALF * rs_density_4x4(:,1) - M_zI * M_HALF * rs_density_4x4(:,4)
    rs_current_density(:,3) = - M_z1 * M_HALF * rs_density_4x4(:,2) + M_z1 * M_HALF * rs_density_4x4(:,3)

  end subroutine transform_rs_densities_to_4x4_rs_densities_backward


  !----------------------------------------------------------
  subroutine calculate_matter_longitudinal_field(gr_mxll, st_mxll, hm_mxll, gr_elec, st_elec, hm_elec, rs_state_matter, geo)
    type(grid_t),                  intent(in)    :: gr_mxll
    type(states_mxll_t),                intent(in)    :: st_mxll
    type(hamiltonian_mxll_t),           intent(in)    :: hm_mxll
    type(grid_t),                  intent(in)    :: gr_elec
    type(states_elec_t),                intent(in)    :: st_elec
    type(hamiltonian_elec_t),           intent(in)    :: hm_elec
    CMPLX,                         intent(inout) :: rs_state_matter(:,:)
    type(geometry_t),    optional, intent(in)    :: geo

    integer            :: idim, ip
    FLOAT              :: dd, width
    CMPLX, allocatable :: tmp_pot_ma_gr(:,:), tmp_pot_mx_gr(:,:), tmp_grad_mx_gr(:,:)

!    !AFE_ALLOCATE(tmp_pot_ma_gr(1:gr_elec%mesh%np_part,1))
    SAFE_ALLOCATE(tmp_pot_mx_gr(1:gr_mxll%mesh%np_part,1))
    SAFE_ALLOCATE(tmp_grad_mx_gr(1:gr_mxll%mesh%np,1:gr_mxll%sb%dim))

    PUSH_SUB(calculate_matter_longitudinal_field)

!    tmp_pot_ma_gr(:,:) = M_z0
!    tmp_pot_ma_gr(1:gr_elec%mesh%np,1) = M_z1 * ( hm_elec%vhartree(1:gr_elec%mesh%np) + hm_elec%ep%vpsl(1:gr_elec%mesh%np) )

!    call zma_mesh_to_mx_mesh(st_mxll, gr_mxll, st_elec, gr_elec, tmp_pot_ma_gr, tmp_pot_mx_gr, 1)

    tmp_pot_mx_gr(:,:) = M_ZERO
    tmp_grad_mx_gr(:,:) = M_ZERO
    call zderivatives_grad(gr_mxll%der, tmp_pot_mx_gr(:,1), tmp_grad_mx_gr(:,:), set_bc = .false.)
    tmp_grad_mx_gr = - tmp_grad_mx_gr

    rs_state_matter = M_z0
    call build_rs_state(real(tmp_grad_mx_gr(1:gr_mxll%mesh%np,:)), aimag(tmp_grad_mx_gr(1:gr_mxll%mesh%np,:)), &
      st_mxll%rs_sign, rs_state_matter(1:gr_mxll%mesh%np,:), gr_mxll%mesh, st_mxll%ep(1:gr_mxll%mesh%np), &
      st_mxll%mu(1:gr_mxll%mesh%np), gr_mxll%mesh%np)

    SAFE_DEALLOCATE_A(tmp_pot_ma_gr)
    SAFE_DEALLOCATE_A(tmp_pot_mx_gr)
    SAFE_DEALLOCATE_A(tmp_grad_mx_gr)

    POP_SUB(calculate_matter_longitudinal_field)
  end subroutine calculate_matter_longitudinal_field


  !----------------------------------------------------------
  subroutine get_vector_pot_and_transverse_field(trans_calc_method, gr_mxll, hm_mxll, st_mxll, tr_mxll, &
    gr, hm, st, tr, poisson_solver, time, field, transverse_field, vector_potential, geo)
    integer,                    intent(in)    :: trans_calc_method
    type(grid_t),               intent(in)    :: gr_mxll
    type(hamiltonian_mxll_t),        intent(in)    :: hm_mxll
    type(states_mxll_t),             intent(in)    :: st_mxll
    type(propagator_mxll_t), intent(in)    :: tr_mxll
    type(grid_t),               intent(in)    :: gr
    type(hamiltonian_elec_t),        intent(in)    :: hm
    type(states_elec_t),             intent(in)    :: st
    type(propagator_t),         intent(in)    :: tr
    type(poisson_t),            intent(in)    :: poisson_solver
    FLOAT,                      intent(in)    :: time
    CMPLX,                      intent(inout) :: field(:,:)
    CMPLX,                      intent(inout) :: transverse_field(:,:)
    FLOAT,                      intent(inout) :: vector_potential(:,:)
    type(geometry_t), optional, intent(in)    :: geo

    integer            :: idim, ip, ip_in, np
    FLOAT              :: abs_r
    CMPLX, allocatable :: tmp_field(:,:), tmp_field_2(:,:)
    logical            :: trans_test, vec_pot_test

    PUSH_SUB(get_vector_pot_and_transverse_field)

    transverse_field = M_z0
    vector_potential = M_ZERO

    np = gr_mxll%mesh%np

    ! !%Variable MaxwellTransverseCalculationTest
    ! !%Type logical
    ! !%Default no
    ! !%Section Hamiltonian
    ! !%Description
    ! !% Description follows
    ! !%End
    ! call parse_variable('MaxwellTransverseCalculationTest', .false., trans_test)

    ! !%Variable MaxwellVectorPotentialTest
    ! !%Type logical
    ! !%Default no
    ! !%Section Hamiltonian
    ! !%Description
    ! !% Description follows
    ! !%End
    ! call parse_variable('MaxwellVectorPotentialTest', .false., vec_pot_test)

    ! if (trans_test) then
    !   do ip=1, maxwell_gr%mesh%np
    !     ! Test 1
    !     field(ip,1) = - sqrt(P_ep/M_TWO) &
    !                 * exp(-(maxwell_gr%mesh%x(ip,1)**2+maxwell_gr%mesh%x(ip,2)**2+maxwell_gr%mesh%x(ip,3)**2)/(M_TWO*M_TWO**2)) &
    !                 * maxwell_gr%mesh%x(ip,2)
    !     field(ip,2) = sqrt(P_ep/M_TWO) &
    !                 * exp(-(maxwell_gr%mesh%x(ip,1)**2+maxwell_gr%mesh%x(ip,2)**2+maxwell_gr%mesh%x(ip,3)**2)/(M_TWO*M_TWO**2)) &
    !                 * maxwell_gr%mesh%x(ip,1)
    !     field(ip,3) = M_z0
    !     ! Test 2
    !     !field(ip,1) = M_z0
    !     !field(ip,2) = M_z0
    !     !!field(ip,3) = sqrt(P_ep/M_TWO) * exp(-maxwell_gr%mesh%x(ip,1)**2/(M_TWO*M_TWO**2))
    !     !field(ip,3) = sqrt(P_ep/M_TWO) * cos(CNST(0.35)*M_pi*maxwell_gr%mesh%x(ip,1))
    !   end do
    !   !call helmholtz_decomposition_transverse_field(maxwell_poisson_solver, maxwell_gr, maxwell_hm, hm, field)
    !   !do ip_in=1, maxwell_hm%maxwell_bc%der_bndry_mask_points_number
    !   !  ip = maxwell_hm%maxwell_bc%der_bndry_mask_points_map(ip_in)
    !   !  field(ip,:) = maxwell_hm%maxwell_bc%der_bndry_mask(ip_in) * field(ip,:)
    !   !end do
    !   maxwell_st%maxwell_rs_state(:,:) = field(:,:)
    ! end if

    if (hm_mxll%ma_mx_coupling) then

      ! if (trans_calc_method == OPTION__MAXWELLTRANSFIELDCALCULATIONMETHOD__TRANS_FIELD_MATTER) then
      !   !AFE_ALLOCATE(tmp_field(1:gr_mxll%mesh%np_part,1:st_mxll%d%dim))
      !   if (hm_mxll%ma_mx_coupling_apply .and. (tr%current_prop_test == 0)) then
      !     call calculate_matter_longitudinal_field(gr_mxll, st_mxll, hm_mxll, gr, st, hm, tmp_field)
      !   end if
      !   transverse_field(1:np,:) = field(1:np,:) - tmp_field(1:np,:)
      !   !AFE_DEALLOCATE_A(tmp_field)

!      else if (trans_calc_method == OPTION__MAXWELLTRANSFIELDCALCULATIONMETHOD__TRANS_FIELD_POISSON) then
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
 
      ! else if (trans_calc_method == OPTION__MAXWELLTRANSFIELDCALCULATIONMETHOD__TRANS_FIELD_POISSON_CORR) then
      !   ! plane waves subtraction
      !   ! longitudinal matter field subtraction to get almost transverse field 
      !   call calculate_matter_longitudinal_field(gr_mxll, st_mxll, hm_mxll, gr, st, hm, transverse_field, geo)
      !   if (tr_mxll%bc_plane_waves .and. hm_mxll%plane_waves_apply) then
      !     transverse_field(1:np,:) = field(1:np,:) - transverse_field(1:np,:) - st_mxll%rs_state_plane_waves(1:np,:)
      !   else
      !     transverse_field(1:np,:) = field(1:np,:) - transverse_field(1:np,:)
      !   end if
      !   ! apply helmholtz decomposition for transverse field
      !   call maxwell_helmholtz_decomposition_trans_field(poisson_solver, gr_mxll, hm_mxll, hm, st_mxll, transverse_field)
      !   ! plane waves addition
      !   if (tr_mxll%bc_plane_waves .and. hm_mxll%plane_waves_apply) then
      !     transverse_field(1:np,:) = transverse_field(1:np,:) + st_mxll%rs_state_plane_waves(1:np,:)
      !   end if

      ! else if (trans_calc_method == OPTION__MAXWELLTRANSFIELDCALCULATIONMETHOD__TRANS_FIELD_POISSON_LONG) then
      !   ! plane wave subtraction
      !   if (tr_mxll%bc_plane_waves .and. hm_mxll%plane_waves_apply) then
      !     transverse_field(1:np,:) = field(1:np,:) - st_mxll%rs_state_plane_waves(1:np,:)
      !   else
      !     transverse_field(1:np,:) = field(1:np,:)
      !   end if
      !   ! apply helmholtz decomposition for longitudinal field
      !   call helmholtz_decomposition_long_field(poisson_solver, gr_mxll, transverse_field)
      !   transverse_field(1:np,:) = field(1:np,:) - transverse_field(1:np,:)
      !   ! plane waves addition
      !   if (tr_mxll%bc_plane_waves .and. hm_mxll%plane_waves_apply) then
      !     transverse_field(1:np,:) = transverse_field(1:np,:) + st_mxll%rs_state_plane_waves(1:np,:)
      !   end if

      ! else if (trans_calc_method == OPTION__MAXWELLTRANSFIELDCALCULATIONMETHOD__TRANS_FIELD_POISSON_LONG_CORR) then
      !   !AFE_ALLOCATE(tmp_field(1:gr_mxll%mesh%np,1:st_mxll%d%dim))
      !   ! longitudinal matter field subtraction to get almost transverse field (1st time)
      !   call calculate_matter_longitudinal_field(gr_mxll, st_mxll, hm_mxll, gr, st, hm, tmp_field, geo)
      !   ! plane waves subtraction
      !   if (tr_mxll%bc_plane_waves) then
      !     tmp_field(1:np,:) = field(1:np,:) - tmp_field(1:np,:) - st_mxll%rs_state_plane_waves(1:np,:)
      !   else
      !     tmp_field(1:np,:) = field(1:np,:) - tmp_field(1:np,:)
      !   end if
      !   ! apply helmholtz decomposition for longitudinal field
      !   call helmholtz_decomposition_long_field(poisson_solver, gr_mxll, tmp_field)
      !   transverse_field(1:np,:) = field(1:np,:) - tmp_field(1:np,:)
      !   ! plane waves addition
      !   if (tr_mxll%bc_plane_waves .and. hm_mxll%plane_waves_apply) then
      !     transverse_field(1:np,:) = transverse_field(1:np,:) + st%rs_state_plane_waves(1:np,:)
      !   end if
      !   !AFE_DEALLOCATE_A(tmp_field)

      ! else if (trans_calc_method == OPTION__MAXWELLTRANSFIELDCALCULATIONMETHOD__TRANS_FIELD_POISSON) then
      !   transverse_field(1:np,:) = field(1:np,:)
      !   ! apply helmholtz decomposition for transverse field
      !   call maxwell_helmholtz_decomposition_trans_field(poisson_solver, gr_mxll, hm_mxll, hm, st_mxll, transverse_field)

      ! end if

    else

      transverse_field(:,:) = field

    end if

    ! if (trans_test) then

    !   !AFE_ALLOCATE(tmp_field(1:maxwell_gr%mesh%np_part,1:maxwell_st%d%dim))
    !   !AFE_ALLOCATE(tmp_field_2(1:maxwell_gr%mesh%np_part,1:maxwell_st%d%dim))

    !   tmp_field = M_z0
    !   tmp_field_2 = M_z0
    !   call zderivatives_div(maxwell_gr%der, transverse_field(:,:), tmp_field(:,1), set_bc = .false.)
    !   maxwell_st%maxwell_test_output(:,1) = sqrt(M_TWO/P_ep) * real(tmp_field(:,1))
    !   maxwell_st%maxwell_test_output(:,2) = sqrt(M_TWO/P_ep) * real(tmp_field(:,2))
    !   maxwell_st%maxwell_test_output(:,3) = sqrt(M_TWO/P_ep) * real(tmp_field(:,3))
    !   tmp_field(1:np,:) = field(1:np,:) - transverse_field(1:np,:)
    !   call zderivatives_curl(maxwell_gr%der, tmp_field(:,:), tmp_field_2(:,:), set_bc = .false.)
    !   maxwell_st%maxwell_test_output(:,:) = real(maxwell_st%maxwell_test_output(:,:)) + M_zI*sqrt(M_TWO/P_ep)*real(tmp_field_2(:,:))

    !   !AFE_DEALLOCATE_A(tmp_field)
    !   !AFE_DEALLOCATE_A(tmp_field_2)

    ! end if

    ! if (vec_pot_test) then
    !   do ip=1, maxwell_gr%mesh%np
    !     field(ip,1) = - M_zI * M_ONE/M_FOUR * &
    !                   exp((-maxwell_gr%mesh%x(ip,1)**2-maxwell_gr%mesh%x(ip,2)**2-maxwell_gr%mesh%x(ip,3)**2)/(M_TWO*M_TWO**2)) * &
    !                   maxwell_gr%mesh%x(ip,2)
    !     field(ip,2) = M_zI * M_ONE/M_FOUR * &
    !                   exp((-maxwell_gr%mesh%x(ip,1)**2-maxwell_gr%mesh%x(ip,2)**2-maxwell_gr%mesh%x(ip,3)**2)/(M_TWO*M_TWO**2)) * &
    !                   maxwell_gr%mesh%x(ip,1)
    !     field(ip,3) = M_z0
    !   end do
    !   maxwell_st%maxwell_test_output(:,:) = field(:,:)
    !   maxwell_st%maxwell_rs_state(:,:) = sqrt(M_ONE/(M_TWO*P_mu)) * field(:,:)
    !   field(:,:) = sqrt(M_ONE/(M_TWO*P_mu)) * field(:,:)
    ! end if

    ! call calculate_vector_potential(maxwell_poisson_solver, maxwell_gr, maxwell_st, maxwell_tr, field, vector_potential)

    ! if (vec_pot_test) then

    !   !AFE_ALLOCATE(tmp_field(1:maxwell_gr%mesh%np_part,1:maxwell_st%d%dim))

    !   tmp_field = M_z0
    !   field = M_zI * vector_potential
    !   call zderivatives_curl(maxwell_gr%der, field(:,:), tmp_field(:,:), set_bc = .false.)

    !   !AFE_DEALLOCATE_A(tmp_field)

    ! end if

    POP_SUB(get_vector_pot_and_transverse_field)

    contains

      subroutine calculate_vector_potential(poisson_solver, gr, st, tr, field, vector_potential)
        type(poisson_t),            intent(in)    :: poisson_solver
        type(grid_t),               intent(in)    :: gr
        type(states_mxll_t),             intent(in)    :: st
        type(propagator_mxll_t), intent(in)    :: tr
        CMPLX,                      intent(in)    :: field(:,:)
        FLOAT,                      intent(inout) :: vector_potential(:,:)

        integer :: idim, wn, ip
        FLOAT   :: k_vector(MAX_DIM), k_vector_abs, e_field(MAX_DIM), b_field(MAX_DIM)
        FLOAT, allocatable :: dtmp(:,:)

        SAFE_ALLOCATE(dtmp(1:gr%mesh%np_part,1:3))

        dtmp = M_ZERO

        call get_magnetic_field_state(field, gr%mesh, st%rs_sign, vector_potential, st%mu, gr%mesh%np_part)
        dtmp = vector_potential
        call dderivatives_curl(gr%der, dtmp, vector_potential, set_bc = .false.)
        do idim=1, st%d%dim
          call dpoisson_solve(poisson_solver, dtmp(:,idim), vector_potential(:,idim), .true.)
        end do
        vector_potential = M_ONE / (M_FOUR * M_PI) * vector_potential

        SAFE_DEALLOCATE_A(dtmp)

      end subroutine calculate_vector_potential

  end subroutine get_vector_pot_and_transverse_field


  ! ---------------------------------------------------------
  subroutine derivatives_boundary_mask(bc, mesh, hm)
    type(bc_mxll_t),  intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    type(hamiltonian_mxll_t), intent(in)    :: hm

    integer :: ip, ip_in, point_info, idim
    FLOAT   :: bounds(2,MAX_DIM), xx(MAX_DIM)
    FLOAT   :: ddv(MAX_DIM), tmp(MAX_DIM), width(MAX_DIM)
    FLOAT, allocatable :: mask(:)

    PUSH_SUB(derivatives_boundary_mask)

    if (hm%bc_zero) then
      bounds(1,:) = ( mesh%idx%nr(2,:) - 2 * mesh%idx%enlarge(:) ) * mesh%spacing(:)
      bounds(2,:) = ( mesh%idx%nr(2,:) -     mesh%idx%enlarge(:) ) * mesh%spacing(:)
    else if (hm%bc_constant) then
      bounds(1,:) = ( mesh%idx%nr(2,:) - 2 * mesh%idx%enlarge(:) ) * mesh%spacing(:)
      bounds(2,:) = ( mesh%idx%nr(2,:) -     mesh%idx%enlarge(:) ) * mesh%spacing(:)
    else if (hm%bc_plane_waves) then
      bounds(1,:) = ( mesh%idx%nr(2,:) - 2 * mesh%idx%enlarge(:) ) * mesh%spacing(:)
      bounds(2,:) = ( mesh%idx%nr(2,:) -     mesh%idx%enlarge(:) ) * mesh%spacing(:)
    end if

    ip_in=0
    do ip=1, mesh%np
      xx(:) = mesh%x(ip,:)
      if ( (abs(xx(1))<=bounds(2,1)) .and. (abs(xx(2))<=bounds(2,2)) .and. (abs(xx(3))<=bounds(2,3)) ) then
        if ( (abs(xx(1))>bounds(1,1)) .or. (abs(xx(2))>bounds(1,2)) .or. (abs(xx(3))>bounds(1,3)) ) then
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
      xx(:) = mesh%x(ip,:)
      if ( (abs(xx(1))<=bounds(2,1)) .and. (abs(xx(2))<=bounds(2,2)) .and. (abs(xx(3))<=bounds(2,3)) ) then
        if ( (abs(xx(1))>bounds(1,1)) .or. (abs(xx(2))>bounds(1,2)) .or. (abs(xx(3))>bounds(1,3)) ) then
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
      ddv(:) = abs(mesh%x(ip,:)) - bounds(1,:)
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

    POP_SUB(derivatives_boundary_mask)
  end subroutine derivatives_boundary_mask


  ! ---------------------------------------------------------
  subroutine inner_and_outer_points_mapping(mesh, st, bounds)
    type(mesh_t),        intent(in)    :: mesh
    type(states_mxll_t),      intent(inout) :: st
    FLOAT,               intent(in)    :: bounds(:,:)

    integer :: ip, ip_in, ip_bd, point_info
    FLOAT   :: xx(3)

    PUSH_SUB(inner_and_outer_points_mapping)

    ! allocate inner and boundary points points map
    ip_in=0
    ip_bd=0
    do ip=1, mesh%np
      xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim)
      if ( (abs(xx(1))<=bounds(2,1)) .and. (abs(xx(2))<=bounds(2,2)) .and. (abs(xx(3))<=bounds(2,3)) ) then
        if ( (abs(xx(1))>bounds(1,1)) .or. (abs(xx(2))>bounds(1,2)) .or. (abs(xx(3))>bounds(1,3)) ) then
          point_info = 1
        else
          point_info = 0
        end if
      else
        point_info = -1
      end if
      if (point_info == 0) then
        ip_in = ip_in + 1
      else
        ip_bd = ip_bd + 1
      end if
    end do
    st%inner_points_number = ip_in
    SAFE_ALLOCATE(st%inner_points_map(1:ip_in))
    st%boundary_points_number = ip_bd
    SAFE_ALLOCATE(st%boundary_points_map(1:ip_bd))

    ! inner and boundary points mapping
    ip_in=0
    ip_bd=0
    do ip=1, mesh%np
      xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim)
      if ( (abs(xx(1))<=bounds(2,1)) .and. (abs(xx(2))<=bounds(2,2)) .and. (abs(xx(3))<=bounds(2,3)) ) then
        if ( (abs(xx(1))>bounds(1,1)) .or. (abs(xx(2))>bounds(1,2)) .or. (abs(xx(3))>bounds(1,3)) ) then
          point_info = 1
        else
          point_info = 0
        end if
      else
        point_info = -1
      end if
      if (point_info == 0) then
        ip_in = ip_in + 1
        st%inner_points_map(ip_in) = ip
      else
        ip_bd = ip_bd + 1
        st%boundary_points_map(ip_bd) = ip
      end if
    end do

    POP_SUB(inner_and_outer_points_mapping)
  end subroutine inner_and_outer_points_mapping


  ! ---------------------------------------------------------
  subroutine surface_grid_points_mapping(mesh, st, bounds)
    type(mesh_t),        intent(in)    :: mesh
    type(states_mxll_t),      intent(inout) :: st
    FLOAT,               intent(in)    :: bounds(:,:)

    integer :: ix, ix_max, iix, iy, iy_max, iiy, iz, iz_max, iiz, idx1, idx2, rankmin, ip_local, ip_global, nn_max
    integer, allocatable :: nn(:,:,:,:)
    FLOAT   :: rr(3), rr_min(3), rr_max(3), delta(3), dmin, vec(2), min_1(3), max_1(3), min_2(3), max_2(3)

    PUSH_SUB(surface_grid_points_mapping)

    st%surface_grid_rows_number(1) = 3
    ix_max  = st%surface_grid_rows_number(1)
    st%surface_grid_rows_number(2) = 3
    iy_max  = st%surface_grid_rows_number(2)
    st%surface_grid_rows_number(3) = 3
    iz_max  = st%surface_grid_rows_number(3)

    delta(1) = M_TWO * abs(bounds(1,1)) / float(ix_max)
    delta(2) = M_TWO * abs(bounds(1,2)) / float(iy_max)
    delta(3) = M_TWO * abs(bounds(1,3)) / float(iz_max)

    st%surface_grid_element(1) = delta(2) * delta(3)
    st%surface_grid_element(2) = delta(1) * delta(3)
    st%surface_grid_element(3) = delta(1) * delta(2)

    SAFE_ALLOCATE(st%surface_grid_center(1:2, 1:st%d%dim, 1:ix_max, 1:iy_max))
    SAFE_ALLOCATE(st%surface_grid_points_number(1:st%d%dim, 1:ix_max, 1:iy_max))
    SAFE_ALLOCATE(nn(1:2,1:3,1:3,1:3))

    st%surface_grid_center(1, 1, :, :) = -bounds(1,1)
    do iy=1, iy_max
      do iz=1, iz_max
        rr(2) = -bounds(1,2) + delta(2)/M_TWO + (iy-1) * delta(2)
        rr(3) = -bounds(1,3) + delta(3)/M_TWO + (iz-1) * delta(3)
        st%surface_grid_center(1, 2, iy, iz) = rr(2)
        st%surface_grid_center(1, 3, iy, iz) = rr(3)
      end do
    end do
    st%surface_grid_center(2, 1, :, :) = bounds(1,1)
    do iy=1, iy_max
      do iz=1, iz_max
        rr(2) = -bounds(1,2) + delta(2)/M_TWO + (iy-1) * delta(2)
        rr(3) = -bounds(1,3) + delta(3)/M_TWO + (iz-1) * delta(3)
        st%surface_grid_center(2, 2, iy, iz) = rr(2)
        st%surface_grid_center(2, 3, iy, iz) = rr(3)
      end do
    end do

    st%surface_grid_center(1, 2, :, :) = -bounds(1,2)
    do ix=1, ix_max
      do iz=1, iz_max
        rr(1) = -bounds(1,1) + delta(1)/M_TWO + (ix-1) * delta(1)
        rr(3) = -bounds(1,3) + delta(3)/M_TWO + (iz-1) * delta(3)
        st%surface_grid_center(1, 1, ix, iz) = rr(1)
        st%surface_grid_center(1, 3, ix, iz) = rr(3)
      end do
    end do
    st%surface_grid_center(2, 2, :, :) = bounds(1,2)
    do ix=1, ix_max
      do iz=1, iz_max
        rr(1) = -bounds(1,2) + delta(1)/M_TWO + (ix-1) * delta(1)
        rr(3) = -bounds(1,3) + delta(3)/M_TWO + (iz-1) * delta(3)
        st%surface_grid_center(2, 1, ix, iz) = rr(1)
        st%surface_grid_center(2, 3, ix, iz) = rr(3)
      end do
    end do

    st%surface_grid_center(1, 3, :, :) = -bounds(1,3)
    do ix=1, ix_max
      do iy=1, iy_max
        rr(1) = -bounds(1,1) + delta(1)/M_TWO + (ix-1) * delta(1)
        rr(2) = -bounds(1,2) + delta(2)/M_TWO + (iy-1) * delta(2)
        st%surface_grid_center(1, 1, ix, iy) = rr(1)
        st%surface_grid_center(1, 2, ix, iy) = rr(2)
      end do
    end do
    st%surface_grid_center(2, 3, :, :) = bounds(1,3)
    do ix=1, ix_max
      do iy=1, iy_max
        rr(1) = -bounds(1,2) + delta(1)/M_TWO + (ix-1) * delta(1)
        rr(2) = -bounds(1,2) + delta(2)/M_TWO + (iy-1) * delta(2)
        st%surface_grid_center(2, 1, ix, iy) = rr(1)
        st%surface_grid_center(2, 2, ix, iy) = rr(2)
      end do
    end do

    st%surface_grid_points_number(:,:,:) = 0

    nn_max = 0

    do iy=1, iy_max
      do iz=1, iz_max
        min_1(iy) = -bounds(1,2) + (iy-1) * delta(2)
        max_1(iy) = -bounds(1,2) + iy * delta(2)
        min_2(iz) = -bounds(1,3) + (iz-1) * delta(3)
        max_2(iz) = -bounds(1,3) + iz * delta(3)
      end do
    end do
    do iiy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
      do iiz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
        vec(1) = iiy * mesh%spacing(2)
        vec(2) = iiz * mesh%spacing(3)
        call get_surface_indices(vec, min_1, max_1, min_2, max_2, idx1, idx2)
        if ((idx1 /= 0) .and. (idx2 /= 0)) then
          st%surface_grid_points_number(1,idx1,idx2) = st%surface_grid_points_number(1,idx1,idx2)+1
          if (nn_max < st%surface_grid_points_number(1,idx1,idx2)) then
            nn_max = st%surface_grid_points_number(1,idx1,idx2)
          end if
        end if
      end do
    end do

    do ix=1, ix_max
      do iz=1, iz_max
        min_1(ix) = -bounds(1,1) + (ix-1) * delta(1)
        max_1(ix) = -bounds(1,1) + ix * delta(1)
        min_2(iz) = -bounds(1,3) + (iz-1) * delta(3)
        max_2(iz) = -bounds(1,3) + iz * delta(3)
      end do
    end do
    do iix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
      do iiz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
        vec(1)     = iix * mesh%spacing(1)
        vec(2)     = iiz * mesh%spacing(3)
        call get_surface_indices(vec, min_1, max_1, min_2, max_2, idx1, idx2)
        if ((idx1 /= 0) .and. (idx2 /= 0)) then
          st%surface_grid_points_number(2,idx1,idx2) = st%surface_grid_points_number(2,idx1,idx2)+1
          if (nn_max < st%surface_grid_points_number(2,idx1,idx2)) then
            nn_max = st%surface_grid_points_number(2,idx1,idx2)
          end if 
        end if
      end do
    end do

    do ix=1, ix_max
      do iy=1, iy_max
        min_1(ix) = -bounds(1,1) + (ix-1) * delta(1)
        max_1(ix) = -bounds(1,1) + ix * delta(1)
        min_2(iy) = -bounds(1,2) + (iy-1) * delta(2)
        max_2(iy) = -bounds(1,2) + iy * delta(2)
      end do
    end do
    do iix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
      do iiy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
        vec(1)     = iix * mesh%spacing(1)
        vec(2)     = iiy * mesh%spacing(2)
        call get_surface_indices(vec, min_1, max_1, min_2, max_2, idx1, idx2)
        if ((idx1 /= 0) .and. (idx2 /= 0)) then
          st%surface_grid_points_number(3,idx1,idx2) = st%surface_grid_points_number(3,idx1,idx2)+1
          if (nn_max < st%surface_grid_points_number(3,idx1,idx2)) then
            nn_max = st%surface_grid_points_number(3,idx1,idx2)
          end if
        end if
      end do
    end do

    SAFE_ALLOCATE(st%surface_grid_points_map(1:2,1:st%d%dim,1:iy_max,1:iz_max,1:nn_max))
    SAFE_ALLOCATE(st%surface_grid_points_map(1:2,1:st%d%dim,1:ix_max,1:iz_max,1:nn_max))
    SAFE_ALLOCATE(st%surface_grid_points_map(1:2,1:st%d%dim,1:ix_max,1:iy_max,1:nn_max))

    nn(:,:,:,:) = 0

    do iy=1, iy_max
      do iz=1, iz_max
        min_1(iy) = -bounds(1,2) + (iy-1) * delta(2)
        max_1(iy) = -bounds(1,2) + iy * delta(2)
        min_2(iz) = -bounds(1,3) + (iz-1) * delta(3)
        max_2(iz) = -bounds(1,3) + iz * delta(3)
      end do
    end do
    do iiy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
      do iiz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
        vec(1) = iiy * mesh%spacing(2)
        vec(2) = iiz * mesh%spacing(3)
        call get_surface_indices(vec, min_1, max_1, min_2, max_2, idx1, idx2)
        if ((idx1 /= 0) .and. (idx2 /= 0)) then
          nn(1, 1, idx1, idx2) = nn(1, 1, idx1, idx2) + 1
          rr(1) = -bounds(1,1)
          rr(2) = iiy * mesh%spacing(2)
          rr(3) = iiz * mesh%spacing(3)
          iix = int(-bounds(1,1)/mesh%spacing(1))
          ip_global = mesh%idx%lxyz_inv(iix,iiy,iiz)
          st%surface_grid_points_map(1, 1, idx1, idx2, nn(1, 1, idx1, idx2)) = ip_global
          nn(2, 1, idx1, idx2) = nn(2, 1, idx1, idx2) + 1
          rr(1) = bounds(1,1)
          rr(2) = iiy * mesh%spacing(2)
          rr(3) = iiz * mesh%spacing(3)
          iix = int(bounds(1,1)/mesh%spacing(1))
          ip_global = mesh%idx%lxyz_inv(iix,iiy,iiz)
          st%surface_grid_points_map(2, 1, idx1, idx2, nn(2, 1, idx1, idx2)) = ip_global
        end if
      end do
    end do

    do ix=1, ix_max
      do iz=1, iz_max
        min_1(ix) = -bounds(1,1) + (ix-1) * delta(1)
        max_1(ix) = -bounds(1,1) + ix * delta(1)
        min_2(iz) = -bounds(1,3) + (iz-1) * delta(3)
        max_2(iz) = -bounds(1,3) + iz * delta(3)
      end do
    end do
    do iix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
      do iiz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
        vec(1) = iix * mesh%spacing(1)
        vec(2) = iiz * mesh%spacing(3)
        call get_surface_indices(vec, min_1, max_1, min_2, max_2, idx1, idx2)
        if ((idx1 /= 0) .and. (idx2 /= 0)) then
          nn(1, 2, idx1, idx2) = nn(1, 2, idx1, idx2) + 1
          rr(1) = iix * mesh%spacing(1)
          rr(2) = -bounds(1,2)
          rr(3) = iiz * mesh%spacing(3)
          iiy = int(-bounds(1,2)/mesh%spacing(2))
          ip_global = mesh%idx%lxyz_inv(iix,iiy,iiz)
          st%surface_grid_points_map(1, 2, idx1, idx2, nn(1, 2, idx1, idx2)) = ip_global
          nn(2, 2, idx1, idx2) = nn(2, 2, idx1, idx2) + 1
          rr(1) = iix * mesh%spacing(1)
          rr(2) = bounds(1,2)
          rr(3) = iiz * mesh%spacing(3)
          iiy = int(bounds(1,2)/mesh%spacing(2))
          ip_global = mesh%idx%lxyz_inv(iix,iiy,iiz)
          st%surface_grid_points_map(2, 2, idx1, idx2, nn(2, 2, idx1, idx2)) = ip_global
        end if
      end do
    end do

    do ix=1, ix_max
      do iy=1, iy_max
        min_1(ix) = -bounds(1,1) + (ix-1) * delta(1)
        max_1(ix) = -bounds(1,1) + ix * delta(1)
        min_2(iy) = -bounds(1,2) + (iy-1) * delta(2)
        max_2(iy) = -bounds(1,2) + iy * delta(2)
      end do
    end do
    do iix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
      do iiy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
        vec(1) = iix * mesh%spacing(1)
        vec(2) = iiy * mesh%spacing(2)
        call get_surface_indices(vec, min_1, max_1, min_2, max_2, idx1, idx2)
        if ((idx1 /= 0) .and. (idx2 /= 0)) then
          nn(1, 3, idx1, idx2) = nn(1, 3, idx1, idx2) + 1
          rr(1) = iix * mesh%spacing(1)
          rr(2) = iiy * mesh%spacing(2)
          rr(3) = -bounds(1,3)
          iiz = int(-bounds(1,3)/mesh%spacing(3))
          ip_global = mesh%idx%lxyz_inv(iix,iiy,iiz)
          st%surface_grid_points_map(1, 3, idx1, idx2, nn(1, 3, idx1, idx2)) = ip_global
          nn(2, 3, idx1, idx2) = nn(2, 3, idx1, idx2) + 1
          rr(1) = iix * mesh%spacing(1)
          rr(2) = iiy * mesh%spacing(2)
          rr(3) = bounds(1,3)
          iiz = int(bounds(1,3)/mesh%spacing(3))
          ip_global = mesh%idx%lxyz_inv(iix,iiy,iiz)
          st%surface_grid_points_map(2, 3, idx1, idx2, nn(2, 3, idx1, idx2)) = ip_global
        end if
      end do
    end do

    SAFE_DEALLOCATE_A(nn)

    POP_SUB(surface_grid_points_mapping)

    contains

      subroutine get_surface_indices(vec, min_1, max_1, min_2, max_2, index_1, index_2)
        FLOAT,   intent(in)  :: vec(:)
        FLOAT,   intent(in)  :: min_1(:)
        FLOAT,   intent(in)  :: max_1(:)
        FLOAT,   intent(in)  :: min_2(:)
        FLOAT,   intent(in)  :: max_2(:)
        integer, intent(out) :: index_1
        integer, intent(out) :: index_2

        if ( ((vec(1) >= min_1(1)) .and. (vec(1) <= max_1(1))) .and. ((vec(2) >= min_2(1)) .and. (vec(2) <= max_2(1))) ) then
          index_1 = 1
          index_2 = 1
        else if ( ((vec(1) >= min_1(2)) .and. (vec(1) <= max_1(2))) .and. ((vec(2) >= min_2(1)) .and. (vec(2) <= max_2(1))) ) then
          index_1 = 2
          index_2 = 1
        else if ( ((vec(1) >= min_1(3)) .and. (vec(1) <= max_1(3))) .and. ((vec(2) >= min_2(1)) .and. (vec(2) <= max_2(1))) ) then
          index_1 = 3
          index_2 = 1
        else if ( ((vec(1) >= min_1(1)) .and. (vec(1) <= max_1(1))) .and. ((vec(2) >= min_2(2)) .and. (vec(2) <= max_2(2))) ) then
          index_1 = 1
          index_2 = 2
        else if ( ((vec(1) >= min_1(2)) .and. (vec(1) <= max_1(2))) .and. ((vec(2) >= min_2(2)) .and. (vec(2) <= max_2(2))) ) then
          index_1 = 2
          index_2 = 2
        else if ( ((vec(1) >= min_1(3)) .and. (vec(1) <= max_1(3))) .and. ((vec(2) >= min_2(2)) .and. (vec(2) <= max_2(2))) ) then
          index_1 = 3
          index_2 = 2
        else if ( ((vec(1) >= min_1(1)) .and. (vec(1) <= max_1(1))) .and. ((vec(2) >= min_2(3)) .and. (vec(2) <= max_2(3))) ) then
          index_1 = 1
          index_2 = 3
        else if ( ((vec(1) >= min_1(2)) .and. (vec(1) <= max_1(2))) .and. ((vec(2) >= min_2(3)) .and. (vec(2) <= max_2(3))) ) then
          index_1 = 2
          index_2 = 3
        else if ( ((vec(1) >= min_1(3)) .and. (vec(1) <= max_1(3))) .and. ((vec(2) >= min_2(3)) .and. (vec(2) <= max_2(3))) ) then
          index_1 = 3
          index_2 = 3
        else
          index_1 = 0
          index_2 = 0
        end if

      end subroutine get_surface_indices

  end subroutine surface_grid_points_mapping


  !----------------------------------------------------------
  subroutine energy_density_calc(gr, st, rs_field, energy_dens, &
                                         e_energy_dens, b_energy_dens, plane_waves_check, &
                                         rs_field_plane_waves, energy_dens_plane_waves)
    type(grid_t),      intent(in)    :: gr
    type(states_mxll_t),    intent(in)    :: st
    CMPLX,             intent(in)    :: rs_field(:,:)
    FLOAT,             intent(inout) :: energy_dens(:)
    FLOAT,             intent(inout) :: e_energy_dens(:)
    FLOAT,             intent(inout) :: b_energy_dens(:)
    logical, optional, intent(in)    :: plane_waves_check
    CMPLX,   optional, intent(in)    :: rs_field_plane_waves(:,:)
    FLOAT,   optional, intent(inout) :: energy_dens_plane_waves(:)

    CMPLX, allocatable :: ztmp(:,:)
    integer            :: idim, ip

    PUSH_SUB(energy_density_calc)

    SAFE_ALLOCATE(ztmp(1:gr%mesh%np_part,1:st%d%dim))

    ztmp(:,:) = rs_field(:,:)

    energy_dens(:) = M_ZERO
    do ip=1, gr%mesh%np
      do idim=1, st%d%dim
        energy_dens(ip) = energy_dens(ip) + conjg(ztmp(ip,idim)) * ztmp(ip,idim)
      end do
    end do

    e_energy_dens(:) = M_ZERO
    do ip=1, gr%mesh%np
      do idim=1, st%d%dim
        e_energy_dens(ip) = e_energy_dens(ip) + real(ztmp(ip,idim))**2
      end do
    end do

    b_energy_dens(:) = M_ZERO
    do ip=1, gr%mesh%np
      do idim=1, st%d%dim
        b_energy_dens(ip) = b_energy_dens(ip) + aimag(ztmp(ip,idim))**2
      end do
    end do

    if (present(rs_field_plane_waves) .and. present(energy_dens_plane_waves) .and. plane_waves_check) then
      ztmp(:,:) = rs_field_plane_waves(:,:)
      energy_dens_plane_waves(:) = M_ZERO
      do ip=1, gr%mesh%np
        do idim=1, st%d%dim
          energy_dens_plane_waves(ip) = energy_dens_plane_waves(ip) + conjg(ztmp(ip,idim)) * ztmp(ip,idim)
        end do
      end do
    end if

    SAFE_DEALLOCATE_A(ztmp)

    POP_SUB(energy_density_calc)
  end subroutine energy_density_calc


  !----------------------------------------------------------
  subroutine energy_mxll_calc(gr, st, hm, rs_field, mx_energy, mx_e_energy, mx_b_energy, &
                                 mx_energy_boundary, rs_field_plane_waves, mx_energy_plane_waves) 
    type(grid_t),        intent(in)  :: gr
    type(states_mxll_t),      intent(in)  :: st
    type(hamiltonian_mxll_t), intent(in)  :: hm
    CMPLX,               intent(in)  :: rs_field(:,:)
    FLOAT,               intent(out) :: mx_energy
    FLOAT,               intent(out) :: mx_e_energy
    FLOAT,               intent(out) :: mx_b_energy
    FLOAT, optional,     intent(out) :: mx_energy_boundary
    CMPLX, optional,     intent(in)  :: rs_field_plane_waves(:,:)
    FLOAT, optional,     intent(out) :: mx_energy_plane_waves

    integer            :: ip, ip_in, ip_bd, idim, dim
    FLOAT              :: dd_energy, dd_e_energy, dd_b_energy
    FLOAT, allocatable :: energy_density(:), energy_density_plane_waves(:), tmp(:), tmp_pw(:)
    FLOAT, allocatable :: e_energy_density(:), tmp_e(:)
    FLOAT, allocatable :: b_energy_density(:), tmp_b(:)
    PUSH_SUB(energy_mxll_calc)

    dim = st%d%dim

    SAFE_ALLOCATE(energy_density(1:gr%mesh%np))
    SAFE_ALLOCATE(energy_density_plane_waves(1:gr%mesh%np))
    SAFE_ALLOCATE(e_energy_density(1:gr%mesh%np))
    SAFE_ALLOCATE(b_energy_density(1:gr%mesh%np))
    SAFE_ALLOCATE(tmp(1:gr%mesh%np))
    SAFE_ALLOCATE(tmp_e(1:gr%mesh%np))
    SAFE_ALLOCATE(tmp_b(1:gr%mesh%np))
    SAFE_ALLOCATE(tmp_pw(1:gr%mesh%np))


    call energy_density_calc(gr, st, rs_field, energy_density, & 
                                     e_energy_density, b_energy_density,       &
                                     hm%plane_waves, &
                                     rs_field_plane_waves, energy_density_plane_waves)

    tmp    = M_ZERO
    tmp_e  = M_ZERO
    tmp_b  = M_ZERO
    tmp_pw = M_ZERO

    do ip_in = 1, st%inner_points_number
      ip =  st%inner_points_map(ip_in)
      tmp(ip)    = energy_density(ip)
      tmp_e(ip)  = e_energy_density(ip)
      tmp_b(ip)  = b_energy_density(ip)
      if (present(rs_field_plane_waves) .and. present(mx_energy_plane_waves))  &
         tmp_pw(ip) = energy_density_plane_waves(ip)
    end do

    mx_energy             = dmf_integrate(gr%mesh, tmp)
    mx_e_energy           = dmf_integrate(gr%mesh, tmp_e)
    mx_b_energy           = dmf_integrate(gr%mesh, tmp_b)
    if (present(rs_field_plane_waves) .and. present(mx_energy_plane_waves))  &
       mx_energy_plane_waves = dmf_integrate(gr%mesh, tmp_pw)

    tmp   = M_ZERO
    tmp_e = M_ZERO
    tmp_b = M_ZERO

    if (present(mx_energy_boundary)) then
      do ip_in = 1, st%boundary_points_number
        ip = st%boundary_points_map(ip_in)
        tmp(ip)   = energy_density(ip)
      end do
      mx_energy_boundary = dmf_integrate(gr%mesh, tmp)
    end if

    SAFE_DEALLOCATE_A(energy_density)
    SAFE_DEALLOCATE_A(energy_density_plane_waves)
    SAFE_DEALLOCATE_A(e_energy_density)
    SAFE_DEALLOCATE_A(b_energy_density)
    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp_e)
    SAFE_DEALLOCATE_A(tmp_b)
    SAFE_DEALLOCATE_A(tmp_pw)

    POP_SUB(energy_mxll_calc)
  end subroutine energy_mxll_calc


  ! ---------------------------------------------------------
  subroutine energy_rate_calc(gr, st, iter, dt, energy_rate, &
                                      delta_energy, energy_via_flux_calc, energy_via_flux_calc_dir)
    type(grid_t),        intent(in)    :: gr
    type(states_mxll_t),      intent(in)    :: st
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt
    FLOAT,               intent(out)   :: energy_rate(:)
    FLOAT,               intent(out)   :: delta_energy(:)
    FLOAT,               intent(out)   :: energy_via_flux_calc(:)
    FLOAT,  optional,    intent(out)   :: energy_via_flux_calc_dir(:,:)

    integer             :: idim, idir, ip_surf, ip_global, ix, ix_max, iy, iy_max, iz, iz_max, ii_max
    FLOAT               :: tmp_sum
    FLOAT,  allocatable :: poynting_vector(:,:), tmp_global(:,:), tmp_surf(:,:,:,:,:)

    PUSH_SUB(energy_rate_calc)

    SAFE_ALLOCATE(poynting_vector(1:gr%mesh%np,1:st%d%dim))
    SAFE_ALLOCATE(tmp_global(1:gr%mesh%np_global,1:st%d%dim))

    call get_poynting_vector(gr, st, st%rs_state, st%rs_sign, &
                                     poynting_vector, ep_field=st%ep, mu_field=st%mu)

    if (gr%mesh%parallel_in_domains) then
      do idim=1, st%d%dim
#if defined(HAVE_MPI)
        call vec_allgather(gr%mesh%vp, tmp_global(:,idim), poynting_vector(:,idim))
        call MPI_Barrier(gr%mesh%vp%comm, mpi_err)
#endif
      end do
    else
      tmp_global(:,:) = poynting_vector(:,:)
    end if

    ix_max = st%surface_grid_rows_number(1)
    iy_max = st%surface_grid_rows_number(2)
    iz_max = st%surface_grid_rows_number(3)
    ii_max = max(ix_max,iy_max,iz_max)

    SAFE_ALLOCATE(tmp_surf(1:2,1:st%d%dim,1:ii_max,1:ii_max,1:st%d%dim))

    tmp_surf = M_ZERO
    tmp_sum  = M_ZERO

    do iy=1, iy_max
      do iz=1, iz_max
        do ip_surf=1, st%surface_grid_points_number(1,iy,iz)
          tmp_surf(1,1,iy,iz,:) = tmp_surf(1,1,iy,iz,:) & 
                                + tmp_global(st%surface_grid_points_map(1,1,iy,iz,ip_surf),:)
          tmp_surf(2,1,iy,iz,:) = tmp_surf(2,1,iy,iz,:) &
                                + tmp_global(st%surface_grid_points_map(2,1,iy,iz,ip_surf),:)
        end do
        tmp_surf(1,1,iy,iz,:) = tmp_surf(1,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
        tmp_surf(2,1,iy,iz,:) = tmp_surf(2,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
      end do
    end do
    do iy=1, iy_max
      do iz=1, iz_max
        tmp_sum = tmp_sum - tmp_surf(1,1,iy,iz,1) * st%surface_grid_element(1)
        tmp_sum = tmp_sum + tmp_surf(2,1,iy,iz,1) * st%surface_grid_element(1)
      end do
    end do

    do ix=1, ix_max
      do iz=1, iz_max
        do ip_surf=1, st%surface_grid_points_number(2,ix,iz)
          tmp_surf(1,2,ix,iz,:) = tmp_surf(1,2,ix,iz,:) &
                                + tmp_global(st%surface_grid_points_map(1,2,ix,iz,ip_surf),:)
          tmp_surf(2,2,ix,iz,:) = tmp_surf(2,2,ix,iz,:) &
                                + tmp_global(st%surface_grid_points_map(2,2,ix,iz,ip_surf),:)
        end do
        tmp_surf(1,2,ix,iz,:) = tmp_surf(1,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
        tmp_surf(2,2,ix,iz,:) = tmp_surf(2,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
      end do
    end do
    do ix=1, ix_max
      do iz=1, iz_max
        tmp_sum = tmp_sum - tmp_surf(1,2,ix,iz,2) * st%surface_grid_element(2)
        tmp_sum = tmp_sum + tmp_surf(2,2,ix,iz,2) * st%surface_grid_element(2)
      end do
    end do

    do ix=1, ix_max
      do iy=1, iy_max
        do ip_surf=1, st%surface_grid_points_number(3,ix,iy)
          tmp_surf(1,3,ix,iy,:) = tmp_surf(1,3,ix,iy,:) &
                                + tmp_global(st%surface_grid_points_map(1,3,ix,iy,ip_surf),:)
          tmp_surf(2,3,ix,iy,:) = tmp_surf(2,3,ix,iy,:) &
                                + tmp_global(st%surface_grid_points_map(2,3,ix,iy,ip_surf),:)
        end do
        tmp_surf(1,3,ix,iy,:) = tmp_surf(1,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
        tmp_surf(2,3,ix,iy,:) = tmp_surf(2,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
      end do
    end do
    do ix=1, ix_max
      do iy=1, iy_max
        tmp_sum = tmp_sum - tmp_surf(1,3,ix,iy,3) * st%surface_grid_element(3)
        tmp_sum = tmp_sum + tmp_surf(2,3,ix,iy,3) * st%surface_grid_element(3)
      end do
    end do

    energy_rate(iter)          = - tmp_sum
    delta_energy(iter)         = energy_rate(iter) * dt
    if (iter > 1) then
      energy_via_flux_calc(iter) = energy_via_flux_calc(iter-1) + delta_energy(iter)
    else if (iter == 1) then
      energy_via_flux_calc(iter) = delta_energy(iter)
    else
      energy_via_flux_calc(iter) = M_ZERO
    end if

    SAFE_DEALLOCATE_A(poynting_vector)
    SAFE_DEALLOCATE_A(tmp_global)

    POP_SUB(energy_rate_calc)
  end subroutine energy_rate_calc


  ! ---------------------------------------------------------
  subroutine energy_rate_calc_plane_waves(gr, st, iter, dt, energy_rate, &
                                                  delta_energy, energy_via_flux_calc, energy_via_flux_calc_dir)
    type(grid_t),        intent(in)    :: gr
    type(states_mxll_t),      intent(in)    :: st
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt
    FLOAT,               intent(out)   :: energy_rate(:)
    FLOAT,               intent(out)   :: delta_energy(:)
    FLOAT,               intent(out)   :: energy_via_flux_calc(:)
    FLOAT,  optional,    intent(out)   :: energy_via_flux_calc_dir(:,:)

    integer             :: idim, idir, ip_surf, ip_global, ix, ix_max, iy, iy_max, iz, iz_max, ii_max
    FLOAT               :: tmp_sum
    FLOAT,  allocatable :: poynting_vector(:,:), tmp_global(:,:), tmp_surf(:,:,:,:,:)

    PUSH_SUB(energy_rate_calc_plane_waves)

    SAFE_ALLOCATE(poynting_vector(1:gr%mesh%np,1:st%d%dim))
    SAFE_ALLOCATE(tmp_global(1:gr%mesh%np_global,1:st%d%dim))

    call get_poynting_vector_plane_waves(gr, st, st%rs_sign, poynting_vector)

    if (gr%mesh%parallel_in_domains) then
      do idim=1, st%d%dim
#if defined(HAVE_MPI)
        call vec_allgather(gr%mesh%vp, tmp_global(:,idim), poynting_vector(:,idim))
        call MPI_Barrier(gr%mesh%vp%comm, mpi_err)
#endif
      end do
    else
      tmp_global(:,:) = poynting_vector(:,:)
    end if

    ix_max = st%surface_grid_rows_number(1)
    iy_max = st%surface_grid_rows_number(2)
    iz_max = st%surface_grid_rows_number(3)
    ii_max = max(ix_max,iy_max,iz_max)

    SAFE_ALLOCATE(tmp_surf(1:2,1:st%d%dim,1:ii_max,1:ii_max,1:st%d%dim))

    tmp_surf = M_ZERO
    tmp_sum  = M_ZERO

    do iy=1, iy_max
      do iz=1, iz_max
        do ip_surf=1, st%surface_grid_points_number(1,iy,iz)
          tmp_surf(1,1,iy,iz,:) = tmp_surf(1,1,iy,iz,:) & 
                                + tmp_global(st%surface_grid_points_map(1,1,iy,iz,ip_surf),:)
          tmp_surf(2,1,iy,iz,:) = tmp_surf(2,1,iy,iz,:) &
                                + tmp_global(st%surface_grid_points_map(2,1,iy,iz,ip_surf),:)
        end do
        tmp_surf(1,1,iy,iz,:) = tmp_surf(1,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
        tmp_surf(2,1,iy,iz,:) = tmp_surf(2,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
      end do
    end do
    do iy=1, iy_max
      do iz=1, iz_max
        tmp_sum = tmp_sum - tmp_surf(1,1,iy,iz,1) * st%surface_grid_element(1)
        tmp_sum = tmp_sum + tmp_surf(2,1,iy,iz,1) * st%surface_grid_element(1)
      end do
    end do

    do ix=1, ix_max
      do iz=1, iz_max
        do ip_surf=1, st%surface_grid_points_number(2,ix,iz)
          tmp_surf(1,2,ix,iz,:) = tmp_surf(1,2,ix,iz,:) &
                                + tmp_global(st%surface_grid_points_map(1,2,ix,iz,ip_surf),:)
          tmp_surf(2,2,ix,iz,:) = tmp_surf(2,2,ix,iz,:) &
                                + tmp_global(st%surface_grid_points_map(2,2,ix,iz,ip_surf),:)
        end do
        tmp_surf(1,2,ix,iz,:) = tmp_surf(1,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
        tmp_surf(2,2,ix,iz,:) = tmp_surf(2,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
      end do
    end do
    do ix=1, ix_max
      do iz=1, iz_max
        tmp_sum = tmp_sum - tmp_surf(1,2,ix,iz,2) * st%surface_grid_element(2)
        tmp_sum = tmp_sum + tmp_surf(2,2,ix,iz,2) * st%surface_grid_element(2)
      end do
    end do

    do ix=1, ix_max
      do iy=1, iy_max
        do ip_surf=1, st%surface_grid_points_number(3,ix,iy)
          tmp_surf(1,3,ix,iy,:) = tmp_surf(1,3,ix,iy,:) &
                                + tmp_global(st%surface_grid_points_map(1,3,ix,iy,ip_surf),:)
          tmp_surf(2,3,ix,iy,:) = tmp_surf(2,3,ix,iy,:) &
                                + tmp_global(st%surface_grid_points_map(2,3,ix,iy,ip_surf),:)
        end do
        tmp_surf(1,3,ix,iy,:) = tmp_surf(1,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
        tmp_surf(2,3,ix,iy,:) = tmp_surf(2,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
      end do
    end do
    do ix=1, ix_max
      do iy=1, iy_max
        tmp_sum = tmp_sum - tmp_surf(1,3,ix,iy,3) * st%surface_grid_element(3)
        tmp_sum = tmp_sum + tmp_surf(2,3,ix,iy,3) * st%surface_grid_element(3)
      end do
    end do

    energy_rate(iter)          = - tmp_sum
    delta_energy(iter)         = energy_rate(iter) * dt
    if (iter > 1) then
      energy_via_flux_calc(iter) = energy_via_flux_calc(iter-1) + delta_energy(iter)
    else if (iter == 1) then
      energy_via_flux_calc(iter) = delta_energy(iter)
    else
      energy_via_flux_calc(iter) = M_ZERO
    end if

    SAFE_DEALLOCATE_A(poynting_vector)
    SAFE_DEALLOCATE_A(tmp_global)

    POP_SUB(energy_rate_calc_plane_waves)
  end subroutine energy_rate_calc_plane_waves


  ! ---------------------------------------------------------
  subroutine poynting_vector_through_box_surfaces(gr, st, iter, dt, poynting_box_surface)
    type(grid_t),        intent(in)    :: gr
    type(states_mxll_t),      intent(in)    :: st
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt
    FLOAT,               intent(out)   :: poynting_box_surface(:,:,:)

    integer             :: idim, idir, ip_surf, ip_global, ix, ix_max, iy, iy_max, iz, iz_max, ii_max
    FLOAT,  allocatable :: poynting_vector(:,:), tmp_global(:,:), tmp_surf(:,:,:,:,:)

    PUSH_SUB(poynting_vector_through_box_surfaces)

    SAFE_ALLOCATE(poynting_vector(1:gr%mesh%np,1:st%d%dim))
    SAFE_ALLOCATE(tmp_global(1:gr%mesh%np_global,1:st%d%dim))

    call get_poynting_vector(gr, st, st%rs_state, st%rs_sign, &
                                     poynting_vector, ep_field=st%ep, mu_field=st%mu)

    if (gr%mesh%parallel_in_domains) then
      do idim=1, st%d%dim
#if defined(HAVE_MPI)
        call vec_allgather(gr%mesh%vp, tmp_global(:,idim), poynting_vector(:,idim))
        call MPI_Barrier(gr%mesh%vp%comm, mpi_err)
#endif
      end do
    else
      tmp_global(:,:) = poynting_vector(:,:)
    end if

    ix_max = st%surface_grid_rows_number(1)
    iy_max = st%surface_grid_rows_number(2)
    iz_max = st%surface_grid_rows_number(3)
    ii_max = max(ix_max,iy_max,iz_max)

    SAFE_ALLOCATE(tmp_surf(1:2,1:st%d%dim,1:ii_max,1:ii_max,1:st%d%dim))

    tmp_surf = M_ZERO
    poynting_box_surface = M_ZERO

    do iy=1, iy_max
      do iz=1, iz_max
        do ip_surf=1, st%surface_grid_points_number(1,iy,iz)
          tmp_surf(1,1,iy,iz,:) = tmp_surf(1,1,iy,iz,:) & 
                                + tmp_global(st%surface_grid_points_map(1,1,iy,iz,ip_surf),:)
          tmp_surf(2,1,iy,iz,:) = tmp_surf(2,1,iy,iz,:) &
                                + tmp_global(st%surface_grid_points_map(2,1,iy,iz,ip_surf),:)
        end do
        tmp_surf(1,1,iy,iz,:) = tmp_surf(1,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
        tmp_surf(2,1,iy,iz,:) = tmp_surf(2,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
      end do
    end do
    do iy=1, iy_max
      do iz=1, iz_max
        poynting_box_surface(1,1,:) = poynting_box_surface(1,1,:) + tmp_surf(1,1,iy,iz,:) * st%surface_grid_element(1)
        poynting_box_surface(2,1,:) = poynting_box_surface(2,1,:) + tmp_surf(2,1,iy,iz,:) * st%surface_grid_element(1)
      end do
    end do

    do ix=1, ix_max
      do iz=1, iz_max
        do ip_surf=1, st%surface_grid_points_number(2,ix,iz)
          tmp_surf(1,2,ix,iz,:) = tmp_surf(1,2,ix,iz,:) &
                                + tmp_global(st%surface_grid_points_map(1,2,ix,iz,ip_surf),:)
          tmp_surf(2,2,ix,iz,:) = tmp_surf(2,2,ix,iz,:) &
                                + tmp_global(st%surface_grid_points_map(2,2,ix,iz,ip_surf),:)
        end do
        tmp_surf(1,2,ix,iz,:) = tmp_surf(1,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
        tmp_surf(2,2,ix,iz,:) = tmp_surf(2,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
      end do
    end do
    do ix=1, ix_max
      do iz=1, iz_max
        poynting_box_surface(1,2,:) = poynting_box_surface(1,2,:) + tmp_surf(1,2,ix,iz,:) * st%surface_grid_element(2)
        poynting_box_surface(2,2,:) = poynting_box_surface(2,2,:) + tmp_surf(2,2,ix,iz,:) * st%surface_grid_element(2)
      end do
    end do

    do ix=1, ix_max
      do iy=1, iy_max
        do ip_surf=1, st%surface_grid_points_number(3,ix,iy)
          tmp_surf(1,3,ix,iy,:) = tmp_surf(1,3,ix,iy,:) &
                                + tmp_global(st%surface_grid_points_map(1,3,ix,iy,ip_surf),:)
          tmp_surf(2,3,ix,iy,:) = tmp_surf(2,3,ix,iy,:) &
                                + tmp_global(st%surface_grid_points_map(2,3,ix,iy,ip_surf),:)
        end do
        tmp_surf(1,3,ix,iy,:) = tmp_surf(1,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
        tmp_surf(2,3,ix,iy,:) = tmp_surf(2,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
      end do
    end do
    do ix=1, ix_max
      do iy=1, iy_max
        poynting_box_surface(1,3,:) = poynting_box_surface(1,3,:) - tmp_surf(1,3,ix,iy,:) * st%surface_grid_element(3)
        poynting_box_surface(2,3,:) = poynting_box_surface(2,3,:) + tmp_surf(2,3,ix,iy,:) * st%surface_grid_element(3)
      end do
    end do

    SAFE_DEALLOCATE_A(poynting_vector)
    SAFE_DEALLOCATE_A(tmp_global)

    POP_SUB(poynting_vector_through_box_surfaces)
  end subroutine poynting_vector_through_box_surfaces


  ! ---------------------------------------------------------
  subroutine poynting_vector_through_box_surfaces_plane_waves(gr, st, iter, dt, poynting_box_surface)
    type(grid_t),        intent(in)    :: gr
    type(states_mxll_t),      intent(in)    :: st
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt
    FLOAT,               intent(out)   :: poynting_box_surface(:,:,:)

    integer             :: idim, idir, ip_surf, ip_global, ix, ix_max, iy, iy_max, iz, iz_max, ii_max
    FLOAT,  allocatable :: poynting_vector(:,:), tmp_global(:,:), tmp_surf(:,:,:,:,:)

    PUSH_SUB(poynting_vector_through_box_surfaces_plane_waves)

    SAFE_ALLOCATE(poynting_vector(1:gr%mesh%np,1:st%d%dim))
    SAFE_ALLOCATE(tmp_global(1:gr%mesh%np_global,1:st%d%dim))

    call get_poynting_vector_plane_waves(gr, st,  st%rs_sign, poynting_vector)

    if (gr%mesh%parallel_in_domains) then
      do idim=1, st%d%dim
#if defined(HAVE_MPI)
        call vec_allgather(gr%mesh%vp, tmp_global(:,idim), poynting_vector(:,idim))
        call MPI_Barrier(gr%mesh%vp%comm, mpi_err)
#endif
      end do
    else
      tmp_global(:,:) = poynting_vector(:,:)
    end if

    ix_max = st%surface_grid_rows_number(1)
    iy_max = st%surface_grid_rows_number(2)
    iz_max = st%surface_grid_rows_number(3)
    ii_max = max(ix_max,iy_max,iz_max)

    SAFE_ALLOCATE(tmp_surf(1:2,1:st%d%dim,1:ii_max,1:ii_max,1:st%d%dim))

    tmp_surf = M_ZERO
    poynting_box_surface = M_ZERO

    do iy=1, iy_max
      do iz=1, iz_max
        do ip_surf=1, st%surface_grid_points_number(1,iy,iz)
          tmp_surf(1,1,iy,iz,:) = tmp_surf(1,1,iy,iz,:) & 
                                + tmp_global(st%surface_grid_points_map(1,1,iy,iz,ip_surf),:)
          tmp_surf(2,1,iy,iz,:) = tmp_surf(2,1,iy,iz,:) &
                                + tmp_global(st%surface_grid_points_map(2,1,iy,iz,ip_surf),:)
        end do
        tmp_surf(1,1,iy,iz,:) = tmp_surf(1,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
        tmp_surf(2,1,iy,iz,:) = tmp_surf(2,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
      end do
    end do
    do iy=1, iy_max
      do iz=1, iz_max
        poynting_box_surface(1,1,:) = poynting_box_surface(1,1,:) + tmp_surf(1,1,iy,iz,:) * st%surface_grid_element(1)
        poynting_box_surface(2,1,:) = poynting_box_surface(2,1,:) + tmp_surf(2,1,iy,iz,:) * st%surface_grid_element(1)
      end do
    end do

    do ix=1, ix_max
      do iz=1, iz_max
        do ip_surf=1, st%surface_grid_points_number(2,ix,iz)
          tmp_surf(1,2,ix,iz,:) = tmp_surf(1,2,ix,iz,:) &
                                + tmp_global(st%surface_grid_points_map(1,2,ix,iz,ip_surf),:)
          tmp_surf(2,2,ix,iz,:) = tmp_surf(2,2,ix,iz,:) &
                                + tmp_global(st%surface_grid_points_map(2,2,ix,iz,ip_surf),:)
        end do
        tmp_surf(1,2,ix,iz,:) = tmp_surf(1,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
        tmp_surf(2,2,ix,iz,:) = tmp_surf(2,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
      end do
    end do
    do ix=1, ix_max
      do iz=1, iz_max
        poynting_box_surface(1,2,:) = poynting_box_surface(1,2,:) + tmp_surf(1,2,ix,iz,:) * st%surface_grid_element(2)
        poynting_box_surface(2,2,:) = poynting_box_surface(2,2,:) + tmp_surf(2,2,ix,iz,:) * st%surface_grid_element(2)
      end do
    end do

    do ix=1, ix_max
      do iy=1, iy_max
        do ip_surf=1, st%surface_grid_points_number(3,ix,iy)
          tmp_surf(1,3,ix,iy,:) = tmp_surf(1,3,ix,iy,:) &
                                + tmp_global(st%surface_grid_points_map(1,3,ix,iy,ip_surf),:)
          tmp_surf(2,3,ix,iy,:) = tmp_surf(2,3,ix,iy,:) &
                                + tmp_global(st%surface_grid_points_map(2,3,ix,iy,ip_surf),:)
        end do
        tmp_surf(1,3,ix,iy,:) = tmp_surf(1,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
        tmp_surf(2,3,ix,iy,:) = tmp_surf(2,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
      end do
    end do
    do ix=1, ix_max
      do iy=1, iy_max
        poynting_box_surface(1,3,:) = poynting_box_surface(1,3,:) - tmp_surf(1,3,ix,iy,:) * st%surface_grid_element(3)
        poynting_box_surface(2,3,:) = poynting_box_surface(2,3,:) + tmp_surf(2,3,ix,iy,:) * st%surface_grid_element(3)
      end do
    end do

    SAFE_DEALLOCATE_A(poynting_vector)
    SAFE_DEALLOCATE_A(tmp_global)

    POP_SUB(poynting_vector_through_box_surfaces_plane_waves)
  end subroutine poynting_vector_through_box_surfaces_plane_waves


  ! ---------------------------------------------------------
  subroutine fields_through_box_surfaces(gr, st, hm, iter, dt, e_field_box_surface, b_field_box_surface)
    type(grid_t),        intent(in)    :: gr
    type(states_mxll_t),      intent(in)    :: st
    type(hamiltonian_mxll_t), intent(in)    :: hm
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt
    FLOAT,               intent(out)   :: e_field_box_surface(:,:,:)
    FLOAT,               intent(out)   :: b_field_box_surface(:,:,:)

    integer             :: idim, idir, ip_surf, ip_global, ix, ix_max, iy, iy_max, iz, iz_max, ii_max, ip_in, ip
    FLOAT,  allocatable :: e_surf(:,:,:,:,:), b_surf(:,:,:,:,:), e_field(:,:), &
      e_field_global(:,:), b_field(:,:), b_field_global(:,:)

    PUSH_SUB(fields_through_box_surfaces)

    SAFE_ALLOCATE(e_field(1:gr%mesh%np,1:st%d%dim))
    SAFE_ALLOCATE(b_field(1:gr%mesh%np,1:st%d%dim))
    SAFE_ALLOCATE(e_field_global(1:gr%mesh%np_global,1:st%d%dim))
    SAFE_ALLOCATE(b_field_global(1:gr%mesh%np_global,1:st%d%dim))

    if (.not. hm%ma_mx_coupling) then
      call get_electric_field_state(st%rs_state_trans+st%rs_state_long, gr%mesh, e_field, st%ep, gr%mesh%np)
      call get_magnetic_field_state(st%rs_state_trans+st%rs_state_long, gr%mesh, st%rs_sign, b_field, &
        st%mu, gr%mesh%np)
    else
      call get_electric_field_state(st%rs_state, gr%mesh, e_field, st%ep, gr%mesh%np)
      call get_magnetic_field_state(st%rs_state, gr%mesh, st%rs_sign, b_field, st%mu, gr%mesh%np)
    end if

    if (gr%mesh%parallel_in_domains) then
      do idim=1, st%d%dim
#if defined(HAVE_MPI)
        call vec_allgather(gr%mesh%vp, e_field_global(:,idim), e_field(:,idim))
        call MPI_Barrier(gr%mesh%vp%comm, mpi_err)
        call vec_allgather(gr%mesh%vp, b_field_global(:,idim), b_field(:,idim))
        call MPI_Barrier(gr%mesh%vp%comm, mpi_err)
#endif
      end do
    else
      e_field_global(:,:) = e_field(:,:)
      b_field_global(:,:) = b_field(:,:)
    end if

    ix_max = st%surface_grid_rows_number(1)
    iy_max = st%surface_grid_rows_number(2)
    iz_max = st%surface_grid_rows_number(3)
    ii_max = max(ix_max,iy_max,iz_max)

    SAFE_ALLOCATE(e_surf(1:2,1:st%d%dim,1:ii_max,1:ii_max,1:st%d%dim))
    SAFE_ALLOCATE(b_surf(1:2,1:st%d%dim,1:ii_max,1:ii_max,1:st%d%dim))

    e_surf = M_ZERO
    b_surf = M_ZERO
    e_field_box_surface = M_ZERO
    b_field_box_surface = M_ZERO

    do iy=1, iy_max
      do iz=1, iz_max
        do ip_surf=1, st%surface_grid_points_number(1,iy,iz)
          e_surf(1,1,iy,iz,:) = e_surf(1,1,iy,iz,:) & 
                                + e_field_global(st%surface_grid_points_map(1,1,iy,iz,ip_surf),:)
          e_surf(2,1,iy,iz,:) = e_surf(2,1,iy,iz,:) &
                                + e_field_global(st%surface_grid_points_map(2,1,iy,iz,ip_surf),:)
          b_surf(1,1,iy,iz,:) = b_surf(1,1,iy,iz,:) & 
                                + b_field_global(st%surface_grid_points_map(1,1,iy,iz,ip_surf),:)
          b_surf(2,1,iy,iz,:) = b_surf(2,1,iy,iz,:) &
                                + b_field_global(st%surface_grid_points_map(2,1,iy,iz,ip_surf),:)
        end do
        e_surf(1,1,iy,iz,:) = e_surf(1,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
        e_surf(2,1,iy,iz,:) = e_surf(2,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
        b_surf(1,1,iy,iz,:) = b_surf(1,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
        b_surf(2,1,iy,iz,:) = b_surf(2,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
      end do
    end do
    do iy=1, iy_max
      do iz=1, iz_max
        e_field_box_surface(1,1,:) = e_field_box_surface(1,1,:) + e_surf(1,1,iy,iz,:) * st%surface_grid_element(1)
        e_field_box_surface(2,1,:) = e_field_box_surface(2,1,:) + e_surf(2,1,iy,iz,:) * st%surface_grid_element(1)
        b_field_box_surface(1,1,:) = b_field_box_surface(1,1,:) + b_surf(1,1,iy,iz,:) * st%surface_grid_element(1)
        b_field_box_surface(2,1,:) = b_field_box_surface(2,1,:) + b_surf(2,1,iy,iz,:) * st%surface_grid_element(1)
      end do
    end do

    do ix=1, ix_max
      do iz=1, iz_max
        do ip_surf=1, st%surface_grid_points_number(2,ix,iz)
          e_surf(1,2,ix,iz,:) = e_surf(1,2,ix,iz,:) &
                              + e_field_global(st%surface_grid_points_map(1,2,ix,iz,ip_surf),:)
          e_surf(2,2,ix,iz,:) = e_surf(2,2,ix,iz,:) &
                              + e_field_global(st%surface_grid_points_map(2,2,ix,iz,ip_surf),:)
          b_surf(1,2,ix,iz,:) = b_surf(1,2,ix,iz,:) &
                              + b_field_global(st%surface_grid_points_map(1,2,ix,iz,ip_surf),:)
          b_surf(2,2,ix,iz,:) = b_surf(2,2,ix,iz,:) &
                              + b_field_global(st%surface_grid_points_map(2,2,ix,iz,ip_surf),:)
        end do
        e_surf(1,2,ix,iz,:) = e_surf(1,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
        e_surf(2,2,ix,iz,:) = e_surf(2,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
        b_surf(1,2,ix,iz,:) = b_surf(1,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
        b_surf(2,2,ix,iz,:) = b_surf(2,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
      end do
    end do
    do ix=1, ix_max
      do iz=1, iz_max
        e_field_box_surface(1,2,:) = e_field_box_surface(1,2,:) + e_surf(1,2,ix,iz,:) * st%surface_grid_element(2)
        e_field_box_surface(2,2,:) = e_field_box_surface(2,2,:) + e_surf(2,2,ix,iz,:) * st%surface_grid_element(2)
        b_field_box_surface(1,2,:) = b_field_box_surface(1,2,:) + b_surf(1,2,ix,iz,:) * st%surface_grid_element(2)
        b_field_box_surface(2,2,:) = b_field_box_surface(2,2,:) + b_surf(2,2,ix,iz,:) * st%surface_grid_element(2)
      end do
    end do

    do ix=1, ix_max
      do iy=1, iy_max
        do ip_surf=1, st%surface_grid_points_number(3,ix,iy)
          e_surf(1,3,ix,iy,:) = e_surf(1,3,ix,iy,:) &
                              + e_field_global(st%surface_grid_points_map(1,3,ix,iy,ip_surf),:)
          e_surf(2,3,ix,iy,:) = e_surf(2,3,ix,iy,:) &
                              + e_field_global(st%surface_grid_points_map(2,3,ix,iy,ip_surf),:)
          b_surf(1,3,ix,iy,:) = b_surf(1,3,ix,iy,:) &
                              + b_field_global(st%surface_grid_points_map(1,3,ix,iy,ip_surf),:)
          b_surf(2,3,ix,iy,:) = b_surf(2,3,ix,iy,:) &
                              + b_field_global(st%surface_grid_points_map(2,3,ix,iy,ip_surf),:)
        end do
        e_surf(1,3,ix,iy,:) = e_surf(1,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
        e_surf(2,3,ix,iy,:) = e_surf(2,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
        b_surf(1,3,ix,iy,:) = b_surf(1,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
        b_surf(2,3,ix,iy,:) = b_surf(2,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
      end do
    end do
    do ix=1, ix_max
      do iy=1, iy_max
        e_field_box_surface(1,3,:) = e_field_box_surface(1,3,:) - e_surf(1,3,ix,iy,:) * st%surface_grid_element(3)
        e_field_box_surface(2,3,:) = e_field_box_surface(2,3,:) + e_surf(2,3,ix,iy,:) * st%surface_grid_element(3)
        b_field_box_surface(1,3,:) = b_field_box_surface(1,3,:) - b_surf(1,3,ix,iy,:) * st%surface_grid_element(3)
        b_field_box_surface(2,3,:) = b_field_box_surface(2,3,:) + b_surf(2,3,ix,iy,:) * st%surface_grid_element(3)
      end do
    end do

    SAFE_DEALLOCATE_A(e_surf)
    SAFE_DEALLOCATE_A(b_surf)
    SAFE_DEALLOCATE_A(e_field)
    SAFE_DEALLOCATE_A(b_field)
    SAFE_DEALLOCATE_A(e_field_global)
    SAFE_DEALLOCATE_A(b_field_global)

    POP_SUB(fields_through_box_surfaces)
  end subroutine fields_through_box_surfaces


  ! ---------------------------------------------------------
  subroutine fields_through_box_surfaces_plane_waves(gr, st, iter, dt, e_field_box_surface, b_field_box_surface)
    type(grid_t),        intent(in)    :: gr
    type(states_mxll_t),      intent(in)    :: st
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt
    FLOAT,               intent(out)   :: e_field_box_surface(:,:,:)
    FLOAT,               intent(out)   :: b_field_box_surface(:,:,:)

    integer             :: idim, idir, ip_surf, ip_global, ix, ix_max, iy, iy_max, iz, iz_max, ii_max
    FLOAT,  allocatable :: e_surf(:,:,:,:,:), b_surf(:,:,:,:,:), e_field(:,:), &
      e_field_global(:,:), b_field(:,:), b_field_global(:,:)

    PUSH_SUB(fields_through_box_surfaces_plane_waves)

    SAFE_ALLOCATE(e_field(1:gr%mesh%np,1:st%d%dim))
    SAFE_ALLOCATE(b_field(1:gr%mesh%np,1:st%d%dim))
    SAFE_ALLOCATE(e_field_global(1:gr%mesh%np_global,1:st%d%dim))
    SAFE_ALLOCATE(b_field_global(1:gr%mesh%np_global,1:st%d%dim))

    call get_electric_field_state(st%rs_state_plane_waves, gr%mesh, e_field, st%ep, gr%mesh%np)
    call get_magnetic_field_state(st%rs_state_plane_waves, gr%mesh, st%rs_sign, b_field, st%mu, gr%mesh%np)

    if (gr%mesh%parallel_in_domains) then
      do idim=1, st%d%dim
#if defined(HAVE_MPI)
        call vec_allgather(gr%mesh%vp, e_field_global(:,idim), e_field(:,idim))
        call MPI_Barrier(gr%mesh%vp%comm, mpi_err)
        call vec_allgather(gr%mesh%vp, b_field_global(:,idim), b_field(:,idim))
        call MPI_Barrier(gr%mesh%vp%comm, mpi_err)
#endif
      end do
    else
      e_field_global(:,:) = e_field(:,:)
      b_field_global(:,:) = b_field(:,:)
    end if

    ix_max = st%surface_grid_rows_number(1)
    iy_max = st%surface_grid_rows_number(2)
    iz_max = st%surface_grid_rows_number(3)
    ii_max = max(ix_max,iy_max,iz_max)

    SAFE_ALLOCATE(e_surf(1:2,1:st%d%dim,1:ii_max,1:ii_max,1:st%d%dim))
    SAFE_ALLOCATE(b_surf(1:2,1:st%d%dim,1:ii_max,1:ii_max,1:st%d%dim))

    e_surf = M_ZERO
    b_surf = M_ZERO
    e_field_box_surface = M_ZERO
    b_field_box_surface = M_ZERO

    do iy=1, iy_max
      do iz=1, iz_max
        do ip_surf=1, st%surface_grid_points_number(1,iy,iz)
          e_surf(1,1,iy,iz,:) = e_surf(1,1,iy,iz,:) & 
                                + e_field_global(st%surface_grid_points_map(1,1,iy,iz,ip_surf),:)
          e_surf(2,1,iy,iz,:) = e_surf(2,1,iy,iz,:) &
                                + e_field_global(st%surface_grid_points_map(2,1,iy,iz,ip_surf),:)
          b_surf(1,1,iy,iz,:) = b_surf(1,1,iy,iz,:) & 
                                + b_field_global(st%surface_grid_points_map(1,1,iy,iz,ip_surf),:)
          b_surf(2,1,iy,iz,:) = b_surf(2,1,iy,iz,:) &
                                + b_field_global(st%surface_grid_points_map(2,1,iy,iz,ip_surf),:)
        end do
        e_surf(1,1,iy,iz,:) = e_surf(1,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
        e_surf(2,1,iy,iz,:) = e_surf(2,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
        b_surf(1,1,iy,iz,:) = b_surf(1,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
        b_surf(2,1,iy,iz,:) = b_surf(2,1,iy,iz,:) / float(st%surface_grid_points_number(1,iy,iz))
      end do
    end do
    do iy=1, iy_max
      do iz=1, iz_max
        e_field_box_surface(1,1,:) = e_field_box_surface(1,1,:) + e_surf(1,1,iy,iz,:) * st%surface_grid_element(1)
        e_field_box_surface(2,1,:) = e_field_box_surface(2,1,:) + e_surf(2,1,iy,iz,:) * st%surface_grid_element(1)
        b_field_box_surface(1,1,:) = b_field_box_surface(1,1,:) + b_surf(1,1,iy,iz,:) * st%surface_grid_element(1)
        b_field_box_surface(2,1,:) = b_field_box_surface(2,1,:) + b_surf(2,1,iy,iz,:) * st%surface_grid_element(1)
      end do
    end do

    do ix=1, ix_max
      do iz=1, iz_max
        do ip_surf=1, st%surface_grid_points_number(2,ix,iz)
          e_surf(1,2,ix,iz,:) = e_surf(1,2,ix,iz,:) &
                              + e_field_global(st%surface_grid_points_map(1,2,ix,iz,ip_surf),:)
          e_surf(2,2,ix,iz,:) = e_surf(2,2,ix,iz,:) &
                              + e_field_global(st%surface_grid_points_map(2,2,ix,iz,ip_surf),:)
          b_surf(1,2,ix,iz,:) = b_surf(1,2,ix,iz,:) &
                              + b_field_global(st%surface_grid_points_map(1,2,ix,iz,ip_surf),:)
          b_surf(2,2,ix,iz,:) = b_surf(2,2,ix,iz,:) &
                              + b_field_global(st%surface_grid_points_map(2,2,ix,iz,ip_surf),:)
        end do
        e_surf(1,2,ix,iz,:) = e_surf(1,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
        e_surf(2,2,ix,iz,:) = e_surf(2,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
        b_surf(1,2,ix,iz,:) = b_surf(1,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
        b_surf(2,2,ix,iz,:) = b_surf(2,2,ix,iz,:) / float(st%surface_grid_points_number(2,ix,iz))
      end do
    end do
    do ix=1, ix_max
      do iz=1, iz_max
        e_field_box_surface(1,2,:) = e_field_box_surface(1,2,:) + e_surf(1,2,ix,iz,:) * st%surface_grid_element(2)
        e_field_box_surface(2,2,:) = e_field_box_surface(2,2,:) + e_surf(2,2,ix,iz,:) * st%surface_grid_element(2)
        b_field_box_surface(1,2,:) = b_field_box_surface(1,2,:) + b_surf(1,2,ix,iz,:) * st%surface_grid_element(2)
        b_field_box_surface(2,2,:) = b_field_box_surface(2,2,:) + b_surf(2,2,ix,iz,:) * st%surface_grid_element(2)
      end do
    end do

    do ix=1, ix_max
      do iy=1, iy_max
        do ip_surf=1, st%surface_grid_points_number(3,ix,iy)
          e_surf(1,3,ix,iy,:) = e_surf(1,3,ix,iy,:) &
                              + e_field_global(st%surface_grid_points_map(1,3,ix,iy,ip_surf),:)
          e_surf(2,3,ix,iy,:) = e_surf(2,3,ix,iy,:) &
                              + e_field_global(st%surface_grid_points_map(2,3,ix,iy,ip_surf),:)
          b_surf(1,3,ix,iy,:) = b_surf(1,3,ix,iy,:) &
                              + b_field_global(st%surface_grid_points_map(1,3,ix,iy,ip_surf),:)
          b_surf(2,3,ix,iy,:) = b_surf(2,3,ix,iy,:) &
                              + b_field_global(st%surface_grid_points_map(2,3,ix,iy,ip_surf),:)
        end do
        e_surf(1,3,ix,iy,:) = e_surf(1,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
        e_surf(2,3,ix,iy,:) = e_surf(2,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
        b_surf(1,3,ix,iy,:) = b_surf(1,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
        b_surf(2,3,ix,iy,:) = b_surf(2,3,ix,iy,:) / float(st%surface_grid_points_number(3,ix,iy))
      end do
    end do
    do ix=1, ix_max
      do iy=1, iy_max
        e_field_box_surface(1,3,:) = e_field_box_surface(1,3,:) - e_surf(1,3,ix,iy,:) * st%surface_grid_element(3)
        e_field_box_surface(2,3,:) = e_field_box_surface(2,3,:) + e_surf(2,3,ix,iy,:) * st%surface_grid_element(3)
        b_field_box_surface(1,3,:) = b_field_box_surface(1,3,:) - b_surf(1,3,ix,iy,:) * st%surface_grid_element(3)
        b_field_box_surface(2,3,:) = b_field_box_surface(2,3,:) + b_surf(2,3,ix,iy,:) * st%surface_grid_element(3)
      end do
    end do

    SAFE_DEALLOCATE_A(e_surf)
    SAFE_DEALLOCATE_A(b_surf)
    SAFE_DEALLOCATE_A(e_field)
    SAFE_DEALLOCATE_A(b_field)
    SAFE_DEALLOCATE_A(e_field_global)
    SAFE_DEALLOCATE_A(b_field_global)

    POP_SUB(fields_through_box_surfaces_plane_waves)
  end subroutine fields_through_box_surfaces_plane_waves


  ! ---------------------------------------------------------
  subroutine mask_absorbing_boundaries(gr, hm, st, tr, time, dt, time_delay, rs_state)
    type(grid_t),               intent(in)    :: gr
    type(hamiltonian_mxll_t),        intent(inout) :: hm
    type(states_mxll_t),             intent(in)    :: st
    type(propagator_mxll_t), intent(in)    :: tr
    FLOAT,                      intent(in)    :: time
    FLOAT,                      intent(in)    :: dt
    FLOAT,                      intent(in)    :: time_delay
    CMPLX,                      intent(inout) :: rs_state(:,:)

    integer            :: ip, ip_in, ip_lim, idim
    FLOAT              :: ka(MAX_DIM), si(MAX_DIM), ee(MAX_DIM), bb(MAX_DIM), ee_aux(MAX_DIM), bb_aux(MAX_DIM), j_aux(MAX_DIM)
    logical            :: mask_check = .false.

    PUSH_SUB(mask_absorbing_boundaries)

    do idim=1, 3
      if (hm%bc%bc_ab_type(idim) == OPTION__MAXWELLABSORBINGBOUNDARIES__MASK) then
        mask_check = .true.
      end if
    end do

    if (mask_check) then
      if (tr%bc_plane_waves .and. hm%plane_waves_apply) then
        call plane_waves_propagation(hm, st, gr, tr, time, dt, time_delay)
        rs_state = rs_state - st%rs_state_plane_waves
        call maxwell_mask(gr, hm, st, rs_state)
        rs_state = rs_state + st%rs_state_plane_waves
      else if (tr%bc_constant .and. hm%spatial_constant_apply) then
        !call constant_at_absorbing_boundaries_calculation(st, hm%bc)
        call constant_boundaries_calculation(tr%bc_constant, hm%bc, hm, st, rs_state)
        do ip_in=1, hm%bc%constant_points_number
          ip = hm%bc%constant_points_map(ip_in)
          rs_state(ip,:) = rs_state(ip,:) - st%rs_state_const(:)
        end do
        call maxwell_mask(gr, hm, st, rs_state)
        do ip_in=1, hm%bc%constant_points_number
          ip = hm%bc%constant_points_map(ip_in)
          rs_state(ip,:) = rs_state(ip,:) + st%rs_state_const(:)
        end do
      else
        call maxwell_mask(gr, hm, st, rs_state)
      end if
    end if

    POP_SUB(mask_absorbing_boundaries)
  end subroutine mask_absorbing_boundaries


  ! ---------------------------------------------------------
  subroutine maxwell_mask(gr, hm, st, rs_state)
    type(grid_t),        intent(in)    :: gr
    type(hamiltonian_mxll_t), intent(in)    :: hm
    type(states_mxll_t),      intent(in)    :: st
    CMPLX,               intent(inout) :: rs_state(:,:)

    integer :: ip, ip_in, idim

    PUSH_SUB(maxwell_mask)

    do idim=1, 3
      if (hm%bc%bc_ab_type(idim) == OPTION__MAXWELLABSORBINGBOUNDARIES__MASK) then
        do ip_in=1, hm%bc%mask_points_number(idim)
          ip = hm%bc%mask_points_map(ip_in,idim)
          rs_state(ip,:) = rs_state(ip,:)*hm%bc%mask(ip_in,idim)
        end do
      end if
    end do

    POP_SUB(maxwell_mask)
  end subroutine maxwell_mask


  ! ---------------------------------------------------------
  subroutine pml_propagation_stage_1(hm, gr, st, tr, ff_rs_state, ff_rs_state_pml)
    type(hamiltonian_mxll_t),        intent(inout) :: hm
    type(grid_t),               intent(in)    :: gr
    type(states_mxll_t),             intent(inout) :: st
    type(propagator_mxll_t), intent(inout) :: tr
    CMPLX,                      intent(in)    :: ff_rs_state(:,:)
    CMPLX,                      intent(inout) :: ff_rs_state_pml(:,:)

    integer            :: ip, ip_in, ff_points, ff_dim
    CMPLX, allocatable :: ff_rs_state_plane_waves(:,:), rs_state_constant(:,:), ff_rs_state_constant(:,:)

    PUSH_SUB(pml_propagation_stage_1)

    if (tr%bc_plane_waves .and. hm%plane_waves_apply) then
      ff_points = size(ff_rs_state(:,1))
      ff_dim    = size(ff_rs_state(1,:))
      SAFE_ALLOCATE(ff_rs_state_plane_waves(1:ff_points,1:ff_dim))
      call transform_rs_state_forward(hm, gr, st, st%rs_state_plane_waves, ff_rs_state_plane_waves)
      ff_rs_state_pml = ff_rs_state - ff_rs_state_plane_waves
      SAFE_DEALLOCATE_A(ff_rs_state_plane_waves)
    else if (tr%bc_constant .and. hm%spatial_constant_apply) then
      ff_dim    = size(ff_rs_state(1,:))
      SAFE_ALLOCATE(rs_state_constant(1,1:3))
      SAFE_ALLOCATE(ff_rs_state_constant(1,1:ff_dim))
      rs_state_constant(1,:) = st%rs_state_const(:)
      call transform_rs_state_forward(hm, gr, st, rs_state_constant, ff_rs_state_constant)
      do ip=1, gr%mesh%np_part
        ff_rs_state_pml(ip,:) = ff_rs_state(ip,:) - ff_rs_state_constant(1,:)
      end do
      SAFE_DEALLOCATE_A(rs_state_constant)
      SAFE_DEALLOCATE_A(ff_rs_state_constant)
    else
      ff_rs_state_pml = ff_rs_state
    end if

    POP_SUB(pml_propagation_stage_1)
  end subroutine pml_propagation_stage_1


  ! ---------------------------------------------------------
  subroutine pml_propagation_stage_2(hm, gr, st, tr, time, dt, time_delay, ff_rs_state_pml, ff_rs_state)
    type(hamiltonian_mxll_t),        intent(inout) :: hm
    type(grid_t),               intent(in)    :: gr
    type(states_mxll_t),             intent(inout) :: st
    type(propagator_mxll_t), intent(inout) :: tr
    FLOAT,                      intent(in)    :: time
    FLOAT,                      intent(in)    :: dt
    FLOAT,                      intent(in)    :: time_delay
    CMPLX,                      intent(inout) :: ff_rs_state_pml(:,:)
    CMPLX,                      intent(inout) :: ff_rs_state(:,:)

    integer            :: ip, ip_in, ff_points, ff_dim
    CMPLX, allocatable :: ff_rs_state_plane_waves(:,:), rs_state_constant(:,:), ff_rs_state_constant(:,:)

    PUSH_SUB(pml_propagation_stage_2)

    if (tr%bc_plane_waves .and. hm%plane_waves_apply) then
      ff_points = size(ff_rs_state(:,1))
      ff_dim    = size(ff_rs_state(1,:))
      hm%cpml_hamiltonian = .true.
      call exponential_mxll_apply(hm, gr, st, tr, time, dt, ff_rs_state_pml)
      hm%cpml_hamiltonian = .false.
      call plane_waves_propagation(hm, st, gr, tr, time, dt, time_delay)
      SAFE_ALLOCATE(ff_rs_state_plane_waves(1:ff_points,1:ff_dim))
      call transform_rs_state_forward(hm, gr, st, st%rs_state_plane_waves, ff_rs_state_plane_waves)
      do ip_in=1, hm%bc%plane_waves_points_number
        ip = hm%bc%plane_waves_points_map(ip_in)
        ff_rs_state(ip,:) = ff_rs_state_pml(ip,:) + ff_rs_state_plane_waves(ip,:)
      end do
      SAFE_DEALLOCATE_A(ff_rs_state_plane_waves)
    else if (tr%bc_constant .and. hm%spatial_constant_apply) then
      ff_dim    = size(ff_rs_state(1,:))
      hm%cpml_hamiltonian = .true.
      call exponential_mxll_apply(hm, gr, st, tr, time, dt, ff_rs_state_pml)
      hm%cpml_hamiltonian = .false.
      SAFE_ALLOCATE(rs_state_constant(1,1:ff_dim))
      SAFE_ALLOCATE(ff_rs_state_constant(1,1:ff_dim))
      rs_state_constant(1,:) = st%rs_state_const(:)
      call transform_rs_state_forward(hm, gr, st, rs_state_constant, ff_rs_state_constant)
      do ip_in=1, hm%bc%constant_points_number
        ip = hm%bc%constant_points_map(ip_in)
        ff_rs_state(ip,:) = ff_rs_state_pml(ip,:) + ff_rs_state_constant(1,:)
      end do
      SAFE_DEALLOCATE_A(rs_state_constant)
      SAFE_DEALLOCATE_A(ff_rs_state_constant)
    end if

    POP_SUB(pml_propagation_stage_2)
  end subroutine pml_propagation_stage_2


  ! ---------------------------------------------------------
  subroutine cpml_conv_function_update(hm, gr, ff_rs_state_pml, ff_dim)
    type(hamiltonian_mxll_t), intent(inout) :: hm
    type(grid_t),        intent(in)    :: gr
    CMPLX,               intent(inout) :: ff_rs_state_pml(:,:)
    integer,             intent(in)    :: ff_dim

    PUSH_SUB(cpml_conv_function_update)

    if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RIEMANN_SILBERSTEIN) then
      call cpml_conv_function_update_via_riemann_silberstein(hm, gr, ff_rs_state_pml, ff_dim)
    else if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__ELECTRIC_MAGNETIC_FIELDS) then
      call cpml_conv_function_update_via_e_b_fields(hm, gr, ff_rs_state_pml, ff_dim)
    end if

    POP_SUB(cpml_conv_function_update)
  end subroutine cpml_conv_function_update


  ! ---------------------------------------------------------
  subroutine cpml_conv_function_update_via_riemann_silberstein(hm, gr, ff_rs_state_pml, ff_dim)
    type(hamiltonian_mxll_t), intent(inout) :: hm
    type(grid_t),        intent(in)    :: gr
    CMPLX,               intent(inout) :: ff_rs_state_pml(:,:)
    integer,             intent(in)    :: ff_dim

    integer :: ip, ip_in, np_part, rs_sign
    CMPLX, allocatable :: tmp_partial(:), tmp_partial_2(:,:)
    CMPLX              :: pml_a(3), pml_b(3), pml_g(3), pml_g_p(3), pml_g_m(3)

    PUSH_SUB(cpml_conv_function_update_via_riemann_silberstein)

    np_part = gr%der%mesh%np_part
    rs_sign = hm%rs_sign

    if (ff_dim == 3) then

      SAFE_ALLOCATE(tmp_partial(np_part))

      ! calculation g(1,2)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,2), tmp_partial(:), 1, set_bc = .false.)
      do ip_in=1, hm%bc%pml_points_number
        ip       = hm%bc%pml_points_map(ip_in)
        pml_a(:) = hm%bc%pml_a(ip_in,:)
        pml_b(:) = hm%bc%pml_b(ip_in,:)
        pml_g(:) = hm%bc%pml_conv_plus(ip_in,1,:)
        pml_g(2) = real(pml_a(1)) * real(tmp_partial(ip)) + real(pml_b(1)) * real(pml_g(2)) + &
                   M_zI * ( aimag(pml_a(1)) * aimag(tmp_partial(ip)) + aimag(pml_b(1)) * aimag(pml_g(2)) )
        hm%bc%pml_conv_plus(ip_in,1,:) = pml_g(:)
      end do
      ! calculation g(2,1)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,1), tmp_partial(:), 2, set_bc = .false.)
      do ip_in=1, hm%bc%pml_points_number
        ip       = hm%bc%pml_points_map(ip_in)
        pml_a(:) = hm%bc%pml_a(ip_in,:)
        pml_b(:) = hm%bc%pml_b(ip_in,:)
        pml_g(:) = hm%bc%pml_conv_plus(ip_in,2,:)
        pml_g(1) = real(pml_a(2)) * real(tmp_partial(ip)) + real(pml_b(2)) * real(pml_g(1)) + &
                   M_zI * ( aimag(pml_a(2)) * aimag(tmp_partial(ip)) + aimag(pml_b(2)) * aimag(pml_g(1)) )
        hm%bc%pml_conv_plus(ip_in,2,:) = pml_g(:)
      end do
      ! calculation g(1,3)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,3), tmp_partial(:), 1, set_bc = .false.)
      do ip_in=1, hm%bc%pml_points_number
        ip       = hm%bc%pml_points_map(ip_in)
        pml_a(:) = hm%bc%pml_a(ip_in,:)
        pml_b(:) = hm%bc%pml_b(ip_in,:)
        pml_g(:) = hm%bc%pml_conv_plus(ip_in,1,:)
        pml_g(3) = real(pml_a(1)) * real(tmp_partial(ip)) + real(pml_b(1)) * real(pml_g(3)) + &
                   M_zI * ( aimag(pml_a(1)) * aimag(tmp_partial(ip)) + aimag(pml_b(1)) * aimag(pml_g(3)) )
        hm%bc%pml_conv_plus(ip_in,1,:) = pml_g(:)
      end do
      ! calculation g(3,1)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,1), tmp_partial(:), 3, set_bc = .false.)
      do ip_in=1, hm%bc%pml_points_number
        ip       = hm%bc%pml_points_map(ip_in)
        pml_a(:) = hm%bc%pml_a(ip_in,:)
        pml_b(:) = hm%bc%pml_b(ip_in,:)
        pml_g(:) = hm%bc%pml_conv_plus(ip_in,3,:)
        pml_g(1) = real(pml_a(3)) * real(tmp_partial(ip)) + real(pml_b(3)) * real(pml_g(1)) + &
                   M_zI * ( aimag(pml_a(3)) * aimag(tmp_partial(ip)) + aimag(pml_b(3)) * aimag(pml_g(1)) )
        hm%bc%pml_conv_plus(ip_in,3,:) = pml_g(:)
      end do
      ! calculation g(2,3)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,3), tmp_partial(:), 2, set_bc = .false.)
      do ip_in=1, hm%bc%pml_points_number
        ip       = hm%bc%pml_points_map(ip_in)
        pml_a(:) = hm%bc%pml_a(ip_in,:)
        pml_b(:) = hm%bc%pml_b(ip_in,:)
        pml_g(:) = hm%bc%pml_conv_plus(ip_in,2,:)
        pml_g(3) = real(pml_a(2)) * real(tmp_partial(ip)) + real(pml_b(2)) * real(pml_g(3)) + &
                   M_zI * ( aimag(pml_a(2)) * aimag(tmp_partial(ip)) + aimag(pml_b(2)) * aimag(pml_g(3)) )
        hm%bc%pml_conv_plus(ip_in,2,:) = pml_g(:)
      end do
      ! calculation g(3,2)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,2), tmp_partial(:), 3, set_bc = .false.)
      do ip_in=1, hm%bc%pml_points_number
        ip       = hm%bc%pml_points_map(ip_in)
        pml_a(:) = hm%bc%pml_a(ip_in,:)
        pml_b(:) = hm%bc%pml_b(ip_in,:)
        pml_g(:) = hm%bc%pml_conv_plus(ip_in,3,:)
        pml_g(2) = real(pml_a(3)) * real(tmp_partial(ip)) + real(pml_b(3)) * real(pml_g(2)) + &
                   M_zI * ( aimag(pml_a(3)) * aimag(tmp_partial(ip)) + aimag(pml_b(3)) * aimag(pml_g(2)) )
        hm%bc%pml_conv_plus(ip_in,3,:) = pml_g(:)
      end do

      SAFE_DEALLOCATE_A(tmp_partial)

    else if (ff_dim == 6) then

      SAFE_ALLOCATE(tmp_partial_2(np_part,2))

      ! calculation g(1,2)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,2), tmp_partial_2(:,1), 1, set_bc = .false.)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,5), tmp_partial_2(:,2), 1, set_bc = .false.)
      do ip_in=1, hm%bc%pml_points_number
        ip         = hm%bc%pml_points_map(ip_in)
        pml_a(:)   = hm%bc%pml_a(ip_in,:)
        pml_b(:)   = hm%bc%pml_b(ip_in,:)
        pml_g_p(:) = hm%bc%pml_conv_plus(ip_in,1,:)
        pml_g_m(:) = hm%bc%pml_conv_minus(ip_in,1,:)
        pml_g_p(2) = real(pml_a(1)) * real(tmp_partial_2(ip,1)) + real(pml_b(1)) * real(pml_g_p(2)) + &
                     M_zI * ( aimag(pml_a(1)) * aimag(tmp_partial_2(ip,1)) + aimag(pml_b(1)) * aimag(pml_g_p(2)) )
        pml_g_m(2) = real(pml_a(1)) * real(tmp_partial_2(ip,2)) + real(pml_b(1)) * real(pml_g_m(2)) + &
                     M_zI * ( aimag(pml_a(1)) * aimag(tmp_partial_2(ip,2)) + aimag(pml_b(1)) * aimag(pml_g_m(2)) )
        hm%bc%pml_conv_plus(ip_in,1,:) = pml_g_p(:)
        hm%bc%pml_conv_minus(ip_in,1,:) = pml_g_m(:)
      end do
      ! calculation g(2,1)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,1), tmp_partial_2(:,1), 2, set_bc = .false.)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,4), tmp_partial_2(:,2), 2, set_bc = .false.)
      do ip_in=1, hm%bc%pml_points_number
        ip         = hm%bc%pml_points_map(ip_in)
        pml_a(:)   = hm%bc%pml_a(ip_in,:)
        pml_b(:)   = hm%bc%pml_b(ip_in,:)
        pml_g_p(:) = hm%bc%pml_conv_plus(ip_in,2,:)
        pml_g_m(:) = hm%bc%pml_conv_minus(ip_in,2,:)
        pml_g_p(1) = real(pml_a(2)) * real(tmp_partial_2(ip,1)) + real(pml_b(2)) * real(pml_g_p(1)) + &
                     M_zI * ( aimag(pml_a(2)) * aimag(tmp_partial_2(ip,1)) + aimag(pml_b(2)) * aimag(pml_g_p(1)) )
        pml_g_m(1) = real(pml_a(2)) * real(tmp_partial_2(ip,2)) + real(pml_b(2)) * real(pml_g_m(1)) + &
                     M_zI * ( aimag(pml_a(2)) * aimag(tmp_partial_2(ip,2)) + aimag(pml_b(2)) * aimag(pml_g_m(1)) )
        hm%bc%pml_conv_plus(ip_in,2,:) = pml_g_p(:)
        hm%bc%pml_conv_minus(ip_in,2,:) = pml_g_m(:)
      end do
      ! calculation g(1,3)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,3), tmp_partial_2(:,1), 1, set_bc = .false.)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,6), tmp_partial_2(:,2), 1, set_bc = .false.)
      do ip_in=1, hm%bc%pml_points_number
        ip         = hm%bc%pml_points_map(ip_in)
        pml_a(:)   = hm%bc%pml_a(ip_in,:)
        pml_b(:)   = hm%bc%pml_b(ip_in,:)
        pml_g_p(:) = hm%bc%pml_conv_plus(ip_in,1,:)
        pml_g_m(:) = hm%bc%pml_conv_minus(ip_in,1,:)
        pml_g_p(3) = real(pml_a(1)) * real(tmp_partial_2(ip,1)) + real(pml_b(1)) * real(pml_g_p(3)) + &
                     M_zI * ( aimag(pml_a(1)) * aimag(tmp_partial_2(ip,1)) + aimag(pml_b(1)) * aimag(pml_g_p(3)) )
        pml_g_m(3) = real(pml_a(1)) * real(tmp_partial_2(ip,2)) + real(pml_b(1)) * real(pml_g_m(3)) + &
                     M_zI * ( aimag(pml_a(1)) * aimag(tmp_partial_2(ip,2)) + aimag(pml_b(1)) * aimag(pml_g_m(3)) )
        hm%bc%pml_conv_plus(ip_in,1,:) = pml_g_p(:)
        hm%bc%pml_conv_minus(ip_in,1,:) = pml_g_m(:)
      end do
      ! calculation g(3,1)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,1), tmp_partial_2(:,1), 3, set_bc = .false.)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,4), tmp_partial_2(:,2), 3, set_bc = .false.)
      do ip_in=1, hm%bc%pml_points_number
        ip         = hm%bc%pml_points_map(ip_in)
        pml_a(:)   = hm%bc%pml_a(ip_in,:)
        pml_b(:)   = hm%bc%pml_b(ip_in,:)
        pml_g_p(:) = hm%bc%pml_conv_plus(ip_in,3,:)
        pml_g_m(:) = hm%bc%pml_conv_minus(ip_in,3,:)
        pml_g_p(1) = real(pml_a(3)) * real(tmp_partial_2(ip,1)) + real(pml_b(3)) * real(pml_g_p(1)) + &
                     M_zI * ( aimag(pml_a(3)) * aimag(tmp_partial_2(ip,1)) + aimag(pml_b(3)) * aimag(pml_g_p(1)) )
        pml_g_m(1) = real(pml_a(3)) * real(tmp_partial_2(ip,2)) + real(pml_b(3)) * real(pml_g_m(1)) + &
                     M_zI * ( aimag(pml_a(3)) * aimag(tmp_partial_2(ip,2)) + aimag(pml_b(3)) * aimag(pml_g_m(1)) )
        hm%bc%pml_conv_plus(ip_in,3,:) = pml_g_p(:)
        hm%bc%pml_conv_minus(ip_in,3,:) = pml_g_m(:)
      end do
      ! calculation g(2,3)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,3), tmp_partial_2(:,1), 2, set_bc = .false.)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,6), tmp_partial_2(:,2), 2, set_bc = .false.)
      do ip_in=1, hm%bc%pml_points_number
        ip         = hm%bc%pml_points_map(ip_in)
        pml_a(:)   = hm%bc%pml_a(ip_in,:)
        pml_b(:)   = hm%bc%pml_b(ip_in,:)
        pml_g_p(:) = hm%bc%pml_conv_plus(ip_in,2,:)
        pml_g_m(:) = hm%bc%pml_conv_minus(ip_in,2,:)
        pml_g_p(3) = real(pml_a(2)) * real(tmp_partial_2(ip,1)) + real(pml_b(2)) * real(pml_g_p(3)) + &
                     M_zI * ( aimag(pml_a(2)) * aimag(tmp_partial_2(ip,1)) + aimag(pml_b(2)) * aimag(pml_g_p(3)) )
        pml_g_m(3) = real(pml_a(2)) * real(tmp_partial_2(ip,2)) + real(pml_b(2)) * real(pml_g_m(3)) + &
                     M_zI * ( aimag(pml_a(2)) * aimag(tmp_partial_2(ip,2)) + aimag(pml_b(2)) * aimag(pml_g_m(3)) )
        hm%bc%pml_conv_plus(ip_in,2,:) = pml_g_p(:)
        hm%bc%pml_conv_minus(ip_in,2,:) = pml_g_m(:)
      end do
      ! calculation g(3,2)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,2), tmp_partial_2(:,1), 3, set_bc = .false.)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,5), tmp_partial_2(:,2), 3, set_bc = .false.)
      do ip_in=1, hm%bc%pml_points_number
        ip         = hm%bc%pml_points_map(ip_in)
        pml_a(:)   = hm%bc%pml_a(ip_in,:)
        pml_b(:)   = hm%bc%pml_b(ip_in,:)
        pml_g_p(:) = hm%bc%pml_conv_plus(ip_in,3,:)
        pml_g_m(:) = hm%bc%pml_conv_minus(ip_in,3,:)
        pml_g_p(2) = real(pml_a(3)) * real(tmp_partial_2(ip,1)) + real(pml_b(3)) * real(pml_g_p(2)) + &
                     M_zI * ( aimag(pml_a(3)) * aimag(tmp_partial_2(ip,1)) + aimag(pml_b(3)) * aimag(pml_g_p(2)) )
        pml_g_m(2) = real(pml_a(3)) * real(tmp_partial_2(ip,2)) + real(pml_b(3)) * real(pml_g_m(2)) + &
                     M_zI * ( aimag(pml_a(3)) * aimag(tmp_partial_2(ip,2)) + aimag(pml_b(3)) * aimag(pml_g_m(2)) )
        hm%bc%pml_conv_plus(ip_in,3,:) = pml_g_p(:)
        hm%bc%pml_conv_minus(ip_in,3,:) = pml_g_m(:)
      end do

      SAFE_DEALLOCATE_A(tmp_partial_2)

    end if

    POP_SUB(cpml_conv_function_update_via_riemann_silberstein)
  end subroutine cpml_conv_function_update_via_riemann_silberstein


  ! ---------------------------------------------------------
  subroutine cpml_conv_function_update_via_e_b_fields(hm, gr, ff_rs_state_pml, ff_dim)
    type(hamiltonian_mxll_t), intent(inout) :: hm
    type(grid_t),        intent(in)    :: gr
    CMPLX,               intent(in)    :: ff_rs_state_pml(:,:)
    integer,             intent(in)    :: ff_dim

    integer :: ip, ip_in, np_part
    FLOAT, allocatable :: tmp_e(:,:), tmp_b(:,:), tmp_partial_e(:), tmp_partial_b(:)
    CMPLX              :: pml_a(3), pml_b(3), pml_g(3)

    PUSH_SUB(cpml_conv_function_update_via_e_b_fields)

    np_part = gr%der%mesh%np_part
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
    do ip_in=1, hm%bc%pml_points_number
      ip       = hm%bc%pml_points_map(ip_in)
      pml_a(:) = hm%bc%pml_a(ip_in,:)
      pml_b(:) = hm%bc%pml_b(ip_in,:)
      pml_g(:) = hm%bc%pml_conv_plus(ip_in,1,:)
      pml_g(2) = real(pml_a(1)) * tmp_partial_e(ip) + real(pml_b(1)) * real(pml_g(2)) + &
                 M_zI * ( aimag(pml_a(1)) * tmp_partial_b(ip) + aimag(pml_b(1)) * aimag(pml_g(2)) )
      hm%bc%pml_conv_plus(ip_in,1,:) = pml_g(:)
    end do

    ! calculation g(2,1)
    call dderivatives_partial(gr%der, tmp_e(:,1), tmp_partial_e(:), 2, set_bc = .false.)
    call dderivatives_partial(gr%der, tmp_b(:,1), tmp_partial_b(:), 2, set_bc = .false.)
    tmp_partial_e(:) = sqrt(P_ep/M_TWO) * tmp_partial_e(:)
    tmp_partial_b(:) = sqrt(M_ONE/(M_TWO*P_mu)) * tmp_partial_b(:)
    do ip_in=1, hm%bc%pml_points_number
      ip       = hm%bc%pml_points_map(ip_in)
      pml_a(:) = hm%bc%pml_a(ip_in,:)
      pml_b(:) = hm%bc%pml_b(ip_in,:)
      pml_g(:) = hm%bc%pml_conv_plus(ip_in,2,:)
      pml_g(1) = real(pml_a(2)) * tmp_partial_e(ip) + real(pml_b(2)) * real(pml_g(1)) + &
                 M_zI * ( aimag(pml_a(2)) * tmp_partial_b(ip) + aimag(pml_b(2)) * aimag(pml_g(1)) )
      hm%bc%pml_conv_plus(ip_in,2,:) = pml_g(:)
    end do

    ! calculation g(1,3)
    call dderivatives_partial(gr%der, tmp_e(:,3), tmp_partial_e(:), 1, set_bc = .false.)
    call dderivatives_partial(gr%der, tmp_b(:,3), tmp_partial_b(:), 1, set_bc = .false.)
    tmp_partial_e(:) = sqrt(P_ep/M_TWO) * tmp_partial_e(:)
    tmp_partial_b(:) = sqrt(M_ONE/(M_TWO*P_mu)) * tmp_partial_b(:)
    do ip_in=1, hm%bc%pml_points_number
      ip       = hm%bc%pml_points_map(ip_in)
      pml_a(:) = hm%bc%pml_a(ip_in,:)
      pml_b(:) = hm%bc%pml_b(ip_in,:)
      pml_g(:) = hm%bc%pml_conv_plus(ip_in,1,:)
      pml_g(3) = real(pml_a(1)) * tmp_partial_e(ip) + real(pml_b(1)) * real(pml_g(3)) + &
                 M_zI * ( aimag(pml_a(1)) * tmp_partial_b(ip) + aimag(pml_b(1)) * aimag(pml_g(3)) )
      hm%bc%pml_conv_plus(ip_in,1,:) = pml_g(:)
    end do

    ! calculation g(3,1)
    call dderivatives_partial(gr%der, tmp_e(:,1), tmp_partial_e(:), 3, set_bc = .false.)
    call dderivatives_partial(gr%der, tmp_b(:,1), tmp_partial_b(:), 3, set_bc = .false.)
    tmp_partial_e(:) = sqrt(P_ep/M_TWO) * tmp_partial_e(:)
    tmp_partial_b(:) = sqrt(M_ONE/(M_TWO*P_mu)) * tmp_partial_b(:)
    do ip_in=1, hm%bc%pml_points_number
      ip       = hm%bc%pml_points_map(ip_in)
      pml_a(:) = hm%bc%pml_a(ip_in,:)
      pml_b(:) = hm%bc%pml_b(ip_in,:)
      pml_g(:) = hm%bc%pml_conv_plus(ip_in,3,:)
      pml_g(1) = real(pml_a(3)) * tmp_partial_e(ip) + real(pml_b(3)) * real(pml_g(1)) + &
                 M_zI * ( aimag(pml_a(3)) * tmp_partial_b(ip) + aimag(pml_b(3)) * aimag(pml_g(1)) )
      hm%bc%pml_conv_plus(ip_in,3,:) = pml_g(:)
    end do

    ! calculation g(2,3)
    call dderivatives_partial(gr%der, tmp_e(:,3), tmp_partial_e(:), 2, set_bc = .false.)
    call dderivatives_partial(gr%der, tmp_b(:,3), tmp_partial_b(:), 2, set_bc = .false.)
    tmp_partial_e(:) = sqrt(P_ep/M_TWO) * tmp_partial_e(:)
    tmp_partial_b(:) = sqrt(M_ONE/(M_TWO*P_mu)) * tmp_partial_b(:)
    do ip_in=1, hm%bc%pml_points_number
      ip       = hm%bc%pml_points_map(ip_in)
      pml_a(:) = hm%bc%pml_a(ip_in,:)
      pml_b(:) = hm%bc%pml_b(ip_in,:)
      pml_g(:) = hm%bc%pml_conv_plus(ip_in,2,:)
      pml_g(3) = real(pml_a(2)) * tmp_partial_e(ip) + real(pml_b(2)) * real(pml_g(3)) + &
                 M_zI * ( aimag(pml_a(2)) * tmp_partial_b(ip) + aimag(pml_b(2)) * aimag(pml_g(3)) )
      hm%bc%pml_conv_plus(ip_in,2,:) = pml_g(:)
    end do

    ! calculation g(3,2)
    call dderivatives_partial(gr%der, tmp_e(:,2), tmp_partial_e(:), 3, set_bc = .false.)
    call dderivatives_partial(gr%der, tmp_b(:,2), tmp_partial_b(:), 3, set_bc = .false.)
    tmp_partial_e(:) = sqrt(P_ep/M_TWO) * tmp_partial_e(:)
    tmp_partial_b(:) = sqrt(M_ONE/(M_TWO*P_mu)) * tmp_partial_b(:)
    do ip_in=1, hm%bc%pml_points_number
      ip       = hm%bc%pml_points_map(ip_in)
      pml_a(:) = hm%bc%pml_a(ip_in,:)
      pml_b(:) = hm%bc%pml_b(ip_in,:)
      pml_g(:) = hm%bc%pml_conv_plus(ip_in,3,:)
      pml_g(2) = real(pml_a(3)) * tmp_partial_e(ip) + real(pml_b(3)) * real(pml_g(2)) + &
                 M_zI * ( aimag(pml_a(3)) * tmp_partial_b(ip) + aimag(pml_b(3)) * aimag(pml_g(2)) )
      hm%bc%pml_conv_plus(ip_in,3,:) = pml_g(:)
    end do

    SAFE_DEALLOCATE_A(tmp_e)
    SAFE_DEALLOCATE_A(tmp_partial_e)
    SAFE_DEALLOCATE_A(tmp_b)
    SAFE_DEALLOCATE_A(tmp_partial_b)

    POP_SUB(cpml_conv_function_update_via_e_b_fields)
  end subroutine cpml_conv_function_update_via_e_b_fields


  ! ---------------------------------------------------------
  subroutine generate_medium_boxes(hm, gr, nr_of_boxes)
    type(hamiltonian_mxll_t), intent(inout) :: hm
    type(grid_t),        intent(in)    :: gr
    integer,             intent(in)    :: nr_of_boxes

    integer :: il, ip, ip_in, ip_in_max, ip_bd, ip_bd_max, ipp, idim, point_info, err
    FLOAT   :: bounds(2,3), xx(MAX_DIM), xxp(MAX_DIM), dd, dd_max, dd_min
    FLOAT, allocatable  :: tmp(:), tmp_grad(:,:)
    character(len=256)  :: string

    PUSH_SUB(generate_medium_boxes)
 
    SAFE_ALLOCATE(tmp(gr%mesh%np_part))
    SAFE_ALLOCATE(tmp_grad(gr%mesh%np_part,1:gr%mesh%sb%dim))

    SAFE_ALLOCATE(hm%medium_box_points_number(nr_of_boxes))
    SAFE_ALLOCATE(hm%medium_box_bdry_number(nr_of_boxes))
    hm%medium_box_number = nr_of_boxes

    ip_in_max = 0
    do il=1, nr_of_boxes
      do idim=1, 3
        bounds(1,idim) = hm%medium_box_center(idim,il) - hm%medium_box_size(idim,il)/M_TWO 
        bounds(2,idim) = hm%medium_box_center(idim,il) + hm%medium_box_size(idim,il)/M_TWO
      end do
      ip_in=0
      ip_bd=0
      do ip=1, gr%mesh%np
        xx(:) = gr%mesh%x(ip,:)
        if (check_point_in_bounds(xx, bounds)) then
          ip_in = ip_in+1
        end if
        if (check_point_on_bounds(xx, bounds)) then
          ip_bd = ip_bd+1
        end if
      end do
      if (ip_in > ip_in_max) ip_in_max = ip_in
      if (ip_bd > ip_bd_max) ip_bd_max = ip_bd
      hm%medium_box_points_number(il) = ip_in
      hm%medium_box_bdry_number(il) = ip_bd
    end do

    dd_max = max(2*gr%mesh%spacing(1),2*gr%mesh%spacing(2),2*gr%mesh%spacing(3))

    SAFE_ALLOCATE(hm%medium_box_points_map(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(hm%medium_box_bdry_map(ip_bd_max,nr_of_boxes))
    SAFE_ALLOCATE(hm%medium_box_aux_ep(ip_in_max,1:3,nr_of_boxes))
    SAFE_ALLOCATE(hm%medium_box_aux_mu(ip_in_max,1:3,nr_of_boxes))
    SAFE_ALLOCATE(hm%medium_box_c(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(hm%medium_box_ep(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(hm%medium_box_mu(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(hm%medium_box_sigma_e(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(hm%medium_box_sigma_m(ip_in_max,nr_of_boxes))

    do il=1, nr_of_boxes
      ip_in=0
      ip_bd=0
      do ip=1, gr%mesh%np
        xx(:) = gr%mesh%x(ip,:)
        if (check_point_in_bounds(xx, bounds)) then
          ip_in = ip_in+1
          hm%medium_box_points_map(ip_in,il) = ip
        end if
        if (check_point_on_bounds(xx, bounds)) then
          ip_bd = ip_bd+1
          hm%medium_box_bdry_map(ip_bd,il) = ip
        end if
      end do
    end do

    do il=1, nr_of_boxes

      do ip_in=1, hm%medium_box_points_number(il)
        ip = hm%medium_box_points_map(ip_in,il)
        if (hm%medium_box_shape(il) == OPTION__MAXWELLMEDIUMBOX__SMOOTH) then
          xx(:) = gr%mesh%x(ip,:)
          dd_min = M_HUGE
          do ip_bd=1, hm%medium_box_bdry_number(il)
            ipp = hm%medium_box_bdry_map(ip_bd,il)
            xxp(:) = gr%mesh%x(ipp,:)
            dd = sqrt((xx(1)-xxp(1))**2+(xx(2)-xxp(2))**2+(xx(3)-xxp(3))**2)
            if (dd < dd_min) dd_min=dd
          end do
          hm%medium_box_ep(ip_in,il) = P_ep &
                                                     + ( ( P_ep * hm%medium_box_ep_factor(il) - P_ep )  & 
                                                         * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min-M_TWO*dd_max)) ) )
          hm%medium_box_mu(ip_in,il) = P_mu & 
                                                     + ( ( P_mu * hm%medium_box_mu_factor(il) - P_mu ) &
                                                     * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min-M_TWO*dd_max)) ) )
          hm%medium_box_c(ip_in,il)  = &
            M_ONE/sqrt(hm%medium_box_ep(ip_in,il)*hm%medium_box_mu(ip_in,il))
          hm%medium_box_sigma_e(ip_in,il) = hm%medium_box_sigma_e_factor(il) &
                                                          * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min-M_TWO*dd_max)) ) 
          hm%medium_box_sigma_m(ip_in,il) = hm%medium_box_sigma_m(ip_in,il) &
                                                          * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min-M_TWO*dd_max)) )
        else if (hm%medium_box_shape(il) == OPTION__MAXWELLMEDIUMBOX__EDGED) then
          hm%medium_box_ep(ip_in,il) = P_ep * hm%medium_box_ep_factor(il)
          hm%medium_box_mu(ip_in,il) = P_mu * hm%medium_box_mu_factor(il)
          hm%medium_box_c(ip_in,il)  = &
            M_ONE/sqrt(hm%medium_box_ep(ip_in,il)*hm%medium_box_mu(ip_in,il))
          hm%medium_box_sigma_e(ip_in,il) = hm%medium_box_sigma_e_factor(il)
          hm%medium_box_sigma_m(ip_in,il) = hm%medium_box_sigma_m(ip_in,il)
        end if
      end do

      tmp = P_ep
      do  ip_in=1, hm%medium_box_points_number(il)
        ip = hm%medium_box_points_map(ip_in,il)
        tmp(ip) =  hm%medium_box_ep(ip_in,il)
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in=1, hm%medium_box_points_number(il)
        ip = hm%medium_box_points_map(ip_in,il)
        hm%medium_box_aux_ep(ip_in,:,il) = &
          tmp_grad(ip,:)/(M_FOUR * hm%medium_box_ep(ip_in,il))
      end do

      tmp = P_mu
      do  ip_in=1, hm%medium_box_points_number(il)
        ip = hm%medium_box_points_map(ip_in,il)
        tmp(ip) =  hm%medium_box_mu(ip_in,il)
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in=1, hm%medium_box_points_number(il)
        ip = hm%medium_box_points_map(ip_in,il)
        hm%medium_box_aux_mu(ip_in,:,il) = &
          tmp_grad(ip,:)/(M_FOUR * hm%medium_box_mu(ip_in,il))
      end do

      ! print information about the medium box -- get from Rene's version in maxwell_propagator.F90

      SAFE_DEALLOCATE_A(tmp)
    end do

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp_grad)

    POP_SUB(generate_medium_boxes)

    contains

      logical pure function check_point_in_bounds(xx, bounds) result (check)
        FLOAT, intent(in) :: xx(:)
        FLOAT, intent(in) :: bounds(:,:)

        check = .false.
        if ( (xx(1)>=bounds(1,1)) .and. (xx(1)<=bounds(2,1)) .and. &
             (xx(2)>=bounds(1,2)) .and. (xx(2)<=bounds(2,2)) .and. &
             (xx(3)>=bounds(1,3)) .and. (xx(3)<=bounds(2,3)) ) then
          check = .true.
        end if

      end function check_point_in_bounds
 
      logical pure function check_point_on_bounds(xx, bounds) result (check)
        FLOAT, intent(in) :: xx(:)
        FLOAT, intent(in) :: bounds(:,:)

        check = .false.
        if ( xx(1)==bounds(1,1) .and. (xx(2)>=bounds(1,2) .and. xx(3)>=bounds(1,3)) &
                                .and. (xx(2)<=bounds(2,2) .and. xx(3)<=bounds(2,3)) .or. &
             xx(2)==bounds(1,2) .and. (xx(1)>=bounds(1,1) .and. xx(3)>=bounds(1,3)) &
                                .and. (xx(1)<=bounds(2,1) .and. xx(3)<=bounds(2,3)) .or. &
             xx(3)==bounds(1,3) .and. (xx(1)>=bounds(1,1) .and. xx(2)>=bounds(1,2)) &
                                .and. (xx(1)<=bounds(2,1) .and. xx(2)<=bounds(2,2)) .or. &
             xx(1)==bounds(2,1) .and. (xx(2)>=bounds(1,2) .and. xx(3)>=bounds(1,3)) &
                                .and. (xx(2)<=bounds(2,2) .and. xx(3)<=bounds(2,3)) .or. &
             xx(2)==bounds(2,2) .and. (xx(1)>=bounds(1,1) .and. xx(3)>=bounds(1,3)) &
                                .and. (xx(1)<=bounds(2,1) .and. xx(3)<=bounds(2,3)) .or. &
             xx(3)==bounds(2,3) .and. (xx(1)>=bounds(1,1) .and. xx(2)>=bounds(1,2)) &
                                .and. (xx(1)<=bounds(2,1) .and. xx(2)<=bounds(2,2)) ) then
          check = .true.
        end if

      end function check_point_on_bounds

      subroutine get_medium_io_function(medium_func, hm, mesh, il, io_func)
        FLOAT,               intent(in)    :: medium_func(:)
        type(hamiltonian_mxll_t), intent(in)    :: hm
        type(mesh_t),        intent(in)    :: mesh
        integer,             intent(in)    :: il
        FLOAT,               intent(inout) :: io_func(:)

        integer :: ip, ip_in

        do ip_in=1, hm%medium_box_points_number(il)
          ip          = hm%medium_box_points_map(ip_in,il)
          io_func(ip) = medium_func(ip_in)
        end do

      end subroutine get_medium_io_function

  end subroutine generate_medium_boxes


  ! ! -------------------------------------------------------
  ! subroutine maxwell_matter_mesh_mapping(gr, gr_mxll)
  !   type(grid_t),         intent(inout) :: gr
  !   type(grid_t),         intent(inout) :: gr_mxll

  !   integer :: ip, idim, ip_ma_local, ip_mx_local, ip_ma_global, ip_mx_global, idx
  !   integer :: ipart, mx_np_part, mx_npart, ma_np_part, ma_npart, ma_comm, mx_comm
  !   integer, allocatable :: ip_ma_local_vec(:), ip_mx_local_vec(:)
  !   FLOAT   :: dmin_ma, dmin_mx, xx(3)
  !   type(mpi_grp_t) :: ma_grp, mx_grp

  !   !USH_SUB(maxwell_matter_mesh_mapping)

  !   ma_np_part = gr%mesh%np_part
  !   ma_npart   = gr%mesh%vp%npart
  !   mx_np_part = gr_mxll%mesh%np_part
  !   mx_npart   = gr_mxll%mesh%vp%npart

  !   !AFE_ALLOCATE(ip_mx_local_vec(gr_mxll%mesh%vp%npart))
  !   !AFE_ALLOCATE(ip_ma_local_vec(gr%mesh%vp%npart))
  !   !AFE_ALLOCATE(gr%mesh%ma_mx_mesh_mapping%global1_to_global2_map(1:gr%mesh%np_part_global))
  !   !AFE_ALLOCATE(gr%mesh%ma_mx_mesh_mapping%local1_overlap(1:gr%mesh%np_part))
  !   !AFE_ALLOCATE(gr%mesh%ma_mx_mesh_mapping%global2_overlap(1:gr%mesh%np_part))
  !   !AFE_ALLOCATE(gr%mesh%ma_mx_mesh_mapping%rank_map(1:gr%mesh%np_part))
  !   !AFE_ALLOCATE(gr%mesh%ma_mx_mesh_mapping%local1_to_local2_map(1:ma_np_part,1:ma_npart))
  !   !AFE_ALLOCATE(gr_mxll%mesh%mx_ma_mesh_mapping%global1_to_global2_map(1:gr_mxll%mesh%np_part_global))
  !   !AFE_ALLOCATE(gr_mxll%mesh%mx_ma_mesh_mapping%local1_overlap(1:gr_mxll%mesh%np_part))
  !   !AFE_ALLOCATE(gr_mxll%mesh%mx_ma_mesh_mapping%global2_overlap(1:gr_mxll%mesh%np_part))
  !   !AFE_ALLOCATE(gr_mxll%mesh%mx_ma_mesh_mapping%rank_map(1:gr_mxll%mesh%np_part))
  !   !AFE_ALLOCATE(gr_mxll%mesh%mx_ma_mesh_mapping%local1_to_local2_map(1:mx_np_part,1:mx_npart))

  !   ! mapping global points from global matter mesh to global Maxwell mesh
  !   do ip_ma_global=1, gr%mesh%np_part_global
  !     if (check_point_overlap(gr%mesh, gr_mxll%mesh, ip_ma_global)) then
  !       ip_mx_global = gr_mxll%mesh%idx%lxyz_inv(gr%mesh%idx%lxyz(ip_ma_global,1), &
  !                                               gr%mesh%idx%lxyz(ip_ma_global,2), &
  !                                               gr%mesh%idx%lxyz(ip_ma_global,3))
  !       gr%mesh%ma_mx_mesh_mapping%global1_to_global2_map(ip_ma_global) = ip_mx_global
  !     else
  !       gr%mesh%ma_mx_mesh_mapping%global1_to_global2_map(ip_ma_global) = 0
  !     end if
  !   end do

  !   ! mapping global points from global Maxwell mesh to global matter mesh
  !   do ip_mx_global=1, gr_mxll%mesh%np_part_global
  !     if (check_point_overlap(gr_mxll%mesh, gr%mesh, ip_mx_global)) then
  !       ip_ma_global = gr%mesh%idx%lxyz_inv(gr_mxll%mesh%idx%lxyz(ip_mx_global,1), &
  !                                           gr_mxll%mesh%idx%lxyz(ip_mx_global,2), &
  !                                           gr_mxll%mesh%idx%lxyz(ip_mx_global,3))
  !       gr_mxll%mesh%mx_ma_mesh_mapping%global1_to_global2_map(ip_mx_global) = ip_ma_global
  !     else
  !       gr_mxll%mesh%mx_ma_mesh_mapping%global1_to_global2_map(ip_mx_global) = 0
  !     end if
  !   end do

  !   ! mapping local points from local matter mesh to global Maxwell mesh
  !   idx=0
  !   do ip_ma_local=1, gr%mesh%np_part
  !     if (gr%mesh%parallel_in_domains) then
  !       ip_ma_global = vec_local2global_part(gr%mesh%vp, ip_ma_local, gr%mesh%vp%partno)
  !     else
  !       ip_ma_global = ip_ma_local
  !     end if
  !     if (check_point_overlap(gr%mesh, gr_mxll%mesh, ip_ma_global)) then
  !       if (gr%mesh%ma_mx_mesh_mapping%global1_to_global2_map(ip_ma_global) /= 0) then
  !         idx = idx+1
  !         ! initializing local1_overlap(idx)
  !         gr%mesh%ma_mx_mesh_mapping%local1_overlap(idx) = ip_ma_local
  !         ! initializing global2_overlap(idx)
  !         gr%mesh%ma_mx_mesh_mapping%global2_overlap(idx) &
  !           = gr%mesh%ma_mx_mesh_mapping%global1_to_global2_map(ip_ma_global)
  !       end if
  !     end if
  !   end do
  !   gr%mesh%ma_mx_mesh_mapping%local1_overlap_number = idx

  !   ! mapping local points from local Maxwell mesh to global matter mesh
  !   idx=0
  !   do ip_mx_local=1, gr_mxll%mesh%np_part
  !     if (gr_mxll%mesh%parallel_in_domains) then
  !       ip_mx_global = vec_local2global_part(gr_mxll%mesh%vp, ip_mx_local, gr_mxll%mesh%vp%partno)
  !     else
  !       ip_mx_global = ip_mx_local
  !     end if
  !     if (check_point_overlap(gr_mxll%mesh, gr%mesh, ip_mx_global)) then
  !       if (gr_mxll%mesh%mx_ma_mesh_mapping%global1_to_global2_map(ip_mx_global) /= 0) then
  !         idx = idx+1
  !         ! initializing local1_overlap(idx)
  !         gr_mxll%mesh%mx_ma_mesh_mapping%local1_overlap(idx) = ip_mx_local
  !         ! initializing global2_overlap(idx)
  !         gr_mxll%mesh%mx_ma_mesh_mapping%global2_overlap(idx) &
  !           = gr_mxll%mesh%mx_ma_mesh_mapping%global1_to_global2_map(ip_mx_global)
  !       end if
  !     end if
  !   end do
  !   gr_mxll%mesh%mx_ma_mesh_mapping%local1_overlap_number = idx

  !   !OP_SUB(maxwell_matter_mesh_mapping)

  !   contains

  !     logical function check_point_overlap(mesh_1, mesh_2, ip_1) result(check)
  !       type(mesh_t), intent(in) :: mesh_1
  !       type(mesh_t), intent(in) :: mesh_2
  !       integer,      intent(in) :: ip_1

  !       if ( (mesh_1%idx%lxyz(ip_1,1) >= mesh_2%idx%nr(1,1)) .and. &
  !            (mesh_1%idx%lxyz(ip_1,1) <= mesh_2%idx%nr(2,1)) .and. &
  !            (mesh_1%idx%lxyz(ip_1,2) >= mesh_2%idx%nr(1,2)) .and. &
  !            (mesh_1%idx%lxyz(ip_1,2) <= mesh_2%idx%nr(2,2)) .and. &
  !            (mesh_1%idx%lxyz(ip_1,3) >= mesh_2%idx%nr(1,3)) .and. &
  !            (mesh_1%idx%lxyz(ip_1,3) <= mesh_2%idx%nr(2,3)) ) then
  !         check = .true.
  !       else
  !         check = .false.
  !       end if

  !     end function check_point_overlap

  !   end subroutine maxwell_matter_mesh_mapping


  ! ! ---------------------------------------------------------
  ! subroutine get_mx_ma_coupling_points(geo, hm, mx_ma_coupling_points)
  !   type(geometry_t),    intent(in)    :: geo
  !   type(hamiltonian_elec_t), intent(in)    :: hm
  !   FLOAT,               intent(inout) :: mx_ma_coupling_points(:,:)

  !   integer :: iatom

  !   !USH_SUB(get_mx_ma_coupling_points)

  !   if (hm%mx_ma_coupling_type == OPTION__MAXWELLCOUPLINGMETHOD__MULTIPOLE_EXPANSION_CENTER_OF_MASS) then
  !     call find_center_of_mass(geo, mx_ma_coupling_points(:,1), .false.)
  !   else if (hm%mx_ma_coupling_type == OPTION__MAXWELLCOUPLINGMETHOD__MULTIPOLE_EXPANSION_IONS) then
  !     do iatom=1, geo%natoms
  !       mx_ma_coupling_points(:,iatom) = geo%atom(iatom)%x(:)
  !     end do
  !   end if

  !   !OP_SUB(get_mx_ma_coupling_points)
  ! end subroutine get_mx_ma_coupling_points


  ! ! ---------------------------------------------------------
  ! subroutine maxwell_grid_points_coupling_points_mapping(gr, mx_ma_coupling_points, mx_ma_coupling_points_number, &
  !                                                        mx_ma_coupling_map)
  !   type(grid_t),     intent(inout) :: gr
  !   FLOAT,            intent(in)    :: mx_ma_coupling_points(:,:)
  !   integer,          intent(in)    :: mx_ma_coupling_points_number
  !   integer,          intent(inout) :: mx_ma_coupling_map(:)

  !   integer :: ip_global, ic, idim, ic_no
  !   FLOAT   :: dd, dd_min, xx(MAX_DIM)

  !   !USH_SUB(maxwell_grid_points_coupling_points_mapping)

  !   do ip_global=1, gr%mesh%np_global
  !     do ic=1, mx_ma_coupling_points_number
  !       dd = M_ZERO
  !       do idim=1, gr%sb%dim
  !         xx(idim) = gr%mesh%idx%lxyz(ip_global,idim) * gr%mesh%spacing(idim)
  !         dd = dd + ( xx(idim) - mx_ma_coupling_points(idim,ic) )**2
  !       end do
  !       dd = sqrt(dd)
  !       if (ic == 1) then
  !         dd_min = dd
  !         ic_no = ic
  !       else
  !         if ( dd < dd_min ) then
  !           dd_min = dd
  !           ic_no = ic
  !         end if
  !       end if
  !     end do
  !     mx_ma_coupling_map(ip_global) = ic_no
  !   end do

  !   !OP_SUB(maxwell_grid_points_coupling_points_mapping)
  ! end subroutine maxwell_grid_points_coupling_points_mapping


  ! ---------------------------------------------------------
  subroutine td_function_mxll_init(st, namespace, gr, hm)
    type(states_mxll_t),      intent(inout) :: st
    type(namespace_t),   intent(in)    :: namespace
    type(grid_t),        intent(in)    :: gr
    type(hamiltonian_mxll_t), intent(in)    :: hm

    type(block_t)        :: blk
    integer              :: il, nlines, idim, ncols, ierr
    FLOAT                :: e_field(MAX_DIM), b_field(MAX_DIM)
    character(len=1024)  :: mxf_expression

    PUSH_SUB(td_function_init)

    !%Variable UserDefinedConstantSpacialMaxwellField
    !%Type block
    !%Section MaxwellStates
    !%Description
    !% Follows
    !%
    !% Example:
    !%
    !% <tt>%UserDefinedConstantSpacialMaxwellFields
    !% <br>&nbsp;&nbsp;   plane_wave_parser      | E_x | E_y | E_z | B_x | B_y | B_z | "tdf_function"
    !% <br>%</tt>
    !%
    !% Description about UserDefinedConstantSpacialMaxwellField follows
    !%
    !%End

    if(parse_block(namespace, 'UserDefinedConstantSpacialMaxwellField', blk) == 0) then
      st%rs_state_const_external = .true.
      nlines = parse_block_n(blk)
      SAFE_ALLOCATE(st%rs_state_const_td_function(nlines))
      SAFE_ALLOCATE(st%rs_state_const_amp(MAX_DIM,nlines))
      ! read all lines
      do il = 1, nlines
        e_field = M_ZERO
        b_field = M_ZERO
        ! Check that number of columns is five or six.
        ncols = parse_block_cols(blk, il - 1)
        if (ncols  /= 7) then
          message(1) = 'Each line in the UserDefinedConstantSpacialMaxwellField block must have'
          message(2) = 'seven columns.'
          call messages_fatal(2)
        end if
        do idim=1, st%d%dim
          call parse_block_float( blk, il - 1, idim-1, e_field(idim))
        end do
        do idim=1, st%d%dim
          call parse_block_float( blk, il - 1, idim+2, b_field(idim))
        end do
        call parse_block_string( blk, il - 1, 6, mxf_expression)
        call build_rs_vector(e_field, b_field, st%rs_sign, st%rs_state_const_amp(:,il))
        call tdf_read(st%rs_state_const_td_function(il), namespace, trim(mxf_expression), ierr)
      end do
    end if

    POP_SUB(td_function_mxll_init)
  end subroutine td_function_mxll_init


  ! ---------------------------------------------------------
  subroutine spatial_constant_calculation(constant_calc, st, gr, hm, time, dt, delay, rs_state)
    logical,             intent(in)    :: constant_calc
    type(states_mxll_t),      intent(in)    :: st
    type(grid_t),        intent(in)    :: gr
    type(hamiltonian_mxll_t), intent(in)    :: hm
    FLOAT,               intent(in)    :: time
    FLOAT,               intent(in)    :: dt
    FLOAT,               intent(in)    :: delay
    CMPLX,               intent(inout) :: rs_state(:,:)

    integer :: ip, ic, icn
    FLOAT   :: tf_old, tf_new 

    PUSH_SUB(spatial_constant_calculation)

    if (hm%spatial_constant_apply) then
      if (constant_calc) then
        icn = size(st%rs_state_const_td_function(:))
        st%rs_state_const(:) = M_z0
        do ic=1, icn
          tf_old = tdf(st%rs_state_const_td_function(ic),time-delay-dt)
          tf_new = tdf(st%rs_state_const_td_function(ic),time-delay)
          forall (ip=1:gr%mesh%np_part)
            rs_state(ip,:) = rs_state(ip,:) + st%rs_state_const_amp(:,ic) * (tf_new-tf_old)
          end forall
          st%rs_state_const(:) = st%rs_state_const(:) + st%rs_state_const_amp(:,ic)
        end do
        st%rs_state_const(:) = st%rs_state_const(:)*tf_new
      end if
    end if

    POP_SUB(spatial_constant_calculation)
  end subroutine spatial_constant_calculation


  ! ---------------------------------------------------------
  subroutine constant_boundaries_calculation(constant_calc, bc, hm, st, rs_state)
    logical,              intent(in)    :: constant_calc
    type(bc_mxll_t),   intent(in)    :: bc
    type(states_mxll_t),       intent(in)    :: st
    type(hamiltonian_mxll_t),  intent(in)    :: hm
    CMPLX,                intent(inout) :: rs_state(:,:)

    integer :: ip_in, ip

    PUSH_SUB(constant_boundaries_calculation)

    if (hm%spatial_constant_apply) then
      if (constant_calc) then
        do ip_in=1, bc%constant_points_number
          ip = bc%constant_points_map(ip_in)
          rs_state(ip,:) = st%rs_state_const(:)
          bc%constant_rs_state(ip_in,:) = st%rs_state_const(:)
        end do
      end if
    end if    

    POP_SUB(constant_boundaries_calculation)
  end subroutine constant_boundaries_calculation


  ! ---------------------------------------------------------
  subroutine mirror_pec_boundaries_calculation(bc, st, rs_state)
    type(bc_mxll_t)                :: bc
    type(states_mxll_t)            :: st
    CMPLX                          :: rs_state(:,:)

    integer                    :: ip, ip_in, idim
    FLOAT                      :: e_field(3), b_field(3)

    PUSH_SUB(mirror_pec_boundaries_calculation)

    do idim = 1, 3
      if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MIRROR_PEC) then
        do ip_in = 1, bc%mirror_points_number(idim)
          ip         = bc%mirror_points_map(ip_in,idim)
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
    type(bc_mxll_t)                :: bc
    type(states_mxll_t)            :: st
    CMPLX                          :: rs_state(:,:)

    integer                    :: ip, ip_in, idim
    FLOAT                      :: e_field(3), b_field(3)

    PUSH_SUB(mirror_pmc_boundaries_calculation)

    do idim=1, 3
      if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MIRROR_PMC) then
        do ip_in=1, bc%mirror_points_number(idim)
          ip         = bc%mirror_points_map(ip_in,idim)
          b_field(:) = M_ZERO
          call get_electric_field_vector(rs_state(ip,:), e_field(:), st%ep(ip))
          call build_rs_vector(e_field(:), b_field(:), st%rs_sign, rs_state(ip,:), st%ep(ip), st%mu(ip))
        end do
      end if
    end do

    POP_SUB(mirror_pmc_boundaries_calculation)
  end subroutine mirror_pmc_boundaries_calculation


  ! ---------------------------------------------------------
  subroutine plane_waves_boundaries_calculation(hm, tr, st, mesh, time, time_delay, rs_state)
    type(hamiltonian_mxll_t)   :: hm
    type(propagator_mxll_t)    :: tr
    type(states_mxll_t)        :: st
    type(mesh_t)               :: mesh
    FLOAT                      :: time
    FLOAT                      :: time_delay
    CMPLX                      :: rs_state(:,:)

    integer                    :: ip, ip_global, ip_in, wn, oam, sam
    FLOAT                      :: dd, x_prop(MAX_DIM), rr, vv(MAX_DIM), k_vector(MAX_DIM), k_vector_abs, nn, rad, phi, test_perp, test_limit
    FLOAT                      :: width, shift, e0(MAX_DIM), e_field(MAX_DIM), dummy(MAX_DIM), b0(MAX_DIM), b_field(MAX_DIM)
    CMPLX                      :: rs_state_add(MAX_DIM), unit_sigma(3), unit_z(3), unit_plus(3), unit_minus(3)
    CMPLX                      :: e0cmplx(MAX_DIM), b0cmplx(MAX_DIM)

    PUSH_SUB(plane_waves_boundaries_calculation)

    test_limit = CNST(1e-9)

    if (hm%plane_waves_apply) then
      do wn = 1, hm%bc%plane_waves_number
        do ip_in = 1, hm%bc%plane_waves_points_number
          ip = hm%bc%plane_waves_points_map(ip_in)
            if (wn == 1) rs_state(ip,:) = M_Z0
            nn           = sqrt(st%ep(ip)/P_ep*st%mu(ip)/P_mu)
            vv(:)        = hm%bc%plane_waves_v_vector(:,wn)
            k_vector(:)  = hm%bc%plane_waves_k_vector(:,wn)
            k_vector_abs = sqrt(sum(k_vector(:)**2))
            x_prop(:)    = mesh%x(ip,:) - vv(:) * (time - time_delay)
            rr           = sqrt(sum(x_prop(:)**2))
            if (hm%bc%plane_waves_modus(wn) == & 
                OPTION__USERDEFINEDMAXWELLINCIDENTWAVES__PLANE_WAVE_MX_FUNCTION) then
              e0(:)      = hm%bc%plane_waves_e_field(:,wn)
              e_field(:) = e0(:) * mxf(hm%bc%plane_waves_mx_function(wn), x_prop(:)) 
            end if
            b_field(1:3) = M_ONE/P_c * M_ONE/k_vector_abs * dcross_product(k_vector,e_field)
            call build_rs_vector(e_field, b_field, st%rs_sign, rs_state_add, st%ep(ip), st%mu(ip))
            rs_state(ip,:) =  rs_state(ip,:) + rs_state_add(:)
          end do
        end do
      else
        do ip_in = 1, hm%bc%plane_waves_points_number
          ip             = hm%bc%plane_waves_points_map(ip_in)
          rs_state(ip,:) = M_z0
        end do
      end if

    POP_SUB(plane_waves_boundaries_calculation)
  end subroutine plane_waves_boundaries_calculation


  ! ---------------------------------------------------------
  subroutine plane_waves_propagation(hm, st, gr, tr, time, dt, time_delay)
    type(hamiltonian_mxll_t)        :: hm
    type(states_mxll_t)             :: st
    type(grid_t)               :: gr
    type(propagator_mxll_t) :: tr
    FLOAT                      :: time
    FLOAT                      :: dt
    FLOAT                      :: time_delay

    integer            :: ip_in, ip, ff_dim
    CMPLX, allocatable :: ff_rs_state(:,:)

    PUSH_SUB(plane_waves_propagation)

    if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then
      ff_dim = 6
    else if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_GAUSS) then
      ff_dim = 4
    else 
      ff_dim = 3
    end if

    SAFE_ALLOCATE(ff_rs_state(1:gr%mesh%np_part,ff_dim))

    call transform_rs_state_forward(hm, gr, st, st%rs_state_plane_waves, ff_rs_state)

    ! Time evolution of RS plane waves state without any coupling with H(inter_time)
    call hamiltonian_mxll_update(hm, time=time)
    !call exponential_apply(tr%te, gr%mesh, hm, ff_rs_state, 1, 1, dt, time)
    call transform_rs_state_backward(hm, gr, st, ff_rs_state, st%rs_state_plane_waves)
    call plane_waves_boundaries_calculation(hm, tr, st, gr%mesh, time+dt, time_delay, st%rs_state_plane_waves)

    SAFE_DEALLOCATE_A(ff_rs_state)

    POP_SUB(plane_waves_propagation)
  end subroutine plane_waves_propagation


  ! ---------------------------------------------------------
  subroutine plane_waves_in_box_calculation(bc, time, gr, st, hm, rs_state)
    type(bc_mxll_t)  , intent(in)    :: bc
    FLOAT,                intent(in)    :: time
    type(grid_t),         intent(in)    :: gr
    type(states_mxll_t),       intent(in)    :: st
    type(hamiltonian_mxll_t),  intent(in)    :: hm
    CMPLX,                intent(inout) :: rs_state(:,:)

    integer                    :: ip, ip_global, ip_in, wn, oam, sam, idim, np
    FLOAT                      :: dd, x_prop(MAX_DIM), rr, vv(MAX_DIM), k_vector(MAX_DIM), k_vector_abs, nn, rad, phi
    FLOAT                      :: width, shift, e0(MAX_DIM), e_field(MAX_DIM), dummy(MAX_DIM), b0(MAX_DIM), b_field(MAX_DIM)
    CMPLX                      :: rs_state_add(MAX_DIM), unit_sigma(3), unit_z(3), unit_plus(3), unit_minus(3)
    CMPLX                      :: e0cmplx(MAX_DIM), b0cmplx(MAX_DIM)

    PUSH_SUB(plane_waves_in_box_calculation)
    
    np            = gr%mesh%np_part
    !rs_state(:,:) = M_z0
    do wn=1, bc%plane_waves_number
      vv(:)        = hm%bc%plane_waves_v_vector(:,wn)
      k_vector(:)  = hm%bc%plane_waves_k_vector(:,wn)
      k_vector_abs = sqrt(sum(k_vector(:)**2))
      e0(:)      = hm%bc%plane_waves_e_field(:,wn)
      if (hm%bc%bessel_beam) then
        do ip = 1, gr%mesh%np
          if (wn==1) rs_state(ip,:) = M_Z0
          nn           = sqrt(st%ep(ip)/P_ep*st%mu(ip)/P_mu)
!          if (hm%bc%plane_waves_mx_function(wn)%mode == MXF_BESSEL_E_FIELD) then
!            e_field(1) = e0(1)*mxf(hm%bc%plane_waves_mx_function(wn), gr%mesh%x(ip,:), 1, gr%mesh%x(ip,:), time)
!            e_field(2) = e0(2)*mxf(hm%bc%plane_waves_mx_function(wn), gr%mesh%x(ip,:), 2, gr%mesh%x(ip,:), time)
!            e_field(3) = e0(3)*mxf(hm%bc%plane_waves_mx_function(wn), gr%mesh%x(ip,:), 3, gr%mesh%x(ip,:), time)
!            b_field(:) = M_ZERO
!          end if
!          if (hm%bc%plane_waves_mx_function(wn)%mode == MXF_BESSEL_B_FIELD) then
!            b_field(1) = e0(1)*k_vector(2)/(P_c*sqrt(k_vector(1)**2+k_vector(2)**2)) &
!                         * mxf(hm%bc%plane_waves_mx_function(wn), gr%mesh%x(ip,:), 1, gr%mesh%x(ip,:), time)
!            b_field(2) = e0(2)*k_vector(2)/(P_c*sqrt(k_vector(1)**2+k_vector(2)**2)) &
!                         * mxf(hm%bc%plane_waves_mx_function(wn), gr%mesh%x(ip,:), 2, gr%mesh%x(ip,:), time)
!            b_field(3) = e0(3)*k_vector(2)/(P_c*sqrt(k_vector(1)**2+k_vector(2)**2)) &
!                         * mxf(hm%bc%plane_waves_mx_function(wn), gr%mesh%x(ip,:), 3, gr%mesh%x(ip,:), time)
!            e_field(:) = M_ZERO
!          end if
          call build_rs_vector(e_field, b_field, st%rs_sign, rs_state_add, st%ep(ip), st%mu(ip))
          rs_state(ip,:) =  rs_state(ip,:) + rs_state_add(:)
        end do
      else
        do ip = 1, gr%mesh%np
          if (wn==1) rs_state(ip,:) = M_Z0
          nn           = sqrt(st%ep(ip)/P_ep*st%mu(ip)/P_mu)
          x_prop       = M_ZERO
          !x_prop(:)    = gr%mesh%x(ip,:) - vv(:) * time
          rr           = sqrt(sum(x_prop(:)**2))
          if (bc%plane_waves_modus(wn) == OPTION__USERDEFINEDMAXWELLINCIDENTWAVES__PLANE_WAVE_PARSER) then
            do idim=1, gr%mesh%sb%dim
               call parse_expression(e_field(idim), dummy(idim), gr%mesh%sb%dim, x_prop, rr, M_ZERO, &
                 bc%plane_waves_e_field_string(idim,wn))
               e_field(idim) = units_to_atomic(units_inp%energy/units_inp%length, e_field(idim))
            end do
          else if (bc%plane_waves_modus(wn) == & 
            OPTION__USERDEFINEDMAXWELLINCIDENTWAVES__PLANE_WAVE_MX_FUNCTION) then
            e0(:)      = bc%plane_waves_e_field(:,wn)
            e_field(:) = e0(:) * mxf(bc%plane_waves_mx_function(wn), x_prop(:)) 
          end if
          b_field = 1/P_c * M_ONE/k_vector_abs * dcross_product(k_vector,e_field)
          call build_rs_vector(e_field, b_field, st%rs_sign, rs_state_add, st%ep(ip), st%mu(ip))
          rs_state(ip,:) = rs_state(ip,:) + rs_state_add(:)
        end do
      end if
    end do

    POP_SUB(plane_waves_in_box_calculation)
  end subroutine plane_waves_in_box_calculation


  ! ---------------------------------------------------------
  subroutine complex_matrix_vector_multiplication(mm, nn, matrix, vector, vector_r)
    integer, intent(in)    :: mm, nn
    CMPLX,   intent(inout) :: matrix(:,:), vector(:)
    CMPLX,   intent(inout) :: vector_r(:)

    CMPLX   :: tmp_sum
    integer :: idx_m, idx_n

    do idx_m=1, mm
      tmp_sum = M_z0
      do idx_n=1, nn
        tmp_sum = tmp_sum + matrix(idx_m,idx_n) * vector(idx_n)
      end do
      vector_r(idx_m) = tmp_sum
    end do

  end subroutine complex_matrix_vector_multiplication


end module propagator_mxll_oct_m

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
  public ::                                  &
    propagator_mxll_init,                    &
    mxll_propagation_step,                   &
    transform_rs_state,                      &
    rs_state_to_cube_map,                    &
    calculate_matter_longitudinal_field,     &
    get_vector_pot_and_transverse_field,     &
    energy_density_calc,                     &
    energy_mxll_calc,                        &
    energy_rate_calc,                        &
    poynting_vector_through_box_surfaces,    &
    poynting_vector_through_box_surfaces_plane_waves, &
    fields_through_box_surfaces,             &
    fields_through_box_surfaces_plane_waves, &
    constant_boundaries_calculation,         &
    maxwell_mask,                            &
    td_function_mxll_init,                   &
    plane_waves_in_box_calculation

  integer, public, parameter ::   &
     RS_TRANS_FORWARD = 1,        &
     RS_TRANS_BACKWARD = 2


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
        select case (hm%bc%bc_type(icol))
          case (MXLL_BC_ZERO)
          string = 'Zero'
          tr%bc_zero = .true.
        case (MXLL_BC_CONSTANT)
          string = 'Constant'
          tr%bc_constant = .true.
          tr%bc_add_ab_region = .true.
          hm%bc_constant = .true.
          hm%bc_add_ab_region = .true.
          SAFE_ALLOCATE(st%rs_state_const(1:st%dim))
          st%rs_state_const = M_z0
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
        end select
        write(message(1),'(a,I1,a,a)') 'Maxwell boundary condition in direction ', icol, ': ', trim(string)
        call messages_info(1)
        if (plane_waves_set .and. .not. (parse_is_defined(namespace, 'MaxwellIncidentWaves')) ) then
          write(message(1),'(a)') 'Input: Maxwell boundary condition option is set to "plane_waves".'
          write(message(2),'(a)') 'Input: User defined Maxwell plane waves have to be defined!'
          call messages_fatal(2)
        end if
      end do

     call parse_block_end(blk)
     call messages_print_stress(stdout, namespace=namespace)

    end if

    !%Variable MaxwellMediumBox
    !%Type block
    !%Section Time-Dependent::Propagation
    !%Description
    !% Defines parameters for the box.
    !%
    !% Example:
    !%
    !% <tt>%MaxwellMediumBox
    !% <br>&nbsp;&nbsp;   center_x | center_y | center_z | x_length | y_length | z_length | epsilon_factor | mu_factor | sigma_e | sigma_m | edged/smooth
    !% <br>%</tt>
    !%
    !% Position of center (three components) and length (three components), followed by permittivity
    !% factor, electric conductivity and magnetic conductivity, and finally type of edge.
    !%
    !%Option edged 1
    !% Medium box edges are steep.
    !%Option smooth 2
    !% Medium box edged and softened.
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
    !% Whether to perform  apriximations to the ETRS propagator.
    !%Option no 0
    !% No approximations.
    !%Option const_steps 1
    !% Use constant current density..
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
    !% not doing any numerical propagation of Maxwells equations.
    !%End
    call parse_variable(namespace, 'MaxwellPlaneWavesInBox', .false., tr%plane_waves_in_box)
    call set_medium_rs_state(st, gr, hm)

    call derivatives_boundary_mask(hm%bc, gr%mesh, hm)

    !tr%te%exp = .true.
    call exponential_init(tr%te, namespace) ! initialize Maxwell propagator

    POP_SUB(propagator_mxll_init)
  end subroutine propagator_mxll_init


  ! ---------------------------------------------------------
  subroutine mxll_propagation_step(hm, namespace, gr, st, tr, rs_state, time, dt)
    type(hamiltonian_mxll_t),        intent(inout) :: hm
    type(namespace_t),           intent(in)    :: namespace
    type(grid_t),               intent(inout) :: gr
    type(states_mxll_t),             intent(inout) :: st
    type(propagator_mxll_t), intent(inout) :: tr
    CMPLX,                      intent(inout) :: rs_state(:,:)
    FLOAT,                      intent(in)    :: time
    FLOAT,                      intent(in)    :: dt

    integer            :: ii, inter_steps, ff_dim, idim
    FLOAT              :: inter_dt, inter_time, delay
    CMPLX, allocatable :: ff_rs_state(:,:)
    CMPLX, allocatable :: ff_rs_state_pml(:,:)
    logical            :: pml_check = .false.

    PUSH_SUB(mxll_propagation_step)



    if (tr%plane_waves_in_box) then
      call plane_waves_in_box_calculation(hm%bc, time+dt, gr, st, hm, rs_state)
      POP_SUB(mxll_propagation_step)
      return
    end if

    do idim=1, 3
      if (hm%bc%bc_ab_type(idim) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML) then
        pml_check = .true.
      end if
    end do

    if (hm%operator == FARADAY_AMPERE_MEDIUM) then
      ff_dim = 6
    else if (hm%operator == FARADAY_AMPERE_GAUSS) then
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

    end if

    do ii=1, inter_steps

      ! intermediate time
      inter_time = time + inter_dt * (ii-1)

      ! transformation of RS state into 3x3 or 4x4 representation
      call transform_rs_state(hm, rs_state, ff_rs_state, RS_TRANS_FORWARD)

      if ((hm%ma_mx_coupling_apply) .or. hm%current_density_ext_flag) then

        message(1) = "Maxwell-matter coupling or external current not implemented yet"
        call messages_fatal(1)

      else

        ! RS state propagation
        call hamiltonian_mxll_update(hm, time=inter_time)
        if (pml_check) then
          call pml_propagation_stage_1(hm, gr, st, tr, ff_rs_state, ff_rs_state_pml)
          hm%cpml_hamiltonian = .true.
        end if
        call exponential_mxll_apply(hm, namespace, gr, st, tr, inter_dt, ff_rs_state)

        if (pml_check) then
          hm%cpml_hamiltonian = .false.
          call pml_propagation_stage_2(hm, namespace, gr, st, tr, inter_time, inter_dt, delay,&
            ff_rs_state_pml, ff_rs_state)
        end if
      end if

      ! PML convolution function update
      if (pml_check) then
        call cpml_conv_function_update(hm, gr, ff_rs_state_pml, ff_dim)
      end if

      ! back tranformation of RS state representation
      call transform_rs_state(hm, rs_state, ff_rs_state, RS_TRANS_BACKWARD)

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
        call plane_waves_boundaries_calculation(hm, st, gr%mesh, inter_time+inter_dt, delay ,rs_state)
      end if

    end do

    SAFE_DEALLOCATE_A(ff_rs_state)

    if (pml_check) then
      SAFE_DEALLOCATE_A(ff_rs_state_pml)
    end if

    POP_SUB(mxll_propagation_step)
  end subroutine mxll_propagation_step


  ! ---------------------------------------------------------
  subroutine exponential_mxll_apply(hm, namespace, gr, st, tr, dt, ff)
    type(hamiltonian_mxll_t),   intent(inout) :: hm
    type(namespace_t),          intent(in)    :: namespace
    type(grid_t),               intent(in)    :: gr
    type(states_mxll_t),        intent(inout) :: st
    type(propagator_mxll_t),    intent(inout) :: tr
    FLOAT,                      intent(in)    :: dt
    CMPLX,                      intent(inout) :: ff(:,:)

    type(batch_t) :: ffbatch
    integer :: istate

    PUSH_SUB(exponential_mxll_apply)

    call zbatch_init(ffbatch, 1, 1, st%dim, gr%mesh%np_part)

    do istate = 1, st%dim
      call batch_set_state(ffbatch, istate, gr%mesh%np_part, ff(:, istate))
    end do

    if (st%pack_states) call ffbatch%do_pack()

    call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, ffbatch, dt)

    do istate = 1, st%dim
      call batch_get_state(ffbatch, istate, gr%mesh%np_part, ff(:, istate))
    end do

    call ffbatch%end()

    POP_SUB(exponential_mxll_apply)
  end subroutine exponential_mxll_apply


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

    do idim=1, st%dim
      if (hm%bc%bc_type(idim) == MXLL_BC_MEDIUM) then
        do ip_in=1, hm%bc%mxmedium%points_number(idim)
          ip = hm%bc%mxmedium%points_map(ip_in,idim)
          st%ep(ip) = hm%bc%mxmedium%ep(ip_in,idim)
          st%mu(ip) = hm%bc%mxmedium%mu(ip_in,idim)
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

    integer :: ipx, ipy, ipz, idx, np(st%dim), idim

    PUSH_SUB(rs_state_to_cube_map)

    do idim=1, st%dim
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

    SAFE_DEALLOCATE_P(st%rs_state_fft_map)
    SAFE_DEALLOCATE_P(st%rs_state_fft_map_inv)

    POP_SUB(rs_state_to_cube_map)
  end subroutine rs_state_to_cube_map


  ! ---------------------------------------------------------
  subroutine transform_rs_state(hm, rs_state, ff_rs_state, sign)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    CMPLX,               intent(inout)      :: rs_state(:,:)
    CMPLX,               intent(inout)      :: ff_rs_state(:,:)
    integer,             intent(in)         :: sign

    ASSERT(sign == RS_TRANS_FORWARD .or. sign == RS_TRANS_BACKWARD)

    if (hm%operator == FARADAY_AMPERE_MEDIUM) then
      message(1) = "Maxwell solver in linear media not yet implemented"
      call messages_fatal(1)
    else if (hm%operator == FARADAY_AMPERE_GAUSS) then
      if (sign == RS_TRANS_FORWARD) then
         call transform_rs_state_to_4x4_rs_state_forward(rs_state, ff_rs_state)
      else
        call transform_rs_state_to_4x4_rs_state_backward(ff_rs_state, rs_state)
      end if
    else
      if (sign == RS_TRANS_FORWARD) then
        ff_rs_state(:, 1:3) = rs_state(:, 1:3)
      else
        rs_state(:, 1:3) = ff_rs_state(:, 1:3)
      end if
    end if

  end subroutine transform_rs_state


  !----------------------------------------------------------
  subroutine transform_rs_state_to_6x6_rs_state_forward(rs_state_3x3_plus, rs_state_3x3_minus, rs_state_6x6)
    CMPLX, intent(in)    :: rs_state_3x3_plus(:,:)
    CMPLX, intent(in)    :: rs_state_3x3_minus(:,:)
    CMPLX, intent(inout) :: rs_state_6x6(:,:)

    integer :: ii

    ! no push_sub, called to frequently
    do ii = 1, 3
      rs_state_6x6(:, ii) = rs_state_3x3_plus(:, ii)
      rs_state_6x6(:, ii+3) = rs_state_3x3_plus(:, ii)
    end do

  end subroutine transform_rs_state_to_6x6_rs_state_forward


  !----------------------------------------------------------
  subroutine transform_rs_state_to_6x6_rs_state_backward(rs_state_6x6, rs_state)
    CMPLX, intent(in)    :: rs_state_6x6(:,:)
    CMPLX, intent(inout) :: rs_state(:,:)

    integer :: ii

    ! no push_sub, called to frequently
    do ii = 1, 3
      rs_state(:, ii) = M_HALF * (rs_state_6x6(:, ii) + conjg(rs_state_6x6(:, ii+3)))
    end do

  end subroutine transform_rs_state_to_6x6_rs_state_backward


  !----------------------------------------------------------
  subroutine transform_rs_densities_to_6x6_rs_densities_forward(rs_charge_density, rs_current_density, rs_density_6x6)
    CMPLX, intent(in)    :: rs_charge_density(:)
    CMPLX, intent(in)    :: rs_current_density(:,:)
    CMPLX, intent(inout) :: rs_density_6x6(:,:)

    integer :: ii

    ! no push_sub, called to frequently
    do ii = 1, 3
      rs_density_6x6(:, ii) = rs_current_density(:, ii)
      rs_density_6x6(:, ii+3) = rs_current_density(:, ii)
    end do

  end subroutine transform_rs_densities_to_6x6_rs_densities_forward


  !----------------------------------------------------------
  subroutine transform_rs_densities_to_6x6_rs_densities_backward(rs_density_6x6, rs_charge_density, rs_current_density)
    CMPLX, intent(in)    :: rs_density_6x6(:,:)
    CMPLX, intent(inout) :: rs_charge_density(:)
    CMPLX, intent(inout) :: rs_current_density(:,:)

    integer :: ii

    ! no push_sub, called to frequently
    do ii = 1, 3
       rs_current_density(:, ii) = M_HALF * real(rs_density_6x6(:, ii) + rs_density_6x6(:, ii+3))
    end do

  end subroutine transform_rs_densities_to_6x6_rs_densities_backward


  !----------------------------------------------------------
  subroutine transform_rs_state_to_4x4_rs_state_forward(rs_state_3x3, rs_state_4x4)
    CMPLX, intent(in)    :: rs_state_3x3(:,:)
    CMPLX, intent(inout) :: rs_state_4x4(:,:)

    ! no push_sub, called to frequently
    rs_state_4x4(:,1) = M_z1 * (-rs_state_3x3(:,1) + rs_state_3x3(:,2))
    rs_state_4x4(:,2) = M_z1 * rs_state_3x3(:,3)
    rs_state_4x4(:,3) = M_z1 * rs_state_3x3(:,3)
    rs_state_4x4(:,4) = M_z1 * (rs_state_3x3(:,1) + rs_state_3x3(:,2))

  end subroutine transform_rs_state_to_4x4_rs_state_forward


  !----------------------------------------------------------
  subroutine transform_rs_state_to_4x4_rs_state_backward(rs_state_4x4, rs_state_3x3)
    CMPLX, intent(in)    :: rs_state_4x4(:,:)
    CMPLX, intent(inout) :: rs_state_3x3(:,:)

    ! no push_sub, called to frequently
    rs_state_3x3(:,1) = M_z1 * M_HALF * (-rs_state_4x4(:,1) + rs_state_4x4(:,4))
    rs_state_3x3(:,2) = M_zI * M_HALF * (-rs_state_4x4(:,1) - rs_state_4x4(:,4))
    rs_state_3x3(:,3) = M_z1 * M_HALF * (rs_state_4x4(:,2) + rs_state_4x4(:,3))

  end subroutine transform_rs_state_to_4x4_rs_state_backward


  !----------------------------------------------------------
  subroutine transform_rs_densities_to_4x4_rs_densities_forward(rs_charge_density, rs_current_density, rs_density_4x4)
    CMPLX, intent(in)    :: rs_charge_density(:)
    CMPLX, intent(in)    :: rs_current_density(:,:)
    CMPLX, intent(inout) :: rs_density_4x4(:,:)

    ! no push_sub, called to frequently
    rs_density_4x4(:,1) = M_z1 * (-rs_current_density(:,1) + rs_current_density(:,2))
    rs_density_4x4(:,2) = M_z1 * (rs_current_density(:,3) - rs_charge_density(:))
    rs_density_4x4(:,3) = M_z1 * (rs_current_density(:,3) + rs_charge_density(:))
    rs_density_4x4(:,4) = M_z1 * (rs_current_density(:,1) + rs_current_density(:,2))

  end subroutine transform_rs_densities_to_4x4_rs_densities_forward


  !----------------------------------------------------------
  subroutine transform_rs_densities_to_4x4_rs_densities_backward(rs_density_4x4, rs_charge_density, rs_current_density)
    CMPLX, intent(in)    :: rs_density_4x4(:,:)
    CMPLX, intent(inout) :: rs_charge_density(:)
    CMPLX, intent(inout) :: rs_current_density(:,:)

    ! no push_sub, called to frequently
    rs_charge_density(:)    = M_z1 * M_HALF * (-rs_density_4x4(:,2) + rs_density_4x4(:,3))
    rs_current_density(:,1) = M_z1 * M_HALF * (-rs_density_4x4(:,1) + rs_density_4x4(:,4))
    rs_current_density(:,2) = M_zI * M_HALF * (-rs_density_4x4(:,1) - rs_density_4x4(:,4))
    rs_current_density(:,3) = M_z1 * M_HALF * (-rs_density_4x4(:,2) + rs_density_4x4(:,3))

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
    call build_rs_state(real(tmp_grad_mx_gr(1:gr_mxll%mesh%np,:)), aimag(tmp_grad_mx_gr(1:gr_mxll%mesh%np,:)), &
      st_mxll%rs_sign, rs_state_matter(1:gr_mxll%mesh%np,:), gr_mxll%mesh, st_mxll%ep(1:gr_mxll%mesh%np), &
      st_mxll%mu(1:gr_mxll%mesh%np), gr_mxll%mesh%np)

    SAFE_DEALLOCATE_A(tmp_pot_mx_gr)
    SAFE_DEALLOCATE_A(tmp_grad_mx_gr)

    POP_SUB(calculate_matter_longitudinal_field)
  end subroutine calculate_matter_longitudinal_field


  !----------------------------------------------------------
  subroutine get_vector_pot_and_transverse_field(trans_calc_method, gr_mxll, hm_mxll, st_mxll, tr_mxll, &
    hm, st, tr, poisson_solver, time, field, transverse_field, vector_potential, geo)
    integer,                    intent(in)    :: trans_calc_method
    type(grid_t),               intent(in)    :: gr_mxll
    type(hamiltonian_mxll_t),        intent(in)    :: hm_mxll
    type(states_mxll_t),             intent(in)    :: st_mxll
    type(propagator_mxll_t), intent(in)    :: tr_mxll
    type(hamiltonian_elec_t),        intent(in)    :: hm
    type(states_elec_t),             intent(in)    :: st
    type(propagator_t),         intent(in)    :: tr
    type(poisson_t),            intent(in)    :: poisson_solver
    FLOAT,                      intent(in)    :: time
    CMPLX,                      intent(inout) :: field(:,:)
    CMPLX,                      intent(inout) :: transverse_field(:,:)
    FLOAT,                      intent(inout) :: vector_potential(:,:)
    type(geometry_t), optional, intent(in)    :: geo

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

    contains

      subroutine calculate_vector_potential(poisson_solver, gr, st, tr, field, vector_potential)
        type(poisson_t),            intent(in)    :: poisson_solver
        type(grid_t),               intent(in)    :: gr
        type(states_mxll_t),             intent(in)    :: st
        type(propagator_mxll_t), intent(in)    :: tr
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

  end subroutine get_vector_pot_and_transverse_field


  ! ---------------------------------------------------------
  subroutine derivatives_boundary_mask(bc, mesh, hm)
    type(bc_mxll_t),  intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    type(hamiltonian_mxll_t), intent(in)    :: hm

    integer :: ip, ip_in, point_info, idim, dim
    FLOAT   :: bounds(2, mesh%sb%dim), xx(mesh%sb%dim)
    FLOAT   :: ddv(mesh%sb%dim), tmp(mesh%sb%dim), width(mesh%sb%dim)
    FLOAT, allocatable :: mask(:)

    PUSH_SUB(derivatives_boundary_mask)

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

    POP_SUB(derivatives_boundary_mask)
  end subroutine derivatives_boundary_mask



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

    SAFE_ALLOCATE(ztmp(1:gr%mesh%np_part,1:st%dim))

    ztmp(:,:) = rs_field(:,:)

    energy_dens(:) = M_ZERO
    do ip = 1, gr%mesh%np
      do idim = 1, st%dim
        energy_dens(ip) = energy_dens(ip) + conjg(ztmp(ip,idim)) * ztmp(ip,idim)
      end do
    end do

    e_energy_dens(:) = M_ZERO
    do ip = 1, gr%mesh%np
      do idim = 1, st%dim
        e_energy_dens(ip) = e_energy_dens(ip) + real(ztmp(ip,idim))**2
      end do
    end do

    b_energy_dens(:) = M_ZERO
    do ip = 1, gr%mesh%np
      do idim = 1, st%dim
        b_energy_dens(ip) = b_energy_dens(ip) + aimag(ztmp(ip,idim))**2
      end do
    end do

    if (present(rs_field_plane_waves) .and. present(energy_dens_plane_waves) .and. plane_waves_check) then
      ztmp(:,:) = rs_field_plane_waves(:,:)
      energy_dens_plane_waves(:) = M_ZERO
      do ip = 1, gr%mesh%np
        do idim = 1, st%dim
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

    integer            :: ip, ip_in
    FLOAT, allocatable :: energy_density(:), energy_density_plane_waves(:), tmp(:), tmp_pw(:)
    FLOAT, allocatable :: e_energy_density(:), tmp_e(:)
    FLOAT, allocatable :: b_energy_density(:), tmp_b(:)
    
    PUSH_SUB(energy_mxll_calc)

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

    integer             :: idim, ip_surf, ix, ix_max, iy, iy_max, iz, iz_max, ii_max
    FLOAT               :: tmp_sum
    FLOAT,  allocatable :: poynting_vector(:,:), tmp_global(:,:), tmp_surf(:,:,:,:,:)

    PUSH_SUB(energy_rate_calc)

    SAFE_ALLOCATE(poynting_vector(1:gr%mesh%np,1:st%dim))
    SAFE_ALLOCATE(tmp_global(1:gr%mesh%np_global,1:st%dim))

    call get_poynting_vector(gr, st, st%rs_state, st%rs_sign, &
                                     poynting_vector, ep_field=st%ep, mu_field=st%mu)

    if (gr%mesh%parallel_in_domains) then
      do idim=1, st%dim
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

    SAFE_ALLOCATE(tmp_surf(1:2,1:st%dim,1:ii_max,1:ii_max,1:st%dim))

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

    integer             :: idim, ip_surf, ix, ix_max, iy, iy_max, iz, iz_max, ii_max
    FLOAT               :: tmp_sum
    FLOAT,  allocatable :: poynting_vector(:,:), tmp_global(:,:), tmp_surf(:,:,:,:,:)

    PUSH_SUB(energy_rate_calc_plane_waves)

    SAFE_ALLOCATE(poynting_vector(1:gr%mesh%np,1:st%dim))
    SAFE_ALLOCATE(tmp_global(1:gr%mesh%np_global,1:st%dim))

    call get_poynting_vector_plane_waves(gr, st, st%rs_sign, poynting_vector)

    if (gr%mesh%parallel_in_domains) then
      do idim=1, st%dim
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

    SAFE_ALLOCATE(tmp_surf(1:2,1:st%dim,1:ii_max,1:ii_max,1:st%dim))

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

    integer             :: idim, ip_surf, ix, ix_max, iy, iy_max, iz, iz_max, ii_max
    FLOAT,  allocatable :: poynting_vector(:,:), tmp_global(:,:), tmp_surf(:,:,:,:,:)

    PUSH_SUB(poynting_vector_through_box_surfaces)

    SAFE_ALLOCATE(poynting_vector(1:gr%mesh%np,1:st%dim))
    SAFE_ALLOCATE(tmp_global(1:gr%mesh%np_global,1:st%dim))

    call get_poynting_vector(gr, st, st%rs_state, st%rs_sign, &
                                     poynting_vector, ep_field=st%ep, mu_field=st%mu)

    if (gr%mesh%parallel_in_domains) then
      do idim=1, st%dim
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

    SAFE_ALLOCATE(tmp_surf(1:2,1:st%dim,1:ii_max,1:ii_max,1:st%dim))

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
    type(states_mxll_t), intent(in)    :: st
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt
    FLOAT,               intent(out)   :: poynting_box_surface(:,:,:)

    integer             :: idim, ip_surf, ix, ix_max, iy, iy_max, iz, iz_max, ii_max
    FLOAT,  allocatable :: poynting_vector(:,:), tmp_global(:,:), tmp_surf(:,:,:,:,:)

    PUSH_SUB(poynting_vector_through_box_surfaces_plane_waves)

    SAFE_ALLOCATE(poynting_vector(1:gr%mesh%np,1:st%dim))
    SAFE_ALLOCATE(tmp_global(1:gr%mesh%np_global,1:st%dim))

    call get_poynting_vector_plane_waves(gr, st,  st%rs_sign, poynting_vector)

    if (gr%mesh%parallel_in_domains) then
      do idim=1, st%dim
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

    SAFE_ALLOCATE(tmp_surf(1:2,1:st%dim,1:ii_max,1:ii_max,1:st%dim))

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

    integer             :: idim, ip_surf, ix, ix_max, iy, iy_max, iz, iz_max, ii_max
    FLOAT,  allocatable :: e_surf(:,:,:,:,:), b_surf(:,:,:,:,:), e_field(:,:), &
      e_field_global(:,:), b_field(:,:), b_field_global(:,:)

    PUSH_SUB(fields_through_box_surfaces)

    SAFE_ALLOCATE(e_field(1:gr%mesh%np,1:st%dim))
    SAFE_ALLOCATE(b_field(1:gr%mesh%np,1:st%dim))
    SAFE_ALLOCATE(e_field_global(1:gr%mesh%np_global,1:st%dim))
    SAFE_ALLOCATE(b_field_global(1:gr%mesh%np_global,1:st%dim))

    if (.not. hm%ma_mx_coupling) then
      ! TODO: change this (trans+long)
      call get_electric_field_state(st%rs_state_trans+st%rs_state_long, gr%mesh, e_field, &
        st%ep, gr%mesh%np)
      call get_magnetic_field_state(st%rs_state_trans+st%rs_state_long, gr%mesh, st%rs_sign, b_field, &
        st%mu, gr%mesh%np)
    else
      call get_electric_field_state(st%rs_state, gr%mesh, e_field, st%ep, gr%mesh%np)
      call get_magnetic_field_state(st%rs_state, gr%mesh, st%rs_sign, b_field, st%mu, gr%mesh%np)
    end if

    if (gr%mesh%parallel_in_domains) then
      do idim=1, st%dim
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

    SAFE_ALLOCATE(e_surf(1:2,1:st%dim,1:ii_max,1:ii_max,1:st%dim))
    SAFE_ALLOCATE(b_surf(1:2,1:st%dim,1:ii_max,1:ii_max,1:st%dim))

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

    integer             :: idim, ip_surf, ix, ix_max, iy, iy_max, iz, iz_max, ii_max
    FLOAT,  allocatable :: e_surf(:,:,:,:,:), b_surf(:,:,:,:,:), e_field(:,:), &
      e_field_global(:,:), b_field(:,:), b_field_global(:,:)

    PUSH_SUB(fields_through_box_surfaces_plane_waves)

    SAFE_ALLOCATE(e_field(1:gr%mesh%np,1:st%dim))
    SAFE_ALLOCATE(b_field(1:gr%mesh%np,1:st%dim))
    SAFE_ALLOCATE(e_field_global(1:gr%mesh%np_global,1:st%dim))
    SAFE_ALLOCATE(b_field_global(1:gr%mesh%np_global,1:st%dim))

    call get_electric_field_state(st%rs_state_plane_waves, gr%mesh, e_field, st%ep, gr%mesh%np)
    call get_magnetic_field_state(st%rs_state_plane_waves, gr%mesh, st%rs_sign, b_field, st%mu, gr%mesh%np)

    if (gr%mesh%parallel_in_domains) then
      do idim=1, st%dim
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

    SAFE_ALLOCATE(e_surf(1:2,1:st%dim,1:ii_max,1:ii_max,1:st%dim))
    SAFE_ALLOCATE(b_surf(1:2,1:st%dim,1:ii_max,1:ii_max,1:st%dim))

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

    integer            :: ip, ip_in, idim
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

    POP_SUB(mask_absorbing_boundaries)
  end subroutine mask_absorbing_boundaries


  ! ---------------------------------------------------------
  subroutine maxwell_mask(hm, rs_state)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    CMPLX,                    intent(inout) :: rs_state(:,:)

    integer :: ip, ip_in, idim

    PUSH_SUB(maxwell_mask)

    do idim = 1, 3
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
    type(hamiltonian_mxll_t),   intent(inout) :: hm
    type(grid_t),               intent(in)    :: gr
    type(states_mxll_t),        intent(inout) :: st
    type(propagator_mxll_t),    intent(inout) :: tr
    CMPLX,                      intent(in)    :: ff_rs_state(:,:)
    CMPLX,                      intent(inout) :: ff_rs_state_pml(:,:)

    integer            :: ip, ff_points, ff_dim
    CMPLX, allocatable :: ff_rs_state_plane_waves(:,:), rs_state_constant(:,:), ff_rs_state_constant(:,:)

    PUSH_SUB(pml_propagation_stage_1)

    if (tr%bc_plane_waves .and. hm%plane_waves_apply) then
      ff_points = size(ff_rs_state(:,1))
      ff_dim    = size(ff_rs_state(1,:))
      SAFE_ALLOCATE(ff_rs_state_plane_waves(1:ff_points,1:ff_dim))
      call transform_rs_state(hm, st%rs_state_plane_waves, ff_rs_state_plane_waves, RS_TRANS_FORWARD)
      ff_rs_state_pml = ff_rs_state - ff_rs_state_plane_waves
      SAFE_DEALLOCATE_A(ff_rs_state_plane_waves)
    else if (tr%bc_constant .and. hm%spatial_constant_apply) then
      ff_dim    = size(ff_rs_state(1,:))
      SAFE_ALLOCATE(rs_state_constant(1,1:3))
      SAFE_ALLOCATE(ff_rs_state_constant(1,1:ff_dim))
      rs_state_constant(1,:) = st%rs_state_const(:)
      call transform_rs_state(hm, rs_state_constant, ff_rs_state_constant, RS_TRANS_FORWARD)
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
  subroutine pml_propagation_stage_2(hm, namespace, gr, st, tr, time, dt, time_delay, ff_rs_state_pml, ff_rs_state)
    type(hamiltonian_mxll_t),   intent(inout) :: hm
    type(namespace_t),          intent(in)    :: namespace
    type(grid_t),               intent(in)    :: gr
    type(states_mxll_t),        intent(inout) :: st
    type(propagator_mxll_t),    intent(inout) :: tr
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
      call exponential_mxll_apply(hm, namespace, gr, st, tr, dt, ff_rs_state_pml)
      hm%cpml_hamiltonian = .false.
      call plane_waves_propagation(hm, st, gr, tr, time, dt, time_delay)
      SAFE_ALLOCATE(ff_rs_state_plane_waves(1:ff_points,1:ff_dim))
      call transform_rs_state(hm, st%rs_state_plane_waves, ff_rs_state_plane_waves, RS_TRANS_FORWARD)
      do ip_in=1, hm%bc%plane_wave%points_number
        ip = hm%bc%plane_wave%points_map(ip_in)
        ff_rs_state(ip,:) = ff_rs_state_pml(ip,:) + ff_rs_state_plane_waves(ip,:)
      end do
      SAFE_DEALLOCATE_A(ff_rs_state_plane_waves)
    else if (tr%bc_constant .and. hm%spatial_constant_apply) then
      ff_dim    = size(ff_rs_state(1,:))
      hm%cpml_hamiltonian = .true.
      call exponential_mxll_apply(hm, namespace, gr, st, tr, dt, ff_rs_state_pml)
      hm%cpml_hamiltonian = .false.
      SAFE_ALLOCATE(rs_state_constant(1,1:ff_dim))
      SAFE_ALLOCATE(ff_rs_state_constant(1,1:ff_dim))
      rs_state_constant(1,:) = st%rs_state_const(:)
      call transform_rs_state(hm, rs_state_constant, ff_rs_state_constant, RS_TRANS_BACKWARD)
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

    if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RS) then
      call cpml_conv_function_update_via_riemann_silberstein(hm, gr, ff_rs_state_pml, ff_dim)
    else if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__EM) then
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
      do ip_in=1, hm%bc%pml%points_number
        ip       = hm%bc%pml%points_map(ip_in)
        pml_a(:) = hm%bc%pml%a(ip_in,:)
        pml_b(:) = hm%bc%pml%b(ip_in,:)
        pml_g(:) = hm%bc%pml%conv_plus(ip_in,1,:)
        pml_g(2) = real(pml_a(1)) * real(tmp_partial(ip)) + real(pml_b(1)) * real(pml_g(2)) + &
                   M_zI * ( aimag(pml_a(1)) * aimag(tmp_partial(ip)) + aimag(pml_b(1)) * aimag(pml_g(2)) )
        hm%bc%pml%conv_plus(ip_in,1,:) = pml_g(:)
      end do
      ! calculation g(2,1)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,1), tmp_partial(:), 2, set_bc = .false.)
      do ip_in=1, hm%bc%pml%points_number
        ip       = hm%bc%pml%points_map(ip_in)
        pml_a(:) = hm%bc%pml%a(ip_in,:)
        pml_b(:) = hm%bc%pml%b(ip_in,:)
        pml_g(:) = hm%bc%pml%conv_plus(ip_in,2,:)
        pml_g(1) = real(pml_a(2)) * real(tmp_partial(ip)) + real(pml_b(2)) * real(pml_g(1)) + &
                   M_zI * ( aimag(pml_a(2)) * aimag(tmp_partial(ip)) + aimag(pml_b(2)) * aimag(pml_g(1)) )
        hm%bc%pml%conv_plus(ip_in,2,:) = pml_g(:)
      end do
      ! calculation g(1,3)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,3), tmp_partial(:), 1, set_bc = .false.)
      do ip_in=1, hm%bc%pml%points_number
        ip       = hm%bc%pml%points_map(ip_in)
        pml_a(:) = hm%bc%pml%a(ip_in,:)
        pml_b(:) = hm%bc%pml%b(ip_in,:)
        pml_g(:) = hm%bc%pml%conv_plus(ip_in,1,:)
        pml_g(3) = real(pml_a(1)) * real(tmp_partial(ip)) + real(pml_b(1)) * real(pml_g(3)) + &
                   M_zI * ( aimag(pml_a(1)) * aimag(tmp_partial(ip)) + aimag(pml_b(1)) * aimag(pml_g(3)) )
        hm%bc%pml%conv_plus(ip_in,1,:) = pml_g(:)
      end do
      ! calculation g(3,1)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,1), tmp_partial(:), 3, set_bc = .false.)
      do ip_in=1, hm%bc%pml%points_number
        ip       = hm%bc%pml%points_map(ip_in)
        pml_a(:) = hm%bc%pml%a(ip_in,:)
        pml_b(:) = hm%bc%pml%b(ip_in,:)
        pml_g(:) = hm%bc%pml%conv_plus(ip_in,3,:)
        pml_g(1) = real(pml_a(3)) * real(tmp_partial(ip)) + real(pml_b(3)) * real(pml_g(1)) + &
                   M_zI * ( aimag(pml_a(3)) * aimag(tmp_partial(ip)) + aimag(pml_b(3)) * aimag(pml_g(1)) )
        hm%bc%pml%conv_plus(ip_in,3,:) = pml_g(:)
      end do
      ! calculation g(2,3)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,3), tmp_partial(:), 2, set_bc = .false.)
      do ip_in=1, hm%bc%pml%points_number
        ip       = hm%bc%pml%points_map(ip_in)
        pml_a(:) = hm%bc%pml%a(ip_in,:)
        pml_b(:) = hm%bc%pml%b(ip_in,:)
        pml_g(:) = hm%bc%pml%conv_plus(ip_in,2,:)
        pml_g(3) = real(pml_a(2)) * real(tmp_partial(ip)) + real(pml_b(2)) * real(pml_g(3)) + &
                   M_zI * ( aimag(pml_a(2)) * aimag(tmp_partial(ip)) + aimag(pml_b(2)) * aimag(pml_g(3)) )
        hm%bc%pml%conv_plus(ip_in,2,:) = pml_g(:)
      end do
      ! calculation g(3,2)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,2), tmp_partial(:), 3, set_bc = .false.)
      do ip_in=1, hm%bc%pml%points_number
        ip       = hm%bc%pml%points_map(ip_in)
        pml_a(:) = hm%bc%pml%a(ip_in,:)
        pml_b(:) = hm%bc%pml%b(ip_in,:)
        pml_g(:) = hm%bc%pml%conv_plus(ip_in,3,:)
        pml_g(2) = real(pml_a(3)) * real(tmp_partial(ip)) + real(pml_b(3)) * real(pml_g(2)) + &
                   M_zI * ( aimag(pml_a(3)) * aimag(tmp_partial(ip)) + aimag(pml_b(3)) * aimag(pml_g(2)) )
        hm%bc%pml%conv_plus(ip_in,3,:) = pml_g(:)
      end do

      SAFE_DEALLOCATE_A(tmp_partial)

    else if (ff_dim == 6) then

      SAFE_ALLOCATE(tmp_partial_2(np_part,2))

      ! calculation g(1,2)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,2), tmp_partial_2(:,1), 1, set_bc = .false.)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,5), tmp_partial_2(:,2), 1, set_bc = .false.)
      do ip_in=1, hm%bc%pml%points_number
        ip         = hm%bc%pml%points_map(ip_in)
        pml_a(:)   = hm%bc%pml%a(ip_in,:)
        pml_b(:)   = hm%bc%pml%b(ip_in,:)
        pml_g_p(:) = hm%bc%pml%conv_plus(ip_in,1,:)
        pml_g_m(:) = hm%bc%pml%conv_minus(ip_in,1,:)
        pml_g_p(2) = real(pml_a(1)) * real(tmp_partial_2(ip,1)) + real(pml_b(1)) * real(pml_g_p(2)) + &
                     M_zI * ( aimag(pml_a(1)) * aimag(tmp_partial_2(ip,1)) + aimag(pml_b(1)) * aimag(pml_g_p(2)) )
        pml_g_m(2) = real(pml_a(1)) * real(tmp_partial_2(ip,2)) + real(pml_b(1)) * real(pml_g_m(2)) + &
                     M_zI * ( aimag(pml_a(1)) * aimag(tmp_partial_2(ip,2)) + aimag(pml_b(1)) * aimag(pml_g_m(2)) )
        hm%bc%pml%conv_plus(ip_in,1,:) = pml_g_p(:)
        hm%bc%pml%conv_minus(ip_in,1,:) = pml_g_m(:)
      end do
      ! calculation g(2,1)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,1), tmp_partial_2(:,1), 2, set_bc = .false.)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,4), tmp_partial_2(:,2), 2, set_bc = .false.)
      do ip_in=1, hm%bc%pml%points_number
        ip         = hm%bc%pml%points_map(ip_in)
        pml_a(:)   = hm%bc%pml%a(ip_in,:)
        pml_b(:)   = hm%bc%pml%b(ip_in,:)
        pml_g_p(:) = hm%bc%pml%conv_plus(ip_in,2,:)
        pml_g_m(:) = hm%bc%pml%conv_minus(ip_in,2,:)
        pml_g_p(1) = real(pml_a(2)) * real(tmp_partial_2(ip,1)) + real(pml_b(2)) * real(pml_g_p(1)) + &
                     M_zI * ( aimag(pml_a(2)) * aimag(tmp_partial_2(ip,1)) + aimag(pml_b(2)) * aimag(pml_g_p(1)) )
        pml_g_m(1) = real(pml_a(2)) * real(tmp_partial_2(ip,2)) + real(pml_b(2)) * real(pml_g_m(1)) + &
                     M_zI * ( aimag(pml_a(2)) * aimag(tmp_partial_2(ip,2)) + aimag(pml_b(2)) * aimag(pml_g_m(1)) )
        hm%bc%pml%conv_plus(ip_in,2,:) = pml_g_p(:)
        hm%bc%pml%conv_minus(ip_in,2,:) = pml_g_m(:)
      end do
      ! calculation g(1,3)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,3), tmp_partial_2(:,1), 1, set_bc = .false.)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,6), tmp_partial_2(:,2), 1, set_bc = .false.)
      do ip_in=1, hm%bc%pml%points_number
        ip         = hm%bc%pml%points_map(ip_in)
        pml_a(:)   = hm%bc%pml%a(ip_in,:)
        pml_b(:)   = hm%bc%pml%b(ip_in,:)
        pml_g_p(:) = hm%bc%pml%conv_plus(ip_in,1,:)
        pml_g_m(:) = hm%bc%pml%conv_minus(ip_in,1,:)
        pml_g_p(3) = real(pml_a(1)) * real(tmp_partial_2(ip,1)) + real(pml_b(1)) * real(pml_g_p(3)) + &
                     M_zI * ( aimag(pml_a(1)) * aimag(tmp_partial_2(ip,1)) + aimag(pml_b(1)) * aimag(pml_g_p(3)) )
        pml_g_m(3) = real(pml_a(1)) * real(tmp_partial_2(ip,2)) + real(pml_b(1)) * real(pml_g_m(3)) + &
                     M_zI * ( aimag(pml_a(1)) * aimag(tmp_partial_2(ip,2)) + aimag(pml_b(1)) * aimag(pml_g_m(3)) )
        hm%bc%pml%conv_plus(ip_in,1,:) = pml_g_p(:)
        hm%bc%pml%conv_minus(ip_in,1,:) = pml_g_m(:)
      end do
      ! calculation g(3,1)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,1), tmp_partial_2(:,1), 3, set_bc = .false.)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,4), tmp_partial_2(:,2), 3, set_bc = .false.)
      do ip_in=1, hm%bc%pml%points_number
        ip         = hm%bc%pml%points_map(ip_in)
        pml_a(:)   = hm%bc%pml%a(ip_in,:)
        pml_b(:)   = hm%bc%pml%b(ip_in,:)
        pml_g_p(:) = hm%bc%pml%conv_plus(ip_in,3,:)
        pml_g_m(:) = hm%bc%pml%conv_minus(ip_in,3,:)
        pml_g_p(1) = real(pml_a(3)) * real(tmp_partial_2(ip,1)) + real(pml_b(3)) * real(pml_g_p(1)) + &
                     M_zI * ( aimag(pml_a(3)) * aimag(tmp_partial_2(ip,1)) + aimag(pml_b(3)) * aimag(pml_g_p(1)) )
        pml_g_m(1) = real(pml_a(3)) * real(tmp_partial_2(ip,2)) + real(pml_b(3)) * real(pml_g_m(1)) + &
                     M_zI * ( aimag(pml_a(3)) * aimag(tmp_partial_2(ip,2)) + aimag(pml_b(3)) * aimag(pml_g_m(1)) )
        hm%bc%pml%conv_plus(ip_in,3,:) = pml_g_p(:)
        hm%bc%pml%conv_minus(ip_in,3,:) = pml_g_m(:)
      end do
      ! calculation g(2,3)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,3), tmp_partial_2(:,1), 2, set_bc = .false.)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,6), tmp_partial_2(:,2), 2, set_bc = .false.)
      do ip_in=1, hm%bc%pml%points_number
        ip         = hm%bc%pml%points_map(ip_in)
        pml_a(:)   = hm%bc%pml%a(ip_in,:)
        pml_b(:)   = hm%bc%pml%b(ip_in,:)
        pml_g_p(:) = hm%bc%pml%conv_plus(ip_in,2,:)
        pml_g_m(:) = hm%bc%pml%conv_minus(ip_in,2,:)
        pml_g_p(3) = real(pml_a(2)) * real(tmp_partial_2(ip,1)) + real(pml_b(2)) * real(pml_g_p(3)) + &
                     M_zI * ( aimag(pml_a(2)) * aimag(tmp_partial_2(ip,1)) + aimag(pml_b(2)) * aimag(pml_g_p(3)) )
        pml_g_m(3) = real(pml_a(2)) * real(tmp_partial_2(ip,2)) + real(pml_b(2)) * real(pml_g_m(3)) + &
                     M_zI * ( aimag(pml_a(2)) * aimag(tmp_partial_2(ip,2)) + aimag(pml_b(2)) * aimag(pml_g_m(3)) )
        hm%bc%pml%conv_plus(ip_in,2,:) = pml_g_p(:)
        hm%bc%pml%conv_minus(ip_in,2,:) = pml_g_m(:)
      end do
      ! calculation g(3,2)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,2), tmp_partial_2(:,1), 3, set_bc = .false.)
      call zderivatives_partial(gr%der, ff_rs_state_pml(:,5), tmp_partial_2(:,2), 3, set_bc = .false.)
      do ip_in=1, hm%bc%pml%points_number
        ip         = hm%bc%pml%points_map(ip_in)
        pml_a(:)   = hm%bc%pml%a(ip_in,:)
        pml_b(:)   = hm%bc%pml%b(ip_in,:)
        pml_g_p(:) = hm%bc%pml%conv_plus(ip_in,3,:)
        pml_g_m(:) = hm%bc%pml%conv_minus(ip_in,3,:)
        pml_g_p(2) = real(pml_a(3)) * real(tmp_partial_2(ip,1)) + real(pml_b(3)) * real(pml_g_p(2)) + &
                     M_zI * ( aimag(pml_a(3)) * aimag(tmp_partial_2(ip,1)) + aimag(pml_b(3)) * aimag(pml_g_p(2)) )
        pml_g_m(2) = real(pml_a(3)) * real(tmp_partial_2(ip,2)) + real(pml_b(3)) * real(pml_g_m(2)) + &
                     M_zI * ( aimag(pml_a(3)) * aimag(tmp_partial_2(ip,2)) + aimag(pml_b(3)) * aimag(pml_g_m(2)) )
        hm%bc%pml%conv_plus(ip_in,3,:) = pml_g_p(:)
        hm%bc%pml%conv_minus(ip_in,3,:) = pml_g_m(:)
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

    POP_SUB(cpml_conv_function_update_via_e_b_fields)
  end subroutine cpml_conv_function_update_via_e_b_fields


  ! ---------------------------------------------------------
  subroutine generate_medium_boxes(hm, gr, nr_of_boxes)
    type(hamiltonian_mxll_t), intent(inout) :: hm
    type(grid_t),        intent(in)    :: gr
    integer,             intent(in)    :: nr_of_boxes

    integer :: il, ip, ip_in, ip_in_max, ip_bd, ip_bd_max, ipp, idim
    FLOAT   :: bounds(2,gr%sb%dim), xx(gr%sb%dim), xxp(gr%sb%dim), dd, dd_max, dd_min
    FLOAT, allocatable  :: tmp(:), tmp_grad(:,:)

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
        ip = hm%medium_box_points_map(ip_in, il)
        tmp(ip) =  hm%medium_box_ep(ip_in, il)
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in=1, hm%medium_box_points_number(il)
        ip = hm%medium_box_points_map(ip_in, il)
        hm%medium_box_aux_ep(ip_in, :, il) = &
          tmp_grad(ip, :)/(M_FOUR * hm%medium_box_ep(ip_in, il))
      end do

      tmp = P_mu
      do  ip_in=1, hm%medium_box_points_number(il)
        ip = hm%medium_box_points_map(ip_in, il)
        tmp(ip) =  hm%medium_box_mu(ip_in, il)
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in=1, hm%medium_box_points_number(il)
        ip = hm%medium_box_points_map(ip_in, il)
        hm%medium_box_aux_mu(ip_in, :, il) = &
          tmp_grad(ip, :)/(M_FOUR * hm%medium_box_mu(ip_in, il))
      end do

      ! print information about the medium box -- get from Renes version in maxwell_propagator.F90

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


  ! ---------------------------------------------------------
  subroutine td_function_mxll_init(st, namespace, hm)
    type(states_mxll_t),      intent(inout) :: st
    type(namespace_t),   intent(in)    :: namespace
    type(hamiltonian_mxll_t), intent(in)    :: hm

    type(block_t)        :: blk
    integer              :: il, nlines, idim, ncols, ierr
    FLOAT                :: e_field(st%dim), b_field(st%dim)
    character(len=1024)  :: mxf_expression

    PUSH_SUB(td_function_init)

    !%Variable UserDefinedConstantSpacialMaxwellField
    !%Type block
    !%Section MaxwellStates
    !%Description
    !% Define parameters of spatially constant field.
    !%
    !% Example:
    !%
    !% <tt>%UserDefinedConstantSpacialMaxwellFields
    !% <br>&nbsp;&nbsp;   plane_wave_parser      | E_x | E_y | E_z | B_x | B_y | B_z | "tdf_function"
    !% <br>%</tt>
    !%
    !% This block defines three components of E field, three components of B field, and reference to
    !% TD function.
    !%
    !%End

    if(parse_block(namespace, 'UserDefinedConstantSpacialMaxwellField', blk) == 0) then
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
          message(1) = 'Each line in the UserDefinedConstantSpacialMaxwellField block must have'
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

    POP_SUB(td_function_mxll_init)
  end subroutine td_function_mxll_init


  ! ---------------------------------------------------------
  subroutine spatial_constant_calculation(constant_calc, st, gr, hm, time, dt, delay, rs_state)
    logical,             intent(in)    :: constant_calc
    type(states_mxll_t),      intent(inout)    :: st
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
    FLOAT                      :: e_field(st%dim), b_field(st%dim)

    PUSH_SUB(mirror_pec_boundaries_calculation)

    do idim = 1, 3
      if (bc%bc_type(idim) == MXLL_BC_MIRROR_PEC) then
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
    FLOAT                      :: e_field(st%dim), b_field(st%dim)

    PUSH_SUB(mirror_pmc_boundaries_calculation)

    do idim=1, 3
      if (bc%bc_type(idim) == MXLL_BC_MIRROR_PMC) then
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
  subroutine plane_waves_boundaries_calculation(hm, st, mesh, time, time_delay, rs_state)
    type(hamiltonian_mxll_t)   :: hm
    type(states_mxll_t)        :: st
    type(mesh_t)               :: mesh
    FLOAT                      :: time
    FLOAT                      :: time_delay
    CMPLX                      :: rs_state(:,:)

    integer                    :: ip, ip_in, wn
    FLOAT                      :: x_prop(mesh%sb%dim), rr, vv(mesh%sb%dim), k_vector(mesh%sb%dim)
    FLOAT                      :: k_vector_abs, nn
    FLOAT                      :: e0(mesh%sb%dim), e_field(mesh%sb%dim), b_field(mesh%sb%dim)
    CMPLX                      :: rs_state_add(mesh%sb%dim)

    PUSH_SUB(plane_waves_boundaries_calculation)

    if (hm%plane_waves_apply) then
      do wn = 1, hm%bc%plane_wave%number
        do ip_in = 1, hm%bc%plane_wave%points_number
          ip = hm%bc%plane_wave%points_map(ip_in)
            if (wn == 1) rs_state(ip,:) = M_Z0
            nn           = sqrt(st%ep(ip)/P_ep*st%mu(ip)/P_mu)
            vv(:)        = hm%bc%plane_wave%v_vector(1:mesh%sb%dim, wn)
            k_vector(:)  = hm%bc%plane_wave%k_vector(1:mesh%sb%dim, wn)
            k_vector_abs = sqrt(sum(k_vector(1:mesh%sb%dim)**2))
            x_prop(1:mesh%sb%dim) = mesh%x(ip,1:mesh%sb%dim) - vv(1:mesh%sb%dim) * (time - time_delay)
            rr           = sqrt(sum(x_prop(1:mesh%sb%dim)**2))
            if (hm%bc%plane_wave%modus(wn) == &
                OPTION__MAXWELLINCIDENTWAVES__PLANE_WAVE_MX_FUNCTION) then
              e0(1:mesh%sb%dim)      = hm%bc%plane_wave%e_field(1:mesh%sb%dim, wn)
              e_field(1:mesh%sb%dim) = e0(1:mesh%sb%dim) * mxf(hm%bc%plane_wave%mx_function(wn), x_prop(1:mesh%sb%dim))
            end if
            b_field(1:3) = M_ONE/P_c * M_ONE/k_vector_abs * dcross_product(k_vector,e_field)
            call build_rs_vector(e_field, b_field, st%rs_sign, rs_state_add, st%ep(ip), st%mu(ip))
            rs_state(ip,:) =  rs_state(ip,:) + rs_state_add(:)
          end do
        end do
      else
        do ip_in = 1, hm%bc%plane_wave%points_number
          ip             = hm%bc%plane_wave%points_map(ip_in)
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

    integer            :: ff_dim
    CMPLX, allocatable :: ff_rs_state(:,:)

    PUSH_SUB(plane_waves_propagation)

    if (hm%operator == FARADAY_AMPERE_MEDIUM) then
      ff_dim = 6
    else if (hm%operator == FARADAY_AMPERE_GAUSS) then
      ff_dim = 4
    else
      ff_dim = 3
    end if

    SAFE_ALLOCATE(ff_rs_state(1:gr%mesh%np_part,ff_dim))

    call transform_rs_state(hm, st%rs_state_plane_waves, ff_rs_state, RS_TRANS_FORWARD)

    ! Time evolution of RS plane waves state without any coupling with H(inter_time)
    call hamiltonian_mxll_update(hm, time=time)
    !call exponential_apply(tr%te, gr%mesh, hm, ff_rs_state, 1, 1, dt, time)
    call transform_rs_state(hm, st%rs_state_plane_waves, ff_rs_state, RS_TRANS_BACKWARD)
    call plane_waves_boundaries_calculation(hm, st, gr%mesh, time+dt, time_delay, st%rs_state_plane_waves)

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

    integer              :: ip, wn, idim, np
    FLOAT                :: x_prop(gr%sb%dim), rr, vv(gr%sb%dim), k_vector(gr%sb%dim), k_vector_abs, nn
    FLOAT                :: e0(gr%sb%dim), e_field(gr%sb%dim), b_field(gr%sb%dim), dummy(gr%sb%dim)
    CMPLX                :: rs_state_add(st%dim)

    PUSH_SUB(plane_waves_in_box_calculation)

    np            = gr%mesh%np_part
    do wn=1, bc%plane_wave%number
      vv(:)        = hm%bc%plane_wave%v_vector(:,wn)
      k_vector(:)  = hm%bc%plane_wave%k_vector(:,wn)
      k_vector_abs = sqrt(sum(k_vector(:)**2))
      e0(:)      = hm%bc%plane_wave%e_field(:,wn)
      do ip = 1, gr%mesh%np
        if (wn==1) rs_state(ip,:) = M_Z0
        nn           = sqrt(st%ep(ip)/P_ep*st%mu(ip)/P_mu)
        x_prop(:)    = gr%mesh%x(ip,:) - vv(:) * time
        rr           = sqrt(sum(x_prop(:)**2))
        if (bc%plane_wave%modus(wn) == OPTION__MAXWELLINCIDENTWAVES__PLANE_WAVE_PARSER) then
          do idim=1, gr%mesh%sb%dim
            call parse_expression(e_field(idim), dummy(idim), gr%mesh%sb%dim, x_prop, rr, M_ZERO, &
                 bc%plane_wave%e_field_string(idim,wn))
            e_field(idim) = units_to_atomic(units_inp%energy/units_inp%length, e_field(idim))
          end do
        else if (bc%plane_wave%modus(wn) == &
             OPTION__MAXWELLINCIDENTWAVES__PLANE_WAVE_MX_FUNCTION) then
           e0(:)      = bc%plane_wave%e_field(:,wn)
           e_field(:) = e0(:) * mxf(bc%plane_wave%mx_function(wn), x_prop(:))
        end if
        b_field(1:3) = M_ONE/(P_c * k_vector_abs) * dcross_product(k_vector, e_field)
        call build_rs_vector(e_field, b_field, st%rs_sign, rs_state_add, st%ep(ip), st%mu(ip))
        rs_state(ip,:) = rs_state(ip,:) + rs_state_add(:)
      end do
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

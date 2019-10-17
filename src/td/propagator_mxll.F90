!! Copyright (C) 2019 F. Bonafe
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
!! along with st program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
#include "global.h"

module propagator_mxll_oct_m
  use grid_oct_m
  use messages_oct_m
  use namespace_oct_m
  use exponential_oct_m
  use propagator_base_oct_m
  use hamiltonian_mxll_oct_m
  use states_mxll_oct_m
  use profiling_oct_m
  use global_oct_m
  use parser_oct_m
  
  implicit none

  public ::                           &
       propagator_mxll_t,             &
       propagator_mxll_init
  
  type, extends(propagator_abst_t) :: propagator_mxll_t
    integer             :: op_method
    logical             :: bc_add_ab_region  = .false.
    logical             :: bc_zero           = .false.
    logical             :: bc_constant       = .false.
    logical             :: bc_mirror_pec     = .false.
    logical             :: bc_mirror_pmc     = .false.
    logical             :: bc_periodic       = .false.
    logical             :: bc_plane_waves    = .false.
    logical             :: bc_medium         = .false.
    integer             :: inter_steps
    FLOAT               :: delay_time
    logical             :: plane_waves_in_box
    integer             :: tr_etrs_approx
  end type propagator_mxll_t

contains  
  
  ! ---------------------------------------------------------
  subroutine propagator_mxll_init(gr, namespace, st, tr, hm)
    type(grid_t),                 intent(in)    :: gr
    type(states_mxll_t),          intent(inout) :: st
    type(hamiltonian_mxll_t),     intent(inout) :: hm
    type(propagator_mxll_t),      intent(inout) :: tr
    type(namespace_t),            intent(in)    :: namespace

    integer :: default_propagator, il, nlines, ncols, icol, idim
    type(block_t) :: blk
    character(len=256) :: string
    logical :: plane_waves_set = .false.

    PUSH_SUB(propagator_mxll_init)

    !%Variable MaxwellTDPropagator
    !%Type integer
    !%Default maxwell_etrs
    !%Section Time-Dependent::Propagation
    !%Description
    !% There are several time-evolution methods for the Maxwell propagation
    !% similar to the methods for the time-evolution of Kohn-Sham orbitals.
    !%Option maxwell_etrs 1
    !% Enforced time-reversal-symmetry propagation (etrs)
    !%End
    default_propagator = OPTION__MAXWELLTDPROPAGATOR__MAXWELL_ETRS
    call parse_variable(namespace, 'MaxwellTDPropagator', default_propagator, tr%method)
    call messages_print_var_option(stdout, 'MaxwellTDPropagator', tr%method)

    !%Variable MaxwellBoundaryConditions
    !%Type block
    !%Section Time-Dependent::Propagation
    !%Description
    !% Follows
    !%
    !% Example:
    !%
    !% <tt>%UserDefinedMaxwellIncidentWaves
    !% <br>&nbsp;&nbsp;   maxwell_zero | maxwell_mirror_pec | maxwell_consant 
    !% <br>%</tt>
    !%
    !% Description follows
    !%
    !%Option maxwell_zero 0
    !% follows ...
    !%Option maxwell_constant 2
    !% follows ...
    !%Option maxwell_mirror_pec 3
    !% follows ...
    !%Option maxwell_mirror_pmc 4
    !% follows ...
    !%Option maxwell_plane_waves 5
    !% follows ...
    !%Option maxwell_periodic 6
    !% follows ...
    !%Option maxwell_medium 7
    !% follows ...
    !%Option maxwell_lossy_layer 8
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
          if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_ZERO) then
             string = 'Zero'
             tr%bc_zero = .true.
          else if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_CONSTANT) then
             string = 'Constant'
             tr%bc_constant = .true.
             tr%bc_add_ab_region = .true.
             hm%bc_constant = .true.
             hm%bc_add_ab_region = .true.
             SAFE_ALLOCATE(st%rs_state_const(1:st%d%dim))
             st%rs_state_const = M_z0
             !          call maxwell_td_function_init(st, gr, hm)
          else if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MIRROR_PEC) then
             string = 'PEC Mirror'
             tr%bc_mirror_pec = .true.
             hm%bc_mirror_pec = .true.
          else if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MIRROR_PMC) then
             string = 'PMC Mirror'
             tr%bc_mirror_pmc = .true.
             hm%bc_mirror_pmc = .true.
          else if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_PERIODIC) then
             string = 'Periodic'
             tr%bc_periodic = .true.
             hm%bc_periodic = .true.
          else if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_PLANE_WAVES) then
             string = 'Plane waves'
             plane_waves_set = .true.
             tr%bc_plane_waves = .true.
             tr%bc_add_ab_region = .true.
             hm%plane_waves = .true.
             hm%bc_plane_waves = .true.
             hm%bc_add_ab_region = .true.
             SAFE_ALLOCATE(st%rs_state_plane_waves(1:gr%mesh%np_part, 1:st%d%dim))
          else if (hm%bc%bc_type(icol) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MEDIUM) then
             string = 'Medium boundary'
          end if
          write(message(1),'(a,I1,a,a)') 'Maxwell boundary condition in direction ', icol, ': ', trim(string)
          call messages_info(1)
          if (plane_waves_set .and. .not. (parse_is_defined(namespace, 'UserDefinedMaxwellIncidentWaves')) ) then
             write(message(1),'(a)') 'Input: Maxwell boundary condition option is set to "maxwell_plane_waves".'
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

       !      call generate_medium_boxes(hm, gr, nlines)

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
    !%Option mxll_test 3
    !% follows
    !%End
    call parse_variable(namespace, 'MaxwellTDETRSApprox', OPTION__MAXWELLTDETRSAPPROX__NO, tr%tr_etrs_approx)
    call messages_print_var_option(stdout, 'MaxwellTDETRSApprox', tr%tr_etrs_approx)

    !%Variable MaxwellTDOperatorMethod
    !%Type integer
    !%Default maxwell_op_fd
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
    !% The self consistency has to be measured against some accuracy 
    !% threshold. This variable controls the value of that Maxwell threshold.
    !%End
    call parse_variable(namespace, 'MaxwellTDSCFThreshold', CNST(1.0e-6), tr%scf_threshold)

    !%Variable MaxwellPlaneWavesInBox
    !%Type logical
    !%Default no
    !%Section States
    !%Description
    !% Follows
    !%End
    call parse_variable(namespace, 'MaxwellPlaneWavesInBox', .false., tr%plane_waves_in_box)

    !    call set_medium_rs_state(st, gr, hm)
    !    call derivatives_boundary_mask(hm%bc, gr%mesh, hm)

    !tr%te%exp_method = .true.  ! For reviewers: I think this is not needed anymore
    call exponential_init(tr%te, namespace) ! initialize Maxwell propagator

    POP_SUB(propagator_mxll_init)
  end subroutine propagator_mxll_init

  
end module propagator_mxll_oct_m

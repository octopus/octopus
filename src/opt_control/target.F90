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
!! $Id: opt_control.F90 2870 2007-04-28 06:26:47Z acastro $

#include "global.h"

module opt_control_target_m
  use datasets_m
  use density_m
  use derivatives_m
  use excited_states_m
  use geometry_m
  use global_m
  use grid_m
  use output_m
  use io_m
  use io_function_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use opt_control_global_m
  use parser_m
  use profiling_m
  use restart_m
  use species_m
  use species_pot_m
  use spectrum_m
  use states_m
  use states_calc_m
  use states_dim_m
  use string_m
  use td_m
  use unit_m
  use unit_system_m
  use types_m
  use varinfo_m

  implicit none

  private
  public :: target_t,                &
            target_get_state,        &
            target_init,             &
            target_end,              &
            target_output,           &
            target_tdcalc,           &
            target_inh,              &
            target_mode,             &
            target_type,             &
            j1_functional,           &
            calc_chi,                &
            target_move_ions,        &
            parse_velocity_target,   &
            current_functional_type

  integer, public, parameter ::       &
    oct_tg_groundstate      = 1,      &
    oct_tg_excited          = 2,      &
    oct_tg_gstransformation = 3,      &
    oct_tg_userdefined      = 4,      &
    oct_tg_density          = 5,      &        
    oct_tg_local            = 6,      &
    oct_tg_td_local         = 7,      &
    oct_tg_exclude_state    = 8,      &
    oct_tg_hhg              = 9,      &
    oct_tg_velocity         = 10,     &
    oct_tg_current          = 11

  integer, public, parameter ::       &
    oct_targetmode_static = 0,        &
    oct_targetmode_td     = 1

  integer, public, parameter ::       &
    oct_no_curr              = 0,     &
    oct_min_curr             = 1,     &
    oct_max_curr             = 2,     &
    oct_max_curr_ring        = 3,     &
    oct_min_curr_td          = 4

  type target_t
    private
    integer :: type
    type(states_t) :: st
    type(excited_states_t) :: est
    FLOAT, pointer :: rho(:) => null()
    FLOAT, pointer :: td_fitness(:) => null()
    character(len=200) :: td_local_target
    character(len=80) :: excluded_states_list
    character(len=4096) :: vel_input_string
    character(len=1024), pointer :: vel_der_array(:,:) => null()
    FLOAT, pointer :: grad_local_pot(:,:,:) => null()
    logical :: move_ions
    integer :: hhg_nks
    integer, pointer :: hhg_k(:) => null()
    FLOAT,   pointer :: hhg_alpha(:) => null()
    FLOAT,   pointer :: hhg_a(:) => null()
    FLOAT   :: hhg_w0
    FLOAT   :: dt
    integer :: curr_functional
    FLOAT   :: curr_weight
    integer :: strt_iter_curr_tg
  end type target_t


  contains


  ! ----------------------------------------------------------------------
  ! This just copies the states_t variable present in target, into st.
  ! ----------------------------------------------------------------------
  subroutine target_get_state(target, st)
    type(target_t), intent(in)    :: target
    type(states_t), intent(inout) :: st

    PUSH_SUB(target_get_state)
    call states_copy(st, target%st)

    POP_SUB(target_get_state)
  end subroutine target_get_state
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  ! The target is initialized, mainly by reading from the inp file.
  ! ----------------------------------------------------------------------
  subroutine target_init(gr, geo, stin, td, w0, target, oct)
    type(grid_t),     intent(in)    :: gr
    type(geometry_t), intent(in)    :: geo
    type(states_t),   intent(inout) :: stin 
    type(td_t),       intent(in)    :: td
    FLOAT,            intent(in)    :: w0
    type(target_t),   intent(inout) :: target
    type(oct_t),      intent(in)    :: oct

    integer             :: ierr, ip, ist, jst, jj, iatom, ib, idim, inst, inik, &
                           id, ik, no_states
    type(block_t)       :: blk
    FLOAT               :: xx(MAX_DIM), rr, psi_re, psi_im
    FLOAT, allocatable  :: vl(:), vl_grad(:,:)
    FLOAT               :: stencil_size, point_dist
    CMPLX, allocatable  :: rotation_matrix(:, :)
    type(states_t)      :: tmp_st
    character(len=1024) :: expression

    PUSH_SUB(target_init)

    !%Variable OCTTargetOperator
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default oct_tg_gstransformation
    !%Description
    !% The variable <tt>OCTTargetOperator</tt> prescribes which kind of target functional is
    !% to be used.
    !%Option oct_tg_groundstate 1 
    !% The target operator is a projection operator on the ground state, <i>i.e.</i> the
    !% objective is to populate the ground state as much as possible.
    !%Option oct_tg_excited 2
    !% The target operator is an "excited state". This means that the target operator
    !% is a linear combination of Slater determinants, each one formed by replacing
    !% in the ground-state Slater determinant one occupied state with one excited
    !% state (<i>i.e.</i> "single excitations"). The description of which excitations are
    !% used, and with which weights, should be given in a file called
    !% <tt>oct-excited-state-target</tt>. This is still in very preliminary, experimental
    !% phase. See the documentation of subroutine <tt>excited_states_init</tt> in the source
    !% code in order to use this feature.
    !%Option oct_tg_gstransformation 3
    !% The target operator is a projection operator on a transformation of the ground-state 
    !% orbitals defined by the block <tt>OCTTargetTransformStates</tt>.
    !%Option oct_tg_userdefined 4
    !% Allows to define target state by using <tt>OCTTargetUserdefined</tt>.
    !%Option oct_tg_density 5
    !% The target operator is a given density, <i>i.e.</i> the final state should have a density
    !% as close as possible as the one given in the input file, either from the variable
    !% <tt>OCTTargetDensityFromState</tt>, or from <tt>OCTTargetDensity</tt>. It can be extended to a
    !% combination with a current functional by setting <tt>OCTCurrentFunctional</tt> and attributing
    !% a value to <tt>OCTCurrentWeight</tt>. 
    !%Option oct_tg_local 6
    !% The target operator is a local operator.
    !%Option oct_tg_td_local 7
    !% The target operator is a time-dependent local operator.
    !%Option oct_tg_exclude_state 8
    !% Target operator is the projection onto the complement of a given state, given by the
    !% block <tt>OCTTargetTransformStates</tt>. This means that the target operator is the unity
    !% operator minus the projector onto that state.
    !%Option oct_tg_hhg 9
    !% The target is the optimization of the HHG yield.
    !%Option oct_tg_velocity 10
    !% The target is a function of the velocities of the nuclei at the end of the influence of
    !% the external field, defined by <tt>OCTVelocityTarget</tt>
    !%Option oct_tg_current 11
    !% The target is exclusively a target in terms of the current. 
    !% If combined with target that involves the density, set variable <tt>OCTTargetOperator</tt>= <tt>OCTTargetDensity</tt> 
    !% and set explicitly <tt>OCTCurrentFunctional</tt>. Only this combination is enabled. All other targets force
    !% <tt>OCTCurrentFunctional</tt>=0.
    !%End
    call parse_integer(datasets_check('OCTTargetOperator'), oct_tg_gstransformation, target%type)
    if(.not.varinfo_valid_option('OCTTargetOperator', target%type)) &
      call input_error('OCTTargetOperator')

    call states_copy(target%st, stin)
    call states_deallocate_wfns(target%st)
    call states_allocate_wfns(target%st, gr%mesh, TYPE_CMPLX)

    nullify(target%td_fitness)

    select case(target%type)
    case(oct_tg_groundstate)
      message(1) =  'Info: Using Ground State for TargetOperator'
      call messages_info(1)
      call restart_read(trim(restart_dir)//GS_DIR, target%st, gr, geo, ierr, exact = .true.)
      
    case(oct_tg_excited) 

      message(1) =  'Info: TargetOperator is a linear combination of Slater determinants.'
      call messages_info(1)

      call states_look (trim(restart_dir)//GS_DIR, gr%mesh%mpi_grp, ip, ip, target%st%nst, ierr)
      target%st%st_start = 1
      target%st%st_end   = target%st%nst

      SAFE_DEALLOCATE_P(target%st%occ)
      SAFE_DEALLOCATE_P(target%st%eigenval)
      SAFE_DEALLOCATE_P(target%st%node)

      SAFE_ALLOCATE(     target%st%occ(1:target%st%nst, 1:target%st%d%nik))
      SAFE_ALLOCATE(target%st%eigenval(1:target%st%nst, 1:target%st%d%nik))
      SAFE_ALLOCATE(    target%st%node(1:target%st%nst))
      if(target%st%d%ispin == SPINORS) then
        SAFE_DEALLOCATE_P(target%st%spin)
        SAFE_ALLOCATE(target%st%spin(1:3, 1:target%st%nst, 1:target%st%d%nik))
      end if
      call states_allocate_wfns(target%st, gr%mesh, TYPE_CMPLX)
      target%st%node(:)  = 0

      call restart_read(trim(restart_dir)//GS_DIR, target%st, gr, geo, ierr, exact = .true.)
      call excited_states_init(target%est, target%st, "oct-excited-state-target") 

    case(oct_tg_exclude_state)

      message(1) =  'Info: The target functional is the exclusion of a number of states defined by'
      message(2) =  '      "OCTExcludedStates".'
      call messages_info(2)
      !%Variable OCTExcludedStates
      !%Type string
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If the target is the exclusion of several targets, ("OCTTargetOperator = oct_exclude_states") 
      !% then you must declare which states are to be excluded, by setting the OCTExcludedStates variable.
      !% It must be a string in "list" format: "1-8", or "2,3,4-9", for example. Be careful to include
      !% in this list only states that have been calculated in a previous "gs" or "unocc" calculation,
      !% or otherwise the error will be silently ignored.
      !%End
      call parse_string(datasets_check('OCTExcludedStates'), "1", target%excluded_states_list)
      call states_deallocate_wfns(target%st)
      call restart_look_and_read(target%st, gr, geo)

    case(oct_tg_gstransformation)  

      message(1) =  'Info: Using Superposition of States for TargetOperator'
      call messages_info(1)

      !%Variable OCTTargetTransformStates
      !%Type block
      !%Default no
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If <tt>OCTTargetOperator = oct_tg_gstransformation</tt>, you must specify a
      !% <tt>OCTTargetTransformStates</tt> block, in order to specify which linear
      !% combination of the states present in <tt>restart/gs</tt> is used to
      !% create the target state.
      !% 
      !% The syntax is the same as the <tt>TransformStates</tt> block.
      !%End
      if(parse_isdef(datasets_check('OCTTargetTransformStates')) .ne. 0) then
        if(parse_block(datasets_check('OCTTargetTransformStates'), blk) == 0) then
          call states_copy(tmp_st, target%st)
          call states_deallocate_wfns(tmp_st)
          call restart_look_and_read(tmp_st, gr, geo)
          SAFE_ALLOCATE(rotation_matrix(1:target%st%nst, 1:tmp_st%nst))
          rotation_matrix = M_z0
          do ist = 1, target%st%nst
            do jst = 1, parse_block_cols(blk, ist - 1)
              call parse_block_cmplx(blk, ist - 1, jst - 1, rotation_matrix(ist, jst))
            end do
          end do

          call states_rotate(gr%mesh, target%st, tmp_st, rotation_matrix)
          SAFE_DEALLOCATE_A(rotation_matrix)
          call states_end(tmp_st)
          call parse_block_end(blk)
          call density_calc(target%st, gr, target%st%rho)
        else
          message(1) = '"OCTTargetTransformStates" has to be specified as block.'
          call messages_info(1)
          call input_error('OCTTargetTransformStates')
        end if
      else
        message(1) = 'Error: if "OCTTargetOperator = oct_tg_superposition", then you must'
        message(2) = 'supply one "OCTTargetTransformStates" block to create the superposition.'
        call messages_info(2)
        call input_error('OCTTargetTransformStates')
      end if

    case(oct_tg_userdefined) 
       message(1) =  'Info: Target is a user-defined state.'
      call messages_info(1)
      
      !%Variable OCTTargetUserdefined
      !%Type block
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% Example:
      !%
      !% <tt>%OCTTargetUserdefined
      !% <br>&nbsp;&nbsp; 1 | 1 | 1 |  "exp(-r^2)*exp(-i*0.2*x)"
      !% <br>%</tt>
      !%  
      !%End
      if(parse_block(datasets_check('OCTTargetUserdefined'), blk) == 0) then
        
        no_states = parse_block_n(blk)
        do ib = 1, no_states
          call parse_block_integer(blk, ib - 1, 0, idim)
          call parse_block_integer(blk, ib - 1, 1, inst)
          call parse_block_integer(blk, ib - 1, 2, inik)

          ! read formula strings and convert to C strings
          do id = 1, target%st%d%dim
            do ist = 1, target%st%nst
              do ik = 1, target%st%d%nik   
                
                ! does the block entry match and is this node responsible?
                if(.not. (id .eq. idim .and. ist .eq. inst .and. ik .eq. inik    &
                  .and. target%st%st_start .le. ist .and. target%st%st_end .ge. ist) ) cycle
                
                ! parse formula string
                call parse_block_string(                            &
                  blk, ib - 1, 3, target%st%user_def_states(id, ist, ik))
                ! convert to C string
                call conv_to_C_string(target%st%user_def_states(id, ist, ik))
                
                do ip = 1, gr%mesh%np
                  xx = gr%mesh%x(ip, :)
                  rr = sqrt(sum(xx(:)**2))
                  
                  ! parse user-defined expressions
                  call parse_expression(psi_re, psi_im, &
                    gr%sb%dim, xx, rr, M_ZERO, target%st%user_def_states(id, ist, ik))
                  ! fill state
                  target%st%zpsi(ip, id, ist, ik) = psi_re + M_zI * psi_im
                end do
                ! normalize orbital
                call zstates_normalize_orbital(gr%mesh, target%st%d%dim, &
                  target%st%zpsi(:,:, ist, ik))
              end do
            end do
          enddo
        end do
        call parse_block_end(blk)
      else
        message(1) = '"OCTTargetUserdefined" has to be specified as block.'
        call messages_fatal(1)
      end if

    case(oct_tg_density) 

      message(1) =  'Info: Target is a density.'
      call messages_info(1)

      !%Variable OCTTargetDensity
      !%Type string
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If <tt>OCTTargetOperator = oct_tg_density</tt>, then one must supply the target density
      !% that should be searched for. This one can do by supplying a string through
      !% the variable <tt>OCTTargetDensity</tt>.
      !%End

      !%Variable OCTTargetDensityFromState
      !%Type block
      !%Default no
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If <tt>OCTTargetOperator = oct_tg_density</tt>, and <tt>OCTLocalTarget = "OCTTargetDensityFromState"</tt>,
      !% you must specify a <tt>OCTTargetDensityState</tt> block, in order to specify which linear
      !% combination of the states present in <tt>restart/gs</tt> is used to
      !% create the target density.
      !%
      !% The syntax is the same as the <tt>TransformStates</tt> block.
      !%End

      if(parse_isdef('OCTTargetDensity').ne.0) then
        SAFE_ALLOCATE(target%rho(1:gr%mesh%np))
        target%rho = M_ZERO
        call parse_string('OCTTargetDensity', "0", expression)


        if(trim(expression) .eq. 'OCTTargetDensityFromState') then

          if(parse_block(datasets_check('OCTTargetDensityFromState'), blk) == 0) then
            call states_copy(tmp_st, target%st)
            SAFE_DEALLOCATE_P(tmp_st%zpsi)
            call restart_look_and_read(tmp_st, gr, geo)
            SAFE_ALLOCATE(rotation_matrix(1:target%st%nst, 1:tmp_st%nst))
            rotation_matrix = M_z0
            do ist = 1, target%st%nst
              do jst = 1, parse_block_cols(blk, ist - 1)
                call parse_block_cmplx(blk, ist - 1, jst - 1, rotation_matrix(ist, jst))
              end do
            end do
            call states_rotate(gr%mesh, target%st, tmp_st, rotation_matrix)
            SAFE_DEALLOCATE_A(rotation_matrix)
            call density_calc(target%st, gr, target%st%rho)
            do ip = 1, gr%mesh%np
              target%rho(ip) = sum(target%st%rho(ip, 1:target%st%d%spin_channels))
            end do
            call states_end(tmp_st)
            call parse_block_end(blk)
          else
            message(1) = '"OCTTargetDensityState" has to be specified as block.'
            call messages_info(1)
            call input_error('OCTTargetDensity')
          end if

        else

          call conv_to_C_string(expression)
          do ip = 1, gr%mesh%np
            call mesh_r(gr%mesh, ip, rr, coords = xx)
            ! parse user-defined expression
            call parse_expression(psi_re, psi_im, gr%sb%dim, xx, rr, M_ZERO, expression)
            target%rho(ip) = psi_re
          end do
          ! Normalize
          rr = dmf_integrate(gr%mesh, target%rho)
          target%rho = (-target%st%val_charge) * target%rho/rr
        end if

      else
        message(1) = 'If OCTTargetOperator = oct_tg_density, then you must give the shape'
        message(2) = 'of this target in variable "OCTTargetDensity".'
        call messages_fatal(2)
      end if
      
      ! For a target combining density and current

      if(parse_isdef('OCTCurrentFunctional') .ne. 0) then
        message(1) =  'Info: Target is also a current.'
        call messages_info(1)

        call parse_integer(datasets_check('OCTCurrentFunctional'), oct_no_curr, target%curr_functional)
        if( target%curr_functional .gt. M_FOUR .or. &
            target%curr_functional .lt. M_ZERO) then
        message(1) = 'This option is not available for a current functional.'
        message(2) = 'Define 0 .le. "OCTCurrentFunctional" .le. 4.'
        call messages_fatal(2)
        end if
        
        call parse_float(datasets_check('OCTCurrentWeight'), M_ZERO, target%curr_weight)
        select case(target%curr_functional)
        case (oct_min_curr, oct_min_curr_td, oct_max_curr) 
          if( target%curr_weight .lt. M_ZERO) then
            message(1) = 'For "OCTCurrentFunctional" = oct_min_curr(_td) or oct_max_curr(_td)'
            message(2) = 'set "OCTCurrentWeight" .ge. 0.0'
            call messages_fatal(2)
          end if
        end select   
        write(message(1), '(a,i3)')   'Info: OCTCurrentFunctional = ', target%curr_functional
        write(message(2), '(a,f8.3)') 'Info: OCTCurrentWeight = ',  target%curr_weight
        call messages_info(2)

        call parse_integer(datasets_check('OCTStartIterCurrTg'), 0, target%strt_iter_curr_tg)
        if (target_mode(target) .eq. oct_targetmode_td) then
          write(message(1), '(a,i3)')   'Info: TargetMode = ', target_mode(target)
          write(message(2), '(a,i8)') 'Info: OCTStartIterCurrTg = ',  target%strt_iter_curr_tg
          call messages_info(2)
          target%dt = td%dt
          SAFE_ALLOCATE(target%td_fitness(0:td%max_iter))
          target%td_fitness = M_ZERO
        else
          target%strt_iter_curr_tg = M_ZERO
        end if
        if (target%strt_iter_curr_tg .lt. M_ZERO) then
          message(1) = 'OCTStartIterCurrTg must be negative'
          call messages_fatal(1)
        end if
        if (target%strt_iter_curr_tg .ge. td%max_iter) then
          message(1) = 'OCTStartIterCurrTg has to be .lt. TDMaximumIter'
          call messages_fatal(1)
        end if

        SAFE_ALLOCATE(stin%current( 1:gr%mesh%np_part, 1:gr%mesh%sb%dim, 1:stin%d%nspin ) )
        stin%current= M_ZERO
      end if


    case(oct_tg_local)
      !%Variable OCTLocalTarget
      !%Type string
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If <tt>OCTTargetOperator = oct_tg_local</tt>, then one must supply a function
      !% that defines the target. This should be done by defining it through a string, using 
      !% the variable <tt>OCTLocalTarget</tt>.
      !%End

      if(parse_isdef('OCTLocalTarget') .ne. 0) then
        SAFE_ALLOCATE(target%rho(1:gr%mesh%np))
        target%rho = M_ZERO
        call parse_string('OCTLocalTarget', "0", expression)
        call conv_to_C_string(expression)
        do ip = 1, gr%mesh%np
          call mesh_r(gr%mesh, ip, rr, coords = xx)
          ! parse user-defined expression
          call parse_expression(psi_re, psi_im, gr%sb%dim, xx, rr, M_ZERO, expression)
          target%rho(ip) = psi_re
        end do
      else
        message(1) = 'If OCTTargetOperator = oct_tg_local, then you must give the shape'
        message(2) = 'of this target in variable "OCTLocalTarget".'
        call messages_fatal(2)
      end if

    case(oct_tg_td_local)
      if(parse_block(datasets_check('OCTTdTarget'),blk)==0) then
        call parse_block_string(blk, 0, 0, target%td_local_target)
        call conv_to_C_string(target%td_local_target)
        SAFE_ALLOCATE(target%rho(1:gr%mesh%np))
      else
        message(1) = 'If OCTTargetMode = oct_targetmode_td, you must suppy a OCTTDTarget block.'
        call messages_fatal(1)
      end if
      target%dt = td%dt
      SAFE_ALLOCATE(target%td_fitness(0:td%max_iter))
      target%td_fitness = M_ZERO
      call target_build_tdlocal(target, gr, M_ZERO)

    case(oct_tg_hhg)
      !%Variable OCTOptimizeHarmonicSpectrum
      !%Type block
      !%Default no
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% WARNING: Experimental
      !%
      !% If <tt>OCTTargetOperator = oct_tg_hhg</tt>, the target is the harmonic emission spectrum.
      !% In that case, you must supply an <tt>OCTOptimizeHarmonicSpectrum</tt> block in the <tt>inp</tt>
      !% file. The target is given, in general, by:
      !%
      !% <math>J_1 = \int_0^\infty d\omega \alpha(\omega) H(\omega)</math>,
      !%
      !% where <math>H(\omega)</math> is the harmonic spectrum generated by the system, and
      !% <math>\alpha(\omega)</math> is some function that determines what exactly we want
      !% to optimize. The role of the <tt>OCTOptimizeHarmonicSpectrum</tt> block is to determine
      !% this <math>\alpha(\omega)</math> function. Currently, this function is defined as:
      !%
      !% <math>\alpha(\omega) = \sum_{L=1}^{M} \frac{\alpha_L}{a_L} \sqcap( (\omega - L\omega_0)/a_L )</math>,
      !%
      !% where <math>omega_0</math> is the carrier frequency. <math>M</math> is
      !% the number of columns in the <tt>OCTOptimizeHarmonicSpectrum</tt> block. The values of <i>L</i> will be listed
      !% in the first row of this block; <math> alpha_L </math> in the second row, and <math>a_L</math> in
      !% the third.
      !% 
      !% Example:
      !%
      !% <tt>%OCTOptimizeHarmonicSpectrum
      !% <br>&nbsp;&nbsp;  7    |  9    | 11
      !% <br>&nbsp;&nbsp; -1    |  1    | -1 
      !% <br>&nbsp;&nbsp;  0.01 |  0.01 |  0.01
      !% <br>%</tt>
      !%
      !%End
      if(parse_isdef(datasets_check('OCTOptimizeHarmonicSpectrum')) .ne. 0) then
        if(parse_block(datasets_check('OCTOptimizeHarmonicSpectrum'), blk) == 0) then
          target%hhg_nks = parse_block_cols(blk, 0)
          SAFE_ALLOCATE(    target%hhg_k(1:target%hhg_nks))
          SAFE_ALLOCATE(target%hhg_alpha(1:target%hhg_nks))
          SAFE_ALLOCATE(    target%hhg_a(1:target%hhg_nks))
          do jj = 1, target%hhg_nks
            call parse_block_integer(blk, 0, jj - 1, target%hhg_k(jj))
            call parse_block_float(blk, 1, jj - 1, target%hhg_alpha(jj))
            call parse_block_float(blk, 2, jj - 1, target%hhg_a(jj))
          end do
        else
          message(1) = '"OCTOptimizeHarmonicSpectrum" has to be specified as a block.'
          call messages_info(1)
          call input_error('OCTOptimizeHarmonicSpectrum')
        end if
      else
        write(message(1), '(a)') 'If "OCTTargetMode = oct_targetmode_hhg", you must supply an'
        write(message(2), '(a)') '"OCTOptimizeHarmonicSpectrum" block.'
        call messages_fatal(2)
      end if

      target%hhg_w0 = w0
      target%dt     = td%dt
      SAFE_ALLOCATE(target%td_fitness(0:td%max_iter))
      target%td_fitness = M_ZERO

    case(oct_tg_velocity)
      !%Variable OCTVelocityTarget
      !%Type block
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If <tt>OCTTargetOperator = oct_tg_velocity</tt>, then one must supply the 
      !% target to optimize in terms of the ionic velocities. This is done by 
      !% supplying a string trough the block <tt>OCTVelocityTarget</tt>.
      !% Each velocity component is supplied by <tt>"v[n_atom,vec_comp]"</tt>,
      !% while "n_atom" is the respective atom number, corresponding to the 
      !% <tt>Coordinates</tt> block and "vec_comp" is the corresponding
      !% vector component of the velocity. The target string can be
      !% supplied by using several lines in the OCTTargetOperator block.
      !% As an example, the following target can be used to maximize the
      !% velocity difference between atom 1 and 2 (in a 3D system):
      !%
      !% <tt>%OCTVelocityTarget</tt>
      !% <tt> "(v[1,1]-v[2,1])^2 + (v[1,2]-v[2,2])^2 + "</tt>
      !% <tt> "(v[1,3]-v[2,3])^2"</tt>
      !% <tt>%</tt>
      !%
      !%End

      !%Variable OCTVelocityDerivatives
      !%Type block
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If <tt>OCTTargetOperator = oct_tg_velocity</tt>, and
      !% <tt>OCTScheme = oct_algorithm_cg</tt> then you must supply 
      !% the target in terms of the ionic velocities AND the derivatives
      !% of the target with respect to the ionic velocity components.
      !% The derivatives are supplied via strings trough the block
      !% <tt>OCTVelocityDerivatives</tt>.
      !% Each velocity component is supplied by <tt>"v[n_atom,vec_comp]"</tt>,
      !% while "n_atom" is the atom number, corresponding to the 
      !% <tt>Coordinates</tt> block and "vec_comp" is the corresponding
      !% vector component of the velocity. The first line of the 
      !% <tt>OCTVelocityDerivatives</tt> block contains the derivatives
      !% with respect to "v[1,*]", the second with respect to "v[2,*]" and so
      !% on. The first column contains all derivatives with respect "v[*,1]",
      !% the second with respect to "v[*,2]" and the third w.r.t. "v[*,3]".
      !% As an example, we show the <tt>OCTVelocityDerivatives</tt> block
      !% corresponding to the target shown in the <tt>OCTVelocityTarget</tt> 
      !% help section:
      !%
      !% <tt>%OCTVelocityDerivatives</tt>
      !% <tt> " 2*(v[1,1]-v[2,1])" | " 2*(v[1,2]-v[2,2])" | " 2*(v[1,3]-v[2,3])" </tt>
      !% <tt> "-2*(v[1,1]-v[2,1])" | "-2*(v[1,2]-v[2,2])" | "-2*(v[1,3]-v[2,3])" </tt>
      !% <tt>%</tt>
      !%
      !%End
       
      !%Variable OCTMoveIons
      !%Type logical
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If <tt>OCTTargetOperator = oct_tg_velocity</tt>, then one must specify
      !% if the ions are assumed to be fixed or if they can move by setting
      !% <tt>OCTMoveIons</tt> to <tt>true</tt> or <tt>false</tt>.
      !%End
       
       if(parse_block(datasets_check('OCTVelocityTarget'),blk)==0) then
          target%vel_input_string = " "
          do jj=0, parse_block_n(blk)-1
             call parse_block_string(blk, jj, 0, expression)
             target%vel_input_string = trim(target%vel_input_string) // trim(expression)
          end do
       else
          message(1) = 'If OCTTargetOperator = oct_tg_velocity, then you must give the shape'
          message(2) = 'of this target in the block "OCTVelocityTarget".'
          call messages_fatal(2)
       end if
       
       if(parse_isdef('OCTMoveIons') .eq. 0) then
          message(1) = 'If OCTTargetOperator = oct_tg_velocity, then you must supply'
          message(2) = 'the variable "OCTMoveIons".'
          call messages_fatal(2)
       else
          call parse_logical('OCTMoveIons', .false., target%move_ions)
       end if
       
       if(oct%algorithm .eq. oct_algorithm_cg) then
          if(parse_block(datasets_check('OCTVelocityDerivatives'),blk)==0) then
             SAFE_ALLOCATE(target%vel_der_array(1:geo%natoms,1:gr%sb%dim))
             do ist=0, geo%natoms-1
                do jst=0, gr%sb%dim-1
                   call parse_block_string(blk, ist, jst, target%vel_der_array(ist+1, jst+1))
                end do
             end do
          else
             message(1) = 'If OCTTargetOperator = oct_tg_velocity, and'
             message(2) = 'OCTScheme = oct_algorithm_cg, then you must define the'
             message(3) = 'blocks "OCTVelocityTarget" AND "OCTVelocityDerivatives"'
             call messages_fatal(3)
          end if
          
          SAFE_ALLOCATE(target%grad_local_pot(1:geo%natoms, 1:gr%mesh%np, 1:gr%sb%dim))
          SAFE_ALLOCATE(vl(1:gr%mesh%np_part))
          SAFE_ALLOCATE(vl_grad(1:gr%mesh%np, 1:gr%sb%dim))
          SAFE_ALLOCATE(target%rho(1:gr%mesh%np))

          ! calculate gradient of each species potential
          do iatom=1, geo%natoms
             vl(:) = M_ZERO
             vl_grad(:,:) = M_ZERO
             call species_get_local(geo%atom(iatom)%spec, gr%mesh, geo%atom(iatom)%x(1:gr%sb%dim), vl, M_ZERO)
             call dderivatives_grad(gr%der, vl, vl_grad)
             forall(ist=1:gr%mesh%np, jst=1:gr%sb%dim)
                target%grad_local_pot(iatom, ist, jst) = vl_grad(ist, jst)
             end forall
          end do
          SAFE_DEALLOCATE_A(vl)
          SAFE_DEALLOCATE_A(vl_grad)
          
          ! fix gradient jump on boundary - PRELIMINARY
          stencil_size = (real(gr%der%order+1)*gr%mesh%spacing(1))**2
          do ist=1, gr%mesh%np
             do jst=gr%mesh%np, gr%mesh%np_part
                point_dist = M_ZERO
                do ip=1, MAX_DIM
                   point_dist = point_dist +(gr%mesh%x(ist,ip)-gr%mesh%x(jst,ip))**2
                end do
                if(point_dist < stencil_size) then
                   target%grad_local_pot(1:geo%natoms, ist, 1:gr%sb%dim) = M_ZERO
                   exit
                end if
             end do
          end do
       end if
       
    case(oct_tg_current)
      message(1) =  'Info: Target is a current.'
      call messages_info(1)
      !%Variable OCTCurrentFunctional
      !%Type integer
      !%Section Calculation Modes::Optimal Control
      !%Default oct_no_curr
      !%Description
      !% The variable <tt>OCTCurrentFunctional</tt> describes which kind of current target functional is
      !% to be used. EXPERIMENTAL!
      !%Option oct_no_curr 0
      !% No current functional is used, no current calculated.
      !%Option oct_min_curr 1
      !% Minimizes current in all space by maximizing J1[j]= -<tt>OCTCurrentWeight</tt>*\int|j(r)|^2 dr. 
      !% Useful in combination with target density in order to obtain stable final target density.
      !% <tt>OCTCurrentWeight</tt> .GE. 0.
      !%Option oct_max_curr 2
      !% Maximizes current in all space by maximizing J1[j]= <tt>OCTCurrentWeight</tt>*\int|j(r)|^2 dr.
      !% Useful, for example, in combination with a target density in order to obtain a high-velocity impact.
      !% <tt>OCTCurrentWeight</tt> .GE. 0.
      !%Option oct_max_curr_ring 3
      !% Maximizes the current of a quantum ring in one direction. The functional maximixes the z projection of the 
      !% outer product between the position \vec{r} and the current \vec{j}: 
      !% J1[j]=<tt>OCTCurrentWeight</tt>*\int (\vec{r} \times \vec{j}) \hat{z} dr. For <tt>OCTCurrentWeight</tt> .GT. 0. the
      !% current flows in counter clockwise direction, while for <tt>OCTCurrentWeight</tt> .LT. 0 the current is clockwise.
      !%Option oct_min_curr_td 4
      !% Minimizes current in time interval [<tt>OCTStartTimeCurrTg</tt>, 
      !% total time = <tt>TDMaximumIter</tt> * <tt>TDTimeStep</tt>]. 
      !% Set <tt>TDEvolutionMethod </tt> = <tt>crank_nicholson</tt>  
      !%End 

      !%Variable OCTCurrentWeight
      !%Type float
      !%Section Calculation Modes::Optimal Control
      !%Default 0.0
      !%Description
      !% In the case of simultaneous optimization of density n and current j, one can tune the importance
      !% of the minimization of the current, as the respective functionals might not provide results on the
      !% same scale of magnitude. J1[n,j]= J1d[n]+ OCTCurrentWeight * J1c[j]. 
      !% For minimal current, keep <tt>OCTCurrentWeight</tt> .GE. 0, in order to have a well defined functional to maximize. 
      !% For maximal current on ring a positive <tt>OCTCurrentWeight</tt> will lead to counter clockwise current, a negative 
      !% one to clockwise current.
      !%End

      !%Variable OCTStartIterCurrTg
      !%Type integer
      !%Section Calculation Modes::Optimal Control
      !%Default 0
      !%Description
      !% Allows for a time dependent target for the current without defining it for the total 
      !% time-interval of the simulation.
      !% Thus it can be switched on at the iteration desired, <tt>OCTStartIterCurrTg</tt> .ge. 0
      !% and  <tt>OCTStartIterCurrTg</tt> .lt. <tt>TDMaximumIter</tt>. 
      !% Tip: If you would like to specify a real time for switching
      !% the functional on rather than the number of steps, just use something
      !% like:
      !% <tt>OCTStartIterCurrTg</tt> = 100.0 / <tt>TDTimeStep</tt>
      !%End
            
      call parse_integer(datasets_check('OCTCurrentFunctional'), oct_no_curr, target%curr_functional)
      if(target%curr_functional .eq. oct_no_curr .or. &
        target%curr_functional .gt. M_THREE) then
        message(1) = 'If you choose the current as a target,'
        message(2) = 'define 1 .le. "OCTCurrentFunctional" .le. 3.'
        call messages_fatal(2)
      end if
      
      call parse_float(datasets_check('OCTCurrentWeight'), M_ZERO, target%curr_weight)
      select case(target%curr_functional)
      case (oct_min_curr, oct_min_curr_td, oct_max_curr) 
        if( target%curr_weight .lt. M_ZERO) then
          message(1) = 'For "OCTCurrentFunctional" = oct_min_curr(_td)'
          message(2) = 'set "OCTCurrentWeight" .ge. 0.0'
          call messages_fatal(2)
        end if
      end select 
       
      write(message(1), '(a,i3)')   'Info: OCTCurrentFunctional = ', target%curr_functional
      write(message(2), '(a,f8.3)') 'Info: OCTCurrentWeight = ',  target%curr_weight
      call messages_info(2)
      
      call parse_integer(datasets_check('OCTStartIterCurrTg'), 0, target%strt_iter_curr_tg)
      if (target_mode(target) .eq. oct_targetmode_td) then
        write(message(1), '(a,i3)')   'Info: TargetMode = ', target_mode(target)
        write(message(2), '(a,i8)') 'Info: OCTStartIterCurrTg = ',  target%strt_iter_curr_tg
        call messages_info(2)
        target%dt = td%dt
        SAFE_ALLOCATE(target%td_fitness(0:td%max_iter))
        target%td_fitness = M_ZERO
      else
        target%strt_iter_curr_tg = M_ZERO
      end if
      if (target%strt_iter_curr_tg .lt. M_ZERO) then
        message(1) = 'OCTStartIterCurrTg must be negative'
        call messages_fatal(1)
      end if
      if (target%strt_iter_curr_tg .ge. td%max_iter) then
        message(1) = 'OCTStartIterCurrTg has to be .lt. TDMaximumIter'
        call messages_fatal(1)
      end if

      SAFE_ALLOCATE(stin%current( 1:gr%mesh%np_part, 1:gr%mesh%sb%dim, 1:stin%d%nspin ) )
      stin%current= M_ZERO

    case default
      write(message(1),'(a)') "Target Operator not properly defined."
      call messages_fatal(1)
    end select

    ! suppress current functional for all other targets other than density and/or current.
    if(target%type .ne. oct_tg_current .and. &
       target%type .ne. oct_tg_density) then
      target%curr_functional = oct_no_curr
      target%curr_weight = M_ZERO
    end if
   
       
       

    POP_SUB(target_init)
  end subroutine target_init
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  subroutine target_end(target, oct)
    type(target_t), intent(inout) :: target
    type(oct_t), intent(in)       :: oct

    PUSH_SUB(target_end)

    call states_end(target%st)
    if(target%type .eq. oct_tg_local .or. &
       target%type .eq. oct_tg_density .or. &
       target%type .eq. oct_tg_td_local) then
      SAFE_DEALLOCATE_P(target%rho)
    end if
    if(target%type .eq. oct_tg_hhg) then
      SAFE_DEALLOCATE_P(target%hhg_k)
      SAFE_DEALLOCATE_P(target%hhg_alpha)
      SAFE_DEALLOCATE_P(target%hhg_a)
    end if
    if(target_mode(target).eq.oct_targetmode_td) then
      SAFE_DEALLOCATE_P(target%td_fitness)
    end if
    if(target_type(target).eq.oct_tg_velocity) then
       if(oct%algorithm .eq. oct_algorithm_cg) then
          SAFE_DEALLOCATE_P(target%vel_der_array)
          SAFE_DEALLOCATE_P(target%grad_local_pot)
          SAFE_DEALLOCATE_P(target%rho)
       end if
    end if
    
    POP_SUB(target_end)
  end subroutine target_end
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  subroutine target_output(target, gr, dir, geo, outp)
    type(target_t), intent(inout) :: target
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(geometry_t),       intent(in)  :: geo
    type(output_t),         intent(in)  :: outp

    integer :: ierr
    PUSH_SUB(target_output)
    
    call loct_mkdir(trim(dir))

    select case(target%type)
    case(oct_tg_local)
      call doutput_function(outp%how, trim(dir), 'local_target', gr%mesh, &
        target%rho, units_out%length**(-gr%sb%dim), ierr, geo = geo)
    case(oct_tg_td_local)
      call target_build_tdlocal(target, gr, M_ZERO)
      call doutput_function(outp%how, trim(dir), 'td_local_target', gr%mesh, &
        target%rho, units_out%length**(-gr%sb%dim), ierr, geo = geo)
    case(oct_tg_density)
      call doutput_function(outp%how, trim(dir), 'density_target', gr%mesh, &
        target%rho, units_out%length**(-gr%sb%dim), ierr, geo = geo)
    case(oct_tg_excited)
      call output_states(target%est%st, gr, geo, trim(dir)//'/st', outp)
      call excited_states_output(target%est, trim(dir))
    case default
      call output_states(target%st, gr, geo, trim(dir), outp)
    end select

    POP_SUB(target_output)
  end subroutine target_output
  ! ----------------------------------------------------------------------


  ! ---------------------------------------------------------
  ! Calculates, at a given point in time marked by the integer
  ! index, the integrand of the target functional:
  ! <Psi(t)|\hat{O}(t)|Psi(t)>.
  ! ---------------------------------------------------------
  subroutine target_tdcalc(target, gr, psi, time)
    type(target_t), intent(inout) :: target
    type(grid_t),   intent(in)    :: gr
    type(states_t), intent(inout) :: psi
    integer,        intent(in)    :: time

    CMPLX, allocatable :: opsi(:, :)
    FLOAT, allocatable :: multipole(:, :)
    integer :: ist, ip, is

    if(target_type(target)  .eq. oct_tg_velocity) return
    if(target_mode(target)  .ne. oct_targetmode_td) return

    PUSH_SUB(target_tdcalc)

    target%td_fitness(time) = M_ZERO

    select case(target%type)
    case(oct_tg_td_local)

      !!!! WARNING Here one should build the time-dependent target.
      select case(psi%d%ispin)
      case(UNPOLARIZED)
        ASSERT(psi%d%nik .eq. 1)
        SAFE_ALLOCATE(opsi(1:gr%mesh%np_part, 1:1))
        opsi = M_z0
        do ist  = psi%st_start, psi%st_end
          do ip = 1, gr%mesh%np
            opsi(ip, 1) = target%rho(ip) * psi%zpsi(ip, 1, ist, 1)
          end do
          target%td_fitness(time) = &
            target%td_fitness(time) + psi%occ(ist, 1) * zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, ist, 1), opsi(:, :))
        end do
        SAFE_DEALLOCATE_A(opsi)
      case(SPIN_POLARIZED)
        message(1) = 'Error in target.target_tdcalc: spin_polarized.'
        call messages_fatal(1)
      case(SPINORS)
        message(1) = 'Error in target.target_tdcalc: spinors.'
        call messages_fatal(1)
      end select

    case(oct_tg_hhg)

      SAFE_ALLOCATE(multipole(1:4, 1:psi%d%nspin))
      do is = 1, psi%d%nspin
        call dmf_multipoles(gr%mesh, psi%rho(:, is), 1, multipole(:, is))
      end do
      target%td_fitness(time) = sum(multipole(2, 1:psi%d%spin_channels))
      SAFE_DEALLOCATE_A(multipole)

    ! case oct_tg_density only active if combined with td current functional
    case(oct_tg_current, oct_tg_density)

      if (time .ge. target%strt_iter_curr_tg) then
        target%td_fitness(time) = jcurr_functional(target, gr, psi)
      end if 

    case default
      message(1) = 'Error in target.target_tdcalc: default.'
      call messages_fatal(1)
    end select

    POP_SUB(target_tdcalc)
  end subroutine target_tdcalc
  ! ----------------------------------------------------------------------


  ! ---------------------------------------------------------------
  ! Calculates the inhomogeneous term that appears in the equation
  ! for chi, and places it into inh.
  ! ---------------------------------------------------------------
  subroutine target_inh(psi, gr, target, time, inh)
    type(states_t),    intent(inout)     :: psi
    type(grid_t),      intent(in)        :: gr
    type(target_t),    intent(inout)     :: target
    FLOAT,             intent(in)        :: time
    type(states_t),    intent(inout)     :: inh
 
    integer :: ik, ist, idim, ip
    
    PUSH_SUB(target_inh)

    select case(target%type)
    case(oct_tg_td_local)
      call target_build_tdlocal(target, gr, time)
      forall(ik = 1:inh%d%nik, ist = inh%st_start:inh%st_end, idim = 1:inh%d%dim, ip = 1:gr%mesh%np)
        inh%zpsi(ip, idim, ist, ik) = -target%rho(ip) * psi%zpsi(ip, idim, ist, ik)
      end forall
      
    case(oct_tg_velocity)
      ! set -(dF/dn)*psi_j as inhomogenous term --> !!! D_KS(r,t) is not implemented yet !!!
      forall(ik = 1:inh%d%nik, ist = inh%st_start:inh%st_end, idim = 1:inh%d%dim, ip = 1:gr%mesh%np)
         inh%zpsi(ip, idim, ist, ik) = -target%rho(ip) * psi%zpsi(ip, idim, ist, ik)
      end forall
   
    case(oct_tg_current, oct_tg_density)
      if (abs(nint(time/target%dt)) .ge. target%strt_iter_curr_tg) then
        inh%zpsi =  -chi_current(target, gr, psi)
      else
        inh%zpsi = M_ZERO
      end if     
  
    end select

    POP_SUB(target_inh)
  end subroutine target_inh
  !----------------------------------------------------------


  !----------------------------------------------------------
  subroutine target_build_tdlocal(target, gr, time)
    type(target_t), intent(inout) :: target
    type(grid_t),   intent(in)    :: gr
    FLOAT,          intent(in)    :: time

    integer :: ip
    FLOAT :: xx(MAX_DIM), rr, re, im

    PUSH_SUB(target_build_tdlocal)

    do ip = 1, gr%mesh%np
      call mesh_r(gr%mesh, ip, rr, coords = xx)
      call parse_expression(re, im, gr%sb%dim, xx, rr, time, target%td_local_target)
      target%rho(ip) = re
    end do

    POP_SUB(target_build_tdlocal)
  end subroutine target_build_tdlocal
  !----------------------------------------------------------


  ! ---------------------------------------------------------
  ! Calculates the J1 functional, i.e.:
  ! <Psi(T)|\hat{O}|Psi(T) in the time-independent
  ! case, or else \int_0^T dt <Psi(t)|\hat{O}(t)|Psi(t) in 
  ! the time-dependent case.
  ! ---------------------------------------------------------
  FLOAT function j1_functional(target, gr, psi, geo) result(j1)
    type(target_t), intent(inout)   :: target
    type(grid_t),   intent(inout)   :: gr
    type(states_t), intent(inout)   :: psi
    type(geometry_t), intent(in), optional :: geo

    integer :: ip, ist, iter, jj, maxiter, ik
    FLOAT :: omega, aa, maxhh, ww, currfunc_tmp
    FLOAT, allocatable :: local_function(:)
    CMPLX, allocatable :: ddipole(:)
    CMPLX, allocatable :: opsi(:, :)
    FLOAT :: f_re, dummy(3)
    character(len=4096) :: inp_string

    PUSH_SUB(j1_functional)

    j1 = M_ZERO
    select case(target%type)
    case(oct_tg_density)

      SAFE_ALLOCATE(local_function(1:gr%mesh%np))
      do ip = 1, gr%mesh%np
        local_function(ip) = - ( sqrt(psi%rho(ip, 1)) - sqrt(target%rho(ip)) )**2
      end do
      j1 = dmf_integrate(gr%mesh, local_function)
      SAFE_DEALLOCATE_A(local_function)

    case(oct_tg_local)

      select case(psi%d%ispin)
      case(UNPOLARIZED)
        ASSERT(psi%d%nik .eq. 1)
        SAFE_ALLOCATE(opsi(1:gr%mesh%np_part, 1:1))
        opsi = M_z0
        j1 = M_ZERO
        do ist = psi%st_start, psi%st_end
          do ip = 1, gr%mesh%np
            opsi(ip, 1) = target%rho(ip) * psi%zpsi(ip, 1, ist, 1)
          end do
          j1 = j1 + psi%occ(ist, 1) * zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, ist, 1), opsi(:, :))
        end do
        SAFE_DEALLOCATE_A(opsi)
      case(SPIN_POLARIZED)
        message(1) = 'Error in target.j1_functional: spin_polarized.'
        call messages_fatal(1)
      case(SPINORS)
        message(1) = 'Error in target.j1_functional: spinors.'
        call messages_fatal(1)
      end select

    case(oct_tg_td_local)
      maxiter = size(target%td_fitness) - 1
      j1 = M_HALF * target%dt * target%td_fitness(0) + & 
           M_HALF * target%dt * target%td_fitness(maxiter) + & 
           target%dt * sum(target%td_fitness(1:maxiter-1))

    case(oct_tg_excited)
      j1 = abs(zstates_mpdotp(gr%mesh, target%est, psi))**2

    case(oct_tg_exclude_state)

      j1 = M_ONE
      do ist = 1, target%st%nst
        if(loct_isinstringlist(ist, target%excluded_states_list)) then
          j1 = j1 - abs(zmf_dotp(gr%mesh, psi%d%dim, &
            target%st%zpsi(:, :, ist, 1), psi%zpsi(:, :, 1, 1)))**2
        end if
      end do

    case(oct_tg_hhg)

      maxiter = size(target%td_fitness) - 1
      SAFE_ALLOCATE(ddipole(0:maxiter))

      ddipole(0) = M_ZERO
      do iter = 1, maxiter - 1
        ddipole(iter) = (target%td_fitness(iter - 1) + target%td_fitness(iter + 1) - &
                      M_TWO * target%td_fitness(iter)) / target%dt**2
      end do
      call interpolate( target%dt*(/ -3, -2, -1 /),   &
                        ddipole(maxiter - 3:maxiter - 1), &
                        M_ZERO, &
                        ddipole(maxiter) )

      call spectrum_hsfunction_init(target%dt, 0, maxiter, maxiter, ddipole)
      do jj = 1, target%hhg_nks
        aa = target%hhg_a(jj) * target%hhg_w0
        ww = target%hhg_k(jj) * target%hhg_w0
        call spectrum_hsfunction_min(ww - aa, ww + aa, ww, omega, maxhh)
        j1 = j1 + target%hhg_alpha(jj) * log(-maxhh)
      end do
      call spectrum_hsfunction_end()

      SAFE_DEALLOCATE_A(ddipole)

    case(oct_tg_velocity)
      f_re = M_ZERO
      dummy(:) = M_ZERO
      inp_string = target%vel_input_string
      call parse_velocity_target(inp_string, geo)
      call conv_to_C_string(inp_string)
      call parse_expression(f_re, dummy(1), 1, dummy(1:3), dummy(1), dummy(1), inp_string)
      j1 = f_re

    case(oct_tg_current)
      ! calculate functional out of this select case(target%type) block,
      ! so it can be combined with other target%type.
      j1 = M_ZERO

    case default
      do ik = 1, psi%d%nik
        do ist = psi%st_start, psi%st_end
          j1 = j1 + psi%occ(ist, ik) * &
            abs(zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, ist, ik), &
                target%st%zpsi(:, :, ist, ik)))**2
        end do
      end do


    end select

    ! current functionals are conveniently combined with others
    if (target%type .eq. oct_tg_current .or. &
        target%curr_functional .ne. oct_no_curr) then
      select case(target_mode(target))
      case(oct_targetmode_static)
        currfunc_tmp = jcurr_functional(target, gr, psi )
      case(oct_targetmode_td)
        maxiter = size(target%td_fitness) - 1
        currfunc_tmp = M_HALF * target%dt * target%td_fitness(target%strt_iter_curr_tg) + & 
                       M_HALF * target%dt * target%td_fitness(maxiter) + & 
                       target%dt * sum(target%td_fitness(target%strt_iter_curr_tg+1:maxiter-1))  
      end select
      if(conf%devel_version) then
        write(message(1), '(6x,a,f12.5)')    " => Other functional   = ", j1
        write(message(2), '(6x,a,f12.5)')    " => Current functional = ", currfunc_tmp
        call messages_info(2)
      end if
      ! accumulating functional values
      j1 = j1 + currfunc_tmp
    end if 

    POP_SUB(j1_functional)
  end function j1_functional
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! calculate |chi(T)> = \hat{O}(T) |psi(T)>
  ! ---------------------------------------------------------
  subroutine calc_chi(target, gr, psi_in, chi_out, geo)
    type(target_t),    intent(inout) :: target
    type(grid_t),      intent(inout) :: gr
    type(states_t),    intent(inout) :: psi_in
    type(states_t),    intent(inout) :: chi_out
    type(geometry_t),  intent(in)    :: geo
    
    CMPLX :: olap, zdet
    CMPLX, allocatable :: cI(:), dI(:), mat(:, :, :), mm(:, :, :, :), mk(:, :), lambda(:, :)
    integer :: ik, ip, ist, jst, idim, no_electrons, ia, ib, n_pairs, nst, kpoints, jj
    character(len=1024) :: temp_string
    FLOAT :: df_dv, dummy(3)

    PUSH_SUB(calc_chi)

    no_electrons = -nint(psi_in%val_charge)

    select case(target%type)

    case(oct_tg_density)

      select case(psi_in%d%ispin)
      case(UNPOLARIZED)

        ASSERT(psi_in%d%nik .eq. 1)

        if(no_electrons .eq. 1) then
          do ip = 1, gr%mesh%np
            chi_out%zpsi(ip, 1, 1, 1) = sqrt(target%rho(ip)) * &
              exp(M_zI * atan2(aimag(psi_in%zpsi(ip, 1, 1, 1)), &
                                real(psi_in%zpsi(ip, 1, 1, 1) )) )
          end do
        else
          do ist = psi_in%st_start, psi_in%st_end
            do ip = 1, gr%mesh%np
              if(psi_in%rho(ip, 1) > CNST(1.0e-8)) then
                chi_out%zpsi(ip, 1, ist, 1) = sqrt(target%rho(ip) / psi_in%rho(ip, 1)) * &
                  psi_in%zpsi(ip, 1, ist, 1)
              else
                chi_out%zpsi(ip, 1, ist, 1) = M_ZERO !sqrt(target%rho(ip))
              end if
            end do
          end do
        end if

      case(SPIN_POLARIZED)
         message(1) = 'Error in target.calc_chi: spin_polarized.'
         call messages_fatal(1)
      case(SPINORS)
         message(1) = 'Error in target.calc_chi: spinors.'
         call messages_fatal(1)
      end select

    case(oct_tg_local)

      select case(psi_in%d%ispin)
      case(UNPOLARIZED)
        ASSERT(psi_in%d%nik .eq. 1)
        do ist = psi_in%st_start, psi_in%st_end
          do ip = 1, gr%mesh%np
            chi_out%zpsi(ip, 1, ist, 1) = target%rho(ip) * psi_in%zpsi(ip, 1, ist, 1)
          end do
        end do
      case(SPIN_POLARIZED)
         message(1) = 'Error in target.calc_chi: spin_polarized.'
         call messages_fatal(1)
      case(SPINORS)
         message(1) = 'Error in target.calc_chi: spinors.'
         call messages_fatal(1)
      end select

    case(oct_tg_td_local)
      !We assume that there is no time-independent operator.
      forall(ik = 1:chi_out%d%nik, ist = chi_out%st_start:chi_out%st_end, idim = 1:chi_out%d%dim, ip = 1:gr%mesh%np)
        chi_out%zpsi(ip, idim, ist, ik) = M_z0
      end forall

    case(oct_tg_excited) 

      n_pairs = target%est%n_pairs
      kpoints = psi_in%d%nik
      nst = psi_in%nst

      SAFE_ALLOCATE(cI(1:n_pairs))
      SAFE_ALLOCATE(dI(1:n_pairs))
      SAFE_ALLOCATE(mat(1:target%est%st%nst, 1:nst, 1:psi_in%d%nik))
      SAFE_ALLOCATE(mm(1:nst, 1:nst, 1:kpoints, 1:n_pairs))
      SAFE_ALLOCATE(mk(1:gr%mesh%np_part, 1:psi_in%d%dim))
      SAFE_ALLOCATE(lambda(1:n_pairs, 1:n_pairs))

      call zstates_matrix(gr%mesh, target%est%st, psi_in, mat)

      do ia = 1, n_pairs
        cI(ia) = target%est%weight(ia)
        call zstates_matrix_swap(mat, target%est%pair(ia))
        mm(1:nst, 1:nst, 1:kpoints, ia) = mat(1:nst, 1:kpoints, 1:kpoints)
        dI(ia) = zstates_mpdotp(gr%mesh, target%est%st, psi_in, mat)
        if(abs(dI(ia)) > CNST(1.0e-12)) then
          do ik = 1, kpoints
            zdet = lalg_inverter(nst, mm(1:nst, 1:nst, ik, ia))
          end do
        end if
        call zstates_matrix_swap(mat, target%est%pair(ia))
      end do

      do ia = 1, n_pairs
        do ib = 1, n_pairs
          lambda(ia, ib) = conjg(cI(ib)) * cI(ia) * conjg(dI(ia)) * dI(ib)
        end do
      end do

      select case(psi_in%d%ispin)
      case(UNPOLARIZED)
        write(message(1), '(a)') 'Internal error in target.calc_chi: unpolarized.'
        call messages_fatal(1)

      case(SPIN_POLARIZED)
        ASSERT(chi_out%d%nik .eq. 2)

        do ik = 1, kpoints
          do ist = chi_out%st_start, chi_out%st_end
            chi_out%zpsi(:, :, ist, ik) = M_z0
            do ia = 1, n_pairs
              if(ik .ne. target%est%pair(ia)%sigma) cycle
              if(abs(dI(ia)) < CNST(1.0e-12)) cycle
              do ib = 1, n_pairs
                if(abs(dI(ib)) < CNST(1.0e-12)) cycle
                mk = M_z0
                do jst = 1, nst
                  if(jst .eq. target%est%pair(ib)%i) jj = target%est%pair(ia)%a
                  mk(:, :) = mk(:, :) + conjg(mm(ist, jst, ik, ib)) * target%est%st%zpsi(:, :, jj, ik)
                end do
                call lalg_axpy(gr%mesh%np_part, psi_in%d%dim, M_z1, lambda(ib, ia) * mk(:, :), chi_out%zpsi(:, :, ist, ik))
              end do
            end do
          end do
        end do
        
      case(SPINORS)
        ASSERT(chi_out%d%nik .eq. 1)

        do ist = chi_out%st_start, chi_out%st_end
          chi_out%zpsi(:, :, ist, 1) = M_z0

          do ia = 1, n_pairs
            if(abs(dI(ia)) < CNST(1.0e-12)) cycle

            do ib = 1, n_pairs
              if(abs(dI(ib)) < CNST(1.0e-12)) cycle

              mk = M_z0
              do jst = 1, nst
                if(jst .eq. target%est%pair(ib)%i) jj = target%est%pair(ia)%a
                mk(:, :) = mk(:, :) + conjg(mm(ist, jst, 1, ib)) * target%est%st%zpsi(:, :, jj, 1)
              end do

              call lalg_axpy(gr%mesh%np_part, 2, M_z1, lambda(ib, ia) * mk(:, :), chi_out%zpsi(:, :, ist, 1))
            end do
          end do
        end do

      end select

      SAFE_DEALLOCATE_A(cI)
      SAFE_DEALLOCATE_A(dI)
      SAFE_DEALLOCATE_A(mat)
      SAFE_DEALLOCATE_A(mm)
      SAFE_DEALLOCATE_A(mk)
      SAFE_DEALLOCATE_A(lambda)

    case(oct_tg_exclude_state)

      chi_out%zpsi(:, :, 1, 1) = psi_in%zpsi(:, :, 1, 1)
      do ist = 1, target%st%nst
        if(loct_isinstringlist(ist, target%excluded_states_list)) then
          olap = zmf_dotp(gr%mesh, psi_in%d%dim, target%st%zpsi(:, :, ist, 1), psi_in%zpsi(:, :, 1, 1))
          chi_out%zpsi(:, :, 1, 1) = chi_out%zpsi(:, :, 1, 1) - olap * target%st%zpsi(:, :, ist, 1)
        end if
      end do

    case(oct_tg_velocity)
      !we have a time-dependent target --> Chi(T)=0
      forall(ip=1:gr%mesh%np, idim=1:chi_out%d%dim, ist=chi_out%st_start:chi_out%st_end, ik=1:chi_out%d%nik)
         chi_out%zpsi(ip, idim, ist, ik) = M_z0
      end forall
      
      !calculate dF/dn, which is the time-independent part of the inhomogenous term for the propagation of Chi
      !NOTE: D_KS(r,t) is not calculated here
      df_dv = M_ZERO
      dummy(:) = M_ZERO
      target%rho(:) = M_ZERO
      do ist=1, geo%natoms
         do jst=1, gr%sb%dim
            temp_string = target%vel_der_array(ist, jst)
            call parse_velocity_target(temp_string, geo)
            call conv_to_C_string(temp_string)
            call parse_expression(df_dv, dummy(1), 1, dummy(1:3), dummy(1), dummy(1), temp_string)
            target%rho(:) = target%rho(:) + df_dv*target%grad_local_pot(ist,:,jst)/species_weight(geo%atom(ist)%spec)
         end do
      end do
    
    case(oct_tg_current)
      ! the chi_out%zpsi is calculated outside this select block by
      ! an accumulating sum in order to combine it with other targets.
      chi_out%zpsi = M_ZERO
  
    case default

      !olap = zstates_mpdotp(gr%mesh, target%st, psi_in)
      do ik = 1, psi_in%d%nik
        do ist = psi_in%st_start, psi_in%st_end
          olap = zmf_dotp(gr%mesh, target%st%zpsi(:, 1, ist, ik), psi_in%zpsi(:, 1, ist, ik))
          chi_out%zpsi(:, :, ist, ik) = olap * target%st%zpsi(:, :, ist, ik)
        end do
      end do

    end select
    

   ! boundary conditions of static current functionals are combined with others.
    ! chi_out%zpsi is accumulated.
    if(target%type .eq. oct_tg_current .or. &
      target%curr_functional .ne. oct_no_curr) then
      if (target_mode(target) .eq. oct_targetmode_static ) then
        chi_out%zpsi = chi_out%zpsi + chi_current(target, gr, psi_in)
      end if
    end if 


    POP_SUB(calc_chi)
  end subroutine calc_chi


  ! ----------------------------------------------------------------------
  integer pure function target_mode(target)
    type(target_t), intent(in) :: target

    select case(target%type)
    case(oct_tg_td_local, oct_tg_hhg, oct_tg_velocity)
      target_mode = oct_targetmode_td
    case default
      target_mode = oct_targetmode_static
    end select

    ! allow specific current functionals to be td
    ! Attention: yet combined with static density target,
    ! the total target is considered td.
    select case(target%curr_functional)
    case(oct_min_curr_td) 
      target_mode = oct_targetmode_td
    end select
  end function target_mode


  ! ----------------------------------------------------------------------
  integer pure function target_type(target)
    type(target_t), intent(in) :: target

    target_type = target%type

  end function target_type
  ! ----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  integer pure function current_functional_type(target)
    type(target_t), intent(in) :: target
    
    current_functional_type = target%curr_functional
 
  end function current_functional_type
  !-----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  logical pure function target_move_ions(target)
    type(target_t), intent(in) :: target
    
    target_move_ions = target%move_ions

  end function target_move_ions
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  ! replaces the "v[:,:]" from inp_string by the corresponding
  ! numbers from geo%atom(:)%v(:) so that the parser can handle inp_string
  ! ----------------------------------------------------------------------
  subroutine parse_velocity_target(inp_string, geo)
    character(len=*), intent(inout)  :: inp_string  ! input string --> OCTVelocityTarget
    type(geometry_t), intent(in)     :: geo         ! velocities of the atoms --> geo%atom(n_atom)%v(coord)
    integer              :: i,m,n_atom,coord,string_length
    character (LEN=100)  :: v_string

    PUSH_SUB(parse_velocity_target)
    
    string_length = len(inp_string)
    do i=1, string_length - 1
       if(inp_string(i:i+1) == "v[") then
          m = 0
          if(inp_string(i+3:i+3) == ",") m = 1
          if(inp_string(i+4:i+4) == ",") m = 2
          if(m == 0) then
             message(1) = "OCTVelocityTarget Input error!"
             message(2) = "Atom number is either larger than 99 or not defined."
             call messages_fatal(2)
          end if
          read(inp_string(i+2:i+1+m),*) n_atom
          read(inp_string(i+3+m:i+3+m),*) coord
          if(coord < 1 .or. coord > 3) then
             message(1) = "OCTVelocityTarget Input error!"
             message(2) = "Vector component is either larger than 3 or smaller than 1."
             call messages_fatal(2)
          end if
          write(v_string,*) geo%atom(n_atom)%v(coord)
          inp_string = inp_string(:i-1) // "(" // trim(v_string) // ")" // inp_string(i+5+m:)
       end if
    end do

    POP_SUB(parse_velocity_target)
  end subroutine parse_velocity_target
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  ! Calculates a current functional that may be combined with
  ! other functionals found in function target.j1_functional.
  ! ----------------------------------------------------------------------
  FLOAT function jcurr_functional(target, gr, psi) result(jcurr)
    type(target_t), intent(in)    :: target
    type(grid_t),   intent(in)    :: gr
    type(states_t), intent(inout) :: psi

    integer :: ip
    FLOAT, allocatable :: semilocal_function(:)
 
    PUSH_SUB(jcurr_functional)
    
    jcurr = M_ZERO    
    ASSERT(psi%d%nik .eq. 1)
    call states_calc_quantities(gr%der, psi, paramagnetic_current=psi%current) 
    SAFE_ALLOCATE(semilocal_function(1:gr%mesh%np))
    semilocal_function = M_ZERO

    select case(target%curr_functional)
    case(oct_no_curr)
      semilocal_function = M_ZERO

    case(oct_min_curr, oct_min_curr_td)
      do ip = 1, gr%mesh%np
        semilocal_function(ip) = - ( sum(psi%current(ip, 1:gr%sb%dim, 1)**2) ) 
      end do
      
    case(oct_max_curr)
      do ip = 1, gr%mesh%np
        semilocal_function(ip) = sum(psi%current(ip, 1:gr%sb%dim, 1)**2)
      end do

    case(oct_max_curr_ring)
      if(gr%sb%dim .ne. M_TWO) then
        write(message(1), '(a)') 'This target only implemented for 2D.'
        call messages_fatal(1)
      end if
      do ip = 1, gr%mesh%np
        ! func = j_y * x - j_x * y 
        semilocal_function (ip) = psi%current(ip, 2, 1) * gr%mesh%spacing(1) * gr%mesh%x(ip,1) -  &
                                    psi%current(ip, 1, 1) * gr%mesh%spacing(2) * gr%mesh%x(ip,2)
      end do
    case default
      message(1) = 'Error in target.jcurr_functional: chosen target does not exist'
      call messages_fatal(1)
    end select

    jcurr = target%curr_weight * dmf_integrate(gr%mesh, semilocal_function)
        
    SAFE_DEALLOCATE_A(semilocal_function)

    POP_SUB(jcurr_functional)
  end function jcurr_functional
  !-----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  ! Calculates current-specific boundary condition
  !-----------------------------------------------------------------------
  function chi_current(target, gr, psi_in) result(chi)
    type(target_t),    intent(in)    :: target
    type(grid_t),      intent(in)    :: gr
    type(states_t),    intent(inout) :: psi_in
    CMPLX                            :: chi( size(psi_in%zpsi,1), 1, size(psi_in%zpsi,3), 1)
        
    CMPLX, allocatable :: grad_psi_in(:,:,:)
    FLOAT, allocatable :: div_curr_psi_in(:,:)   
    integer :: ip, ist

    PUSH_SUB(chi_current)

    chi = M_ZERO    

    SAFE_ALLOCATE(grad_psi_in(1:gr%der%mesh%np_part, 1:gr%der%mesh%sb%dim, 1))

    if(target_mode(target) .eq. oct_targetmode_td ) then 
      call states_calc_quantities(gr%der, psi_in, paramagnetic_current=psi_in%current) 
    end if

    select case(target%curr_functional)
    case(oct_no_curr)
      chi(ip, 1, ist, 1) = M_ZERO

    case(oct_min_curr, oct_min_curr_td, oct_max_curr)
      SAFE_ALLOCATE(div_curr_psi_in(1:gr%der%mesh%np_part,1))
      call dderivatives_div(gr%der, psi_in%current(:,:,1), div_curr_psi_in(1:gr%der%mesh%np,1)) 
        
      do ist = psi_in%st_start, psi_in%st_end
        call zderivatives_grad(gr%der, psi_in%zpsi(:,1, ist, 1), grad_psi_in(:,:,1))
            
        do ip = 1, gr%mesh%np 
          chi(ip, 1, ist, 1) =  M_zI * target%curr_weight * &
               ( M_TWO * sum(psi_in%current(ip, 1:gr%sb%dim, 1) * grad_psi_in(ip, 1:gr%sb%dim, 1))+ &
               div_curr_psi_in(ip,1) * psi_in%zpsi(ip, 1, ist, 1) )
        end do
      end do
      SAFE_DEALLOCATE_A(div_curr_psi_in)

    case(oct_max_curr_ring)
      do ist = psi_in%st_start, psi_in%st_end

        call zderivatives_grad(gr%der, psi_in%zpsi(:,1, ist, 1), grad_psi_in(:,:,1))
            
        do ip = 1, gr%mesh%np 
          chi(ip, 1, ist, 1) =  M_zI * target%curr_weight * &
               ( grad_psi_in(ip, 1, 1) * gr%mesh%spacing(2) * gr%mesh%x(ip,2) - &
               grad_psi_in(ip, 2, 1) * gr%mesh%spacing(1) * gr%mesh%x(ip,1)   ) 
        end do
      end do

    case default
      message(1) = 'Error in target.chi_current: chosen target does not exist'
      call messages_fatal(1)
    end select
        
    SAFE_DEALLOCATE_A(grad_psi_in)       

    POP_SUB(chi_current)
  end function chi_current
  ! ----------------------------------------------------------------------




end module opt_control_target_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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
!! $Id: target.F90 2870 2007-04-28 06:26:47Z acastro $

#include "global.h"

module opt_control_target_m
  use datasets_m
  use density_m
  use derivatives_m
  use epot_m
  use excited_states_m
  use fft_m
  use forces_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use opt_control_global_m
  use output_m
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
  use td_calc_m
  use td_m
  use types_m
  use unit_m
  use unit_system_m
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
            target_j1,               &
            target_chi,              &
            target_move_ions,        &
            target_curr_functional,  &
            target_init_propagation


  integer, public, parameter ::       &
    oct_tg_groundstate      = 1,      &
    oct_tg_excited          = 2,      &
    oct_tg_gstransformation = 3,      &
    oct_tg_userdefined      = 4,      &
    oct_tg_jdensity         = 5,      &        
    oct_tg_local            = 6,      &
    oct_tg_td_local         = 7,      &
    oct_tg_exclude_state    = 8,      &
    oct_tg_hhg              = 9,      &
    oct_tg_velocity         = 10,     &
    oct_tg_hhgnew           = 12


  integer, public, parameter ::       &
    oct_targetmode_static = 0,        &
    oct_targetmode_td     = 1

  integer, public, parameter ::       &
    oct_no_curr              = 0,     &
    oct_curr_square          = 1,     &
    oct_max_curr_ring        = 2,     &
    oct_curr_square_td       = 3

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
    FLOAT   :: density_weight
    FLOAT   :: curr_weight
    integer :: strt_iter_curr_tg
    FLOAT, pointer :: spatial_curr_wgt(:) => null()
    character(len=1000) :: plateau_string
    CMPLX, pointer :: acc(:, :)
    CMPLX, pointer :: vel(:, :)
    CMPLX, pointer :: gvec(:, :)
    FLOAT, pointer :: alpha(:)
    type(fft_t) :: fft_handler
  end type target_t


contains

  ! ---------------------------------------------------------
  !> This routine performs all the things that must be initialized
  !! prior to a forward evolution, regarding the target. Right now
  !! some of those initizalizations are not done here, and should
  !! be moved.
  subroutine target_init_propagation(tg)
    type(target_t), intent(inout)    :: tg
    PUSH_SUB(target_init_propagation)

    select case(tg%type)
    case(oct_tg_hhgnew)
      tg%vel = M_z0
      tg%gvec = M_z0
      tg%acc = M_z0
    end select

    POP_SUB(target_init_propagation)
  end subroutine target_init_propagation


  ! ----------------------------------------------------------------------
  !> This just copies the states_t variable present in target, into st.
  subroutine target_get_state(tg, st)
    type(target_t), intent(in)    :: tg
    type(states_t), intent(inout) :: st

    PUSH_SUB(target_get_state)
    call states_copy(st, tg%st)

    POP_SUB(target_get_state)
  end subroutine target_get_state
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> The target is initialized, mainly by reading from the inp file.
  subroutine target_init(gr, geo, stin, td, w0, tg, oct, ep)
    type(grid_t),     intent(in)    :: gr
    type(geometry_t), intent(in)    :: geo
    type(states_t),   intent(inout) :: stin 
    type(td_t),       intent(in)    :: td
    FLOAT,            intent(in)    :: w0
    type(target_t),   intent(inout) :: tg
    type(oct_t),      intent(in)    :: oct
    type(epot_t),     intent(inout) :: ep

    integer             :: ierr, ip, ist, jst, jj, iatom, ib, idim, inst, inik, &
                           id, ik, no_states, no_constraint, no_ptpair, cstr_dim(MAX_DIM), iunit, &
                           nn(3), optimize_parity(3)
    logical             :: optimize(3)
    type(block_t)       :: blk
    FLOAT               :: xx(MAX_DIM), rr, psi_re, psi_im, xstart, xend, dw, ww
    FLOAT, allocatable  :: vl(:), vl_grad(:,:), xp(:), tmp_box(:,:)
    FLOAT               :: fact
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
    !%Option oct_tg_jdensity 5
    !% EXPERIMENTAL: 
    !%Option oct_tg_local 6
    !% The target operator is a local operator.
    !%Option oct_tg_td_local 7
    !% The target operator is a time-dependent local operator.
    !%Option oct_tg_exclude_state 8
    !% Target operator is the projection onto the complement of a given state, given by the
    !% block <tt>OCTTargetTransformStates</tt>. This means that the target operator is the unity
    !% operator minus the projector onto that state.
    !%Option oct_tg_hhg 9
    !% The target is the optimization of the HHG yield. You must supply the OCTOptimizeHarmonicSpectrum
    !% block, and it attempts to optimize te maximum of the spectrum around each harmonic peak. You may
    !% use only one of the gradient-less optimization schemes.
    !%Option oct_tg_velocity 10
    !% The target is a function of the velocities of the nuclei at the end of the influence of
    !% the external field, defined by <tt>OCTVelocityTarget</tt>
    !%Option oct_tg_hhgnew 12
    !% EXPERIMENTAL: The  target is the optimization of the HHG yield. You must supply the
    !% OCTHarmonicWeigth string. It attempts to optimized the integral of the harmonic spectrum multiplied
    !% by some user defined weight function.
    !%End
    call parse_integer(datasets_check('OCTTargetOperator'), oct_tg_gstransformation, tg%type)
    if(.not.varinfo_valid_option('OCTTargetOperator', tg%type)) &
      call input_error('OCTTargetOperator')

    call states_copy(tg%st, stin)
    call states_deallocate_wfns(tg%st)
    call states_allocate_wfns(tg%st, gr%mesh, TYPE_CMPLX)

    nullify(tg%td_fitness)


    ! WARNING This should probably go somewhere else
    tg%curr_functional = oct_no_curr
    tg%curr_weight = M_ZERO
    tg%strt_iter_curr_tg = 0

    select case(tg%type)
    case(oct_tg_groundstate)
      message(1) =  'Info: Using Ground State for TargetOperator'
      call messages_info(1)
      call restart_read(trim(restart_dir)//GS_DIR, tg%st, gr, ierr, exact = .true.)
      
    case(oct_tg_excited) 

      message(1) =  'Info: TargetOperator is a linear combination of Slater determinants.'
      call messages_info(1)

      call states_look (trim(restart_dir)//GS_DIR, gr%mesh%mpi_grp, ip, ip, tg%st%nst, ierr)
      tg%st%st_start = 1
      tg%st%st_end   = tg%st%nst

      SAFE_DEALLOCATE_P(tg%st%occ)
      SAFE_DEALLOCATE_P(tg%st%eigenval)
      SAFE_DEALLOCATE_P(tg%st%node)

      SAFE_ALLOCATE(     tg%st%occ(1:tg%st%nst, 1:tg%st%d%nik))
      SAFE_ALLOCATE(tg%st%eigenval(1:tg%st%nst, 1:tg%st%d%nik))
      SAFE_ALLOCATE(    tg%st%node(1:tg%st%nst))
      if(tg%st%d%ispin == SPINORS) then
        SAFE_DEALLOCATE_P(tg%st%spin)
        SAFE_ALLOCATE(tg%st%spin(1:3, 1:tg%st%nst, 1:tg%st%d%nik))
      end if
      call states_allocate_wfns(tg%st, gr%mesh, TYPE_CMPLX)
      tg%st%node(:)  = 0

      call restart_read(trim(restart_dir)//GS_DIR, tg%st, gr, ierr, exact = .true.)
      call excited_states_init(tg%est, tg%st, "oct-excited-state-target") 

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
      call parse_string(datasets_check('OCTExcludedStates'), "1", tg%excluded_states_list)
      call states_deallocate_wfns(tg%st)
      call restart_look_and_read(tg%st, gr)

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
      if(parse_isdef(datasets_check('OCTTargetTransformStates')) /= 0) then
        if(parse_block(datasets_check('OCTTargetTransformStates'), blk) == 0) then
          call states_copy(tmp_st, tg%st)
          call states_deallocate_wfns(tmp_st)
          call restart_look_and_read(tmp_st, gr)
          SAFE_ALLOCATE(rotation_matrix(1:tg%st%nst, 1:tmp_st%nst))
          rotation_matrix = M_z0
          do ist = 1, tg%st%nst
            do jst = 1, parse_block_cols(blk, ist - 1)
              call parse_block_cmplx(blk, ist - 1, jst - 1, rotation_matrix(ist, jst))
            end do
          end do

          call states_rotate(gr%mesh, tg%st, tmp_st, rotation_matrix)
          SAFE_DEALLOCATE_A(rotation_matrix)
          call states_end(tmp_st)
          call parse_block_end(blk)
          call density_calc(tg%st, gr, tg%st%rho)
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
      !% Define a target state. Syntax follows the one of the <tt>UserDefinedStates</tt> block.
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
          do id = 1, tg%st%d%dim
            do ist = 1, tg%st%nst
              do ik = 1, tg%st%d%nik   
                
                ! does the block entry match and is this node responsible?
                if(.not. (id  ==  idim .and. ist  ==  inst .and. ik  ==  inik    &
                  .and. tg%st%st_start  <=  ist .and. tg%st%st_end >= ist) ) cycle
                
                ! parse formula string
                call parse_block_string(                            &
                  blk, ib - 1, 3, tg%st%user_def_states(id, ist, ik))
                ! convert to C string
                call conv_to_C_string(tg%st%user_def_states(id, ist, ik))
                
                do ip = 1, gr%mesh%np
                  xx = gr%mesh%x(ip, :)
                  rr = sqrt(sum(xx(:)**2))
                  
                  ! parse user-defined expressions
                  call parse_expression(psi_re, psi_im, &
                    gr%sb%dim, xx, rr, M_ZERO, tg%st%user_def_states(id, ist, ik))
                  ! fill state
                  tg%st%zpsi(ip, id, ist, ik) = psi_re + M_zI * psi_im
                end do
                ! normalize orbital
                call zstates_normalize_orbital(gr%mesh, tg%st%d%dim, &
                  tg%st%zpsi(:,:, ist, ik))
              end do
            end do
          enddo
        end do
        call parse_block_end(blk)
        call density_calc(tg%st, gr, tg%st%rho)
      else
        message(1) = '"OCTTargetUserdefined" has to be specified as block.'
        call messages_fatal(1)
      end if

    case(oct_tg_jdensity)
      call target_init_density(gr, tg, stin, td)

    case(oct_tg_local)
      !%Variable OCTLocalTarget
      !%Type string
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If <tt>OCTTargetOperator = oct_tg_local</tt>, then one must supply a function
      !% that defines the target. This should be done by defining it through a string, using 
      !% the variable <tt>OCTLocalTarget</tt>.
      !%End

      if(parse_isdef('OCTLocalTarget') /= 0) then
        SAFE_ALLOCATE(tg%rho(1:gr%mesh%np))
        tg%rho = M_ZERO
        call parse_string(datasets_check('OCTLocalTarget'), "0", expression)
        call conv_to_C_string(expression)
        do ip = 1, gr%mesh%np
          call mesh_r(gr%mesh, ip, rr, coords = xx)
          ! parse user-defined expression
          call parse_expression(psi_re, psi_im, gr%sb%dim, xx, rr, M_ZERO, expression)
          tg%rho(ip) = psi_re
        end do
      else
        message(1) = 'If OCTTargetOperator = oct_tg_local, then you must give the shape'
        message(2) = 'of this target in variable "OCTLocalTarget".'
        call messages_fatal(2)
      end if

    case(oct_tg_td_local)
      if(parse_block(datasets_check('OCTTdTarget'),blk)==0) then
        call parse_block_string(blk, 0, 0, tg%td_local_target)
        call conv_to_C_string(tg%td_local_target)
        SAFE_ALLOCATE(tg%rho(1:gr%mesh%np))
        call parse_block_end(blk)
      else
        message(1) = 'If OCTTargetMode = oct_targetmode_td, you must suppy a OCTTDTarget block.'
        call messages_fatal(1)
      end if
      tg%dt = td%dt
      SAFE_ALLOCATE(tg%td_fitness(0:td%max_iter))
      tg%td_fitness = M_ZERO
      call target_build_tdlocal(tg, gr, M_ZERO)

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
      if(parse_isdef(datasets_check('OCTOptimizeHarmonicSpectrum')) /= 0) then
        if(parse_block(datasets_check('OCTOptimizeHarmonicSpectrum'), blk) == 0) then
          tg%hhg_nks = parse_block_cols(blk, 0)
          SAFE_ALLOCATE(    tg%hhg_k(1:tg%hhg_nks))
          SAFE_ALLOCATE(tg%hhg_alpha(1:tg%hhg_nks))
          SAFE_ALLOCATE(    tg%hhg_a(1:tg%hhg_nks))
          do jj = 1, tg%hhg_nks
            call parse_block_integer(blk, 0, jj - 1, tg%hhg_k(jj))
            call parse_block_float(blk, 1, jj - 1, tg%hhg_alpha(jj))
            call parse_block_float(blk, 2, jj - 1, tg%hhg_a(jj))
          end do
          call parse_block_end(blk)
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

      tg%hhg_w0 = w0
      tg%dt     = td%dt
      SAFE_ALLOCATE(tg%td_fitness(0:td%max_iter))
      tg%td_fitness = M_ZERO

    case(oct_tg_hhgnew)
       
      if(parse_isdef('OCTMoveIons')  ==  0) then
         message(1) = 'If OCTTargetOperator = oct_tg_hhgnew, then you must supply'
         message(2) = 'the variable "OCTMoveIons".'
         call messages_fatal(2)
      else
         call parse_logical(datasets_check('OCTMoveIons'), .false., tg%move_ions)
      end if
      
       ! We allocate many things that are perhaps not necessary if we use a direct optimization scheme.


      SAFE_ALLOCATE(tg%vel(td%max_iter+1, MAX_DIM))
      SAFE_ALLOCATE(tg%acc(td%max_iter+1, MAX_DIM))
      SAFE_ALLOCATE(tg%gvec(td%max_iter+1, MAX_DIM))
      SAFE_ALLOCATE(tg%alpha(td%max_iter))
          
      ! The following is a temporary hack, that assumes only one atom at the origin of coordinates.
      if(geo%natoms > 1) then
        message(1) = 'If "OCTTargetOperator = oct_tg_hhgnew", then you can only have one atom.'
        call messages_fatal(1)
      end if

      SAFE_ALLOCATE(tg%grad_local_pot(1:geo%natoms, 1:gr%mesh%np, 1:gr%sb%dim))
      SAFE_ALLOCATE(vl(1:gr%mesh%np_part))
      SAFE_ALLOCATE(vl_grad(1:gr%mesh%np, 1:gr%sb%dim))
      SAFE_ALLOCATE(tg%rho(1:gr%mesh%np))

      vl(:) = M_ZERO
      vl_grad(:,:) = M_ZERO
      call epot_local_potential(ep, gr%der, gr%dgrid, geo, 1, vl)
      call dderivatives_grad(gr%der, vl, vl_grad)
      forall(ist=1:gr%mesh%np, jst=1:gr%sb%dim)
        tg%grad_local_pot(1, ist, jst) = vl_grad(ist, jst)
      end forall

      ! Note that the calculation of the gradient of the potential
      ! is wrong at the borders of the box, since it assumes zero boundary
      ! conditions. The best way to solve this problems is to define the 
      ! target making use of the definition of the forces based on the gradient
      ! of the density, rather than on the gradient of the potential.
          
      !%Variable OCTHarmonicWeight
      !%Type string
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% EXPERIMENTAL: If "OCTTargetOperator = oct_tg_plateau", then the function to optimize is the integral of the
      !% harmonic spectrum H(w), weighted with a function f(w) that is defined as a string here. For example, if 
      !% you set OCTHarmonicWeight  = "step(w-1)", the function to optimize is the integral of step(w-1)*H(w) or, i.e.
      !% the integral of H(w) from one to infinity. In practice, it is better if you also set an upper limit, i.e.
      !% for example f(w) = step(w-1)*step(2-w).
      !%End
      call parse_string(datasets_check('OCTHarmonicWeight'), "1", tg%plateau_string)
      tg%dt = td%dt
      SAFE_ALLOCATE(tg%td_fitness(0:td%max_iter))
      tg%td_fitness = M_ZERO

      iunit = io_open('.alpha', action = 'write')
      dw = (M_TWO * M_PI) / (td%max_iter * tg%dt)
      do jj = 0, td%max_iter - 1
        ww = jj * dw
        call parse_expression(psi_re, psi_im, "w", ww, tg%plateau_string)
        tg%alpha(jj+1) = psi_re
        write(iunit, *) ww, psi_re
      end do
      call io_close(iunit)

      nn(1:3) = (/ td%max_iter, 1, 1 /)
      optimize(1:3) = .false.
      optimize_parity(1:3) = -1
      call fft_init(tg%fft_handler, nn(1:3), 1, FFT_COMPLEX, FFTLIB_FFTW, optimize, optimize_parity)


    case(oct_tg_velocity)
      !%Variable OCTVelocityTarget
      !%Type block
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If <tt>OCTTargetOperator = oct_tg_velocity</tt>, then one must supply the 
      !% target to optimize in terms of the ionic velocities. This is done by 
      !% supplying a string through the block <tt>OCTVelocityTarget</tt>.
      !% Each velocity component is supplied by <tt>"v[n_atom,vec_comp]"</tt>,
      !% while "n_atom" is the respective atom number, corresponding to the 
      !% <tt>Coordinates</tt> block and "vec_comp" is the corresponding
      !% vector component of the velocity. The target string can be
      !% supplied by using several lines in the OCTTargetOperator block.
      !% As an example, the following target can be used to maximize the
      !% velocity difference between atom 1 and 2 (in a 3D system):
      !%
      !% <tt>%OCTVelocityTarget
      !% <br> "(v[1,1]-v[2,1])^2 + (v[1,2]-v[2,2])^2 + "
      !% <br> "(v[1,3]-v[2,3])^2"
      !% <br>%</tt>
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
      !% The derivatives are supplied via strings through the block
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
      !% <tt>%OCTVelocityDerivatives
      !% <br> " 2*(v[1,1]-v[2,1])" | " 2*(v[1,2]-v[2,2])" | " 2*(v[1,3]-v[2,3])"
      !% <br> "-2*(v[1,1]-v[2,1])" | "-2*(v[1,2]-v[2,2])" | "-2*(v[1,3]-v[2,3])"
      !% <br>%</tt>
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
          tg%vel_input_string = " "
          do jj=0, parse_block_n(blk)-1
             call parse_block_string(blk, jj, 0, expression)
             tg%vel_input_string = trim(tg%vel_input_string) // trim(expression)
          end do
          call parse_block_end(blk)
       else
          message(1) = 'If OCTTargetOperator = oct_tg_velocity, then you must give the shape'
          message(2) = 'of this target in the block "OCTVelocityTarget".'
          call messages_fatal(2)
       end if
       
       if(parse_isdef('OCTMoveIons')  ==  0) then
          message(1) = 'If OCTTargetOperator = oct_tg_velocity, then you must supply'
          message(2) = 'the variable "OCTMoveIons".'
          call messages_fatal(2)
       else
          call parse_logical(datasets_check('OCTMoveIons'), .false., tg%move_ions)
       end if
       
       if(oct%algorithm  ==  oct_algorithm_cg) then
          if(parse_block(datasets_check('OCTVelocityDerivatives'),blk)==0) then
             SAFE_ALLOCATE(tg%vel_der_array(1:geo%natoms,1:gr%sb%dim))
             do ist=0, geo%natoms-1
                do jst=0, gr%sb%dim-1
                   call parse_block_string(blk, ist, jst, tg%vel_der_array(ist+1, jst+1))
                end do
             end do
             call parse_block_end(blk)
          else
             message(1) = 'If OCTTargetOperator = oct_tg_velocity, and'
             message(2) = 'OCTScheme = oct_algorithm_cg, then you must define the'
             message(3) = 'blocks "OCTVelocityTarget" AND "OCTVelocityDerivatives"'
             call messages_fatal(3)
          end if
          
          SAFE_ALLOCATE(tg%grad_local_pot(1:geo%natoms, 1:gr%mesh%np, 1:gr%sb%dim))
          SAFE_ALLOCATE(vl(1:gr%mesh%np_part))
          SAFE_ALLOCATE(vl_grad(1:gr%mesh%np, 1:gr%sb%dim))
          SAFE_ALLOCATE(tg%rho(1:gr%mesh%np))

          ! calculate gradient of each species potential
          do iatom=1, geo%natoms
             vl(:) = M_ZERO
             vl_grad(:,:) = M_ZERO
             call epot_local_potential(ep, gr%der, gr%dgrid, geo, iatom, vl)
             call dderivatives_grad(gr%der, vl, vl_grad)
             forall(ist=1:gr%mesh%np, jst=1:gr%sb%dim)
                tg%grad_local_pot(iatom, ist, jst) = vl_grad(ist, jst)
             end forall
          end do
          SAFE_DEALLOCATE_A(vl)
          SAFE_DEALLOCATE_A(vl_grad)

          ! Note that the calculation of the gradient of the potential
          ! is wrong at the borders of the box, since it assumes zero boundary
          ! conditions. The best way to solve this problems is to define the 
          ! target making use of the definition of the forces based on the gradient
          ! of the density, rather than on the gradient of the potential.
          
       end if

       tg%dt = td%dt
       SAFE_ALLOCATE(tg%td_fitness(0:td%max_iter))
       tg%td_fitness = M_ZERO
       
    case default
      write(message(1),'(a)') "Target Operator not properly defined."
      call messages_fatal(1)
    end select

    POP_SUB(target_init)
  end subroutine target_init
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  subroutine target_end(tg, oct)
    type(target_t), intent(inout) :: tg
    type(oct_t), intent(in)       :: oct

    PUSH_SUB(target_end)

    call states_end(tg%st)
    if(tg%type  ==  oct_tg_local .or. &
       tg%type  ==  oct_tg_jdensity .or. &
       tg%type  ==  oct_tg_td_local) then
      SAFE_DEALLOCATE_P(tg%rho)
    end if
    if(tg%type  ==  oct_tg_hhg) then
      SAFE_DEALLOCATE_P(tg%hhg_k)
      SAFE_DEALLOCATE_P(tg%hhg_alpha)
      SAFE_DEALLOCATE_P(tg%hhg_a)
    end if
    if(target_mode(tg) == oct_targetmode_td) then
      SAFE_DEALLOCATE_P(tg%td_fitness)
    end if
    if(target_type(tg) == oct_tg_velocity) then
       if(oct%algorithm  ==  oct_algorithm_cg) then
          SAFE_DEALLOCATE_P(tg%vel_der_array)
          SAFE_DEALLOCATE_P(tg%grad_local_pot)
          SAFE_DEALLOCATE_P(tg%rho)
       end if
    end if
    if(target_type(tg) == oct_tg_hhgnew) then
       if(oct%algorithm  ==  oct_algorithm_cg) then
          SAFE_DEALLOCATE_P(tg%grad_local_pot)
          SAFE_DEALLOCATE_P(tg%rho)
          SAFE_DEALLOCATE_P(tg%vel)
          SAFE_DEALLOCATE_P(tg%acc)
          SAFE_DEALLOCATE_P(tg%gvec)
          SAFE_DEALLOCATE_P(tg%alpha)
          call fft_end(tg%fft_handler)
       end if
    end if
    if(tg%type  ==  oct_tg_jdensity) then
      SAFE_DEALLOCATE_P(tg%spatial_curr_wgt)
    end if
    POP_SUB(target_end)
  end subroutine target_end
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  subroutine target_output(tg, gr, dir, geo, outp)
    type(target_t), intent(inout) :: tg
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(geometry_t),       intent(in)  :: geo
    type(output_t),         intent(in)  :: outp

    integer :: ierr
    PUSH_SUB(target_output)
    
    call loct_mkdir(trim(dir))

    select case(tg%type)
    case(oct_tg_local)
      if(outp%how /= 0) then
        call dio_function_output(outp%how, trim(dir), 'local_target', gr%mesh, &
          tg%rho, units_out%length**(-gr%sb%dim), ierr, geo = geo)
      end if
    case(oct_tg_td_local)
      call target_build_tdlocal(tg, gr, M_ZERO)
      if(outp%how /= 0) then
        call dio_function_output(outp%how, trim(dir), 'td_local_target', gr%mesh, &
          tg%rho, units_out%length**(-gr%sb%dim), ierr, geo = geo)
      end if
    case(oct_tg_jdensity)
      if(outp%how /= 0) then
        if(tg%density_weight > M_ZERO) then
          call dio_function_output(outp%how, trim(dir), 'density_target', gr%mesh, &
            tg%rho, units_out%length**(-gr%sb%dim), ierr, geo = geo)
        end if
      end if
    case(oct_tg_excited)
      call output_states(tg%est%st, gr, geo, trim(dir)//'/st', outp)
      call excited_states_output(tg%est, trim(dir))
    case default
      call output_states(tg%st, gr, geo, trim(dir), outp)
    end select

    POP_SUB(target_output)
  end subroutine target_output
  ! ----------------------------------------------------------------------


  ! ---------------------------------------------------------
  !> Calculates, at a given point in time marked by the integer
  !! index, the integrand of the target functional:
  !! <Psi(t)|\hat{O}(t)|Psi(t)>.
  subroutine target_tdcalc(tg, hm, gr, geo, psi, time, max_time)
    type(target_t),      intent(inout) :: tg
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(inout) :: geo
    type(states_t),      intent(inout) :: psi
    integer,             intent(in)    :: time
    integer,             intent(in)    :: max_time

    CMPLX, allocatable :: opsi(:, :)
    integer :: ist, ip, ia, iw
    FLOAT :: acc(MAX_DIM), dt, dw
    integer :: iatom, idim, ik

    if(target_mode(tg)  /= oct_targetmode_td) return

    PUSH_SUB(target_tdcalc)

    tg%td_fitness(time) = M_ZERO

    select case(tg%type)

    case(oct_tg_hhgnew)

      ! If the ions move, the tg is computed in the propagation routine.
      if(.not.target_move_ions(tg)) then

        SAFE_ALLOCATE(opsi(1:gr%mesh%np_part, 1:1))

        opsi = M_z0
        ! WARNING This does not work for spinors.
        ! The following is a temporary hack. It assumes only one atom at the origin.
        acc = M_ZERO
        do ik = 1, psi%d%nik
          do ist = 1, psi%nst
            do idim = 1, gr%sb%dim
              opsi(1:gr%mesh%np, 1) = tg%grad_local_pot(1, 1:gr%mesh%np, idim) * psi%zpsi(1:gr%mesh%np, 1, ist, ik)
              acc(idim) = acc(idim) + real( psi%occ(ist, ik) * &
                  zmf_dotp(gr%mesh, psi%d%dim, opsi, psi%zpsi(:, :, ist, ik)), REAL_PRECISION )
              tg%acc(time+1, idim) = tg%acc(time+1, idim) + psi%occ(ist, ik) * &
                  zmf_dotp(gr%mesh, psi%d%dim, opsi, psi%zpsi(:, :, ist, ik))
            end do
          end do
        end do

        SAFE_DEALLOCATE_A(opsi)
      end if

      dt = tg%dt
      dw = (M_TWO * M_PI/(max_time * tg%dt))
      if(time  ==  max_time) then
        tg%acc(1, 1:gr%sb%dim) = M_HALF * (tg%acc(1, 1:gr%sb%dim) + tg%acc(max_time+1, 1:gr%sb%dim))
        do ia = 1, gr%sb%dim
          call zfft_forward1(tg%fft_handler, tg%acc(1:max_time, ia), tg%vel(1:max_time, ia))
        end do
        tg%vel = tg%vel * tg%dt
        do iw = 1, max_time
          ! We add the one-half dt term because when doing the propagation we want the value at interpolated times.
          tg%acc(iw, 1:gr%sb%dim) = tg%vel(iw, 1:gr%sb%dim) * tg%alpha(iw) * exp(M_zI * (iw-1) * dw * M_HALF * dt)
        end do
        do ia = 1, gr%sb%dim
          call zfft_backward1(tg%fft_handler, tg%acc(1:max_time, ia), tg%gvec(1:max_time, ia))
        end do
        tg%gvec(max_time + 1, 1:gr%sb%dim) = tg%gvec(1, 1:gr%sb%dim)
        tg%gvec = tg%gvec * (M_TWO * M_PI/ tg%dt)
      end if


    case(oct_tg_velocity)

      ! If the ions move, the target is computed in the propagation routine.
      if(.not.target_move_ions(tg)) then

      SAFE_ALLOCATE(opsi(1:gr%mesh%np_part, 1:1))
      opsi = M_z0
      ! WARNING This does not work for spinors.
      do iatom = 1, geo%natoms
        geo%atom(iatom)%f(1:gr%sb%dim) = hm%ep%fii(1:gr%sb%dim, iatom)
        do ik = 1, psi%d%nik
          do ist = 1, psi%nst
            do idim = 1, gr%sb%dim
              opsi(1:gr%mesh%np, 1) = tg%grad_local_pot(iatom, 1:gr%mesh%np, idim) * psi%zpsi(1:gr%mesh%np, 1, ist, ik)
              geo%atom(iatom)%f(idim) = geo%atom(iatom)%f(idim) + real(psi%occ(ist, ik) * &
                zmf_dotp(gr%mesh, psi%d%dim, opsi, psi%zpsi(:, :, ist, ik)), REAL_PRECISION)
            end do
          end do
        end do
      end do
      SAFE_DEALLOCATE_A(opsi)

      end if

      dt = tg%dt
      if( (time  ==  0) .or. (time  ==  max_time) ) dt = tg%dt * M_HALF
      do iatom = 1, geo%natoms
         geo%atom(iatom)%v(1:MAX_DIM) = geo%atom(iatom)%v(1:MAX_DIM) + &
           geo%atom(iatom)%f(1:MAX_DIM) * dt / species_weight(geo%atom(iatom)%spec)
      end do

    case(oct_tg_td_local)

      !!!! WARNING Here one should build the time-dependent target.
      select case(psi%d%ispin)
      case(UNPOLARIZED)
        ASSERT(psi%d%nik  ==  1)
        SAFE_ALLOCATE(opsi(1:gr%mesh%np_part, 1:1))
        opsi = M_z0
        do ist  = psi%st_start, psi%st_end
          do ip = 1, gr%mesh%np
            opsi(ip, 1) = tg%rho(ip) * psi%zpsi(ip, 1, ist, 1)
          end do
          tg%td_fitness(time) = &
            tg%td_fitness(time) + psi%occ(ist, 1) * &
              real(zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, ist, 1), opsi(:, :)), REAL_PRECISION)
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

     call td_calc_tacc(gr, geo, psi, hm, acc, time*tg%dt)
     tg%td_fitness(time) = acc(1)

    case(oct_tg_jdensity)

      if (time >= tg%strt_iter_curr_tg) then
        tg%td_fitness(time) = jcurr_functional(tg, gr, psi)
      end if 

    case default
      message(1) = 'Error in target.target_tdcalc: default.'
      call messages_fatal(1)
    end select

    POP_SUB(target_tdcalc)
  end subroutine target_tdcalc
  ! ----------------------------------------------------------------------


  ! ---------------------------------------------------------------
  !> Calculates the inhomogeneous term that appears in the equation
  !! for chi, and places it into inh.
  subroutine target_inh(psi, gr, tg, time, inh, iter)
    type(states_t),    intent(inout)     :: psi
    type(grid_t),      intent(in)        :: gr
    type(target_t),    intent(inout)     :: tg
    FLOAT,             intent(in)        :: time
    type(states_t),    intent(inout)     :: inh
    integer,           intent(in)        :: iter
 
    integer :: ik, ist, ip, idim
    CMPLX :: gvec(MAX_DIM)

    PUSH_SUB(target_inh)

    select case(tg%type)
    case(oct_tg_td_local)
      call target_build_tdlocal(tg, gr, time)
      forall(ik = 1:inh%d%nik, ist = inh%st_start:inh%st_end, idim = 1:inh%d%dim, ip = 1:gr%mesh%np)
        inh%zpsi(ip, idim, ist, ik) = - psi%occ(ist, ik) * tg%rho(ip) * psi%zpsi(ip, idim, ist, ik)
      end forall

    case(oct_tg_hhgnew)
      gvec(1:gr%sb%dim) = real(tg%gvec(iter+1, 1:gr%sb%dim), REAL_PRECISION)
      forall(ik = 1:inh%d%nik, ist = inh%st_start:inh%st_end, idim = 1:inh%d%dim, ip = 1:gr%mesh%np)
        inh%zpsi(ip, idim, ist, ik) = &
           - psi%occ(ist, ik) * M_TWO * sum(tg%grad_local_pot(1, ip, 1:gr%sb%dim) * gvec(1:gr%sb%dim)) * &
           psi%zpsi(ip, idim, ist, ik)
      end forall

    case(oct_tg_velocity)
      forall(ik = 1:inh%d%nik, ist = inh%st_start:inh%st_end, idim = 1:inh%d%dim, ip = 1:gr%mesh%np)
         inh%zpsi(ip, idim, ist, ik) = - psi%occ(ist, ik) * tg%rho(ip) * psi%zpsi(ip, idim, ist, ik)
      end forall
   
    case(oct_tg_jdensity)
      if (abs(nint(time/tg%dt)) >= tg%strt_iter_curr_tg) then
        inh%zpsi =  -chi_current(tg, gr, psi)
      else
        inh%zpsi = M_ZERO
      end if     

    case default
      write(message(1),'(a)') 'Internal error in target_inh'
      call messages_fatal(1)
  
    end select

    POP_SUB(target_inh)
  end subroutine target_inh
  !----------------------------------------------------------


  !----------------------------------------------------------
  subroutine target_build_tdlocal(tg, gr, time)
    type(target_t), intent(inout) :: tg
    type(grid_t),   intent(in)    :: gr
    FLOAT,          intent(in)    :: time

    integer :: ip
    FLOAT :: xx(MAX_DIM), rr, re, im

    PUSH_SUB(target_build_tdlocal)

    do ip = 1, gr%mesh%np
      call mesh_r(gr%mesh, ip, rr, coords = xx)
      call parse_expression(re, im, gr%sb%dim, xx, rr, time, tg%td_local_target)
      tg%rho(ip) = re
    end do

    POP_SUB(target_build_tdlocal)
  end subroutine target_build_tdlocal
  !----------------------------------------------------------


  ! ---------------------------------------------------------
  !> Calculates the J1 functional, i.e.:
  !! <Psi(T)|\hat{O}|Psi(T) in the time-independent
  !! case, or else \int_0^T dt <Psi(t)|\hat{O}(t)|Psi(t) in 
  !! the time-dependent case.
  FLOAT function target_j1(tg, gr, psi, geo) result(j1)
    type(target_t), intent(inout)   :: tg
    type(grid_t),   intent(inout)   :: gr
    type(states_t), intent(inout)   :: psi
    type(geometry_t), intent(in), optional :: geo

    integer :: is, ip, ist, jj, maxiter, ik, i
    FLOAT :: omega, aa, maxhh, ww, currfunc_tmp, dw
    FLOAT, allocatable :: local_function(:)
    CMPLX, allocatable :: ddipole(:)
    FLOAT :: f_re, dummy(3)
    character(len=4096) :: inp_string

    PUSH_SUB(target_j1)

    j1 = M_ZERO
    select case(tg%type)
    case(oct_tg_jdensity)
      j1 = target_j1_density(gr, tg, psi)

    case(oct_tg_local)
      j1 = M_ZERO
      do is = 1, psi%d%spin_channels
        j1 = j1 + dmf_dotp(gr%mesh, tg%rho, psi%rho(:, is))
      end do

    case(oct_tg_td_local)
      maxiter = size(tg%td_fitness) - 1
      j1 = M_HALF * tg%dt * tg%td_fitness(0) + & 
           M_HALF * tg%dt * tg%td_fitness(maxiter) + & 
           tg%dt * sum(tg%td_fitness(1:maxiter-1))

    case(oct_tg_excited)
      j1 = abs(zstates_mpdotp(gr%mesh, tg%est, psi))**2

    case(oct_tg_exclude_state)

      j1 = M_ONE
      do ist = 1, tg%st%nst
        if(loct_isinstringlist(ist, tg%excluded_states_list)) then
          j1 = j1 - abs(zmf_dotp(gr%mesh, psi%d%dim, &
            tg%st%zpsi(:, :, ist, 1), psi%zpsi(:, :, 1, 1)))**2
        end if
      end do

    case(oct_tg_hhg)

      maxiter = size(tg%td_fitness) - 1
      SAFE_ALLOCATE(ddipole(0:maxiter))
      ddipole = M_z0
      ddipole = tg%td_fitness

      call spectrum_hsfunction_init(tg%dt, 0, maxiter, maxiter, ddipole)
      do jj = 1, tg%hhg_nks
        aa = tg%hhg_a(jj) * tg%hhg_w0
        ww = tg%hhg_k(jj) * tg%hhg_w0
        call spectrum_hsfunction_min(ww - aa, ww + aa, omega, maxhh)
        j1 = j1 + tg%hhg_alpha(jj) * log(-maxhh)
      end do
      call spectrum_hsfunction_end()

      SAFE_DEALLOCATE_A(ddipole)

    case(oct_tg_hhgnew)
      maxiter = size(tg%td_fitness) - 1
      dw = (M_TWO * M_PI) / (maxiter * tg%dt)
      j1 = M_ZERO
      do i = 0, maxiter - 1
        ww = i * dw
        j1 = j1 + dw * tg%alpha(i+1) * sum(abs(tg%vel(i+1, 1:gr%sb%dim))**2)
      end do

    case(oct_tg_velocity)
      f_re = M_ZERO
      dummy(:) = M_ZERO
      inp_string = tg%vel_input_string
      call target_parse_velocity(inp_string, geo)
      call conv_to_C_string(inp_string)
      call parse_expression(f_re, dummy(1), 1, dummy(1:3), dummy(1), dummy(1), inp_string)
      j1 = f_re

    case default
      do ik = 1, psi%d%nik
        do ist = psi%st_start, psi%st_end
          j1 = j1 + psi%occ(ist, ik) * &
            abs(zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, ist, ik), &
                tg%st%zpsi(:, :, ist, ik)))**2
        end do
      end do
    end select

    POP_SUB(target_j1)
  end function target_j1
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Calculate |chi(T)> = \hat{O}(T) |psi(T)>
  subroutine target_chi(tg, gr, psi_in, chi_out, geo)
    type(target_t),    intent(inout) :: tg
    type(grid_t),      intent(inout) :: gr
    type(states_t),    intent(inout) :: psi_in
    type(states_t),    intent(inout) :: chi_out
    type(geometry_t),  intent(in)    :: geo
    
    CMPLX :: olap, zdet
    CMPLX, allocatable :: cI(:), dI(:), mat(:, :, :), mm(:, :, :, :), mk(:, :), lambda(:, :)
    integer :: ik, ip, ist, jst, idim, no_electrons, ia, ib, n_pairs, nst, kpoints, jj
    character(len=1024) :: temp_string
    FLOAT :: df_dv, dummy(3)


    PUSH_SUB(target_chi)

    no_electrons = -nint(psi_in%val_charge)

    select case(tg%type)

    case(oct_tg_jdensity)
      call target_chi_density(tg, gr, psi_in, chi_out)

    case(oct_tg_local)
      do ik = 1, psi_in%d%nik
        do idim = 1, psi_in%d%dim
          do ist = psi_in%st_start, psi_in%st_end
            do ip = 1, gr%mesh%np
              chi_out%zpsi(ip, idim, ist, ik) = psi_in%occ(ist, ik) * tg%rho(ip) * psi_in%zpsi(ip, idim, ist, ik)
            end do
          end do
        end do
      end do

    case(oct_tg_td_local)
      !We assume that there is no time-independent operator.
      forall(ik = 1:chi_out%d%nik, ist = chi_out%st_start:chi_out%st_end, idim = 1:chi_out%d%dim, ip = 1:gr%mesh%np)
        chi_out%zpsi(ip, idim, ist, ik) = M_z0
      end forall

    case(oct_tg_excited) 

      n_pairs = tg%est%n_pairs
      kpoints = psi_in%d%nik
      nst = psi_in%nst

      SAFE_ALLOCATE(cI(1:n_pairs))
      SAFE_ALLOCATE(dI(1:n_pairs))
      SAFE_ALLOCATE(mat(1:tg%est%st%nst, 1:nst, 1:psi_in%d%nik))
      SAFE_ALLOCATE(mm(1:nst, 1:nst, 1:kpoints, 1:n_pairs))
      SAFE_ALLOCATE(mk(1:gr%mesh%np_part, 1:psi_in%d%dim))
      SAFE_ALLOCATE(lambda(1:n_pairs, 1:n_pairs))

      call zstates_matrix(gr%mesh, tg%est%st, psi_in, mat)

      do ia = 1, n_pairs
        cI(ia) = tg%est%weight(ia)
        call zstates_matrix_swap(mat, tg%est%pair(ia))
        mm(1:nst, 1:nst, 1:kpoints, ia) = mat(1:nst, 1:kpoints, 1:kpoints)
        dI(ia) = zstates_mpdotp(gr%mesh, tg%est%st, psi_in, mat)
        if(abs(dI(ia)) > CNST(1.0e-12)) then
          do ik = 1, kpoints
            zdet = lalg_inverter(nst, mm(1:nst, 1:nst, ik, ia))
          end do
        end if
        call zstates_matrix_swap(mat, tg%est%pair(ia))
      end do

      do ia = 1, n_pairs
        do ib = 1, n_pairs
          lambda(ia, ib) = conjg(cI(ib)) * cI(ia) * conjg(dI(ia)) * dI(ib)
        end do
      end do

      select case(psi_in%d%ispin)
      case(UNPOLARIZED)
        write(message(1), '(a)') 'Internal error in target.target_chi: unpolarized.'
        call messages_fatal(1)

      case(SPIN_POLARIZED)
        ASSERT(chi_out%d%nik  ==  2)

        do ik = 1, kpoints
          do ist = chi_out%st_start, chi_out%st_end
            chi_out%zpsi(:, :, ist, ik) = M_z0
            do ia = 1, n_pairs
              if(ik /= tg%est%pair(ia)%sigma) cycle
              if(abs(dI(ia)) < CNST(1.0e-12)) cycle
              do ib = 1, n_pairs
                if(abs(dI(ib)) < CNST(1.0e-12)) cycle
                mk = M_z0
                do jst = 1, nst
                  if(jst  ==  tg%est%pair(ib)%i) jj = tg%est%pair(ia)%a
                  mk(:, :) = mk(:, :) + conjg(mm(ist, jst, ik, ib)) * tg%est%st%zpsi(:, :, jj, ik)
                end do
                call lalg_axpy(gr%mesh%np_part, psi_in%d%dim, M_z1, lambda(ib, ia) * mk(:, :), chi_out%zpsi(:, :, ist, ik))
              end do
            end do
          end do
        end do
        
      case(SPINORS)
        ASSERT(chi_out%d%nik  ==  1)

        do ist = chi_out%st_start, chi_out%st_end
          chi_out%zpsi(:, :, ist, 1) = M_z0

          do ia = 1, n_pairs
            if(abs(dI(ia)) < CNST(1.0e-12)) cycle

            do ib = 1, n_pairs
              if(abs(dI(ib)) < CNST(1.0e-12)) cycle

              mk = M_z0
              do jst = 1, nst
                if(jst  ==  tg%est%pair(ib)%i) jj = tg%est%pair(ia)%a
                mk(:, :) = mk(:, :) + conjg(mm(ist, jst, 1, ib)) * tg%est%st%zpsi(:, :, jj, 1)
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
      do ist = 1, tg%st%nst
        if(loct_isinstringlist(ist, tg%excluded_states_list)) then
          olap = zmf_dotp(gr%mesh, psi_in%d%dim, tg%st%zpsi(:, :, ist, 1), psi_in%zpsi(:, :, 1, 1))
          chi_out%zpsi(:, :, 1, 1) = chi_out%zpsi(:, :, 1, 1) - olap * tg%st%zpsi(:, :, ist, 1)
        end if
      end do

    case(oct_tg_hhgnew)
      !we have a time-dependent target --> Chi(T)=0
      forall(ip=1:gr%mesh%np, idim=1:chi_out%d%dim, ist=chi_out%st_start:chi_out%st_end, ik=1:chi_out%d%nik)
         chi_out%zpsi(ip, idim, ist, ik) = M_z0
      end forall

    case(oct_tg_velocity)
      !we have a time-dependent target --> Chi(T)=0
      forall(ip=1:gr%mesh%np, idim=1:chi_out%d%dim, ist=chi_out%st_start:chi_out%st_end, ik=1:chi_out%d%nik)
         chi_out%zpsi(ip, idim, ist, ik) = M_z0
      end forall
      
      !calculate dF/dn, which is the time-independent part of the inhomogenous term for the propagation of Chi
      !NOTE: D_KS(r,t) is not calculated here
      df_dv = M_ZERO
      dummy(:) = M_ZERO
      tg%rho(:) = M_ZERO
      do ist=1, geo%natoms
         do jst=1, gr%sb%dim
            temp_string = tg%vel_der_array(ist, jst)
            call target_parse_velocity(temp_string, geo)
            call conv_to_C_string(temp_string)
            call parse_expression(df_dv, dummy(1), 1, dummy(1:3), dummy(1), dummy(1), temp_string)
            tg%rho(:) = tg%rho(:) + df_dv*tg%grad_local_pot(ist,:,jst)/species_weight(geo%atom(ist)%spec)
         end do
      end do
    
    case default

      !olap = zstates_mpdotp(gr%mesh, tg%st, psi_in)
      do ik = 1, psi_in%d%nik
        do ist = psi_in%st_start, psi_in%st_end
          olap = zmf_dotp(gr%mesh, tg%st%zpsi(:, 1, ist, ik), psi_in%zpsi(:, 1, ist, ik))
          chi_out%zpsi(:, :, ist, ik) = olap * tg%st%zpsi(:, :, ist, ik)
        end do
      end do

    end select

    POP_SUB(target_chi)
  end subroutine target_chi


  ! ----------------------------------------------------------------------
  integer pure function target_mode(tg)
    type(target_t), intent(in) :: tg

    select case(tg%type)
    case(oct_tg_td_local, oct_tg_hhg, oct_tg_velocity, oct_tg_hhgnew)
      target_mode = oct_targetmode_td
    case default
      target_mode = oct_targetmode_static
    end select

    ! allow specific current functionals to be td
    ! Attention: yet combined with static density target,
    ! the total target is considered td.
    select case(tg%curr_functional)
    case(oct_curr_square_td) 
      target_mode = oct_targetmode_td
    end select

  end function target_mode


  ! ----------------------------------------------------------------------
  integer pure function target_type(tg)
    type(target_t), intent(in) :: tg

    target_type = tg%type

  end function target_type
  ! ----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  integer pure function target_curr_functional(tg)
    type(target_t), intent(in) :: tg
   
    target_curr_functional = tg%curr_functional
 
  end function target_curr_functional
  !-----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  logical pure function target_move_ions(tg)
    type(target_t), intent(in) :: tg
    
    target_move_ions = tg%move_ions

  end function target_move_ions
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  logical pure function is_spatial_curr_wgt(tg)
    type(target_t), intent(in) :: tg

    is_spatial_curr_wgt = associated(tg%spatial_curr_wgt)
   
  end function is_spatial_curr_wgt
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  ! replaces the "v[:,:]" from inp_string by the corresponding
  ! numbers from geo%atom(:)%v(:) so that the parser can handle inp_string
  ! ----------------------------------------------------------------------
  subroutine target_parse_velocity(inp_string, geo)
    character(len=*), intent(inout)  :: inp_string  ! input string --> OCTVelocityTarget
    type(geometry_t), intent(in)     :: geo         ! velocities of the atoms --> geo%atom(n_atom)%v(coord)
    integer              :: i,m,n_atom,coord,string_length
    character (LEN=100)  :: v_string

    PUSH_SUB(target_parse_velocity)
    
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

    POP_SUB(target_parse_velocity)
  end subroutine target_parse_velocity
  ! ----------------------------------------------------------------------


#include "target_density_inc.F90"

end module opt_control_target_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

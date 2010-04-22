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
  use excited_states_m
  use geometry_m
  use global_m
  use grid_m
  use h_sys_output_m
  use io_m
  use io_function_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use parser_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use profiling_m
  use restart_m
  use spectrum_m
  use states_m
  use states_calc_m
  use states_dim_m
  use string_m
  use td_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public :: target_t,         &
            target_get_state, &
            target_init,      &
            target_end,       &
            target_output,    &
            target_tdcalc,    &
            target_inh,       &
            target_mode,      &
            target_type,      &
            j1_functional,    &
            calc_chi,         &
            target_move_ions, &
            parse_velocity_target

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
    oct_tg_velocity         = 10

  integer, public, parameter ::       &
    oct_targetmode_static = 0,        &
    oct_targetmode_td     = 1

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
    logical :: move_ions
    integer :: hhg_nks
    integer, pointer :: hhg_k(:) => null()
    FLOAT,   pointer :: hhg_alpha(:) => null()
    FLOAT,   pointer :: hhg_a(:) => null()
    FLOAT   :: hhg_w0
    FLOAT   :: dt
  end type target_t


  contains


  ! ----------------------------------------------------------------------
  ! This just copies the states_t variable present in target, into st.
  ! ----------------------------------------------------------------------
  subroutine target_get_state(target, st)
    type(target_t), intent(in)    :: target
    type(states_t), intent(inout) :: st

    call push_sub('target.target_get_state')
    call states_copy(st, target%st)

    call pop_sub('target.target_get_state')
  end subroutine target_get_state
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  ! The target is initialized, mainly by reading from the inp file.
  ! ----------------------------------------------------------------------
  subroutine target_init(gr, geo, stin, td, w0, target)
    type(grid_t),     intent(in)    :: gr
    type(geometry_t), intent(in)    :: geo
    type(states_t),   intent(in)    :: stin
    type(td_t),       intent(in)    :: td
    FLOAT,            intent(in)    :: w0
    type(target_t),   intent(inout) :: target

    integer             :: ierr, ip, ist, jst, jj
    type(block_t)       :: blk
    FLOAT               :: xx(MAX_DIM), rr, psi_re, psi_im
    CMPLX, allocatable  :: rotation_matrix(:, :)
    type(states_t)      :: tmp_st
    character(len=1024) :: expression

    call push_sub('target.target_init')

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
    !% WARNING: not implemented; reserved for use in future releases.
    !%Option oct_tg_density 5
    !% The target operator is a given density, <i>i.e.</i> the final state should have a density
    !% as close as possible as the one given in the input file, either from the variable
    !% <tt>OCTTargetDensityFromState</tt>, or from <tt>OCTTargetDensity</tt>.
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
    !%End
    call parse_integer(datasets_check('OCTTargetOperator'), oct_tg_gstransformation, target%type)
    if(.not.varinfo_valid_option('OCTTargetOperator', target%type)) &
      call input_error('OCTTargetOperator')

    call states_copy(target%st, stin)
    call states_deallocate_wfns(target%st)
    call states_allocate_wfns(target%st, gr%mesh, M_CMPLX)

    nullify(target%td_fitness)

    select case(target%type)
    case(oct_tg_groundstate)
      message(1) =  'Info: Using Ground State for TargetOperator'
      call write_info(1)
      call restart_read(trim(restart_dir)//GS_DIR, target%st, gr, geo, ierr, exact = .true.)
      
    case(oct_tg_excited) 

      message(1) =  'Info: TargetOperator is a linear combination of Slater determinants.'
      call write_info(1)

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
      call states_allocate_wfns(target%st, gr%mesh, M_CMPLX)
      target%st%node(:)  = 0

      call restart_read(trim(restart_dir)//GS_DIR, target%st, gr, geo, ierr, exact = .true.)
      call excited_states_init(target%est, target%st, "oct-excited-state-target") 

    case(oct_tg_exclude_state)

      message(1) =  'Info: The target functional is the exclusion of a number of states defined by'
      message(2) =  '      "OCTExcludedStates".'
      call write_info(2)
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
      call write_info(1)

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
          call states_end(tmp_st)
          call parse_block_end(blk)
          call states_calc_dens(target%st, gr)
        else
          message(1) = '"OCTTargetTransformStates" has to be specified as block.'
          call write_info(1)
          call input_error('OCTTargetTransformStates')
        end if
      else
        message(1) = 'Error: if "OCTTargetOperator = oct_tg_superposition", then you must'
        message(2) = 'supply one "OCTTargetTransformStates" block to create the superposition.'
        call write_info(2)
        call input_error('OCTTargetTransformStates')
      end if

    case(oct_tg_userdefined) 
      message(1) =  'Error: Option oct_tg_userdefined is disabled in this version.'
      call write_fatal(1)

    case(oct_tg_density) 

      message(1) =  'Info: Target is a density.'
      call write_info(1)

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
            call states_calc_dens(target%st, gr)
            do ip = 1, gr%mesh%np
              target%rho(ip) = sum(target%st%rho(ip, 1:target%st%d%spin_channels))
            end do
            call states_end(tmp_st)
            call parse_block_end(blk)
          else
            message(1) = '"OCTTargetDensityState" has to be specified as block.'
            call write_info(1)
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
        call write_fatal(2)
      end if

    case(oct_tg_local)
      !%Variable OCTLocalTarget
      !%Type string
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If <tt>OCTTargetOperator = oct_tg_local</tt>, then one must supply the target density
      !% that should be searched for. One can do this by supplying a string through
      !% the variable <tt>OCTLocalTarget</tt>.
      !%End

      if(parse_isdef('OCTTargetLocal') .ne. 0) then
        SAFE_ALLOCATE(target%rho(1:gr%mesh%np))
        target%rho = M_ZERO
        call parse_string('OCTTargetLocal', "0", expression)
        call conv_to_C_string(expression)
        do ip = 1, gr%mesh%np
          call mesh_r(gr%mesh, ip, rr, coords = xx)
          ! parse user-defined expression
          call parse_expression(psi_re, psi_im, gr%sb%dim, xx, rr, M_ZERO, expression)
          target%rho(ip) = psi_re
        end do
      else
        message(1) = 'If OCTTargetOperator = oct_tg_local, then you must give the shape'
        message(2) = 'of this target in variable "OCTTargetLocal".'
        call write_fatal(2)
      end if

    case(oct_tg_td_local)
      if(parse_block(datasets_check('OCTTdTarget'),blk)==0) then
        call parse_block_string(blk, 0, 0, target%td_local_target)
        call conv_to_C_string(target%td_local_target)
        SAFE_ALLOCATE(target%rho(1:gr%mesh%np))
      else
        message(1) = 'If OCTTargetMode = oct_targetmode_td, you must suppy a OCTTDTarget block.'
        call write_fatal(1)
      end if
      target%dt = td%dt
      SAFE_ALLOCATE(target%td_fitness(0:td%max_iter))
      target%td_fitness = M_ZERO


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
          call write_info(1)
          call input_error('OCTOptimizeHarmonicSpectrum')
        end if
      else
        write(message(1), '(a)') 'If "OCTTargetMode = oct_targetmode_hhg", you must supply an'
        write(message(2), '(a)') '"OCTOptimizeHarmonicSpectrum" block.'
        call write_fatal(2)
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
          call write_fatal(2)
       end if
       
       if(parse_isdef('OCTMoveIons') .eq. 0) then
          message(1) = 'If OCTTargetOperator = oct_tg_velocity, then you must supply'
          message(2) = 'the variable "OCTMoveIons".'
          call write_fatal(2)
       else
          call parse_logical('OCTMoveIons', .false., target%move_ions)
       end if

    case default
      write(message(1),'(a)') "Target Operator not properly defined."
      call write_fatal(1)
    end select

    call pop_sub('target.target_init')
  end subroutine target_init
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  subroutine target_end(target)
    type(target_t), intent(inout) :: target

    call push_sub('target.target_end')

    call states_end(target%st)
    if(target%type .eq. oct_tg_local .or. &
       target%type .eq. oct_tg_density .or. &
       target%type .eq. oct_tg_td_local) then
      SAFE_DEALLOCATE_P(target%rho)
      nullify(target%rho)
    end if
    if(target%type .eq. oct_tg_hhg) then
      SAFE_DEALLOCATE_P(target%hhg_k); nullify(target%hhg_k)
      SAFE_DEALLOCATE_P(target%hhg_alpha); nullify(target%hhg_alpha)
      SAFE_DEALLOCATE_P(target%hhg_a); nullify(target%hhg_a)
    end if
    if(target_mode(target).eq.oct_targetmode_td) then
      SAFE_DEALLOCATE_P(target%td_fitness); nullify(target%td_fitness)
    end if

    call pop_sub('target.target_end')
  end subroutine target_end
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  subroutine target_output(target, gr, dir, geo, outp)
    type(target_t), intent(inout) :: target
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(geometry_t),       intent(in)  :: geo
    type(h_sys_output_t),   intent(in)  :: outp

    integer :: ierr
    call push_sub('target.target_output')
    
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
      call h_sys_output_states(target%est%st, gr, geo, trim(dir)//'/st', outp)
      call excited_states_output(target%est, trim(dir))
    case default
      call h_sys_output_states(target%st, gr, geo, trim(dir), outp)
    end select

    call pop_sub('target.target_output')
  end subroutine target_output
  ! ----------------------------------------------------------------------


  ! ---------------------------------------------------------
  ! Calculates, at a given point in time marked by the integer
  ! index, the integrand of the target functional:
  ! <Psi(t)|\hat{O}(t)|Psi(t)>.
  ! ---------------------------------------------------------
  subroutine target_tdcalc(target, gr, psi, time)
    type(target_t), intent(inout) :: target
    type(grid_t),      intent(in) :: gr
    type(states_t),    intent(in) :: psi
    integer,           intent(in) :: time

    CMPLX, allocatable :: opsi(:, :)
    FLOAT, allocatable :: multipole(:, :)
    integer :: ist, ip, is

    if(target_mode(target)  .ne. oct_targetmode_td) return

    call push_sub('target.target_tdcalc')

    target%td_fitness(time) = M_ZERO

    select case(target%type)
    case(oct_tg_td_local)

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
            target%td_fitness(time) + zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, ist, 1), opsi(:, :))
        end do
        SAFE_DEALLOCATE_A(opsi)
      case(SPIN_POLARIZED); stop 'Error'
      case(SPINORS);        stop 'Error'
      end select

    case(oct_tg_hhg)

      SAFE_ALLOCATE(multipole(1:4, 1:psi%d%nspin))
      do is = 1, psi%d%nspin
        call dmf_multipoles(gr%mesh, psi%rho(:, is), 1, multipole(:, is))
      end do
      target%td_fitness(time) = sum(multipole(2, 1:psi%d%spin_channels))
      SAFE_DEALLOCATE_A(multipole)

    case default
      stop 'Error at target.target_tdcalc.'
    end select

    call pop_sub('target.target_tdcalc')
  end subroutine target_tdcalc
  ! ----------------------------------------------------------------------


  ! ---------------------------------------------------------------
  ! Calculates the inhomogeneous term that appears in the equation
  ! for chi, and places it into inh.
  ! ---------------------------------------------------------------
  subroutine target_inh(psi, gr, target, time, inh)
    type(states_t),    intent(in)        :: psi
    type(grid_t),      intent(in)        :: gr
    type(target_t),    intent(inout)     :: target
    FLOAT,             intent(in)        :: time
    type(states_t),    intent(inout)     :: inh
 
    integer :: ik, ist, idim, ip
    
    call push_sub('target.target_inh')

    select case(target%type)
    case(oct_tg_td_local)
      call target_build_tdlocal(target, gr, time)
      forall(ik = 1:inh%d%nik, ist = inh%st_start:inh%st_end, idim = 1:inh%d%dim, ip = 1:gr%mesh%np)
        inh%zpsi(ip, idim, ist, ik) = -target%rho(ip) * psi%zpsi(ip, idim, ist, ik)
      end forall
      
    end select

    call pop_sub('target.target_inh')
  end subroutine target_inh
  !----------------------------------------------------------


  !----------------------------------------------------------
  subroutine target_build_tdlocal(target, gr, time)
    type(target_t), intent(inout) :: target
    type(grid_t),   intent(in)    :: gr
    FLOAT,          intent(in)    :: time

    integer :: ip
    FLOAT :: xx(MAX_DIM), rr, re, im

    call push_sub('target.target_build_tdlocal')

    do ip = 1, gr%mesh%np
      call mesh_r(gr%mesh, ip, rr, coords = xx)
      call parse_expression(re, im, gr%sb%dim, xx, rr, time, target%td_local_target)
      target%rho(ip) = re
    end do

    call pop_sub('target.target_build_tdlocal')
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
    FLOAT :: omega, aa, maxhh, ww
    FLOAT, allocatable :: local_function(:)
    CMPLX, allocatable :: ddipole(:)
    CMPLX, allocatable :: opsi(:, :)
    FLOAT :: f_re, dummy(3)
    character(len=4096) :: inp_string

    call push_sub('target.j1_functional')

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
          j1 = j1 + zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, ist, 1), opsi(:, :))
        end do
        SAFE_DEALLOCATE_A(opsi)
      case(SPIN_POLARIZED); stop 'Error'
      case(SPINORS);        stop 'Error'
      end select

    case(oct_tg_td_local)
      j1 = sum(target%td_fitness) * target%dt

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

    case default
      do ik = 1, psi%d%nik
        do ist = psi%st_start, psi%st_end
          j1 = j1 + psi%occ(ist, ik) * &
            abs(zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, ist, ik), &
              target%st%zpsi(:, :, ist, ik)))**2
        end do
      end do

    end select

    call pop_sub('target.j1_functional')
  end function j1_functional
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! calculate |chi(T)> = \hat{O}(T) |psi(T)>
  ! ---------------------------------------------------------
  subroutine calc_chi(target, gr, psi_in, chi_out)
    type(target_t),    intent(inout) :: target
    type(grid_t),      intent(inout) :: gr
    type(states_t),    intent(inout) :: psi_in
    type(states_t),    intent(inout) :: chi_out
    
    CMPLX :: olap, zdet
    CMPLX, allocatable :: cI(:), dI(:), mat(:, :, :), mm(:, :, :, :), mk(:, :), lambda(:, :)
    integer :: ik, ip, ist, jst, idim, no_electrons, ia, ib, n_pairs, nst, kpoints, jj

    call push_sub('target.calc_chi')

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
            do ip = 1, gr%mesh%np_part
              if(psi_in%rho(ip, 1) > CNST(1.0e-8)) then
                chi_out%zpsi(ip, 1, ist, 1) = sqrt(target%rho(ip) / psi_in%rho(ip, 1)) * &
                  psi_in%zpsi(ip, 1, ist, 1)
              else
                chi_out%zpsi(ip, 1, ist, 1) = M_ZERO !sqrt(target%rho(ip))
              end if
            end do
          end do
        end if

      case(SPIN_POLARIZED); stop 'Error'
      case(SPINORS);        stop 'Error'
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
      case(SPIN_POLARIZED); stop 'Error'
      case(SPINORS);        stop 'Error'
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
        write(message(1), '(a)') 'Internal error in target.calc_chi.'
        call write_fatal(2)

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

    case default

      !olap = zstates_mpdotp(gr%mesh, target%st, psi_in)
      do ik = 1, psi_in%d%nik
        do ist = psi_in%st_start, psi_in%st_end
          olap = zmf_dotp(gr%mesh, target%st%zpsi(:, 1, ist, ik), psi_in%zpsi(:, 1, ist, ik))
          chi_out%zpsi(:, :, ist, ik) = olap * target%st%zpsi(:, :, ist, ik)
        end do
      end do

    end select

    call pop_sub('target.calc_chi')
  end subroutine calc_chi


  ! ----------------------------------------------------------------------
  integer pure function target_mode(target)
    type(target_t), intent(in) :: target

    select case(target%type)
    case(oct_tg_td_local, oct_tg_hhg)
      target_mode = oct_targetmode_td
    case(oct_tg_velocity)
      target_mode = oct_tg_velocity
    case default
      target_mode = oct_targetmode_static
    end select

  end function target_mode


  ! ----------------------------------------------------------------------
  integer pure function target_type(target)
    type(target_t), intent(in) :: target

    target_type = target%type

  end function target_type
  ! ----------------------------------------------------------------------


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
    CHARACTER (LEN=100)  :: v_string
    
    string_length = len(inp_string)
    do i=1, string_length 
       if(inp_string(i:i+1) == "v[") then
          m = 0
          if(inp_string(i+3:i+3) == ",") m = 1
          if(inp_string(i+4:i+4) == ",") m = 2
          if(m == 0) then
             message(1) = "OCTVelocityTarget Input error!"
             message(2) = "Atom number is either larger than 99 or not defined."
             call write_fatal(2)
          end if
          read(inp_string(i+2:i+1+m),*) n_atom
          read(inp_string(i+3+m:i+3+m),*) coord
          if(coord < 1 .or. coord > 3) then
             message(1) = "OCTVelocityTarget Input error!"
             message(2) = "Vector component is either larger than 3 or smaller than 1."
             call write_fatal(2)
          end if
          write(v_string,*) geo%atom(n_atom)%v(coord)
          inp_string = inp_string(:i-1) // "(" // trim(v_string) // ")" // inp_string(i+5+m:)
       end if
    end do
    
  end subroutine parse_velocity_target
  ! ----------------------------------------------------------------------


end module opt_control_target_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

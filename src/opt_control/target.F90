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
  use varinfo_m
  use messages_m
  use loct_parser_m
  use loct_m
  use lalg_basic_m
  use lalg_adv_m
  use io_m
  use global_m
  use string_m
  use states_m
  use states_dim_m
  use states_calc_m
  use excited_states_m
  use grid_m
  use h_sys_output_m
  use io_function_m
  use geometry_m
  use mesh_m
  use mesh_function_m
  use profiling_m
  use restart_m
  use td_m
  use spectrum_m
  use math_m

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
            calc_chi

  integer, public, parameter ::       &
    oct_tg_groundstate      = 1,      &
    oct_tg_excited          = 2,      &
    oct_tg_gstransformation = 3,      &
    oct_tg_userdefined      = 4,      &
    oct_tg_density          = 5,      &        
    oct_tg_local            = 6,      &
    oct_tg_td_local         = 7,      &
    oct_tg_exclude_state    = 8,      &
    oct_tg_hhg              = 9

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
    integer :: excluded_states
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
    call states_copy(st, target%st)
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

    integer           :: ierr, ip, ist, jst, j
    type(block_t)         :: blk
    FLOAT             :: x(MAX_DIM), r, psi_re, psi_im
    CMPLX, allocatable :: rotation_matrix(:, :)
    type(states_t)    :: tmp_st
    character(len=1024) :: expression

    call push_sub('target.target_init')

    !%Variable OCTTargetOperator
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default 3
    !%Description
    !% The variable OCTTargetOperator prescribes which kind of target functional is
    !% to be used.
    !%
    !% The possible arguments are:
    !%
    !%Option oct_tg_groundstate 1 
    !% The target operator is a projection operator on the ground state, i.e. the
    !% objective is to populate the ground state as much as possible.
    !%Option oct_tg_excited 2
    !% The target operator is an "excited state". This means that the target operator
    !% is a linear combination of Slater determinants, formed each one by swapping,
    !% from the ground state Slater determinant, one occupied state by one excited
    !% state (i.e. "single excitations"). The description of which excitations are
    !% used, and with which weight, should be given in a file called
    !% "oct-excited-state-target". This is still in very preliminary, experimental
    !% phase. See the documentation of subroutine "excited_states_init" in the source
    !% code in order to use this feature.
    !%Option oct_tg_gstransformation 3
    !% The target operator is a projection operator on a transformation of the ground state 
    !% orbitals defined by the block "OCTTargetTransformStates".
    !%Option oct_tg_userdefined 4
    !% WARNING: not implemented; reserved for use in future releases.
    !%Option oct_tg_density 5
    !% The target operator is a given density, i.e. the final state should have a density
    !% as close as possible as the one given in the input file, either from the variable
    !% "OCTTargetDensityFromState", or from "OCTTargetDensity".
    !%Option oct_tg_local 6
    !% The target operator is a local operator.
    !%Option oct_tg_td_local 7
    !% The target operator is a time-dependent local operator.
    !%Option oct_tg_exclude_state 8
    !% Target operator is the projection onto the complement of a given state, given by the
    !% block OCTTargetTransformStates. This means that the target operator is the unity
    !% operator minus the projector onto that state.
    !%Option oct_tg_hhg 9
    !% The target is the optimization of the HHG yield.
    !%End
    call loct_parse_int(datasets_check('OCTTargetOperator'), oct_tg_gstransformation, target%type)
    if(.not.varinfo_valid_option('OCTTargetOperator', target%type)) &
      call input_error('OCTTargetOperator')

    call states_copy(target%st, stin)
    call states_allocate_wfns(target%st, gr%mesh, M_CMPLX)

    nullify(target%td_fitness)

    select case(target%type)
    case(oct_tg_groundstate)
      message(1) =  'Info: Using Ground State for TargetOperator'
      call write_info(1)
      call restart_read(trim(restart_dir)//'gs', target%st, gr, geo, ierr)
      if(ierr.ne.0) then
        write(message(1),'(a)') 'Could not read ground-state wavefunctions from '//trim(restart_dir)//'gs.'
        call write_fatal(1)
      end if
      
    case(oct_tg_excited) 

      message(1) =  'Info: TargetOperator is a linear combination of Slater determinants.'
      call write_info(1)

      call states_look (trim(restart_dir)//'gs', gr%mesh%mpi_grp, ip, ip, target%st%nst, ierr)
      target%st%st_start = 1
      target%st%st_end   = target%st%nst

      SAFE_DEALLOCATE_P(target%st%occ)
      SAFE_DEALLOCATE_P(target%st%eigenval)
      SAFE_DEALLOCATE_P(target%st%momentum)
      SAFE_DEALLOCATE_P(target%st%node)

      ALLOCATE(target%st%occ(target%st%nst, target%st%d%nik), target%st%nst*target%st%d%nik)
      ALLOCATE(target%st%eigenval(target%st%nst, target%st%d%nik), target%st%nst*target%st%d%nik)
      ALLOCATE(target%st%momentum(3,target%st%nst, target%st%d%nik), 3*target%st%nst*target%st%d%nik)
      ALLOCATE(target%st%node(target%st%nst), target%st%nst)
      if(target%st%d%ispin == SPINORS) then
        SAFE_DEALLOCATE_P(target%st%spin)
        ALLOCATE(target%st%spin(3, target%st%nst, target%st%d%nik), target%st%nst*target%st%d%nik*3)
      end if
      call states_allocate_wfns(target%st, gr%mesh, M_CMPLX)
      target%st%node(:)  = 0

      call restart_read(trim(restart_dir)//'gs', target%st, gr, geo, ierr)
      if(ierr.ne.0) then
        write(message(1),'(a)') 'Could not read ground-state wavefunctions from '//trim(restart_dir)//'gs.'
        call write_fatal(1)
      end if
      call excited_states_init(target%est, target%st, "oct-excited-state-target") 

    case(oct_tg_exclude_state)

      message(1) =  'Info: The target functional is the exclusion of a number of states defined by'
      message(2) =  '      "OCTExcludeStates".'
      call write_info(2)
      !%Variable OCTExcludeStates
      !%Type integer
      !%Default 1
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% WARNING: Experimental
      !%End
      call loct_parse_int(datasets_check('OCTExcludeStates'), 1, target%excluded_states)
      call restart_look_and_read(target%st, gr, geo)

    case(oct_tg_gstransformation)  

      message(1) =  'Info: Using Superposition of States for TargetOperator'
      call write_info(1)

      !%Variable OCTTargetTransformStates
      !%Type block
      !%Default no
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If OCTTargetOperator = oct_tg_gstransformation, you must specify one
      !% OCTTargetTransformStates block, in order to specify which linear
      !% combination of the states present in "restart/gs" is used to
      !% create the target state.
      !% 
      !% The syntax is equivalent to the one used for the TransformStates
      !% block.
      !%End
      if(loct_parse_isdef(datasets_check('OCTTargetTransformStates')).ne.0) then
        if(loct_parse_block(datasets_check('OCTTargetTransformStates'), blk) == 0) then
          call states_copy(tmp_st, target%st)
          SAFE_DEALLOCATE_P(tmp_st%zpsi)
          call restart_look_and_read(tmp_st, gr, geo)
          ALLOCATE(rotation_matrix(target%st%nst, tmp_st%nst), target%st%nst*tmp_st%nst)
          rotation_matrix = M_z0
          do ist = 1, target%st%nst
            do jst = 1, loct_parse_block_cols(blk, ist-1)
              call loct_parse_block_cmplx(blk, ist-1, jst-1, rotation_matrix(ist, jst))
            end do
          end do
          call rotate_states(gr%mesh, target%st, tmp_st, rotation_matrix)
          SAFE_DEALLOCATE_A(rotation_matrix)
          call states_end(tmp_st)
          call loct_parse_block_end(blk)
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
      message(1) =  'Error: Option oct_tg_userdefined is disabled in this version'
      call write_fatal(1)

    case(oct_tg_density) 

      message(1) =  'Info: Target is a density.'
      call write_info(1)

      !%Variable OCTTargetDensity
      !%Type string
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If OCTTargetOperator = oct_tg_density, then one must supply the target density
      !% that should be searched for. This one can do by supplying a string through
      !% the variable OCTLocalTarget.
      !%End


      !%Variable OCTTargetDensityFromState
      !%Type block
      !%Default no
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If OCTTargetOperator = oct_tg_density, and OCTLocalTarget = "OCTTargetDensityFromState",
      !% you must specify one OCTTargetDensityState block, in order to specify which linear
      !% combination of the states present in "restart/gs" is used to
      !% create the target density.
      !%
      !% The syntax is equivalent to the one used for the TransformStates
      !% block.
      !%End

      if(loct_parse_isdef('OCTTargetDensity').ne.0) then
        ALLOCATE(target%rho(gr%mesh%np), gr%mesh%np)
        target%rho = M_ZERO
        call loct_parse_string('OCTTargetDensity', "0", expression)


        if(trim(expression).eq.'OCTTargetDensityFromState') then

          if(loct_parse_block(datasets_check('OCTTargetDensityFromState'), blk) == 0) then
            tmp_st = target%st
            SAFE_DEALLOCATE_P(tmp_st%zpsi)
            call restart_look_and_read(tmp_st, gr, geo)
            ALLOCATE(rotation_matrix(target%st%nst, tmp_st%nst), target%st%nst*tmp_st%nst)
            rotation_matrix = M_z0
            do ist = 1, target%st%nst
              do jst = 1, loct_parse_block_cols(blk, ist-1)
                call loct_parse_block_cmplx(blk, ist-1, jst-1, rotation_matrix(ist, jst))
              end do
            end do
            call rotate_states(gr%mesh, target%st, tmp_st, rotation_matrix)
            SAFE_DEALLOCATE_A(rotation_matrix)
            call states_calc_dens(target%st, gr%mesh%np_part)
            do ip = 1, gr%mesh%np
              target%rho(ip) = sum(target%st%rho(ip, 1:target%st%d%spin_channels))
            end do
            call states_end(tmp_st)
            call loct_parse_block_end(blk)
          else
            message(1) = '"OCTTargetDensityState" has to be specified as block.'
            call write_info(1)
            call input_error('OCTTargetDensity')
          end if

        else

          call conv_to_C_string(expression)
          do ip = 1, gr%mesh%np
            call mesh_r(gr%mesh, ip, r, x = x)
            ! parse user defined expression
            call loct_parse_expression(psi_re, psi_im, gr%sb%dim, x, r, M_ZERO, expression)
            target%rho(ip) = psi_re
          end do
          ! Normalize
          r = dmf_integrate(gr%mesh, target%rho)
          target%rho = (-target%st%val_charge) * target%rho/r
        end if

      else
        message(1) = 'If OCTTargetOperator = oct_tg_density, then you must give the shape'
        message(2) = 'of this target in variable "OCTTargetDensity"'
        call write_fatal(2)
      end if

    case(oct_tg_local)
      !%Variable OCTTargetLocal
      !%Type string
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If OCTTargetOperator = oct_tg_local, then one must supply the target density
      !% that should be searched for. This one can do by supplying a string through
      !% the variable OCTLocalTarget.
      !%End

      if(loct_parse_isdef('OCTTargetLocal').ne.0) then
        ALLOCATE(target%rho(gr%mesh%np), gr%mesh%np)
        target%rho = M_ZERO
        call loct_parse_string('OCTTargetLocal', "0", expression)
        call conv_to_C_string(expression)
        do ip = 1, gr%mesh%np
          call mesh_r(gr%mesh, ip, r, x = x)
          ! parse user defined expression
          call loct_parse_expression(psi_re, psi_im, gr%sb%dim, x, r, M_ZERO, expression)
          target%rho(ip) = psi_re
        end do
      else
        message(1) = 'If OCTTargetOperator = oct_tg_local, then you must give the shape'
        message(2) = 'of this target in variable "OCTTargetLocal"'
        call write_fatal(2)
      end if

    case(oct_tg_td_local)
      if(loct_parse_block(datasets_check('OCTTdTarget'),blk)==0) then
        call loct_parse_block_string(blk, 0, 0, target%td_local_target)
        call conv_to_C_string(target%td_local_target)
        ALLOCATE(target%rho(gr%mesh%np), gr%mesh%np)
      else
        message(1) = 'If OCTTargetMode = oct_targetmode_td, you must suppy a OCTTDTarget block'
        call write_fatal(1)
      end if
      target%dt     = td%dt
      ALLOCATE(target%td_fitness(0:td%max_iter), td%max_iter+1)


    case(oct_tg_hhg)
      !%Variable OCTOptimizeHarmonicSpectrum
      !%Type block
      !%Default no
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% WARNING: Experimental
      !%
      !% If "OCTTargetOperator = oct_tg_hhg", the target is the harmonic emission spectrum.
      !% In that case, you must supply a "OCTOptimizeHarmonicSpectrum" block in the "inp"
      !% file. The target is given, in general, by:
      !%
      !% <math>J_1 = \int_0^\infty d\omega \alpha(\omega) H(\omega)</math>,
      !%
      !% where <math>H(\omega)</math> is the harmonic spectrum generated by the system, and
      !% <math>\alpha(\omega)</math> is some function that determines what exactly we want
      !% to optimize. The role of the "OCTOptimizeHarmonicSpectrum" block is to determine
      !% this <math>\alpha(\omega)</math> function. Currently, this function is defined as:
      !%
      !% <math>\alpha(\omega) = \sum_{L=1}^{M} \frac{\alpha_L}{a_L} \sqcap( (\omega - L\omega_0)/a_L )</math>,
      !%
      !% where <math>omega_0</math> is the carrier frequency. <math>M</math> is
      !% the number of columns in the "OCTOptimizeHarmonicSpectrum". The values of $L$ will be listed
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
      if(loct_parse_isdef(datasets_check('OCTOptimizeHarmonicSpectrum')).ne.0) then
        if(loct_parse_block(datasets_check('OCTOptimizeHarmonicSpectrum'), blk) == 0) then
          target%hhg_nks = loct_parse_block_cols(blk, 0)
          ALLOCATE(target%hhg_k(target%hhg_nks), target%hhg_nks)
          ALLOCATE(target%hhg_alpha(target%hhg_nks), target%hhg_nks)
          ALLOCATE(target%hhg_a(target%hhg_nks), target%hhg_nks)
          do j = 1, target%hhg_nks
            call loct_parse_block_int(blk, 0, j-1, target%hhg_k(j))
            call loct_parse_block_float(blk, 1, j-1, target%hhg_alpha(j))
            call loct_parse_block_float(blk, 2, j-1, target%hhg_a(j))
          end do
        else
          message(1) = '"OCTOptimizeHarmonicSpectrum" has to be specified as block.'
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
      ALLOCATE(target%td_fitness(0:td%max_iter), td%max_iter+1)
      target%td_fitness = M_ZERO

    case default
      write(message(1),'(a)') "Target Operator not properly defined."
      call write_fatal(1)
    end select

    call pop_sub()
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

    call pop_sub()
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
      call doutput_function(outp%how, trim(dir), 'local_target', gr%mesh, gr%sb, &
        target%rho, M_ONE, ierr, geo = geo)
    case(oct_tg_td_local)
      call tdtarget_build_tdlocal(target, gr, M_ZERO)
      call doutput_function(outp%how, trim(dir), 'td_local_target', gr%mesh, gr%sb, &
        target%rho, M_ONE, ierr, geo = geo)
    case(oct_tg_density)
      call doutput_function(outp%how, trim(dir), 'density_target', gr%mesh, gr%sb, &
        target%rho, M_ONE, ierr, geo = geo)
    case(oct_tg_excited)
      call h_sys_output_states(target%est%st, gr, geo, trim(dir)//'/st', outp)
      call excited_states_output(target%est, trim(dir))
    case default
      call h_sys_output_states(target%st, gr, geo, trim(dir), outp)
    end select

    call pop_sub()
  end subroutine target_output
  ! ----------------------------------------------------------------------


  ! ---------------------------------------------------------
  ! Calculates, at a given point in time marked by the integer
  ! index "i", the integrand of the target functional:
  ! <Psi(t)|\hat{O}(t)|Psi(t)>.
  ! This subroutine that the operator O is already updated
  ! ---------------------------------------------------------
  subroutine target_tdcalc(target, gr, psi, i)
    type(target_t), intent(inout) :: target
    type(grid_t),      intent(in) :: gr
    type(states_t),    intent(in) :: psi
    integer,           intent(in) :: i
    CMPLX, allocatable :: opsi(:, :)
    FLOAT, allocatable :: multipole(:, :)
    integer :: p, j, is

    if(target_mode(target) .ne. oct_targetmode_td) return

    call push_sub('target.target_tdcalc')

    target%td_fitness(i) = M_ZERO

    select case(target%type)
    case(oct_tg_td_local)

      select case(psi%d%ispin)
      case(UNPOLARIZED)
        ASSERT(psi%d%nik.eq.1)
        ALLOCATE(opsi(gr%mesh%np_part, 1), gr%mesh%np_part)
        opsi = M_z0
        do p  = psi%st_start, psi%st_end
          do j = 1, gr%mesh%np
            opsi(j, 1) = target%rho(j) * psi%zpsi(j, 1, p, 1)
          end do
          target%td_fitness(i) = &
            target%td_fitness(i) + zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, p, 1), opsi(:, :))
        end do
        SAFE_DEALLOCATE_A(opsi)
      case(SPIN_POLARIZED); stop 'Error'
      case(SPINORS);        stop 'Error'
      end select

    case(oct_tg_hhg)

      ALLOCATE(multipole(4, psi%d%nspin), 4*psi%d%nspin)
      do is = 1, psi%d%nspin
        call dmf_multipoles(gr%mesh, psi%rho(:, is), 1, multipole(:, is))
      end do
      target%td_fitness(i) = sum(multipole(2, 1:psi%d%spin_channels))
      SAFE_DEALLOCATE_A(multipole)

    case default
      stop 'Error at calc_tdfitness'
    end select

    call pop_sub()
  end subroutine target_tdcalc
  ! ----------------------------------------------------------------------


  ! ---------------------------------------------------------------
  ! Calculates the inhomogeneous term that appears in the equation
  ! for chi, and places it into inh.
  ! ---------------------------------------------------------------
  subroutine target_inh(psi, gr, target, t, inh)
    type(states_t),    intent(in)        :: psi
    type(grid_t),      intent(in)        :: gr
    type(target_t),    intent(inout)     :: target
    FLOAT,             intent(in)        :: t
    type(states_t),    intent(inout)     :: inh
 
    integer :: ik, ist, idim, i
    
    call push_sub('target.target_inh')

    select case(target%type)
    case(oct_tg_td_local)
      call tdtarget_build_tdlocal(target, gr, t)
      do ik = 1, inh%d%nik
        do ist = inh%st_start, inh%st_end
          do idim = 1, inh%d%dim
            do i = 1, gr%mesh%np
              inh%zpsi(i, idim, ist, ik) = M_zI * target%rho(i) * psi%zpsi(i, idim, ist, ik)
            end do
          end do
        end do
      end do
      
    end select

    call pop_sub()
  end subroutine target_inh
  !----------------------------------------------------------


  !----------------------------------------------------------
  subroutine tdtarget_build_tdlocal(target, gr, t)
    type(target_t), intent(inout) :: target
    type(grid_t),      intent(in) :: gr
    FLOAT, intent(in)             :: t
    integer :: i
    FLOAT :: xx(MAX_DIM), r, re, im
    call push_sub('target.target_build_tdlocal')

    do i = 1, gr%mesh%np
      call mesh_r(gr%mesh, i, r, x = xx)
      call loct_parse_expression(re, im, gr%sb%dim, xx, r, t, target%td_local_target)
      target%rho(i) = re
    end do

    call pop_sub()
  end subroutine tdtarget_build_tdlocal
  !----------------------------------------------------------


  ! ---------------------------------------------------------
  ! Calculates the J1 functional, i.e.:
  ! <Psi(T)|\hat{O}|Psi(T) in the time-independent
  ! case, or else \int_0^T dt <Psi(t)|\hat{O}(t)|Psi(t) in 
  ! the time-dependent case.
  ! ---------------------------------------------------------
  FLOAT function j1_functional(target, gr, psi) result(j1)
    type(target_t), intent(inout)   :: target
    type(grid_t),   intent(inout)   :: gr
    type(states_t), intent(inout)   :: psi

    integer :: i, p, j, maxiter
    FLOAT :: omega, a, maxhh, w
    FLOAT, allocatable :: local_function(:)
    CMPLX, allocatable :: ddipole(:)
    CMPLX, allocatable :: opsi(:, :)

    call push_sub('target.j1_functional')

    j1 = M_ZERO
    select case(target%type)
    case(oct_tg_density)

      ALLOCATE(local_function(gr%mesh%np), gr%mesh%np)
      do i = 1, gr%mesh%np
        local_function(i) = - ( sqrt(psi%rho(i, 1)) - sqrt(target%rho(i)) )**2
      end do
      j1 = dmf_integrate(gr%mesh, local_function)
      SAFE_DEALLOCATE_A(local_function)

    case(oct_tg_local)

      select case(psi%d%ispin)
      case(UNPOLARIZED)
        ASSERT(psi%d%nik.eq.1)
        ALLOCATE(opsi(gr%mesh%np_part, 1), gr%mesh%np_part)
        opsi = M_z0
        j1 = M_ZERO
        do p  = psi%st_start, psi%st_end
          do j = 1, gr%mesh%np
            opsi(j, 1) = target%rho(j) * psi%zpsi(j, 1, p, 1)
          end do
          j1 = j1 + zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, p, 1), opsi(:, :))
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
      do i = 1, target%excluded_states
        j1 = j1 - abs(zmf_dotp(gr%mesh, psi%d%dim, &
          target%st%zpsi(:, :, i, 1), psi%zpsi(:, :, 1, 1)))**2
      end do

    case(oct_tg_hhg)

      maxiter = size(target%td_fitness) - 1
      ALLOCATE(ddipole(0:maxiter), maxiter+1)

      ddipole(0) = M_ZERO
      do i = 1, maxiter - 1
        ddipole(i) = (target%td_fitness(i-1)+target%td_fitness(i+1)- &
                      M_TWO*target%td_fitness(i))/target%dt**2
      end do
      call interpolate( target%dt*(/ -3, -2, -1 /),   &
                        ddipole(maxiter-3:maxiter-1), &
                        M_ZERO, &
                        ddipole(maxiter) )

      call spectrum_hsfunction_init(target%dt, 0, maxiter, maxiter, ddipole)
      do j = 1, target%hhg_nks
        a = target%hhg_a(j) * target%hhg_w0
        w = target%hhg_k(j)*target%hhg_w0
        call spectrum_hsfunction_min(w -a, w +a, w, omega, maxhh)
        j1 = j1 + target%hhg_alpha(j) * log(-maxhh)
      end do
      call spectrum_hsfunction_end()

      SAFE_DEALLOCATE_A(ddipole)

    case default
      j1 = abs(zstates_mpdotp(gr%mesh, psi, target%st))**2

    end select

    call pop_sub()
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
    
    CMPLX   :: olap, zdet
    CMPLX, allocatable :: cI(:), dI(:), mat(:, :, :), mm(:, :, :, :), mk(:, :), lambda(:, :)
    integer :: ik, p, dim, k, j, no_electrons, ia, ib, n_pairs, nst, kpoints, jj

    call push_sub('target.calc_chi')

    no_electrons = -nint(psi_in%val_charge)

    select case(target%type)

    case(oct_tg_density)

      select case(psi_in%d%ispin)
      case(UNPOLARIZED)

        ASSERT(psi_in%d%nik.eq.1)

        if(no_electrons .eq. 1) then
          do j = 1, gr%mesh%np
            chi_out%zpsi(j, 1, 1, 1) = sqrt(target%rho(j)) * &
              exp( M_z1 * atan2(aimag(psi_in%zpsi(j, 1, 1, 1)), &
                                real(psi_in%zpsi(j, 1, 1, 1)  )) )
          end do
        else
          do p  = psi_in%st_start, psi_in%st_end
            do j = 1, gr%mesh%np_part
              if(psi_in%rho(j, 1) > CNST(1.0e-8)) then
                chi_out%zpsi(j, 1, p, 1) = sqrt(target%rho(j)/psi_in%rho(j, 1)) * &
                  psi_in%zpsi(j, 1, p, 1)
              else
                chi_out%zpsi(j, 1, p, 1) = M_ZERO!sqrt(target%rho(j))
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
        ASSERT(psi_in%d%nik.eq.1)
        do p  = psi_in%st_start, psi_in%st_end
          do j = 1, gr%mesh%np
            chi_out%zpsi(j, 1, p, 1) = target%rho(j) * psi_in%zpsi(j, 1, p, 1)
          end do
        end do
      case(SPIN_POLARIZED); stop 'Error'
      case(SPINORS);        stop 'Error'
      end select

    case(oct_tg_td_local)
      !We assume that there is no time-independent operator.
      forall(ik = 1:chi_out%d%nik, p = chi_out%st_start:chi_out%st_end, &
             dim = 1:chi_out%d%dim, j = 1:gr%mesh%np)
        chi_out%zpsi(j, dim, p, ik) = M_z0
      end forall

    case(oct_tg_excited) 

      n_pairs = target%est%n_pairs
      kpoints = psi_in%d%nik
      nst = psi_in%nst

      ALLOCATE(cI(n_pairs), n_pairs)
      ALLOCATE(dI(n_pairs), n_pairs)
      ALLOCATE(mat(target%est%st%nst, nst, psi_in%d%nik), target%est%st%nst*nst*psi_in%d%nik)
      ALLOCATE(mm(nst, nst, kpoints, n_pairs), nst*nst*kpoints*n_pairs)
      ALLOCATE(mk(gr%mesh%np_part, psi_in%d%dim), gr%mesh%np_part * psi_in%d%dim)
      ALLOCATE(lambda(n_pairs, n_pairs), n_pairs*n_pairs)

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
        write(message(1), '(a)') 'Internal error in aux.calc_chi'
        call write_fatal(2)

      case(SPIN_POLARIZED)
        ASSERT(chi_out%d%nik .eq. 2)

        do ik = 1, kpoints
          do k = chi_out%st_start, chi_out%st_end
            chi_out%zpsi(:, :, k, ik) = M_z0
            do ia = 1, n_pairs
              if(ik .ne. target%est%pair(ia)%sigma) cycle
              if(abs(dI(ia)) < CNST(1.0e-12)) cycle
              do ib = 1, n_pairs
                if(abs(dI(ib)) < CNST(1.0e-12)) cycle
                mk = M_z0
                do j = 1, nst
                  if(j .eq. target%est%pair(ib)%i) jj = target%est%pair(ia)%a
                  mk(:, :) = mk(:, :) + conjg(mm(k, j, ik, ib)) * target%est%st%zpsi(:, :, jj, ik)
                end do
                call lalg_axpy(gr%mesh%np_part, psi_in%d%dim, M_z1, lambda(ib, ia)*mk(:, :), chi_out%zpsi(:, :, k, ik))
              end do
            end do
          end do
        end do
        
      case(SPINORS)
        ASSERT(chi_out%d%nik .eq. 1)

        do k = chi_out%st_start, chi_out%st_end
          chi_out%zpsi(:, :, k, 1) = M_z0

          do ia = 1, n_pairs
            if(abs(dI(ia)) < CNST(1.0e-12)) cycle

            do ib = 1, n_pairs
              if(abs(dI(ib)) < CNST(1.0e-12)) cycle

              mk = M_z0
              do j = 1, nst
                if(j .eq. target%est%pair(ib)%i) jj = target%est%pair(ia)%a
                mk(:, :) = mk(:, :) + conjg(mm(k, j, 1, ib)) * target%est%st%zpsi(:, :, jj, 1)
              end do

              call lalg_axpy(gr%mesh%np_part, 2, M_z1, lambda(ib, ia)*mk(:, :), chi_out%zpsi(:, :, k, 1))
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
      do p = 1, target%excluded_states
        olap = zmf_dotp(gr%mesh, psi_in%d%dim, target%st%zpsi(:, :, p, 1), psi_in%zpsi(:, :, 1, 1))
        chi_out%zpsi(:, :, 1, 1) = chi_out%zpsi(:, :, 1, 1) - olap*target%st%zpsi(:, :, p, 1)
      end do

    case default

      olap = zstates_mpdotp(gr%mesh, target%st, psi_in)
      do ik = 1, psi_in%d%nik
        do p  = psi_in%st_start, psi_in%st_end
          chi_out%zpsi(:, :, p, ik) = olap*target%st%zpsi(:, :, p, ik)
        end do
      end do

    end select

    call pop_sub()
  end subroutine calc_chi
  ! ---------------------------------------------------------


  ! ----------------------------------------------------------------------
  integer pure function target_mode(target)
    type(target_t), intent(in) :: target
    select case(target%type)
    case(oct_tg_td_local, oct_tg_hhg)
      target_mode = oct_targetmode_td
    case default
      target_mode = oct_targetmode_static
    end select
  end function target_mode
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  integer pure function target_type(target)
    type(target_t), intent(in) :: target
    target_type = target%type
  end function target_type
  ! ----------------------------------------------------------------------


end module opt_control_target_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

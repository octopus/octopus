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
  use lib_oct_parser_m
  use lib_oct_m
  use io_m
  use global_m
  use string_m
  use states_m
  use excited_states_m
  use grid_m
  use output_m
  use geometry_m
  use states_output_m
  use mesh_m
  use mesh_function_m
  use restart_m
  use timedep_m
  use opt_control_constants_m

  implicit none

  private
  public :: target_t,       &
            target_init,    &
            target_end,     &
            target_output,  &
            calc_tdfitness, &
            calc_inh

  type target_t
    integer :: type
    integer :: mode
    type(states_t) :: st
    type(excited_states_t) :: est
    FLOAT, pointer :: rho(:)
    FLOAT, pointer :: td_fitness(:)
    character(len=200) :: td_local_target
    integer :: excluded_states
  end type target_t

  contains


  ! ----------------------------------------------------------------------
  subroutine target_init(gr, geo, stin, td, target)
    type(grid_t), intent(in)      :: gr
    type(geometry_t), intent(in)  :: geo
    type(states_t), intent(in)    :: stin
    type(td_t), intent(in)        :: td
    type(target_t), intent(inout) :: target

    integer           :: ierr, ip, ist, jst
    C_POINTER         :: blk
    FLOAT             :: x(MAX_DIM), r, psi_re, psi_im
    CMPLX, allocatable :: rotation_matrix(:, :)
    type(states_t)    :: tmp_st
    character(len=1024) :: expression

    call push_sub('target.target_init')

    !%Variable OCTTargetOperator
    !%Type integer
    !%Section Optimal Control
    !%Default 3
    !%Description
    !% The variable OCTTargetOperator prescribes which kind of target functional is
    !% to be used.
    !%
    !% The possible arguments are:
    !%
    !%Option oct_tg_groundstate 1 
    !% Targetoperator is a projection operator on the ground state
    !%Option oct_tg_excited 2
    !% The target operator
    !%Option oct_tg_gstransformation 3
    !% Targetoperator is a projection operator on a transformation of the ground state 
    !% orbitals defined by the block OCTTargetTransformStates
    !%Option oct_tg_userdefined 4
    !% Targetoperator is a projection operator on a user defined state
    !%Option oct_tg_density 5
    !% Targetoperator is a given density.
    !%Option oct_tg_local 6
    !% Target operator is a local operator.
    !%Option oct_tg_td_local 7
    !% Target operator is a time-dependent local operator
    !%Option oct_tg_exclude_state 8
    !% Target operator is the projection onto the complement of a given state, given by the
    !% block OCTTargetTransformStates. This means that the target operator is the unity
    !% operator minus the projector onto that state.
    !%End
    call loct_parse_int(check_inp('OCTTargetOperator'), oct_tg_gstransformation, target%type)
    if(.not.varinfo_valid_option('OCTTargetOperator', target%type)) call input_error('OCTTargetOperator')

    call states_copy(target%st, stin)

    select case(target%type)
    case(oct_tg_groundstate)
      message(1) =  'Info: Using Ground State for TargetOperator'
      call write_info(1)
      target%mode = oct_targetmode_static
      call restart_read(trim(tmpdir)//'gs', target%st, gr, geo, ierr)
      if(ierr.ne.0) then
        write(message(1),'(a)') 'Could not read ground-state wavefunctions from '//trim(tmpdir)//'gs.'
        call write_fatal(1)
      end if

      
    case(oct_tg_excited) 

      message(1) =  'Info: TargetOperator is a linear combination of Slater determinants.'
      call write_info(1)
      target%mode = oct_targetmode_static

      call states_look (trim(tmpdir)//'gs', gr%m, ip, ip, target%st%nst, ierr)
      target%st%st_start = 1
      target%st%st_end   = target%st%nst
      deallocate(target%st%occ, target%st%eigenval, target%st%momentum, target%st%node)
      ALLOCATE(target%st%occ(target%st%nst, target%st%d%nik), target%st%nst*target%st%d%nik)
      ALLOCATE(target%st%eigenval(target%st%nst, target%st%d%nik), target%st%nst*target%st%d%nik)
      ALLOCATE(target%st%momentum(3,target%st%nst, target%st%d%nik), 3*target%st%nst*target%st%d%nik)
      ALLOCATE(target%st%node(target%st%nst), target%st%nst)
      if(target%st%d%ispin == SPINORS) then
        deallocate(target%st%spin)
        ALLOCATE(target%st%spin(3, target%st%nst, target%st%d%nik), target%st%nst*target%st%d%nik*3)
      end if
      call states_allocate_wfns(target%st, gr%m, M_CMPLX)
      target%st%node(:)  = 0

      call restart_read(trim(tmpdir)//'gs', target%st, gr, geo, ierr)
      if(ierr.ne.0) then
        write(message(1),'(a)') 'Could not read ground-state wavefunctions from '//trim(tmpdir)//'gs.'
        call write_fatal(1)
      end if
      call excited_states_init(target%est, target%st, "oct-excited-state-target") 

    case(oct_tg_exclude_state)

      message(1) =  'Info: The target functional is the exclusion of a number of states defined by'
      message(2) =  '      "OCTExcludeStates".'
      call write_info(2)
      target%mode = oct_targetmode_static
      !%Variable OCTExcludeStates
      !%Type integer
      !%Default 1
      !%Section Optimal Control
      !%Description
      !% WARNING: Experimental
      !%End
      call loct_parse_int(check_inp('OCTExcludeStates'), 1, target%excluded_states)
      call restart_look_and_read("tmp", target%st, gr, geo, ierr)

    case(oct_tg_gstransformation)  

      message(1) =  'Info: Using Superposition of States for TargetOperator'
      call write_info(1)
      target%mode = oct_targetmode_static

      !%Variable OCTTargetTransformStates
      !%Type block
      !%Default no
      !%Section Optimal Control
      !%Description
      !% If OCTTargetOperator = oct_tg_gstransformation, you must specify one
      !% OCTTargetTransformStates block, in order to specify which linear
      !% combination of the states present in "restart/gs" is used to
      !% create the target state.
      !% 
      !% The syntax is equivalent to the one used for the TransformStates
      !% block.
      !%End
      if(loct_parse_isdef(check_inp('OCTTargetTransformStates')).ne.0) then
        if(loct_parse_block(check_inp('OCTTargetTransformStates'), blk) == 0) then
          tmp_st = target%st
          deallocate(tmp_st%zpsi)
          call restart_look_and_read("tmp", tmp_st, gr, geo, ierr)
          ALLOCATE(rotation_matrix(target%st%nst, tmp_st%nst), target%st%nst*tmp_st%nst)
          rotation_matrix = M_z0
          do ist = 1, target%st%nst
            do jst = 1, loct_parse_block_cols(blk, ist-1)
              call loct_parse_block_cmplx(blk, ist-1, jst-1, rotation_matrix(ist, jst))
            end do
          end do
          call rotate_states(gr%m, target%st, tmp_st, rotation_matrix)
          deallocate(rotation_matrix)
          call states_end(tmp_st)
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
      target%mode = oct_targetmode_static

      !%Variable OCTTargetDensity
      !%Type string
      !%Section Optimal Control
      !%Description
      !% If OCTTargetOperator = oct_tg_local, then one must supply the target density
      !% that should be searched for. This one can do by supplying a string through
      !% the variable OCTLocalTarget.
      !%End


      !%Variable OCTTargetDensityFromState
      !%Type block
      !%Default no
      !%Section Optimal Control
      !%Description
      !% If OCTTargetOperator = oct_tg_local, and OCTLocalTarget = "OCTTargetDensityFromState",
      !% you must specify one OCTTargetDensityState block, in order to specify which linear
      !% combination of the states present in "restart/gs" is used to
      !% create the target density.
      !%
      !% The syntax is equivalent to the one used for the TransformStates
      !% block.
      !%End

      if(loct_parse_isdef('OCTTargetDensity').ne.0) then
        ALLOCATE(target%rho(NP), NP)
        target%rho = M_ZERO
        call loct_parse_string('OCTTargetDensity', "0", expression)


        if(trim(expression).eq.'OCTTargetDensityFromState') then

          if(loct_parse_block(check_inp('OCTTargetDensityFromState'), blk) == 0) then
            tmp_st = target%st
            deallocate(tmp_st%zpsi)
            call restart_look_and_read("tmp", tmp_st, gr, geo, ierr)
            ALLOCATE(rotation_matrix(target%st%nst, tmp_st%nst), target%st%nst*tmp_st%nst)
            rotation_matrix = M_z0
            do ist = 1, target%st%nst
              do jst = 1, loct_parse_block_cols(blk, ist-1)
                call loct_parse_block_cmplx(blk, ist-1, jst-1, rotation_matrix(ist, jst))
              end do
            end do
            call rotate_states(gr%m, target%st, tmp_st, rotation_matrix)
            deallocate(rotation_matrix)
            call states_calc_dens(target%st, NP_PART, target%st%rho)
            do ip = 1, NP
              target%rho(ip) = sum(target%st%rho(ip, 1:target%st%d%spin_channels))
            end do
            call states_end(tmp_st)
          else
            message(1) = '"OCTTargetDensityState" has to be specified as block.'
            call write_info(1)
            call input_error('OCTTargetDensity')
          end if

        else

          call conv_to_C_string(expression)
          do ip = 1, NP
            call mesh_r(gr%m, ip, r, x = x)
            ! parse user defined expression
            call loct_parse_expression(psi_re, psi_im, &
              x(1), x(2), x(3), r, M_ZERO, expression)
            target%rho(ip) = psi_re
          end do
          ! Normalize
          r = dmf_integrate(gr%m, target%rho)
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
      !%Section Optimal Control
      !%Description
      !% If OCTTargetOperator = oct_tg_local, then one must supply the target density
      !% that should be searched for. This one can do by supplying a string through
      !% the variable OCTLocalTarget.
      !%End
      target%mode = oct_targetmode_static

      if(loct_parse_isdef('OCTTargetLocal').ne.0) then
        ALLOCATE(target%rho(NP), NP)
        target%rho = M_ZERO
        call loct_parse_string('OCTTargetLocal', "0", expression)
        call conv_to_C_string(expression)
        do ip = 1, NP
          call mesh_r(gr%m, ip, r, x = x)
          ! parse user defined expression
          call loct_parse_expression(psi_re, psi_im, &
            x(1), x(2), x(3), r, M_ZERO, expression)
          target%rho(ip) = psi_re
        end do
      else
        message(1) = 'If OCTTargetOperator = oct_tg_local, then you must give the shape'
        message(2) = 'of this target in variable "OCTTargetLocal"'
        call write_fatal(2)
      end if

    case(oct_tg_td_local)
      target%mode = oct_targetmode_td
      call tdtarget_init(target, gr, td)

    case default
      write(message(1),'(a)') "Target Operator not properly defined."
      call write_fatal(1)
    end select


    call pop_sub()
  end subroutine target_init


  ! ----------------------------------------------------------------------
  subroutine target_end(target)
    type(target_t), intent(inout) :: target
    integer :: i

    call push_sub('target.target_end')

    call states_end(target%st)
    if(target%type .eq. oct_tg_local .or. &
       target%type .eq. oct_tg_density .or. &
       target%type .eq. oct_tg_td_local) then
      deallocate(target%rho)
      nullify(target%rho)
    end if
    if(associated(target%td_fitness)) then
      deallocate(target%td_fitness); nullify(target%td_fitness)
    end if

    call pop_sub()
  end subroutine target_end


  ! ----------------------------------------------------------------------
  subroutine target_output(target, gr, dir, outp)
    type(target_t), intent(inout) :: target
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(output_t),   intent(in)  :: outp

    integer :: ierr
    call push_sub('target.target_output')

    select case(target%type)
    case(oct_tg_local)
      call doutput_function(outp%how, trim(dir), 'local_target', gr%m, gr%sb, &
        target%rho, M_ONE, ierr)
    case(oct_tg_td_local)
      call tdtarget_build_tdlocal(target, gr, M_ZERO)
      call doutput_function(outp%how, trim(dir), 'td_local_target', gr%m, gr%sb, &
        target%rho, M_ONE, ierr)
    case(oct_tg_density)
      call doutput_function(outp%how, trim(dir), 'density_target', gr%m, gr%sb, &
        target%rho, M_ONE, ierr)
    case(oct_tg_excited)
      call excited_states_output(target%est, gr, trim(dir), outp)
    case default
      call states_output(target%st, gr, trim(dir), outp)
    end select

    call pop_sub()
  end subroutine target_output


  ! ---------------------------------------------------------
  ! Calculates, at a given point in time marked by the integer
  ! index "i", the integrand of the target functional:
  ! <Psi(t)|\hat{O}(t)|Psi(t)>.
  ! This subroutine that the operator O is already updated
  ! ---------------------------------------------------------
  subroutine calc_tdfitness(target, gr, psi, i)
    type(target_t), intent(inout) :: target
    type(grid_t),      intent(in) :: gr
    type(states_t),    intent(in) :: psi
    integer,           intent(in) :: i
    CMPLX, allocatable :: opsi(:, :)
    integer :: p, j

    call push_sub('target.calc_tdfitness')

    target%td_fitness(i) = M_ZERO

    select case(target%type)
    case(oct_tg_td_local)

      select case(psi%d%ispin)
      case(UNPOLARIZED)
        ASSERT(psi%d%nik.eq.1)
        ALLOCATE(opsi(NP_PART, 1), NP_PART)
        opsi = M_z0
        do p  = psi%st_start, psi%st_end
          do j = 1, NP
            opsi(j, 1) = target%rho(j) * psi%zpsi(j, 1, p, 1)
          end do
          target%td_fitness(i) = target%td_fitness(i) + zstates_dotp(gr%m, psi%d%dim, psi%zpsi(:, :, p, 1), opsi(:, :))
        end do
        deallocate(opsi)
      case(SPIN_POLARIZED); stop 'Error'
      case(SPINORS);        stop 'Error'
      end select


    case default
      stop 'Error at calc_tdfitness'
    end select


    call pop_sub()
  end subroutine calc_tdfitness


  ! ---------------------------------------------------------------
  ! Calculates the inhomogeneous term that appears in the equation
  ! for chi, and places it into inh.
  ! ---------------------------------------------------------------
  subroutine calc_inh(psi, gr, target, t, inh)
    type(states_t),    intent(in)        :: psi
    type(grid_t),      intent(in)        :: gr
    type(target_t),    intent(inout)     :: target
    FLOAT,             intent(in)        :: t
    type(states_t),    intent(inout)     :: inh
 
    integer :: ik, ist, idim, i
    
    call push_sub('target.calc_inh')

    select case(target%type)
    case(oct_tg_td_local)
      call tdtarget_build_tdlocal(target, gr, t)
      do ik = 1, inh%d%nik
        do ist = inh%st_start, inh%st_end
          do idim = 1, inh%d%dim
            do i = 1, NP
              inh%zpsi(i, idim, ist, ik) = -M_zI * target%rho(i) * psi%zpsi(i, idim, ist, ik)
            end do
          end do
        end do
      end do
      
    end select

    call pop_sub()
  end subroutine calc_inh


  !----------------------------------------------------------
  ! 
  !----------------------------------------------------------
  subroutine tdtarget_init(target, gr, td)
    type(target_t), intent(inout) :: target
    type(grid_t),      intent(in) :: gr
    type(td_t),        intent(in) :: td

    C_POINTER         :: blk

    call push_sub('target.tdtarget_init')

    !%Variable OCTTdTarget
    !%Type block
    !%Section Optimal Control
    !%Description
    !% octopus features also time-dependent targets, i.e., one wants to optimize a laser 
    !% field that achieves a predefined time-dependent target. An example, could be the 
    !% evolution of occupation numbers in time. A time-dependent target consists of two 
    !% parts, i.e., the operator itself (a projection or a local operator) and its 
    !% time-dependence. 
    !%End
    if(loct_parse_block(check_inp('OCTTdTarget'),blk)==0) then
      select case(target%type)
      case(oct_tg_td_local)
        call loct_parse_block_string(blk, 0, 0, target%td_local_target)
        call conv_to_C_string(target%td_local_target)
        ALLOCATE(target%rho(NP), NP)
      case default
        stop 'ERROR in tdtarget_init'
      end select
    else
      message(1) = 'If OCTTargetMode = oct_targetmode_td, you must suppy a OCTTDTarget block'
      call write_fatal(1)
    end if

    if(target%mode .ne. oct_targetmode_td ) then
      nullify(target%td_fitness)
    else
      ALLOCATE(target%td_fitness(0:td%max_iter), td%max_iter+1)
    end if

    call pop_sub()
  end subroutine tdtarget_init
  !----------------------------------------------------------


  !----------------------------------------------------------
  subroutine tdtarget_build_tdlocal(target, gr, t)
    type(target_t), intent(inout) :: target
    type(grid_t),      intent(in) :: gr
    FLOAT, intent(in)             :: t
    integer :: i
    FLOAT :: xx(MAX_DIM), r, re, im
    call push_sub('target.target_build_tdlocal')

    do i = 1, NP
      call mesh_r(gr%m, i, r, x = xx)
      call loct_parse_expression(re, im, xx(1), xx(2), xx(3), r, t, target%td_local_target)
      target%rho(i) = re
    end do

    call pop_sub()
  end subroutine tdtarget_build_tdlocal
  !----------------------------------------------------------




end module opt_control_target_m

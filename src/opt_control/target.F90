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
  use grid_m
  use output_m
  use geometry_m
  use states_output_m
  use mesh_m
  use mesh_function_m
  use restart_m
  use opt_control_constants_m

  implicit none

  private
  public :: target_t,    &
            target_init, &
            target_end,  &
            target_output

  type target_t
    integer :: totype
    type(states_t) :: st
    FLOAT, pointer :: rho(:)
  end type target_t

  contains

  ! ----------------------------------------------------------------------
  subroutine target_output(gr, dir, outp, target)
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(output_t),   intent(in)  :: outp
    type(target_t), intent(inout) :: target

    integer :: ierr
    call push_sub('target.target_output')

    select case(target%totype)
    case(oct_tg_local)
      call doutput_function(outp%how, trim(dir), 'local_target', gr%m, gr%sb, &
        target%rho, M_ONE, ierr, is_tmp = .false.)
    case default
      call states_output(target%st, gr, trim(dir), outp)
    end select

    call pop_sub()
  end subroutine target_output

  ! ----------------------------------------------------------------------
  subroutine target_init(gr, geo, stin, target)
    type(grid_t), intent(in)      :: gr
    type(geometry_t), intent(in)  :: geo
    type(states_t), intent(in)    :: stin
    type(target_t), intent(inout) :: target

    integer           :: i, kpoints, dim, nst, no_c, state, ierr, isize, ik, ib, &
                         no_states, kk, p, ip, no_blk, ist, jst
    C_POINTER         :: blk
    FLOAT             :: x(MAX_DIM), r, psi_re, psi_im
    CMPLX             :: c_weight
    CMPLX, allocatable :: rotation_matrix(:, :)
    type(states_t)    :: tmp_st
    character(len=1024) :: expression
    character(len=10) :: fname

    call push_sub('target.target_init')

    !%Variable OCTTargetOperator
    !%Type integer
    !%Section Optimal Control
    !%Default 2
    !%Description
    !% The string OCTTargetOperator describes the initial state of the quantum system
    !% Possible arguments are:
    !%Option oct_tg_groundstate 1 
    !% Targetoperator is a projection operator on the ground state
    !%Option oct_tg_excited 2
    !% Currently not in use.
    !%Option oct_tg_gstransformation 3
    !% Targetoperator is a projection operator on a transformation of the ground state 
    !% orbitals defined by the block OCTTargetTransformStates
    !%Option oct_tg_userdefined 4
    !% Targetoperator is a projection operator on a user defined state
    !%Option oct_tg_local 5
    !% Targetoperator is a local operator. Specify the shape and position within the block OCTLocalTarget.
    !%Option oct_tg_td_local 6
    !% Target operator is time-dependent, please specify block OCTTdTarget
    !%End
    call loct_parse_int(check_inp('OCTTargetOperator'),oct_tg_excited, target%totype)
    if(.not.varinfo_valid_option('OCTTargetOperator', target%totype)) call input_error('OCTTargetOperator')    

    call states_copy(target%st, stin)

    select case(target%totype)
    case(oct_tg_groundstate)
      message(1) =  'Info: Using Ground State for TargetOperator'
      call write_info(1)
      call restart_read(trim(tmpdir)//'gs', target%st, gr, geo, ierr)
      if(ierr.ne.0) then
        write(message(1),'(a)') 'Could not read ground-state wavefunctions from '//trim(tmpdir)//'gs.'
        call write_fatal(1)
      end if
      
    case(oct_tg_excited) 
      message(1) = 'Error: using an excited state as the target state for an '
      message(2) = 'optimal control run is not possible yet.'
      message(3) = 'Try using "OCTInitialState = oct_is_transformation" instead.'
      call write_fatal(3)

    case(oct_tg_gstransformation)  

      message(1) =  'Info: Using Superposition of States for TargetOperator'
      call write_info(1)

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
      message(1) =  'Info: Using userdefined state for Targetoperator'
      call write_info(1)
      !%Variable OCTInitialUserdefined
      !%Type block
      !%Section Optimal Control
      !%Description
      !% 
      !% Example:
      !%
      !% <tt>%UserDefinedStates
      !% <br>&nbsp;&nbsp; 1 | 1 | 1 |  "exp(-r^2)*exp(-i*0.2*x)"
      !% <br>%</tt>
      !%  
      !%End

    case(oct_tg_local) 

      message(1) =  'Info: Using Local Target'
      call write_info(1)

      !%Variable OCTLocalTarget
      !%Type string
      !%Section Optimal Control
      !%Description
      !% If OCTTargetOperator = oct_tg_local, then one must supply the target density
      !% that should be searched for. This one can do by supplying a string through
      !% the variable OCTLocalTarget.
      !%End
      if(loct_parse_isdef('OCTLocalTarget').ne.0) then
        ALLOCATE(target%rho(NP), NP)
        target%rho = M_ZERO
        call loct_parse_string('OCTLocalTarget', "0", expression)
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
        !call lalg_scal(NP, M_ONE/r, target%rho)
        target%rho = target%rho/r
      else
        message(1) = 'If OCTTargetOperator = oct_tg_local, then you must give the shape'
        message(2) = 'of this target in variable "OCTLocalTarget"'
        call write_fatal(2)
      end if
      
    case(oct_tg_td_local)
      target%st%zpsi(:,1,1,1) = M_z0

    case default
      write(message(1),'(a)') "Target Operator not properly defined."
      call write_fatal(1)
    end select
    
    call pop_sub()
  end subroutine target_init


  subroutine target_end(target)
    type(target_t), intent(inout) :: target
    call push_sub('target.target_end')

    call states_end(target%st)
    if(target%totype.eq.oct_tg_local) then
      deallocate(target%rho)
      nullify(target%rho)
    end if

    call pop_sub()
  end subroutine target_end


end module opt_control_target_m

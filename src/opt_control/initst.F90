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

module opt_control_initst_m
  use datasets_m
  use varinfo_m
  use messages_m
  use loct_parser_m
  use global_m
  use string_m
  use states_m
  use states_calc_m
  use grid_m
  use geometry_m
  use profiling_m
  use restart_m

  implicit none

  private
  public :: initial_state_init

  integer, parameter ::  &
    oct_is_groundstate      = 1,      &
    oct_is_excited          = 2,      &
    oct_is_gstransformation = 3,      &
    oct_is_userdefined      = 4         


  contains


  ! ---------------------------------------------------------
  subroutine initial_state_init(gr, geo, st, initial_state)
    type(grid_t), intent(in)      :: gr
    type(geometry_t), intent(in)  :: geo
    type(states_t), intent(in)    :: st
    type(states_t), intent(inout) :: initial_state

    integer           :: ist, jst, ik, ib, idim, inst, inik, id, is, ip, ierr, no_states, istype
    type(block_t)     :: blk
    type(states_t)    :: tmp_st 
    FLOAT             :: x(MAX_DIM), r, psi_re, psi_im
    CMPLX, allocatable :: rotation_matrix(:, :)
    call push_sub('defstates.initial_state_init')

    call states_copy(initial_state, st)
    call states_allocate_wfns(initial_state, gr%mesh, M_CMPLX)

    !%Variable OCTInitialState
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default 1
    !%Description
    !% The string OCTInitialState describes the initial state of the quantum system
    !% Possible arguments are:
    !%Option oct_is_groundstate 1
    !% start in the ground state 
    !%Option oct_is_excited 2
    !% Currently not in use.
    !%Option oct_is_gstransformation 3
    !% start in a transformation of the ground-state orbitals, as defined in the
    !% block OCTInitialTransformStates
    !%Option oct_is_userdefined 4
    !% start in a userdefined state 
    !%End
    call loct_parse_int(datasets_check('OCTInitialState'), oct_is_groundstate, istype)
    if(.not.varinfo_valid_option('OCTInitialState', istype)) call input_error('OCTInitialState')    

    select case(istype)
    case(oct_is_groundstate) 
      message(1) =  'Info: Using Ground State for InitialState'
      call write_info(1)
      call restart_read(trim(restart_dir)//GS_DIR, initial_state, gr, geo, ierr)

    case(oct_is_excited)  
      message(1) = 'Error: using an excited state as the starting state for an '
      message(2) = 'optimal control run is not possible yet.'
      message(3) = 'Try using "OCTInitialState = oct_is_transformation" instead.'
      call write_fatal(3)

    case(oct_is_gstransformation)   
      message(1) =  'Info: Using Superposition of States for InitialState'
      call write_info(1)


      !%Variable OCTInitialTransformStates
      !%Type block
      !%Default no
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% If OCTInitialState = oct_is_gstransformation, you must specify one
      !% OCTInitialTransformStates block, in order to specify which linear
      !% combination of the states present in "restart/gs" is used to
      !% create the initial state.
      !% 
      !% The syntax is equivalent to the one used for the TransformStates
      !% block.
      !%End
      if(loct_parse_isdef(datasets_check('OCTInitialTransformStates')).ne.0) then
        if(loct_parse_block(datasets_check('OCTInitialTransformStates'), blk) == 0) then
          call states_copy(tmp_st, initial_state)
          SAFE_DEALLOCATE_P(tmp_st%zpsi)
          call restart_look_and_read(tmp_st, gr, geo)
          SAFE_ALLOCATE(rotation_matrix(1:initial_state%nst, 1:tmp_st%nst))
          rotation_matrix = M_z0
          do ist = 1, initial_state%nst
            do jst = 1, loct_parse_block_cols(blk, ist-1)
              call loct_parse_block_cmplx(blk, ist-1, jst-1, rotation_matrix(ist, jst))
            end do
          end do
          call rotate_states(gr%mesh, initial_state, tmp_st, rotation_matrix)
          SAFE_DEALLOCATE_A(rotation_matrix)
          call states_end(tmp_st)
        else
          message(1) = '"OCTInitialTransformStates" has to be specified as block.'
          call write_info(1)
          call input_error('OCTInitialTransformStates')
        end if
      else
        message(1) = 'Error: if "OCTInitialState = oct_is_gstransformation", then you must'
        message(2) = 'supply one "OCTInitialTransformStates" block to define the transformation.'
        call write_info(2)
        call input_error('OCTInitialTransformStates')
      end if


    case(oct_is_userdefined) 
      message(1) =  'Info: Building userdefined InitialState'
      call write_info(1)
      
      !%Variable OCTInitialUserdefined
      !%Type block
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% 
      !% Example:
      !%
      !% <tt>%UserDefinedStates
      !% <br>&nbsp;&nbsp; 1 | 1 | 1 |  "exp(-r^2)*exp(-i*0.2*x)"
      !% <br>%</tt>
      !%  
      !%End
      if(loct_parse_block(datasets_check('OCTInitialUserdefined'),blk)==0) then
        
        no_states = loct_parse_block_n(blk)
        do ib = 1, no_states
          write(6,*) 'HELLO USERDEF' 
          call loct_parse_block_int(blk, ib-1, 0, idim)
          call loct_parse_block_int(blk, ib-1, 1, inst)
          call loct_parse_block_int(blk, ib-1, 2, inik)
          write(6,*) ' DEBUG: ', idim,inst,inik
          ! read formula strings and convert to C strings
          do id = 1, initial_state%d%dim
            do is = 1, initial_state%nst
              do ik = 1, initial_state%d%nik   
                
                ! does the block entry match and is this node responsible?
                if(.not.(id.eq.idim .and. is.eq.inst .and. ik.eq.inik    &
                  .and. initial_state%st_start.le.is .and. initial_state%st_end.ge.is) ) cycle
                
                ! parse formula string
                call loct_parse_block_string(                            &
                  blk, ib-1, 3, initial_state%user_def_states(id, is, ik))
                ! convert to C string
                call conv_to_C_string(initial_state%user_def_states(id, is, ik))
                
                do ip = 1, gr%mesh%np
                  x = gr%mesh%x(ip, :)
                  r = sqrt(sum(x(:)**2))
                  
                  ! parse user defined expressions
                  call loct_parse_expression(psi_re, psi_im, gr%sb%dim, x, r, M_ZERO, initial_state%user_def_states(id, is, ik))
                  ! fill state
                  !write(6,*) psi_re, psi_im
                  initial_state%zpsi(ip, id, is, ik) = psi_re + M_zI*psi_im
                end do
                ! normalize orbital
                call zstates_normalize_orbital(gr%mesh, initial_state%d%dim, initial_state%zpsi(:,:, is, ik))
              end do
            end do
          enddo
        end do
        call loct_parse_block_end(blk)
      else
        message(1) = '"UserDefinedStates" has to be specified as block.'
        call write_fatal(1)
      end if
      
    case default
      write(message(1),'(a)') "No valid initial state defined."
      write(message(2),'(a)') "Choosing the ground state."
      call write_info(2)
    end select

    call states_calc_dens(initial_state, gr%mesh%np_part)
    
    call pop_sub()
  end subroutine initial_state_init

end module opt_control_initst_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:


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


  ! ---------------------------------------------------------
  subroutine def_istate(gr, geo, initial_state)
    type(grid_t), intent(in)      :: gr
    type(geometry_t), intent(in)  :: geo
    type(states_t), intent(inout) :: initial_state

    integer           :: i, kk, no_c, state, no_blk, kpoints, nst, dim, ist, jst
    C_POINTER         :: blk
    integer           :: p, ik, ib, idim, inst, inik
    integer           :: id ,is, ip, ierr, no_states, isize
    character(len=10) :: fname
    type(states_t)    :: tmp_st 
    FLOAT             :: x(MAX_DIM), r, psi_re, psi_im
    CMPLX             :: c_weight
    CMPLX, allocatable :: rotation_matrix(:, :)
    integer           :: istype
    
    call push_sub('opt_control.def_istate')


    !%Variable OCTInitialState
    !%Type integer
    !%Section Optimal Control
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
    call loct_parse_int(check_inp('OCTInitialState'), oct_is_groundstate, istype)
    if(.not.varinfo_valid_option('OCTInitialState', istype)) call input_error('OCTInitialState')    

    select case(istype)
    case(oct_is_groundstate) 
      message(1) =  'Info: Using Ground State for InitialState'
      call write_info(1)
      call restart_read(trim(tmpdir)//'gs', initial_state, gr, geo, ierr)

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
      !%Section OptimalControl
      !%Description
      !% If OCTInitialState = oct_is_gstransformation, you must specify one
      !% OCTInitialTransformStates block, in order to specify which linear
      !% combination of the states present in "restart/gs" is used to
      !% create the initial state.
      !% 
      !% The syntax is equivalent to the one used for the TransformStates
      !% block.
      !%End
      if(loct_parse_isdef(check_inp('OCTInitialTransformStates')).ne.0) then
        if(loct_parse_block(check_inp('OCTInitialTransformStates'), blk) == 0) then
          tmp_st = initial_state
          deallocate(tmp_st%zpsi)
          call restart_look_and_read("tmp", tmp_st, gr, geo, ierr)
          ALLOCATE(rotation_matrix(initial_state%nst, tmp_st%nst), initial_state%nst*tmp_st%nst)
          rotation_matrix = M_z0
          do ist = 1, initial_state%nst
            do jst = 1, loct_parse_block_cols(blk, ist-1)
              call loct_parse_block_cmplx(blk, ist-1, jst-1, rotation_matrix(ist, jst))
            end do
          end do
          call rotate_states(gr%m, initial_state, tmp_st, rotation_matrix)
          deallocate(rotation_matrix)
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
      if(loct_parse_block(check_inp('OCTInitialUserdefined'),blk)==0) then
        
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
                
                do ip = 1, gr%m%np
                  x = gr%m%x(ip, :)
                  r = sqrt(sum(x(:)**2))
                  
                  ! parse user defined expressions
                  call loct_parse_expression(psi_re, psi_im,             &
                    x(1), x(2), x(3), r, M_ZERO, initial_state%user_def_states(id, is, ik))
                  ! fill state
                  !write(6,*) psi_re, psi_im
                  initial_state%zpsi(ip, id, is, ik) = psi_re + M_zI*psi_im
                end do
                ! normalize orbital
                call zstates_normalize_orbital(gr%m, initial_state%d%dim, initial_state%zpsi(:,:, is, ik))
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
    
    call pop_sub()
  end subroutine def_istate


  ! ----------------------------------------------------------------------
  subroutine def_toperator(oct, gr, geo, target_state)
    type(oct_t), intent(in)       :: oct
    type(grid_t), intent(in)      :: gr
    type(geometry_t), intent(in)  :: geo
    type(states_t), intent(inout) :: target_state

    integer           :: i, kpoints, dim, nst, no_c, state, ierr, isize, ik, ib, &
                         no_states, kk, p, ip, no_blk, ist, jst
    C_POINTER         :: blk
    FLOAT             :: x(MAX_DIM), r, psi_re, psi_im
    CMPLX             :: c_weight
    CMPLX, allocatable :: rotation_matrix(:, :)
    type(states_t)    :: tmp_st
    character(len=1024) :: expression
    character(len=10) :: fname


    call push_sub('opt_control.def_toperator')

    select case(oct%totype)
    case(oct_tg_groundstate)
      message(1) =  'Info: Using Ground State for TargetOperator'
      call write_info(1)
      call restart_read(trim(tmpdir)//'gs', target_state, gr, geo, ierr)
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
      !%Section OptimalControl
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
          tmp_st = target_state
          deallocate(tmp_st%zpsi)
          call restart_look_and_read("tmp", tmp_st, gr, geo, ierr)
          ALLOCATE(rotation_matrix(target_state%nst, tmp_st%nst), target_state%nst*tmp_st%nst)
          rotation_matrix = M_z0
          do ist = 1, target_state%nst
            do jst = 1, loct_parse_block_cols(blk, ist-1)
              call loct_parse_block_cmplx(blk, ist-1, jst-1, rotation_matrix(ist, jst))
            end do
          end do
          call rotate_states(gr%m, target_state, tmp_st, rotation_matrix)
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
      !%Type block
      !%Section Optimal Control
      !%Description
      !% This block defines the shape and position of the local target operator. 
      !% It is possible to define several local operators which are then summed up.
      !% The syntax of each line is, then:
      !%
      !% <tt>%OCTLocalTarget
      !% <br>&nbsp;&nbsp; "exp(-r^2)*exp(-i*0.2*x)"
      !% <br>%</tt>
      !%  
      !% Example
      !%
      !% <tt>%OCTLocalTarget
      !% <br>&nbsp;&nbsp; "exp(-r^2)*exp(-i*0.2*x)"
      !% <br>%</tt>
      !%
      !%End
      ! parse input
      
      if(loct_parse_block(check_inp('OCTLocalTarget'),blk)==0) then
        no_states = loct_parse_block_n(blk)
        target_state%zpsi(:, 1, 1, 1) = m_z0

        do ib = 1, no_states
          ! parse formula string
          call loct_parse_block_string(blk, ib-1, 0, expression)
          ! convert to C string
          call conv_to_C_string(expression)
          
          do ip = 1, gr%m%np
            x = gr%m%x(ip, :)
            r = sqrt(sum(x(:)**2))
            
            ! parse user defined expressions
            call loct_parse_expression(psi_re, psi_im, &
              x(1), x(2), x(3), r, M_ZERO, expression)
            
            ! fill state
            target_state%zpsi(ip, 1, 1, 1) = target_state%zpsi(ip, 1, 1, 1) &
              + psi_re + M_zI*psi_im
          end do
          
          ! normalize orbital
          call zstates_normalize_orbital(gr%m, target_state%d%dim, &
            target_state%zpsi(:,:, 1, 1))
        end do
        call loct_parse_block_end(blk)
      else
        message(1) = '"OCTLocalTarget" has to be specified as block.'
        call write_fatal(1)
      end if

    case(oct_tg_td_local)
      target_state%zpsi(:,1,1,1) = M_z0

    case default
      write(message(1),'(a)') "Target Operator not properly defined."
      call write_fatal(1)
    end select
    
    call pop_sub()
  end subroutine def_toperator

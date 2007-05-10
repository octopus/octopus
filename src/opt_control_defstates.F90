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
  subroutine def_istate(oct, gr, geo, initial_state)
    type(oct_t), intent(in)       :: oct
    type(grid_t), intent(in)      :: gr
    type(geometry_t), intent(in)  :: geo
    type(states_t), intent(inout) :: initial_state

    integer           :: i, kk, no_c, state, no_blk, kpoints, nst, dim
    C_POINTER         :: blk
    integer           :: p, ik, ib, idim, inst, inik
    integer           :: id ,is, ip, ierr, no_states, isize
    character(len=10) :: fname
    type(states_t)    :: tmp_st 
    FLOAT             :: x(MAX_DIM), r, psi_re, psi_im
    CMPLX             :: c_weight
    
    call push_sub('opt_control.def_istate')

    select case(oct%istype)
    case(oct_is_groundstate) 
      message(1) =  'Info: Using Ground State for InitialState'
      call write_info(1)
      call restart_read(trim(tmpdir)//'restart_gs', initial_state, gr, geo, ierr)

    case(oct_is_excited)  
      message(1) =  'Info: Using Excited State for InitialState'
      call write_info(1)

      !%Variable OCTTargetStateNumber
      !%Type integer
      !%Section Optimal Control
      !%Default 2
      !%Description
      !% Specify the target state, ordered by energy
      !%End
      call loct_parse_int(check_inp('OCTInitialStateNumber'), 2, state)

      !TODO: The following lines of code do not look too clear, and will probably break easily.
      !They should also be isolated and taken away, since they are repeated in several places.
      tmp_st = initial_state
      call restart_look(trim(tmpdir)//'restart_gs', gr%m, kpoints, dim, nst, ierr)
      tmp_st%nst    = nst
      tmp_st%st_end = nst
      deallocate(tmp_st%eigenval)
      deallocate(tmp_st%occ)
      call states_allocate_wfns(tmp_st, gr%m)
      ALLOCATE(tmp_st%eigenval(tmp_st%nst, tmp_st%d%nik), tmp_st%nst*tmp_st%d%nik)
      ALLOCATE(tmp_st%occ(tmp_st%nst, tmp_st%d%nik), tmp_st%nst*tmp_st%d%nik)
      if(tmp_st%d%ispin == SPINORS) then
        ALLOCATE(tmp_st%spin(3, tmp_st%nst, tmp_st%d%nik), tmp_st%nst*tmp_st%d%nik*3)
        tmp_st%spin = M_ZERO
      end if
      tmp_st%eigenval = huge(REAL_PRECISION)
      tmp_st%occ      = M_ZERO
      call restart_read(trim(tmpdir)//'restart_gs', tmp_st, gr, geo, ierr)
      
      initial_state%zpsi(:, :, 1, 1) = tmp_st%zpsi(:, :, state, 1)

    case(oct_is_superposition)   
      message(1) =  'Info: Using Superposition of States for InitialState'
      call write_info(1)

      !%Variable OCTInitialSuperposition
      !%Type block
      !%Section Optimal Control
      !%Description
      !% The syntax of each line is, then:
      !%
      !% <tt>%OCTInitialSuperposition
      !% <br>&nbsp;&nbsp;state1 | weight1 | state2 | weight2 | ... 
      !% <br>%</tt>
      !%
      !% Note that weight can be complex, to produce current carrying superpositions.
      !% Example:
      !%
      !% <tt>%OCTInitialSuperposition
      !% <br>&nbsp;&nbsp;1 | 1 | 2 | -1 | ... 
      !% <br>%</tt>
      !%
      !%End 
      if(loct_parse_block(check_inp('OCTInitialSuperposition'),blk)==0) then

        tmp_st = initial_state
        call restart_look(trim(tmpdir)//'restart_gs', gr%m, kpoints, dim, nst, ierr)
        tmp_st%nst    = nst
        tmp_st%st_end = nst
        deallocate(tmp_st%eigenval)
        deallocate(tmp_st%occ)
        call states_allocate_wfns(tmp_st, gr%m)
        ALLOCATE(tmp_st%eigenval(tmp_st%nst, tmp_st%d%nik), tmp_st%nst*tmp_st%d%nik)
        ALLOCATE(tmp_st%occ(tmp_st%nst, tmp_st%d%nik), tmp_st%nst*tmp_st%d%nik)
        if(tmp_st%d%ispin == SPINORS) then
          ALLOCATE(tmp_st%spin(3, tmp_st%nst, tmp_st%d%nik), 3*tmp_st%nst*tmp_st%d%nik)
          tmp_st%spin = M_ZERO
        end if
        tmp_st%eigenval = huge(REAL_PRECISION)
        tmp_st%occ      = M_ZERO
        call restart_read(trim(tmpdir)//'restart_gs', tmp_st, gr, geo, ierr)


        no_blk = loct_parse_block_n(blk)
        ! TODO: One should check that the states requested are actually present in tmp_st   
        do i=1, no_blk
          ! The structure of the block is:
          ! domain | function_type | center | width | weight 
          no_c = loct_parse_block_cols(blk, i-1)
          initial_state%zpsi(:,:,i,1) = M_z0
          do kk=1, no_c, 2
            call loct_parse_block_int(blk, i-1, kk-1, state)
            call loct_parse_block_cmplx(blk, i-1, kk, c_weight)
            initial_state%zpsi(:,:,1,1) = initial_state%zpsi(:,:,1,1) + c_weight * tmp_st%zpsi(:,:,state,1)
          enddo
        end do
        call states_end(tmp_st)

        ! normalize state
        do ik = 1, initial_state%d%nik
          do p  = initial_state%st_start, initial_state%st_end
            call zstates_normalize_orbital(gr%m, initial_state%d%dim, &
              initial_state%zpsi(:,:, p, ik))
          enddo
        enddo

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
                         no_states, kk, p, ip, no_blk
    C_POINTER         :: blk
    FLOAT             :: x(MAX_DIM), r, psi_re, psi_im
    CMPLX             :: c_weight      
    type(states_t)    :: tmp_st
    character(len=1024) :: expression
    character(len=10) :: fname


    call push_sub('opt_control.def_toperator')

    select case(oct%totype)
    case(oct_tg_groundstate)
      message(1) =  'Info: Using Ground State for TargetOperator'
      call write_info(1)
      call restart_read(trim(tmpdir)//'restart_gs', target_state, gr, geo, ierr)
      
    case(oct_tg_excited) 
      message(1) =  'Info: Using Excited State for TargetOperator'
      call write_info(1)

      !%Variable OCTTargetStateNumber
      !%Type integer
      !%Section Optimal Control
      !%Default 2
      !%Description
      !% Specify the target state, ordered by energy.
      !%End
      call loct_parse_int(check_inp('OCTTargetStateNumber'), 2, state)

      !TODO: The following lines of code do not look too clear, and will probably break easily.
      tmp_st = target_state
      call restart_look(trim(tmpdir)//'restart_gs', gr%m, kpoints, dim, nst, ierr)
      tmp_st%nst    = nst
      tmp_st%st_end = nst
      deallocate(tmp_st%eigenval)
      deallocate(tmp_st%occ)
      call states_allocate_wfns(tmp_st, gr%m)
      ALLOCATE(tmp_st%eigenval(tmp_st%nst, tmp_st%d%nik), tmp_st%nst*tmp_st%d%nik)
      ALLOCATE(tmp_st%occ(tmp_st%nst, tmp_st%d%nik), tmp_st%nst*tmp_st%d%nik)
      if(tmp_st%d%ispin == SPINORS) then
        ALLOCATE(tmp_st%spin(3, tmp_st%nst, tmp_st%d%nik), tmp_st%nst*tmp_st%d%nik*3)
        tmp_st%spin = M_ZERO
      end if
      tmp_st%eigenval = huge(REAL_PRECISION)
      tmp_st%occ      = M_ZERO
      call restart_read(trim(tmpdir)//'restart_gs', tmp_st, gr, geo, ierr)

      target_state%zpsi(:, :, 1, 1) = tmp_st%zpsi(:, :, state, 1)
      call states_end(tmp_st)

    case(oct_tg_superposition)  
      message(1) =  'Info: Using Superposition of States for TargetOperator'
      call write_info(1)

      !%Variable OCTTargetSuperposition
      !%Type block
      !%Section Optimal Control
      !%Description
      !% The syntax of each line is, then:
      !%
      !% <tt>%OCTTargetSuperposition
      !% <br>&nbsp;&nbsp;state1 | weight1 | state2 | weight2 | ... 
      !% <br>%</tt>
      !% 
      !% Example:
      !%
      !% <tt>%OCTTargetSuperposition
      !% <br>&nbsp;&nbsp;1 | 1 | 2 | -1 | ... 
      !% <br>%</tt>
      !%
      !%End
      if(loct_parse_block(check_inp('OCTTargetSuperposition'),blk)==0) then

        tmp_st = target_state
        call restart_look(trim(tmpdir)//'restart_gs', gr%m, kpoints, dim, nst, ierr)
        tmp_st%nst    = nst
        tmp_st%st_end = nst
        deallocate(tmp_st%eigenval)
        deallocate(tmp_st%occ)
        call states_allocate_wfns(tmp_st, gr%m)
        ALLOCATE(tmp_st%eigenval(tmp_st%nst, tmp_st%d%nik), tmp_st%nst*tmp_st%d%nik)
        ALLOCATE(tmp_st%occ(tmp_st%nst, tmp_st%d%nik), tmp_st%nst*tmp_st%d%nik)
        if(tmp_st%d%ispin == SPINORS) then
          ALLOCATE(tmp_st%spin(3, tmp_st%nst, tmp_st%d%nik), tmp_st%nst*tmp_st%d%nik*3)
          tmp_st%spin = M_ZERO
        end if
        tmp_st%eigenval = huge(REAL_PRECISION)
        tmp_st%occ      = M_ZERO
        call restart_read(trim(tmpdir)//'restart_gs', tmp_st, gr, geo, ierr)

        no_blk = loct_parse_block_n(blk)
        ! TODO: for each orbital            
        do i=1, no_blk
          ! The structure of the block is:
          ! domain | function_type | center | width | weight 
          no_c = loct_parse_block_cols(blk, i-1)
          target_state%zpsi(:,:,i,1) = M_z0
          do kk=1, no_c, 2
            call loct_parse_block_int(blk, i-1, kk-1, state)
            call loct_parse_block_cmplx(blk, i-1, kk, c_weight)
            target_state%zpsi(:, :, 1, 1) = target_state%zpsi(:, :, 1, 1) + c_weight * tmp_st%zpsi(:, :, state, 1)
          enddo
        end do
        call states_end(tmp_st)

        ! normalize state
        do ik = 1, target_state%d%nik
          do p  = target_state%st_start, target_state%st_end
            call zstates_normalize_orbital(gr%m, target_state%d%dim, target_state%zpsi(:,:, p, ik))
          enddo
        enddo

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

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
  subroutine penalty_init(penalty, oct, gr, steps, ctr_iter_max, dt)
    type(oct_penalty_t), intent(inout) :: penalty
    type(oct_t), intent(in)  :: oct
    type(grid_t), pointer    :: gr
    integer,      intent(in) :: steps
    integer,      intent(in) :: ctr_iter_max
    FLOAT,        intent(in) :: dt

    C_POINTER                :: blk
    integer                  :: no_lines, no_col, i
    type(filter_t),  pointer :: tdp(:)

    call push_sub('opt_control_penalty.penalty_init')

    !%Variable OCTPenalty
    !%Type float
    !%Section Optimal Control
    !%Default 1.0
    !%Description
    !% The variable specificies the value of the penalty factor for the 
    !% integrated field strength (fluence). Large value - small fluence. The value is 
    !% always positive. A transient shape can be specified using the block OCTLaserEnvelope.
    !% In this case OCTPenalty is multiplied with time-dependent function. 
    !% The value depends on the coupling between the states. A good start might be a 
    !% value from 0.1 (strong fields) to 10 (weak fields). 
    !%End
    call loct_parse_float(check_inp('OCTPenalty'), M_ONE, penalty%penalty)
    ! penalty array for fixed fluence run 
    ! the array is only interesting for the development of new algorithms

    ALLOCATE(penalty%a_penalty(0:ctr_iter_max+1), ctr_iter_max+2)
    penalty%a_penalty = penalty%penalty
    ALLOCATE(penalty%tdpenalty(NDIM, 0:2*steps), NDIM*(2*steps+1))
    penalty%tdpenalty = M_ONE

    penalty%mode_tdpenalty = .false.

    !%Variable OCTLaserEnvelope
    !%Type block
    !%Section Optimal Control
    !%Description
    !% Often a predefined time-dependent envelope on the laser field is required. 
    !% This can be achieved by making the penalty factor time-dependent. 
    !% Here, you may specify the required time dependent envelope.
    !% The code allows for more than one enevelope, all of them will be added together and can weighted against each other. 
    !% It is possible to choose different envelopes for different polarization directions. 
    !% The block is similar to the time filter.
    !%End
        
    if (loct_parse_block(check_inp('OCTLaserEnvelope'), blk)==0) then
      no_lines = loct_parse_block_n(blk)
      ALLOCATE(tdp(no_lines), no_lines)
      penalty%mode_tdpenalty = .true.
      
      ! The structure of the block is:
      ! polarization | function | weight 
      do i = 1, no_lines
        no_col = loct_parse_block_cols(blk, i-1)
        call loct_parse_block_cmplx(blk, i-1, 0, tdp(i)%fpol(1)) 
        call loct_parse_block_cmplx(blk, i-1, 1, tdp(i)%fpol(2))
        call loct_parse_block_cmplx(blk, i-1, 2, tdp(i)%fpol(3))
        ! parse formula string
        call loct_parse_block_string(                            &
          blk, i-1, 3, tdp(i)%expression)  
        ! convert to C string
        call conv_to_C_string(tdp(i)%expression)
        call loct_parse_block_float(blk, i-1, 4, tdp(i)%weight)
        !call loct_parse_block_int(blk, i-1, 1, tdp(i)%ftype)
        !call loct_parse_block_float(blk, i-1, 2, tdp(i)%center)
        !call loct_parse_block_float(blk, i-1, 3, tdp(i)%width)
        
        tdp(i)%domain = filter_time
        
        ALLOCATE(tdp(i)%numerical(NDIM,0:2*steps), NDIM*(2*steps+1))
        call build_filter(gr, tdp(i), steps, dt)
      end do

      do i = 1, no_lines      
        tdp(i)%weight = tdp(i)%weight/SUM(tdp(1:no_lines)%weight)
      end do

      penalty%tdpenalty = M_ZERO
      do i=1,no_lines
        
        penalty%tdpenalty =  penalty%tdpenalty &
          + tdp(i)%weight * real(tdp(i)%numerical,REAL_PRECISION)
      end do

      penalty%tdpenalty = M_ONE / (penalty%tdpenalty + CNST(0.0000001) )
      
      ! all we want to know is tdpenalty
      do i = 1, no_lines
        if(associated(tdp(i)%numerical)) then
          deallocate(tdp(i)%numerical) 
          nullify(tdp(i)%numerical)
        end if
      end do
      call loct_parse_block_end(blk)
    end if

    penalty%tdpenalty = penalty%tdpenalty * penalty%penalty

    call pop_sub()
  end subroutine penalty_init

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

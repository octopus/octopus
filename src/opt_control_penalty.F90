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
  subroutine penalty_init(penalty, oct, par, steps, ctr_iter_max, dt)
    type(oct_penalty_t), intent(inout)         :: penalty
    type(oct_t), intent(in)                    :: oct
    type(oct_control_parameters_t), intent(in) :: par
    integer,      intent(in)                   :: steps
    integer,      intent(in)                   :: ctr_iter_max
    FLOAT,        intent(in)                   :: dt

    character(len=1024)      :: expression
    FLOAT                    :: weight, t, f_re, f_im
    C_POINTER                :: blk
    integer                  :: no_lines, i, j
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

    !!! WARNING This should be done only if there are td penalties required?
    ALLOCATE(penalty%td_penalty(par%no_parameters), par%no_parameters)
    do i = 1, par%no_parameters
      call tdf_init_numerical(penalty%td_penalty(i), 2*steps, dt*M_HALF, initval = M_z1)
    end do

    penalty%mode_tdpenalty = .false.

    !%Variable OCTLaserEnvelope
    !%Type block
    !%Section Optimal Control
    !%Description
    !% Often a predefined time-dependent envelope on the control parameter is desired. 
    !% This can be achieved by making the penalty factor time-dependent. 
    !% Here, you may specify the required time dependent envelope.
    !%
    !% It is possible to choose different envelopes for different control parameters.
    !% There should be one line for each control parameter.
    !%End
        
    if (loct_parse_block(check_inp('OCTLaserEnvelope'), blk)==0) then
      no_lines = loct_parse_block_n(blk)
      if(no_lines .ne.par%no_parameters) call input_error('OCTLaserEnvelope')

      penalty%mode_tdpenalty = .true.
      
      ! The structure of the block is:
      ! function | weight 
      do i = 1, no_lines
        call loct_parse_block_string(blk, i-1, 0, expression)  
        call conv_to_C_string(expression)
        call loct_parse_block_float(blk, i-1, 1, weight)
        do j = 1, 2*steps+1
          t = (j-1)*M_HALF*dt
          call loct_parse_expression(f_re, f_im, "t", t, expression)
          call tdf_set_numerical(penalty%td_penalty(i), j, M_ONE /(f_re + CNST(1.0e-7))  )
        end do
      end do

      call loct_parse_block_end(blk)
    end if

    do i = 1, par%no_parameters
      call tdf_scalar_multiply(penalty%penalty, penalty%td_penalty(i))
    end do

    call pop_sub()
  end subroutine penalty_init

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_classical(ions, namespace, tg, td, oct)
    type(ions_t),      intent(in)    :: ions
    type(namespace_t), intent(in)    :: namespace
    type(target_t),    intent(inout) :: tg
    type(td_t),        intent(in)    :: td
    type(oct_t),       intent(in)    :: oct

    integer             :: jj, ist, jst
    type(block_t)       :: blk
    character(len=1024) :: expression
    PUSH_SUB(target_init_classical)

    tg%move_ions = ion_dynamics_ions_move(td%ions_dyn)
    tg%dt = td%dt
    ASSERT(tg%move_ions)

    !%Variable OCTClassicalTarget
    !%Type block
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% If <tt>OCTTargetOperator = oct_tg_classical</tt>, the you must supply this block.
    !% It should contain a string (e.g. "(q[1,1]-q[1,2])*p[2,1]") with a mathematical
    !% expression in terms of two arrays, q and p, that represent the position and momenta
    !% of the classical variables. The first index runs through the various classical particles,
    !% and the second index runs through the spatial dimensions.
    !%
    !% In principle, the block only contains one entry (string). However, if the expression is very
    !% long, you can split it into various lines (one column each) that will be concatenated.
    !%
    !% The QOCT algorithm will attempt to maximize this expression, at the end of the propagation.
    !%End
    if(parse_block(namespace, 'OCTClassicalTarget', blk)==0) then
      tg%classical_input_string = " "
      do jj=0, parse_block_n(blk)-1
        call parse_block_string(blk, jj, 0, expression)
        tg%classical_input_string = trim(tg%classical_input_string) // trim(expression)
      end do
      call parse_block_end(blk)
    else
      message(1) = 'If OCTTargetOperator = oct_tg_classical, then you must give the shape'
      message(2) = 'of this target in the block "OCTClassicalTarget".'
      call messages_fatal(2)
    end if

    !%Variable OCTMomentumDerivatives
    !%Type block
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% This block should contain the derivatives of the expression given in  
    !% <tt>OCTClassicalTarget</tt> with respect to the p array components.
    !% Each line corresponds to a different classical particle, whereas the 
    !% columns correspond to each spatial dimension: the (i,j) block component
    !% corresponds with the derivative wrt p[i,j].
    !%End
    if( parse_block(namespace, 'OCTMomentumDerivatives', blk)==0   ) then
      SAFE_ALLOCATE(tg%mom_der_array(1:ions%natoms,1:ions%space%dim))
      do ist=0, ions%natoms-1
        do jst=0, ions%space%dim-1
          call parse_block_string(blk, ist, jst, tg%mom_der_array(ist+1, jst+1))
        end do
      end do
      call parse_block_end(blk)
    elseif(oct%algorithm  ==  OPTION__OCTSCHEME__OCT_CG .or. oct%algorithm == OPTION__OCTSCHEME__OCT_BFGS) then
      message(1) = 'If "OCTTargetOperator = oct_classical" and "OCTScheme = oct_cg" or'
      message(2) = '"OCTScheme = oct_bfgs", then you must define the blocks "OCTClassicalTarget",' 
      message(3) = '"OCTPositionDerivatives" AND "OCTMomentumDerivatives"'
      call messages_fatal(3)
    end if

    !%Variable OCTPositionDerivatives
    !%Type block
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% This block should contain the derivatives of the expression given in  
    !% <tt>OCTClassicalTarget</tt> with respect to the q array components.
    !% Each line corresponds to a different classical particle, whereas the 
    !% columns correspond to each spatial dimension: the (i,j) block component
    !% corresponds with the derivative wrt q[i,j].
    !%End
    if( parse_block(namespace, 'OCTPositionDerivatives', blk)==0 ) then
      SAFE_ALLOCATE(tg%pos_der_array(1:ions%natoms,1:ions%space%dim))
      do ist=0, ions%natoms-1
        do jst=0, ions%space%dim-1
          call parse_block_string(blk, ist, jst, tg%pos_der_array(ist+1, jst+1))
        end do
      end do
      call parse_block_end(blk)
    elseif(oct%algorithm  ==  OPTION__OCTSCHEME__OCT_CG .or. oct%algorithm == OPTION__OCTSCHEME__OCT_BFGS) then
      message(1) = 'If "OCTTargetOperator = oct_tg_classical" and "OCTScheme = oct_cg" or'
      message(2) = '"OCTScheme = oct_bfgs", then you must define the blocks "OCTClassicalTarget",'
      message(3) = '"OCTPositionDerivatives" AND "OCTMomentumDerivatives"'
      call messages_fatal(3)
    end if

    POP_SUB(target_init_classical)
  end subroutine target_init_classical
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_classical(tg)
    type(target_t),   intent(inout) :: tg

    PUSH_SUB(target_end_classical)

    SAFE_DEALLOCATE_A(tg%pos_der_array)
    SAFE_DEALLOCATE_A(tg%mom_der_array)

    POP_SUB(target_end_classical)
  end subroutine target_end_classical
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  subroutine target_output_classical
    PUSH_SUB(target_output_classical)
    POP_SUB(target_output_classical)
  end subroutine target_output_classical
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_classical(tg, qcpsi) result(j1)
    type(target_t),            intent(inout) :: tg
    type(opt_control_state_t), intent(in)    :: qcpsi

    FLOAT, pointer :: q(:, :), p(:, :)
    FLOAT :: dummy(3)
    character(len=4096) :: inp_string

    PUSH_SUB(target_j1_classical)

    q => opt_control_point_q(qcpsi)
    p => opt_control_point_p(qcpsi)

    inp_string = tg%classical_input_string
    call parse_array(inp_string, q, 'q')
    call parse_array(inp_string, p, 'p')
    call conv_to_C_string(inp_string)
    call parse_expression(j1, dummy(1), 1, dummy(1:3), dummy(1), dummy(1), inp_string)

    nullify(q)
    POP_SUB(target_j1_classical)
  end function target_j1_classical
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_classical(tg, qcpsi, qcchi, ions)
    type(target_t),            intent(inout) :: tg
    type(opt_control_state_t), intent(inout) :: qcpsi
    type(opt_control_state_t), intent(inout) :: qcchi
    type(ions_t),     intent(in)    :: ions

    integer :: ist, jst, ib, iqn
    character(len=1024) :: temp_string
    FLOAT :: df_dv, dummy(3)
    FLOAT, pointer :: q(:, :), p(:, :), tq(:, :), tp(:, :)
    type(states_elec_t), pointer :: chi
    PUSH_SUB(target_chi_classical)
      
    tq => opt_control_point_q(qcchi)
    tp => opt_control_point_p(qcchi)
    q => opt_control_point_q(qcpsi)
    p => opt_control_point_p(qcpsi)
    tq = M_ZERO
    tp = M_ZERO

    do ist = 1, ions%natoms
      do jst=1, ions%space%dim
        temp_string = tg%mom_der_array(ist, jst)
        call parse_array(temp_string, p, 'p')
        call parse_array(temp_string, q, 'q')
        call conv_to_C_string(temp_string)
        call parse_expression(df_dv, dummy(1), 1, dummy(1:3), dummy(1), dummy(1), temp_string)
        tq(ist, jst) = -df_dv
      end do
    end do

    do ist = 1, ions%natoms
      do jst=1, ions%space%dim
        temp_string = tg%pos_der_array(ist, jst)
        call parse_array(temp_string, p, 'p')
        call parse_array(temp_string, q, 'q')
        call conv_to_C_string(temp_string)
        call parse_expression(df_dv, dummy(1), 1, dummy(1:3), dummy(1), dummy(1), temp_string)
        tp(ist, jst) = df_dv
      end do
    end do

    chi => opt_control_point_qs(qcchi)

    !We assume that there is no time-independent operator.
    do iqn = chi%d%kpt%start, chi%d%kpt%end
      do ib = chi%group%block_start, chi%group%block_end
        call batch_set_zero(chi%group%psib(ib, iqn))
      end do
    end do

    nullify(tq)
    nullify(tp)
    nullify(q)
    nullify(p)
    nullify(chi)
    POP_SUB(target_chi_classical)
  end subroutine target_chi_classical
  ! ----------------------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

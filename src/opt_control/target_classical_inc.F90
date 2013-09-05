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
!! $Id: target_velocity_inc.F90 $


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_classical(geo, tg, td)
    type(geometry_t), intent(in)    :: geo
    type(target_t),   intent(inout) :: tg
    type(td_t),       intent(in)    :: td

    integer             :: jj, ist, jst
    type(block_t)       :: blk
    character(len=1024) :: expression
    PUSH_SUB(target_init_classical)

    tg%move_ions = ion_dynamics_ions_move(td%ions)
    tg%dt = td%dt
    ASSERT(tg%move_ions)

    !%Variable OCTClassicalTarget
    !%Type block
    !%Section Calculation Modes::Optimal Control
    !%Description
    !%
    !%End

    !%Variable OCTClassicalDerivatives
    !%Type block
    !%Section Calculation Modes::Optimal Control
    !%Description
    !%
    !%End

    if(parse_block(datasets_check('OCTClassicalTarget'),blk)==0) then
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

    if( parse_block(datasets_check('OCTMomentumDerivatives'),blk)==0   ) then
      SAFE_ALLOCATE(tg%mom_der_array(1:geo%natoms,1:geo%space%dim))
      do ist=0, geo%natoms-1
        do jst=0, geo%space%dim-1
          call parse_block_string(blk, ist, jst, tg%mom_der_array(ist+1, jst+1))
        end do
      end do
      call parse_block_end(blk)
    else
      message(1) = 'If OCTTargetOperator = oct_tg_classical, then you must define the'
      message(2) = 'blocks "OCTClassicalTarget", "OCTPositionDerivatives" AND "OCTVelocityDerivatives"'
      call messages_fatal(2)
    end if

    if( parse_block(datasets_check('OCTPositionDerivatives'),blk)==0 ) then
      SAFE_ALLOCATE(tg%pos_der_array(1:geo%natoms,1:geo%space%dim))
      do ist=0, geo%natoms-1
        do jst=0, geo%space%dim-1
          call parse_block_string(blk, ist, jst, tg%pos_der_array(ist+1, jst+1))
        end do
      end do
      call parse_block_end(blk)
    else
      message(1) = 'If OCTTargetOperator = oct_tg_classical, then you must define the'
      message(2) = 'blocks "OCTClassicalTarget", "OCTPositionDerivatives" AND "OCTVelocityDerivatives"'
      call messages_fatal(2)
    end if

    POP_SUB(target_init_classical)
  end subroutine target_init_classical
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_classical(tg)
    type(target_t),   intent(inout) :: tg
    PUSH_SUB(target_end_classical)
    SAFE_DEALLOCATE_P(tg%pos_der_array)
    SAFE_DEALLOCATE_P(tg%mom_der_array)
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
    !j1 = q(1, 1)

    nullify(q)
    POP_SUB(target_j1_classical)
  end function target_j1_classical
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_classical(tg, qcpsi, qcchi, geo)
    type(target_t),            intent(inout) :: tg
    type(opt_control_state_t), intent(inout) :: qcpsi
    type(opt_control_state_t), intent(inout) :: qcchi
    type(geometry_t), intent(in)    :: geo

    integer :: ist, jst
    character(len=1024) :: temp_string
    FLOAT :: df_dv, dummy(3)
    FLOAT, pointer :: q(:, :), p(:, :), tq(:, :), tp(:, :)
    PUSH_SUB(target_chi_classical)
      
    tq => opt_control_point_q(qcchi)
    tp => opt_control_point_p(qcchi)
    q => opt_control_point_q(qcpsi)
    p => opt_control_point_p(qcpsi)
    tq = M_ZERO
    tp = M_ZERO

    do ist=1, geo%natoms
      do jst=1, geo%space%dim
        temp_string = tg%mom_der_array(ist, jst)
        call parse_array(temp_string, p, 'p')
        call parse_array(temp_string, q, 'q')
        call conv_to_C_string(temp_string)
        call parse_expression(df_dv, dummy(1), 1, dummy(1:3), dummy(1), dummy(1), temp_string)
        tq(ist, jst) = -df_dv
      end do
    end do

    do ist=1, geo%natoms
      do jst=1, geo%space%dim
        temp_string = tg%pos_der_array(ist, jst)
        call parse_array(temp_string, p, 'p')
        call parse_array(temp_string, q, 'q')
        call conv_to_C_string(temp_string)
        call parse_expression(df_dv, dummy(1), 1, dummy(1:3), dummy(1), dummy(1), temp_string)
        tp(ist, jst) = df_dv
      end do
    end do

    nullify(tq)
    nullify(tp)
    nullify(q)
    nullify(p)
    POP_SUB(target_chi_classical)
  end subroutine target_chi_classical
  ! ----------------------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

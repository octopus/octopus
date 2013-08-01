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
  subroutine target_init_classical(tg, td)
    type(target_t),   intent(inout) :: tg
    type(td_t),       intent(in)    :: td
    PUSH_SUB(target_init_classical)

    tg%move_ions = ion_dynamics_ions_move(td%ions)
    tg%dt = td%dt

    POP_SUB(target_init_classical)
  end subroutine target_init_classical
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_classical
    PUSH_SUB(target_end_classical)
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

    FLOAT, pointer :: q(:, :)
    PUSH_SUB(target_j1_classical)

    q => opt_control_point_q(qcpsi)
    j1 = q(1, 1)

    nullify(q)
    POP_SUB(target_j1_classical)
  end function target_j1_classical
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_classical(tg, qcchi)
    type(target_t),            intent(inout) :: tg
    type(opt_control_state_t), intent(inout) :: qcchi

    FLOAT, pointer :: q(:, :), p(:, :)
    PUSH_SUB(target_chi_classical)
      
    q => opt_control_point_q(qcchi)
    p => opt_control_point_p(qcchi)

    q(1, 1) = M_ZERO
    p(1, 1) = M_ONE

    nullify(q)
    nullify(p)
    POP_SUB(target_chi_classical)
  end subroutine target_chi_classical
  ! ----------------------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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
!! $Id: target_tdlocal_inc.F90 $


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_tdlocal(gr, tg, td)
    type(grid_t),     intent(in)    :: gr
    type(target_t),   intent(inout) :: tg
    type(td_t),       intent(in)    :: td

    type(block_t)       :: blk
    PUSH_SUB(target_init_tdlocal)

    if(parse_block(datasets_check('OCTTdTarget'), blk)==0) then
      call parse_block_string(blk, 0, 0, tg%td_local_target)
      call conv_to_C_string(tg%td_local_target)
      SAFE_ALLOCATE(tg%rho(1:gr%mesh%np))
      call parse_block_end(blk)
    else
      message(1) = 'If OCTTargetMode = oct_targetmode_td, you must suppy a OCTTDTarget block.'
      call messages_fatal(1)
    end if
    tg%dt = td%dt
    SAFE_ALLOCATE(tg%td_fitness(0:td%max_iter))
    tg%td_fitness = M_ZERO
    call target_build_tdlocal(tg, gr, M_ZERO)

    POP_SUB(target_init_tdlocal)
  end subroutine target_init_tdlocal


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_tdlocal(tg) result(j1)
    type(target_t),   intent(inout) :: tg

    integer :: maxiter
    PUSH_SUB(target_j1_tdlocal)

    maxiter = size(tg%td_fitness) - 1
    j1 = M_HALF * tg%dt * tg%td_fitness(0) + & 
         M_HALF * tg%dt * tg%td_fitness(maxiter) + & 
         tg%dt * sum(tg%td_fitness(1:maxiter-1))


    POP_SUB(target_j1_tdlocal)
  end function target_j1_tdlocal


  !----------------------------------------------------------
  subroutine target_build_tdlocal(tg, gr, time)
    type(target_t), intent(inout) :: tg
    type(grid_t),   intent(in)    :: gr
    FLOAT,          intent(in)    :: time

    integer :: ip
    FLOAT :: xx(MAX_DIM), rr, re, im

    PUSH_SUB(target_build_tdlocal)

    do ip = 1, gr%mesh%np
      call mesh_r(gr%mesh, ip, rr, coords = xx)
      call parse_expression(re, im, gr%sb%dim, xx, rr, time, tg%td_local_target)
      tg%rho(ip) = re
    end do

    POP_SUB(target_build_tdlocal)
  end subroutine target_build_tdlocal
  !----------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

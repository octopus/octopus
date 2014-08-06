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
!! $Id$


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_spin(tg)
    type(target_t),   intent(inout) :: tg

    type(block_t)       :: blk
    integer :: jst
    CMPLX :: alpha(3)

    PUSH_SUB(target_init_spin)


    message(1) =  'Info: Using a spin target'
    call messages_info(1)

    !%Variable OCTTargetSpin
    !%Type block
    !%Default no
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% (EXPERIMENTAL)
    !%
    !%End
    if(parse_isdef(datasets_check('OCTTargetSpin')) /= 0) then
      if(parse_block(datasets_check('OCTTargetSpin'), blk) == 0) then
        alpha = M_z0
        do jst = 1, parse_block_cols(blk, 0)
          call parse_block_cmplx(blk, 0, jst - 1, alpha(jst))
        end do
        call parse_block_end(blk)
        
        alpha = alpha/sqrt(dot_product(alpha, alpha))

        tg%spin_matrix(1, 1) = alpha(3)
        tg%spin_matrix(2, 2) = -alpha(3)
        tg%spin_matrix(1, 2) = alpha(1) - M_zI * alpha(2)
        tg%spin_matrix(2, 1) = alpha(1) + M_zI * alpha(2)
        
      else
        message(1) = '"OCTTargetSpin" has to be specified as block.'
        call messages_info(1)
        call input_error('OCTTargetSpin')
      end if

    else
      message(1) = 'Error: if "OCTTargetOperator = oct_tg_spin", then you must'
      message(2) = 'supply one "OCTTargetSpin" block.'
      call messages_info(2)
      call input_error('OCTTargetSpin')
    end if


    POP_SUB(target_init_spin)
  end subroutine target_init_spin



  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_spin(tg, gr, psi) result(j1)
    type(target_t),   intent(inout) :: tg
    type(grid_t),     intent(inout) :: gr
    type(states_t),   intent(inout) :: psi
    
    integer :: i, j

    PUSH_SUB(target_j1_spin)

    j1 = M_ZERO
    do i = 1, 2
      do j = 1, 2
        j1 = j1 + tg%spin_matrix(i,j) * zmf_dotp(gr%mesh, psi%zpsi(:, i, 1, 1), psi%zpsi(:, j, 1, 1))
      end do
    end do
    

    POP_SUB(target_j1_spin)
  end function target_j1_spin



  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_spin(tg, gr, psi_in, chi_out)
    type(target_t),    intent(inout) :: tg
    type(grid_t),      intent(inout) :: gr
    type(states_t),    intent(inout) :: psi_in
    type(states_t),    intent(inout) :: chi_out
    
    integer :: i, j

    PUSH_SUB(target_chi_spin)
    
    chi_out%zpsi(:, :, 1, 1) = M_ZERO
    do i = 1, 2
      do j = 1, 2
        chi_out%zpsi(:, i, 1, 1) = chi_out%zpsi(:, i, 1, 1) + tg%spin_matrix(i,j) * psi_in%zpsi(:, j, 1, 1)
      end do
    end do

    POP_SUB(target_chi_spin)
  end subroutine target_chi_spin


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

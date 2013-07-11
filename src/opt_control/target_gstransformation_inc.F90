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
!! $Id: target_gstransformation_inc.F90 $


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_gstransformation(gr, tg)
    type(grid_t),     intent(in)    :: gr
    type(target_t),   intent(inout) :: tg

    type(block_t)       :: blk
    integer :: ist, jst
    CMPLX, allocatable  :: rotation_matrix(:, :)
    type(states_t)      :: tmp_st
    PUSH_SUB(target_init_gstransformation)

    message(1) =  'Info: Using Superposition of States for TargetOperator'
    call messages_info(1)

    !%Variable OCTTargetTransformStates
    !%Type block
    !%Default no
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% If <tt>OCTTargetOperator = oct_tg_gstransformation</tt>, you must specify a
    !% <tt>OCTTargetTransformStates</tt> block, in order to specify which linear
    !% combination of the states present in <tt>restart/gs</tt> is used to
    !% create the target state.
    !% 
    !% The syntax is the same as the <tt>TransformStates</tt> block.
    !%End
    if(parse_isdef(datasets_check('OCTTargetTransformStates')) /= 0) then
      if(parse_block(datasets_check('OCTTargetTransformStates'), blk) == 0) then
        call states_copy(tmp_st, tg%st)
        call states_deallocate_wfns(tmp_st)
        call restart_look_and_read(tmp_st, gr)
        SAFE_ALLOCATE(rotation_matrix(1:tg%st%nst, 1:tmp_st%nst))
        rotation_matrix = M_z0
        do ist = 1, tg%st%nst
          do jst = 1, parse_block_cols(blk, ist - 1)
            call parse_block_cmplx(blk, ist - 1, jst - 1, rotation_matrix(ist, jst))
          end do
        end do

        call states_rotate(gr%mesh, tg%st, tmp_st, rotation_matrix)
        SAFE_DEALLOCATE_A(rotation_matrix)
        call states_end(tmp_st)
        call parse_block_end(blk)
        call density_calc(tg%st, gr, tg%st%rho)
      else
        message(1) = '"OCTTargetTransformStates" has to be specified as block.'
        call messages_info(1)
        call input_error('OCTTargetTransformStates')
      end if
    else
      message(1) = 'Error: if "OCTTargetOperator = oct_tg_superposition", then you must'
      message(2) = 'supply one "OCTTargetTransformStates" block to create the superposition.'
      call messages_info(2)
      call input_error('OCTTargetTransformStates')
    end if

    POP_SUB(target_init_gstransformation)
  end subroutine target_init_gstransformation


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_gstransformation()
    PUSH_SUB(target_init_gstransformation)


    POP_SUB(target_init_gstransformation)
  end subroutine target_end_gstransformation


  ! ----------------------------------------------------------------------
  subroutine target_output_gstransformation(tg, gr, dir, geo, outp)
    type(target_t), intent(inout) :: tg
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(geometry_t),       intent(in)  :: geo
    type(output_t),         intent(in)  :: outp
    PUSH_SUB(target_output_gstransformation)
    
    call loct_mkdir(trim(dir))
    call output_states(tg%st, gr, geo, trim(dir), outp)

    POP_SUB(target_output_gstransformation)
  end subroutine target_output_gstransformation
  ! ----------------------------------------------------------------------



  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_gstransformation(tg, gr, psi) result(j1)
    type(target_t),   intent(inout) :: tg
    type(grid_t),     intent(inout) :: gr
    type(states_t),   intent(inout) :: psi

    integer :: ik, ist
    PUSH_SUB(target_j1_gstransformation)

    j1 = M_ZERO
    do ik = 1, psi%d%nik
      do ist = psi%st_start, psi%st_end
        j1 = j1 + psi%occ(ist, ik) * &
          abs(zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, ist, ik), &
              tg%st%zpsi(:, :, ist, ik)))**2
      end do
    end do

    POP_SUB(target_j1_gstransformation)
  end function target_j1_gstransformation


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_gstransformation(tg, gr, psi_in, chi_out)
    type(target_t),    intent(inout) :: tg
    type(grid_t),      intent(inout) :: gr
    type(states_t),    intent(inout) :: psi_in
    type(states_t),    intent(inout) :: chi_out

    integer :: ik, ist
    CMPLX :: olap
    PUSH_SUB(target_chi_gstransformation)

    do ik = 1, psi_in%d%nik
      do ist = psi_in%st_start, psi_in%st_end
        olap = zmf_dotp(gr%mesh, tg%st%zpsi(:, 1, ist, ik), psi_in%zpsi(:, 1, ist, ik))
        chi_out%zpsi(:, :, ist, ik) = olap * tg%st%zpsi(:, :, ist, ik)
      end do
    end do

    POP_SUB(target_chi_gstransformation)
  end subroutine target_chi_gstransformation

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

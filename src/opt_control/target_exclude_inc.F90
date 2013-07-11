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
!! $Id: target_exclude_inc.F90 $


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_exclude(gr, tg)
    type(grid_t),     intent(in)    :: gr
    type(target_t),   intent(inout) :: tg

    PUSH_SUB(target_init_exclude)

    message(1) =  'Info: The target functional is the exclusion of a number of states defined by'
    message(2) =  '      "OCTExcludedStates".'
    call messages_info(2)
    !%Variable OCTExcludedStates
    !%Type string
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% If the target is the exclusion of several targets, ("OCTTargetOperator = oct_exclude_states") 
    !% then you must declare which states are to be excluded, by setting the OCTExcludedStates variable.
    !% It must be a string in "list" format: "1-8", or "2,3,4-9", for example. Be careful to include
    !% in this list only states that have been calculated in a previous "gs" or "unocc" calculation,
    !% or otherwise the error will be silently ignored.
    !%End
    call parse_string(datasets_check('OCTExcludedStates'), "1", tg%excluded_states_list)
    call states_deallocate_wfns(tg%st)
    call restart_look_and_read(tg%st, gr)

    POP_SUB(target_init_exclude)
  end subroutine target_init_exclude


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_exclude()
    PUSH_SUB(target_end_exclude)

    POP_SUB(target_end_exclude)
  end subroutine target_end_exclude


  ! ----------------------------------------------------------------------
  subroutine target_output_exclude(tg, gr, dir, geo, outp)
    type(target_t), intent(inout) :: tg
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(geometry_t),       intent(in)  :: geo
    type(output_t),         intent(in)  :: outp

    PUSH_SUB(target_output_exclude)
    
    call loct_mkdir(trim(dir))
    call output_states(tg%st, gr, geo, trim(dir), outp)

    POP_SUB(target_output_exclude)
  end subroutine target_output_exclude
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_exclude(gr, tg, psi) result(j1)
    type(grid_t),     intent(inout) :: gr
    type(target_t),   intent(inout) :: tg
    type(states_t),   intent(inout) :: psi

    integer :: ist
    PUSH_SUB(target_j1_local)

    j1 = M_ONE
    do ist = 1, tg%st%nst
      if(loct_isinstringlist(ist, tg%excluded_states_list)) then
        j1 = j1 - abs(zmf_dotp(gr%mesh, psi%d%dim, &
          tg%st%zpsi(:, :, ist, 1), psi%zpsi(:, :, 1, 1)))**2
      end if
    end do

    POP_SUB(target_j1_exclude)
  end function target_j1_exclude


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_exclude(tg, gr, psi_in, chi_out)
    type(target_t),    intent(inout) :: tg
    type(grid_t),      intent(inout) :: gr
    type(states_t),    intent(inout) :: psi_in
    type(states_t),    intent(inout) :: chi_out

    integer :: ist
    CMPLX :: olap
    PUSH_SUB(target_chi_exclude)

    chi_out%zpsi(:, :, 1, 1) = psi_in%zpsi(:, :, 1, 1)
    do ist = 1, tg%st%nst
      if(loct_isinstringlist(ist, tg%excluded_states_list)) then
        olap = zmf_dotp(gr%mesh, psi_in%d%dim, tg%st%zpsi(:, :, ist, 1), psi_in%zpsi(:, :, 1, 1))
        chi_out%zpsi(:, :, 1, 1) = chi_out%zpsi(:, :, 1, 1) - olap * tg%st%zpsi(:, :, ist, 1)
      end if
    end do

    POP_SUB(target_chi_exclude)
  end subroutine target_chi_exclude

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

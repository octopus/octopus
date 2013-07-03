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

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

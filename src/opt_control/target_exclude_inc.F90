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
  subroutine target_init_exclude(mesh, namespace, space, tg, td, restart, kpoints)
    type(mesh_t),      intent(in)    :: mesh
    type(namespace_t), intent(in)    :: namespace
    type(space_t),     intent(in)    :: space
    type(target_t),    intent(inout) :: tg
    type(td_t),        intent(in)    :: td
    type(restart_t),   intent(inout) :: restart
    type(kpoints_t),   intent(in)    :: kpoints

    PUSH_SUB(target_init_exclude)

    tg%move_ions = ion_dynamics_ions_move(td%ions_dyn)
    tg%dt = td%dt

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
    call parse_variable(namespace, 'OCTExcludedStates', '1', tg%excluded_states_list)
    call states_elec_deallocate_wfns(tg%st)

    call states_elec_look_and_load(restart, namespace, space, tg%st, mesh, kpoints)

    POP_SUB(target_init_exclude)
  end subroutine target_init_exclude


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_exclude()
    PUSH_SUB(target_end_exclude)

    POP_SUB(target_end_exclude)
  end subroutine target_end_exclude


  ! ----------------------------------------------------------------------
  subroutine target_output_exclude(tg, namespace, space, gr, dir, ions, hm, outp)
    type(target_t),      intent(in) :: tg
    type(namespace_t),   intent(in) :: namespace
    type(space_t),       intent(in) :: space
    type(grid_t),        intent(in) :: gr
    character(len=*),    intent(in) :: dir
    type(ions_t),        intent(in) :: ions
    type(hamiltonian_elec_t), intent(in) :: hm
    type(output_t),      intent(in) :: outp

    PUSH_SUB(target_output_exclude)
    
    call io_mkdir(trim(dir), namespace)
    call output_states(outp, namespace, space, trim(dir), tg%st, gr, ions, hm, -1)

    POP_SUB(target_output_exclude)
  end subroutine target_output_exclude
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_exclude(gr, tg, psi) result(j1)
    type(grid_t),   intent(in) :: gr
    type(target_t), intent(in) :: tg
    type(states_elec_t), intent(in) :: psi

    integer :: ist
    CMPLX, allocatable :: zpsi1(:, :), zpsi(:, :)
    
    PUSH_SUB(target_j1_local)

    SAFE_ALLOCATE(zpsi(1:gr%mesh%np, 1:tg%st%d%dim))
    SAFE_ALLOCATE(zpsi1(1:gr%mesh%np, 1:tg%st%d%dim))

    call states_elec_get_state(psi, gr%mesh, 1, 1, zpsi1)
    
    j1 = M_ONE
    do ist = 1, tg%st%nst
      if(loct_isinstringlist(ist, tg%excluded_states_list)) then
        call states_elec_get_state(tg%st, gr%mesh, ist, 1, zpsi)
        j1 = j1 - abs(zmf_dotp(gr%mesh, psi%d%dim, zpsi, zpsi1))**2
      end if
    end do

    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(zpsi1)

    POP_SUB(target_j1_exclude)
  end function target_j1_exclude


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_exclude(tg, gr, psi_in, chi_out)
    type(target_t),    intent(in)    :: tg
    type(grid_t),      intent(in)    :: gr
    type(states_elec_t),    intent(in)    :: psi_in
    type(states_elec_t),    intent(inout) :: chi_out

    integer :: ist, ib
    CMPLX :: olap
    CMPLX, allocatable :: zpsi(:, :), zst(:, :), zchi(:, :)
    PUSH_SUB(target_chi_exclude)

    do ib = chi_out%group%block_start, chi_out%group%block_end 
      call psi_in%group%psib(ib, 1)%copy_data_to(gr%mesh%np, chi_out%group%psib(ib, 1))
    end do

    SAFE_ALLOCATE(zpsi(1:gr%mesh%np, 1:tg%st%d%dim))
    SAFE_ALLOCATE(zst(1:gr%mesh%np, 1:tg%st%d%dim))
    SAFE_ALLOCATE(zchi(1:gr%mesh%np, 1:tg%st%d%dim))

    call states_elec_get_state(chi_out, gr%mesh, 1, 1, zchi)
    
    do ist = 1, tg%st%nst
      if(loct_isinstringlist(ist, tg%excluded_states_list)) then
        call states_elec_get_state(psi_in, gr%mesh, ist, 1, zpsi)
        call states_elec_get_state(tg%st, gr%mesh, ist, 1, zst)
        olap = zmf_dotp(gr%mesh, psi_in%d%dim, zst, zpsi)
        zchi(1:gr%mesh%np, 1:tg%st%d%dim) = zchi(1:gr%mesh%np, 1:tg%st%d%dim) - olap*zst(1:gr%mesh%np, 1:tg%st%d%dim)
      end if
    end do

    call states_elec_set_state(chi_out, gr%mesh, 1, 1, zchi)

    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(zst)
    SAFE_DEALLOCATE_A(zchi)
    
    POP_SUB(target_chi_exclude)
  end subroutine target_chi_exclude

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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
  subroutine target_init_gstransformation(gr, namespace, space, tg, td, restart, kpoints)
    type(grid_t),      intent(in)    :: gr
    type(namespace_t), intent(in)    :: namespace
    type(space_t),     intent(in)    :: space
    type(target_t),    intent(inout) :: tg
    type(td_t),        intent(in)    :: td
    type(restart_t),   intent(inout) :: restart
    type(kpoints_t),   intent(in)    :: kpoints

    PUSH_SUB(target_init_gstransformation)

    message(1) =  'Info: Using Superposition of States for TargetOperator'
    call messages_info(1)

    tg%move_ions = ion_dynamics_ions_move(td%ions_dyn)
    tg%dt = td%dt

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
    call states_elec_transform(tg%st, namespace, space, restart, gr%mesh, kpoints, prefix = "OCTTarget")

    if(.not. parse_is_defined(namespace, 'OCTTargetTransformStates')) then
      message(1) = 'If "OCTTargetOperator = oct_tg_superposition", then you must'
      message(2) = 'supply an "OCTTargetTransformStates" block to create the superposition.'
      call messages_fatal(2)
    end if
    call density_calc(tg%st, gr, tg%st%rho)
    
    POP_SUB(target_init_gstransformation)
  end subroutine target_init_gstransformation


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_gstransformation()
    PUSH_SUB(target_init_gstransformation)


    POP_SUB(target_init_gstransformation)
  end subroutine target_end_gstransformation


  ! ----------------------------------------------------------------------
  subroutine target_output_gstransformation(tg, namespace, space, gr, dir, ions, hm, outp)
    type(target_t),      intent(in) :: tg
    type(namespace_t),   intent(in) :: namespace
    type(space_t),       intent(in) :: space
    type(grid_t),        intent(in) :: gr
    character(len=*),    intent(in) :: dir
    type(ions_t),        intent(in) :: ions
    type(hamiltonian_elec_t), intent(in) :: hm
    type(output_t),      intent(in) :: outp
    PUSH_SUB(target_output_gstransformation)
    
    call io_mkdir(trim(dir), namespace)
    call output_states(outp, namespace, space, trim(dir), tg%st, gr, ions, hm, -1)

    POP_SUB(target_output_gstransformation)
  end subroutine target_output_gstransformation
  ! ----------------------------------------------------------------------



  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_gstransformation(tg, gr, psi) result(j1)
    type(target_t),      intent(in) :: tg
    type(grid_t),        intent(in) :: gr
    type(states_elec_t), intent(in) :: psi

    integer :: ik, ist

    CMPLX, allocatable :: zpsi(:, :), zst(:, :)
    
    PUSH_SUB(target_j1_gstransformation)

    SAFE_ALLOCATE(zpsi(1:gr%mesh%np, 1:tg%st%d%dim))
    SAFE_ALLOCATE(zst(1:gr%mesh%np, 1:tg%st%d%dim))
    
    j1 = M_ZERO
    do ik = 1, psi%d%nik
      do ist = psi%st_start, psi%st_end

        call states_elec_get_state(psi, gr%mesh, ist, ik, zpsi)
        call states_elec_get_state(tg%st, gr%mesh, ist, ik, zst)
        
        j1 = j1 + abs(zmf_dotp(gr%mesh, psi%d%dim, zpsi, zst))**2
        
      end do
    end do

    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(zst)
    
    POP_SUB(target_j1_gstransformation)
  end function target_j1_gstransformation


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_gstransformation(tg, gr, psi_in, chi_out)
    type(target_t),      intent(in)    :: tg
    type(grid_t),        intent(in)    :: gr
    type(states_elec_t), intent(in)    :: psi_in
    type(states_elec_t), intent(inout) :: chi_out

    integer :: ik, ist
    CMPLX :: olap
    CMPLX, allocatable :: zpsi(:, :), zst(:, :), zchi(:, :)
    
    PUSH_SUB(target_chi_gstransformation)

    SAFE_ALLOCATE(zpsi(1:gr%mesh%np, 1:tg%st%d%dim))
    SAFE_ALLOCATE(zst(1:gr%mesh%np, 1:tg%st%d%dim))
    SAFE_ALLOCATE(zchi(1:gr%mesh%np, 1:tg%st%d%dim))
    
    do ik = 1, psi_in%d%nik
      do ist = psi_in%st_start, psi_in%st_end

        call states_elec_get_state(psi_in, gr%mesh, ist, ik, zpsi)
        call states_elec_get_state(tg%st, gr%mesh, ist, ik, zst)
        
        olap = zmf_dotp(gr%mesh, zst(:, 1), zpsi(:, 1))
        zchi(1:gr%mesh%np, 1:tg%st%d%dim) = olap*zst(1:gr%mesh%np, 1:tg%st%d%dim)
        
        call states_elec_set_state(chi_out, gr%mesh, ist, ik, zchi)

      end do
    end do

    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(zst)
    SAFE_DEALLOCATE_A(zchi)
    
    POP_SUB(target_chi_gstransformation)
  end subroutine target_chi_gstransformation

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

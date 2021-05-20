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
  subroutine target_init_local(gr, namespace, tg, td)
    type(grid_t),      intent(in)    :: gr
    type(namespace_t), intent(in)    :: namespace
    type(target_t),    intent(inout) :: tg
    type(td_t),        intent(in)    :: td

    integer             :: ip
    FLOAT               :: xx(1:gr%sb%dim), rr, psi_re, psi_im
    character(len=1024) :: expression
    PUSH_SUB(target_init_local)

    tg%move_ions = ion_dynamics_ions_move(td%ions_dyn)
    tg%dt = td%dt

    !%Variable OCTLocalTarget
    !%Type string
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% If <tt>OCTTargetOperator = oct_tg_local</tt>, then one must supply a function
    !% that defines the target. This should be done by defining it through a string, using 
    !% the variable <tt>OCTLocalTarget</tt>.
    !%End
    if(parse_is_defined(namespace, 'OCTLocalTarget')) then
      SAFE_ALLOCATE(tg%rho(1:gr%mesh%np))
      tg%rho = M_ZERO
      call parse_variable(namespace, 'OCTLocalTarget', "0", expression)
      call conv_to_C_string(expression)
      do ip = 1, gr%mesh%np
        call mesh_r(gr%mesh, ip, rr, coords = xx)
        ! parse user-defined expression
        call parse_expression(psi_re, psi_im, gr%sb%dim, xx, rr, M_ZERO, expression)
        tg%rho(ip) = psi_re
      end do
    else
      message(1) = 'If OCTTargetOperator = oct_tg_local, then you must give the shape'
      message(2) = 'of this target in variable "OCTLocalTarget".'
      call messages_fatal(2)
    end if

    POP_SUB(target_init_local)
  end subroutine target_init_local


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_local(tg)
    type(target_t),   intent(inout) :: tg

    PUSH_SUB(target_end_local)

    SAFE_DEALLOCATE_A(tg%rho)

    POP_SUB(target_end_local)
  end subroutine target_end_local


  ! ----------------------------------------------------------------------
  subroutine target_output_local(tg, namespace, space, mesh, dir, ions, outp)
    type(target_t),    intent(in) :: tg
    type(namespace_t), intent(in) :: namespace
    type(space_t),     intent(in) :: space
    type(mesh_t),      intent(in) :: mesh
    character(len=*),  intent(in) :: dir
    type(ions_t),      intent(in) :: ions
    type(output_t),    intent(in) :: outp

    integer :: ierr
    PUSH_SUB(target_output_local)
    
    call io_mkdir(trim(dir), namespace)
    call dio_function_output(outp%how(0), trim(dir), 'local_target', namespace, space, mesh, &
      tg%rho, units_out%length**(-space%dim), ierr, ions = ions)


    POP_SUB(target_output_local)
  end subroutine target_output_local
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_local(mesh, tg, psi) result(j1)
    type(mesh_t),        intent(in) :: mesh
    type(target_t),      intent(in) :: tg
    type(states_elec_t), intent(in) :: psi

    integer :: is
    PUSH_SUB(target_j1_local)

    j1 = M_ZERO
    do is = 1, psi%d%spin_channels
      j1 = j1 + dmf_dotp(mesh, tg%rho, psi%rho(:, is))
    end do

    POP_SUB(target_j1_local)
  end function target_j1_local


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_local(tg, mesh, psi_in, chi_out)
    type(target_t),      intent(in)    :: tg
    type(mesh_t),        intent(in)    :: mesh
    type(states_elec_t), intent(in)    :: psi_in
    type(states_elec_t), intent(inout) :: chi_out

    integer :: ik, idim, ist, ip
    CMPLX, allocatable :: zpsi(:, :)
    
    PUSH_SUB(target_chi_local)

    SAFE_ALLOCATE(zpsi(1:mesh%np, 1:psi_in%d%dim))
    
    do ik = 1, psi_in%d%nik
      do idim = 1, psi_in%d%dim
        do ist = psi_in%st_start, psi_in%st_end
          call states_elec_get_state(psi_in, mesh, ist, ik, zpsi)
          do ip = 1, mesh%np
            zpsi(ip, idim) = psi_in%occ(ist, ik)*tg%rho(ip)*zpsi(ip, idim)
          end do
          call states_elec_set_state(chi_out, mesh, ist, ik, zpsi)
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(zpsi)

    POP_SUB(target_chi_local)
  end subroutine target_chi_local



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

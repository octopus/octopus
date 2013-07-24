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
!! $Id: opt_control.F90 2875 2007-04-30 16:54:15Z acastro $


  ! ---------------------------------------------------------
  subroutine oct_finalcheck(sys, hm, td)
    type(system_t), intent(inout)      :: sys
    type(hamiltonian_t), intent(inout) :: hm
    type(td_t), intent(inout)          :: td

    type(states_t) :: psi
    type(opt_control_state_t) :: qcpsi    
    FLOAT :: j1, jfunctional, fluence, j2

    type(controlfunction_t), pointer :: par

    PUSH_SUB(oct_finalcheck)

    if(.not.oct%oct_double_check) then
      POP_SUB(oct_finalcheck)
      return
    end if

    call oct_iterator_bestpar(par, iterator)

    call opt_control_state_copy(qcpsi, initial_st)
    call propagate_forward(sys, hm, td, par, oct_target, qcpsi, write_iter = .true.)
    call opt_control_get_qs(psi, qcpsi)

    j1 = target_j1(oct_target, sys%gr, qcpsi)
    call opt_control_state_end(qcpsi)

    fluence = controlfunction_fluence(par)
    j2 = controlfunction_j2(par)
    jfunctional = j1 + j2

    write(message(1), '(a)') 'Final propagation with the best field'
    call messages_print_stress(stdout, trim(message(1)))
    write(message(1), '(6x,a,f12.5)') " => J1       = ", j1
    write(message(2), '(6x,a,f12.5)') " => J        = ", jfunctional
    write(message(3), '(6x,a,f12.5)') " => J2       = ", j2
    write(message(4), '(6x,a,f12.5)') " => Fluence  = ", fluence
    call messages_info(4)
    call messages_print_stress(stdout)

    call output_states(psi, sys%gr, sys%geo, OCT_DIR//'final', sys%outp)

    nullify(par)
    call states_end(psi)
    POP_SUB(oct_finalcheck)
  end subroutine oct_finalcheck

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

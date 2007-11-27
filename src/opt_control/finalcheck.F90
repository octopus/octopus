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
!! $Id: opt_control.F90 2875 2007-04-30 16:54:15Z acastro $


  ! ---------------------------------------------------------
  subroutine oct_finalcheck(oct, initial_st, target, par, sys, h, td)
    type(oct_t), intent(in)            :: oct
    type(states_t), intent(in)         :: initial_st
    type(target_t), intent(inout)      :: target
    type(system_t), intent(inout)      :: sys
    type(oct_control_parameters_t), intent(inout) :: par
    type(hamiltonian_t), intent(inout) :: h
    type(td_t), intent(inout)          :: td

    type(states_t) :: psi
    FLOAT :: j1, jfunctional, fluence, j2

    call push_sub('finalcheck.oct_finalcheck')

    if(.not.oct%oct_double_check) then
      call pop_sub(); return
    end if

    psi = initial_st
    
    call parameters_to_h(par, h%ep)
    call propagate_forward(sys, h, td, target, psi, write_iter = .true.)

    j1 = j1_functional(sys%gr, sys%geo, h%ep, psi, target)
    fluence = laser_fluence(par)
    j2 = j2_functional(par)
    jfunctional = j1 + j2

    write(message(1), '(a)') 'Final propagation with the best field'
    call messages_print_stress(stdout, trim(message(1)))
    write(message(1), '(6x,a,f12.5)') " => J1       = ", j1
    write(message(2), '(6x,a,f12.5)') " => J        = ", jfunctional
    write(message(3), '(6x,a,f12.5)') " => J2       = ", j2
    write(message(4), '(6x,a,f12.5)') " => Fluence  = ", fluence
    call write_info(4)
    call messages_print_stress(stdout)

    call states_end(psi)
    call pop_sub()
  end subroutine oct_finalcheck

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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


  ! ---------------------------------------------------------
  !> Propagator specifically designed for the QOCT+TDDFT problem
  subroutine td_qoct_tddft_propagator(hm, xc, gr, st, tr, t, dt, ions, geo)
    type(hamiltonian_t), intent(inout) :: hm
    type(xc_t),          intent(in)    :: xc
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
    type(propagator_t),  intent(inout) :: tr
    FLOAT,               intent(in)    :: t, dt
    type(ion_dynamics_t),            intent(inout) :: ions
    type(geometry_t),                intent(inout) :: geo

    type(ion_state_t) :: ions_state
    PUSH_SUB(td_qoct_tddft_propagator)
    
    if( (hm%theory_level /= INDEPENDENT_PARTICLES) .and. &
        (.not.hamiltonian_oct_exchange(hm)) ) then
      call vksinterp_interpolate(tr%vksold, 2, .false., gr%mesh%np, st%d%nspin, hm%vhxc, hm%imvhxc, t, dt, t-dt/M_TWO)
    end if

    !move the ions to time 'time - dt/2'
    if(ion_dynamics_ions_move(ions)) then
      call ion_dynamics_save_state(ions, geo, ions_state)
      call ion_dynamics_propagate(ions, gr%sb, geo, t - dt/M_TWO, M_HALF*dt)
      call hamiltonian_epot_generate(hm, gr, geo, st, time = t - dt/M_TWO)
    end if

    call hamiltonian_update(hm, gr%mesh, time = t-dt/M_TWO)
    call exponential_apply_all(tr%te, gr%der, hm, xc, st, dt, t - dt/M_TWO)

    !restore to time 'time - dt'
    if(ion_dynamics_ions_move(ions)) call ion_dynamics_restore_state(ions, geo, ions_state)

    POP_SUB(td_qoct_tddft_propagator)
  end subroutine td_qoct_tddft_propagator
  ! ---------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

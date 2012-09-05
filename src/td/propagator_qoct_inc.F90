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
!! $Id: propagator.F90 5406 2009-05-16 18:17:47Z xavier $


  ! ---------------------------------------------------------
  !> Propagator specifically designed for the QOCT+TDDFT problem
  subroutine td_qoct_tddft_propagator(hm, gr, st, tr, t, dt)!, gauge_force, ions, geo)
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
    type(propagator_t),  intent(inout) :: tr
    FLOAT,               intent(in)    :: t, dt

    PUSH_SUB(td_qoct_tddft_propagator)
    
    if( (hm%theory_level .ne. INDEPENDENT_PARTICLES) .and. &
        (.not.hamiltonian_oct_exchange(hm)) ) then
      call interpolate( (/t, t-dt/), tr%v_old(:, :, 0:1), t-dt/M_TWO, hm%vhxc(:, :))
    end if

    call hamiltonian_update(hm, gr%mesh, time = t-dt/M_TWO)
    call exponential_apply_all(tr%te, gr%der, hm, st, dt, t - dt/M_TWO)

    POP_SUB(td_qoct_tddft_propagator)
  end subroutine td_qoct_tddft_propagator
  ! ---------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

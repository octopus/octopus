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

#include "global.h"

module propagator_qoct_oct_m
  use density_oct_m
  use exponential_oct_m
  use grid_oct_m
  use geometry_oct_m
  use global_oct_m
  use hamiltonian_oct_m
  use ion_dynamics_oct_m
  use lda_u_oct_m
  use messages_oct_m
  use oct_exchange_oct_m
  use potential_interpolation_oct_m
  use propagator_base_oct_m
  use states_oct_m
  use xc_oct_m

  implicit none

  private

  public ::                    &
    td_qoct_tddft_propagator

contains

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
      (.not. oct_exchange_enabled(hm%oct_exchange)) ) then
      !TODO: This does not support complex scaling
      if(hm%family_is_mgga_with_exc) then
        call potential_interpolation_interpolate(tr%vksold, 2, t, dt, t-dt/M_TWO, hm%vhxc, vtau = hm%vtau)
      else
        call potential_interpolation_interpolate(tr%vksold, 2, t, dt, t-dt/M_TWO, hm%vhxc)
      end if
    end if

    !move the ions to time 'time - dt/2'
    if(ion_dynamics_ions_move(ions)) then
      call ion_dynamics_save_state(ions, geo, ions_state)
      call ion_dynamics_propagate(ions, gr%sb, geo, t - dt/M_TWO, M_HALF*dt)
      call hamiltonian_epot_generate(hm, gr, geo, st, time = t - dt/M_TWO)
    end if

    call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = t-dt/M_TWO)
    call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy )
    call exponential_apply_all(tr%te, gr%der, hm, xc, st, dt)

    call density_calc(st, gr, st%rho)

    !restore to time 'time - dt'
    if(ion_dynamics_ions_move(ions)) then
      call ion_dynamics_restore_state(ions, geo, ions_state)
    end if

    POP_SUB(td_qoct_tddft_propagator)
  end subroutine td_qoct_tddft_propagator
  ! ---------------------------------------------------------

end module propagator_qoct_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

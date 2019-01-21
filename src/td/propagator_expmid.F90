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

module propagator_expmid_oct_m
  use density_oct_m
  use exponential_oct_m
  use gauge_field_oct_m
  use grid_oct_m
  use geometry_oct_m
  use global_oct_m
  use hamiltonian_oct_m
  use ion_dynamics_oct_m
  use lda_u_oct_m
  use math_oct_m
  use messages_oct_m
  use potential_interpolation_oct_m
  use propagator_base_oct_m
  use states_oct_m
  use xc_oct_m

  implicit none

  private

  public ::                       &
    exponential_midpoint

contains
  
  ! ---------------------------------------------------------
  !> Exponential midpoint
  subroutine exponential_midpoint(hm, gr, st, tr, time, dt, ionic_scale, ions, geo, move_ions)
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    type(propagator_t),  target,     intent(inout) :: tr
    FLOAT,                           intent(in)    :: time
    FLOAT,                           intent(in)    :: dt
    FLOAT,                           intent(in)    :: ionic_scale
    type(ion_dynamics_t),            intent(inout) :: ions
    type(geometry_t),                intent(inout) :: geo
    logical,                         intent(in)    :: move_ions

    integer :: ib, ik
    type(ion_state_t) :: ions_state
    FLOAT :: vecpot(1:MAX_DIM), vecpot_vel(1:MAX_DIM)
    CMPLX :: zt, zdt

    PUSH_SUB(propagator_dt.exponential_midpoint)

    ! the half step of this propagator screws with the gauge field kick
    ASSERT(hm%ep%gfield%with_gauge_field .eqv. .false.)

    vecpot(:)     = M_ZERO
    vecpot_vel(:) = M_ZERO
    
    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
        call potential_interpolation_interpolate(tr%vksold, 3, &
          time, dt, time - dt/M_TWO, hm%vhxc)
    end if

    !move the ions to time 'time - dt/2'
    if(move_ions .and.  ion_dynamics_ions_move(ions)) then
      call ion_dynamics_save_state(ions, geo, ions_state)
      call ion_dynamics_propagate(ions, gr%sb, geo, time - dt/M_TWO, ionic_scale*CNST(0.5)*dt)
      call hamiltonian_epot_generate(hm, gr, geo, st, time = time - dt/M_TWO)
    end if

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_get_vec_pot(hm%ep%gfield, vecpot)
      call gauge_field_get_vec_pot_vel(hm%ep%gfield, vecpot_vel)
      call gauge_field_propagate(hm%ep%gfield, M_HALF*dt, time)
    end if
    call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = time - M_HALF*dt)
    !We update the occupation matrices
    call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy )
    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end

        call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, dt)

      end do
    end do

    !restore to time 'time - dt'
    if(move_ions .and. ion_dynamics_ions_move(ions)) then
      call ion_dynamics_restore_state(ions, geo, ions_state)
    end if

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_set_vec_pot(hm%ep%gfield, vecpot)
      call gauge_field_set_vec_pot_vel(hm%ep%gfield, vecpot_vel)
      call hamiltonian_update(hm, gr%mesh, gr%der%boundaries)
    end if

    call density_calc(st, gr, st%rho)

    POP_SUB(propagator_dt.exponential_midpoint)
  end subroutine exponential_midpoint
! ---------------------------------------------------------

end module propagator_expmid_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End: 

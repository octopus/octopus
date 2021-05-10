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
  use grid_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use ion_dynamics_oct_m
  use ions_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use potential_interpolation_oct_m
  use propagator_base_oct_m
  use propagation_ops_elec_oct_m
  use space_oct_m
  use states_elec_oct_m
  use xc_oct_m

  implicit none

  private

  public ::                       &
    exponential_midpoint

contains
  
  ! ---------------------------------------------------------
  !> Exponential midpoint
  subroutine exponential_midpoint(hm, namespace, space, gr, st, tr, time, dt, ionic_scale, ions_dyn, ions, move_ions)
    type(hamiltonian_elec_t), target, intent(inout) :: hm
    type(namespace_t),                intent(in)    :: namespace
    type(space_t),                    intent(in)    :: space
    type(grid_t),             target, intent(inout) :: gr
    type(states_elec_t),      target, intent(inout) :: st
    type(propagator_base_t),  target, intent(inout) :: tr
    FLOAT,                            intent(in)    :: time
    FLOAT,                            intent(in)    :: dt
    FLOAT,                            intent(in)    :: ionic_scale
    type(ion_dynamics_t),             intent(inout) :: ions_dyn
    type(ions_t),                     intent(inout) :: ions
    logical,                          intent(in)    :: move_ions

    PUSH_SUB(propagator_dt.exponential_midpoint)

    ! the half step of this propagator screws with the gauge field kick
    ASSERT(hm%ep%gfield%with_gauge_field .eqv. .false.)

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      if(family_is_mgga_with_exc(hm%xc)) then
        call potential_interpolation_interpolate(tr%vksold, 3, &
          time, dt, time - dt/M_TWO, hm%vhxc, vtau = hm%vtau)
      else
        call potential_interpolation_interpolate(tr%vksold, 3, &
          time, dt, time - dt/M_TWO, hm%vhxc)
      end if
    end if

    !move the ions to time 'time - dt/2'
    call propagation_ops_elec_move_ions(tr%propagation_ops_elec, gr, hm, st, namespace, space, ions_dyn, ions, &
            time - M_HALF*dt, ionic_scale*M_HALF*dt, save_pos = .true., move_ions = move_ions)

    call propagation_ops_elec_propagate_gauge_field(tr%propagation_ops_elec, namespace, hm, M_HALF*dt, time, save_gf = .true.)

    call propagation_ops_elec_update_hamiltonian(namespace, space, st, gr%mesh, hm, time - dt*M_HALF)

    call propagation_ops_elec_fuse_density_exp_apply(tr%te, namespace, st, gr, hm, dt)

    !restore to time 'time - dt'
    call propagation_ops_elec_restore_ions(tr%propagation_ops_elec, ions_dyn, ions, move_ions = move_ions)

    call propagation_ops_elec_restore_gauge_field(tr%propagation_ops_elec, namespace, space, hm, gr%mesh)

    POP_SUB(propagator_dt.exponential_midpoint)
  end subroutine exponential_midpoint
! ---------------------------------------------------------

end module propagator_expmid_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End: 

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

module propagator_magnus_oct_m
  use batch_oct_m
  use density_oct_m
  use exponential_oct_m
  use gauge_field_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use ion_dynamics_oct_m
  use ions_oct_m
  use lasers_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use potential_interpolation_oct_m
  use profiling_oct_m
  use propagation_ops_elec_oct_m
  use propagator_base_oct_m
  use propagator_rk_oct_m
  use space_oct_m
  use states_elec_oct_m
  use v_ks_oct_m
  use xc_oct_m

  implicit none

  private

  public ::                    &
    td_magnus,                 &
    td_cfmagnus4

contains
  
  ! ---------------------------------------------------------
  !> Magnus propagator
  subroutine td_magnus(hm, gr, st, tr, namespace, time, dt)
    type(hamiltonian_elec_t), target, intent(inout) :: hm
    type(grid_t),             target, intent(inout) :: gr
    type(states_elec_t),      target, intent(inout) :: st
    type(propagator_base_t),  target, intent(inout) :: tr
    type(namespace_t),                intent(in)    :: namespace
    FLOAT,                            intent(in)    :: time
    FLOAT,                            intent(in)    :: dt

    integer :: j, is, i
    FLOAT :: atime(2)
    FLOAT, allocatable :: vaux(:, :, :), pot(:)

    PUSH_SUB(propagator_dt.td_magnus)

    ASSERT(.not. family_is_mgga_with_exc(hm%xc))

    SAFE_ALLOCATE(vaux(1:gr%mesh%np, 1:st%d%nspin, 1:2))

    atime(1) = (M_HALF-sqrt(M_THREE)/CNST(6.0))*dt
    atime(2) = (M_HALF+sqrt(M_THREE)/CNST(6.0))*dt

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      do j = 1, 2
        call potential_interpolation_interpolate(tr%vksold, 3, time, dt, time - dt + atime(j), vaux(:, :, j))
      end do
    else
      vaux = M_ZERO
    end if

    do j = 1, 2
      ! WARNING: This should be carefully tested, and extended to allow for velocity-gauge laser fields.
      do i = 1, hm%ext_lasers%no_lasers
        select case(laser_kind(hm%ext_lasers%lasers(i)))
        case(E_FIELD_ELECTRIC)
          SAFE_ALLOCATE(pot(1:gr%mesh%np))
          pot = M_ZERO
          call laser_potential(hm%ext_lasers%lasers(i), gr%mesh, pot, time - dt + atime(j))
          do is = 1, st%d%nspin
            vaux(:, is, j) = vaux(:, is, j) + pot(:)
          end do
          SAFE_DEALLOCATE_A(pot)
        case(E_FIELD_MAGNETIC, E_FIELD_VECTOR_POTENTIAL)
          write(message(1),'(a)') 'The Magnus propagator cannot be used with magnetic fields, or'
          write(message(2),'(a)') 'with an electric field described in the velocity gauge.'
          call messages_fatal(2, namespace=namespace)
        end select
      end do
    end do

    tr%vmagnus(:, :, 2)  = M_HALF*(vaux(:, :, 1) + vaux(:, :, 2))
    tr%vmagnus(:, :, 1) = (sqrt(M_THREE)/CNST(12.0))*dt*(vaux(:, :, 2) - vaux(:, :, 1))

    call propagation_ops_elec_fuse_density_exp_apply(tr%te, namespace, st, gr, hm, dt, vmagnus = tr%vmagnus)

    SAFE_DEALLOCATE_A(vaux)
    POP_SUB(propagator_dt.td_magnus)
  end subroutine td_magnus


  ! ---------------------------------------------------------
  !> Commutator-free Magnus propagator of order 4.
  subroutine td_cfmagnus4(ks, namespace, space, hm, gr, st, tr, time, dt, ions_dyn, ions, iter)
    type(v_ks_t),             target, intent(inout) :: ks
    type(namespace_t),                intent(in)    :: namespace
    type(space_t),                    intent(in)    :: space
    type(hamiltonian_elec_t), target, intent(inout) :: hm
    type(grid_t),             target, intent(inout) :: gr
    type(states_elec_t),      target, intent(inout) :: st
    type(propagator_base_t),  target, intent(inout) :: tr
    FLOAT,                            intent(in)    :: time
    FLOAT,                            intent(in)    :: dt
    type(ion_dynamics_t),             intent(inout) :: ions_dyn
    type(ions_t),                     intent(inout) :: ions
    integer,                          intent(in)    :: iter

    FLOAT :: alpha1, alpha2, c1, c2, t1, t2
    FLOAT, allocatable :: vhxc1(:, :), vhxc2(:, :)

    if(ion_dynamics_ions_move(ions_dyn) .or. gauge_field_is_applied(hm%ep%gfield)) then
      message(1) = "The commutator-free Magnus expansion cannot be used with moving ions or gauge fields"
      call messages_fatal(1, namespace=namespace)
    end if

    PUSH_SUB(propagator_dt.td_cfmagnus4)

    if(iter < 4) then
      call td_explicit_runge_kutta4(ks, namespace, space, hm, gr, st, time, dt, ions_dyn, ions)
      POP_SUB(propagator_dt.td_cfmagnus4)
      return
    end if


    alpha1 = (M_THREE - M_TWO * sqrt(M_THREE))/CNST(12.0)
    alpha2 = (M_THREE + M_TWO * sqrt(M_THREE))/CNST(12.0)
    c1 = M_HALF - sqrt(M_THREE)/CNST(6.0)
    c2 = M_HALF + sqrt(M_THREE)/CNST(6.0)
    t1 = time - dt + c1*dt
    t2 = time - dt + c2*dt

    SAFE_ALLOCATE(vhxc1(1:gr%mesh%np, 1:st%d%nspin))
    SAFE_ALLOCATE(vhxc2(1:gr%mesh%np, 1:st%d%nspin))

    call potential_interpolation_interpolate(tr%vksold, 4, time, dt, t1, vhxc1)
    call potential_interpolation_interpolate(tr%vksold, 4, time, dt, t2, vhxc2)

    hm%vhxc = M_TWO * (alpha2 * vhxc1 + alpha1 * vhxc2)
    call hamiltonian_elec_update2(hm, gr%mesh, space, (/ t1, t2 /), (/ M_TWO * alpha2, M_TWO * alpha1/) )
    ! propagate by dt/2 
    call propagation_ops_elec_exp_apply(tr%te, namespace, st, gr%mesh, hm, M_HALF*dt)

    hm%vhxc = M_TWO * (alpha1 * vhxc1 + alpha2 * vhxc2)
    call hamiltonian_elec_update2(hm, gr%mesh, space, (/ t1, t2 /), (/ M_TWO * alpha1, M_TWO * alpha2/) )
    ! propagate by dt/2
    !TODO: fuse this with density calc
    call propagation_ops_elec_exp_apply(tr%te, namespace, st, gr%mesh, hm, M_HALF*dt)

    call density_calc(st, gr, st%rho)

    SAFE_DEALLOCATE_A(vhxc1)
    SAFE_DEALLOCATE_A(vhxc2)
    POP_SUB(propagator_dt.td_cfmagnus4)
  end subroutine td_cfmagnus4

end module propagator_magnus_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End: 

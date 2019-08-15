!! Copyright (C) 2019 N. Tancogne-Dejean
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

module propagation_ops_elec_oct_m
  use batch_oct_m
  use density_oct_m  
  use exponential_oct_m
  use geometry_oct_m
  use gauge_field_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use ion_dynamics_oct_m
  use lda_u_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use potential_interpolation_oct_m
  use profiling_oct_m
  use states_elec_oct_m
  use varinfo_oct_m
  use propagation_ops_abst_oct_m
  use xc_oct_m

  implicit none

  private
  public ::                             &
    propagation_ops_elec_t,                      &
    propagation_ops_elec_update_hamiltonian,     &
    propagation_ops_elec_exp_apply,              &
    propagation_ops_elec_fuse_density_exp_apply, &
    propagation_ops_elec_move_ions,              &
    propagation_ops_elec_restore_ions,           &
    propagation_ops_elec_propagate_gauge_field,  &
    propagation_ops_elec_restore_gauge_field,    &
    propagation_ops_elec_interpolate_get

  type, extends(propagation_ops_abst_t) :: propagation_ops_elec_t
    private

    type(ion_state_t) :: ions_state
    FLOAT :: vecpot(1:MAX_DIM), vecpot_vel(1:MAX_DIM)

  contains

    procedure :: init => propagation_ops_elec_init
    procedure :: end => propagation_ops_elec_end
  end type propagation_ops_elec_t


contains

  subroutine propagation_ops_elec_init(wo)
    class(propagation_ops_elec_t),  intent(inout) :: wo

    PUSH_SUB(propagation_ops_elec_init)

    POP_SUB(propagation_ops_elec_init)
  end subroutine

  subroutine propagation_ops_elec_end(wo)
    class(propagation_ops_elec_t),  intent(inout) :: wo

    PUSH_SUB(propagation_ops_elec_end)

    POP_SUB(propagation_ops_elec_end)
  end subroutine

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_update_hamiltonian(namespace, st, mesh, hm, time)
    type(namespace_t),        intent(in)    :: namespace
    type(states_elec_t),      intent(inout) :: st
    type(mesh_t),             intent(in)    :: mesh
    type(hamiltonian_elec_t), intent(inout) :: hm
    FLOAT,                    intent(in)    :: time

    type(profile_t), save :: prof

    PUSH_SUB(propagation_ops_elec_update_hamiltonian)

    call profiling_in(prof, 'ELEC_UPDATE_H')

    call hamiltonian_elec_update(hm, mesh, namespace, time = time)
    call lda_u_update_occ_matrices(hm%lda_u, mesh, st, hm%hm_base, hm%energy)

    call profiling_out(prof)

    POP_SUB(propagation_ops_elec_update_hamiltonian)

  end subroutine propagation_ops_elec_update_hamiltonian


  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_move_ions(wo, gr, hm, st, namespace, ions, geo, time, ion_time, save_pos, move_ions)
    class(propagation_ops_elec_t), intent(inout) :: wo
    type(grid_t),                  intent(in)    :: gr
    type(hamiltonian_elec_t),      intent(inout) :: hm
    type(states_elec_t),           intent(inout) :: st
    type(namespace_t),             intent(in)    :: namespace
    type(ion_dynamics_t),          intent(inout) :: ions
    type(geometry_t),              intent(inout) :: geo
    FLOAT,                         intent(in)    :: time
    FLOAT,                         intent(in)    :: ion_time
    logical, optional,             intent(in)    :: save_pos
    logical, optional,             intent(in)    :: move_ions

    type(profile_t), save :: prof

    PUSH_SUB(propagation_ops_elec_move_ions)    

    call profiling_in(prof, 'ELEC_MOVE_IONS')
   
    if (ion_dynamics_ions_move(ions) .and. optional_default(move_ions, .true.)) then
      if (optional_default(save_pos, .false.)) then
        call ion_dynamics_save_state(ions, geo, wo%ions_state)
      end if
      call ion_dynamics_propagate(ions, gr%sb, geo, time, ion_time)
      call hamiltonian_elec_epot_generate(hm, namespace, gr, geo, st, time = time)
    end if

    POP_SUB(propagation_ops_elec_move_ions)

    call profiling_out(prof)

  end subroutine propagation_ops_elec_move_ions

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_restore_ions(wo, ions, geo, move_ions)
    class(propagation_ops_elec_t),    intent(inout) :: wo
    type(ion_dynamics_t),    intent(inout) :: ions
    type(geometry_t),        intent(inout) :: geo
    logical, optional,       intent(in)    :: move_ions

    type(profile_t), save :: prof

    PUSH_SUB(propagation_ops_elec_restore_ions)    

    call profiling_in(prof, 'ELEC_RESTORE_IONS')

    if (ion_dynamics_ions_move(ions) .and. optional_default(move_ions, .true.)) then
      call ion_dynamics_restore_state(ions, geo, wo%ions_state)
    end if

    call profiling_out(prof)

    POP_SUB(propagation_ops_elec_retore_ions)

  end subroutine propagation_ops_elec_restore_ions

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_propagate_gauge_field(wo, hm, dt, time, save_gf)
    class(propagation_ops_elec_t), intent(inout) :: wo
    type(hamiltonian_elec_t),      intent(inout) :: hm
    FLOAT,                         intent(in)    :: dt
    FLOAT,                         intent(in)    :: time
    logical,  optional,            intent(in)    :: save_gf

    type(profile_t), save :: prof

    PUSH_SUB(propagation_ops_elec_propagate_gauge_field)

    call profiling_in(prof, 'ELEC_MOVE_GAUGE')

    if(gauge_field_is_applied(hm%ep%gfield)) then
      if(optional_default(save_gf, .false.)) then
        call gauge_field_get_vec_pot(hm%ep%gfield, wo%vecpot)
        call gauge_field_get_vec_pot_vel(hm%ep%gfield, wo%vecpot_vel)
      end if
      call gauge_field_propagate(hm%ep%gfield, dt, time)
    end if

    POP_SUB(propagation_ops_elec_propagate_gauge_field)

    call profiling_out(prof)

  end subroutine propagation_ops_elec_propagate_gauge_field

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_restore_gauge_field(wo, namespace, hm, mesh)
    class(propagation_ops_elec_t), intent(in)    :: wo
    type(namespace_t),             intent(in)    :: namespace
    type(hamiltonian_elec_t),      intent(inout) :: hm
    type(mesh_t),                  intent(in)    :: mesh

    type(profile_t), save :: prof

    PUSH_SUB(propagation_ops_elec_restore_gauge_field)

    call profiling_in(prof, 'ELEC_RESTORE_GAUGE')

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_set_vec_pot(hm%ep%gfield, wo%vecpot)
      call gauge_field_set_vec_pot_vel(hm%ep%gfield, wo%vecpot_vel)
      call hamiltonian_elec_update(hm, mesh, namespace)
    end if

    call profiling_out(prof)

    POP_SUB(propagation_ops_elec_retore_gauge_field)

  end subroutine propagation_ops_elec_restore_gauge_field

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_exp_apply(te, st, mesh, hm, dt)
    type(exponential_t),      intent(inout) :: te 
    type(states_elec_t),      intent(inout) :: st
    type(mesh_t),             intent(in)    :: mesh
    type(hamiltonian_elec_t), intent(inout) :: hm
    FLOAT,                    intent(in)    :: dt

    integer :: ik, ib
    type(profile_t), save :: prof

    PUSH_SUB(propagation_ops_elec_exp_apply)

    call profiling_in(prof, 'ELEC_EXP_APPLY')

    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call exponential_apply_batch(te, mesh, hm, st%group%psib(ib, ik), ik, dt)
      end do
    end do

    call profiling_out(prof)

    POP_SUB(propagation_ops_elec_exp_apply)

  end subroutine propagation_ops_elec_exp_apply

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_fuse_density_exp_apply(te, st, gr, hm, dt, dt2)
    type(exponential_t),      intent(inout) :: te
    type(states_elec_t),      intent(inout) :: st
    type(grid_t),             intent(inout) :: gr
    type(hamiltonian_elec_t), intent(inout) :: hm
    FLOAT,                    intent(in)    :: dt
    FLOAT, optional,          intent(in)    :: dt2

    integer :: ik, ib
    type(batch_t) :: zpsib_dt
    type(density_calc_t) :: dens_calc
    type(profile_t), save :: prof

    PUSH_SUB(propagation_ops_elec_fuse_density_exp_apply)

    call profiling_in(prof, 'ELEC_FUSE_DENS_EXP_APPLY')

    call density_calc_init(dens_calc, st, gr, st%rho)

    if(present(dt2)) then
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end

          call batch_copy(st%group%psib(ib, ik), zpsib_dt)
          if(batch_is_packed(st%group%psib(ib, ik))) call batch_pack(zpsib_dt, copy = .false.)

          !propagate the state dt/2 and dt, simultaneously, with H(time - dt)
          call exponential_apply_batch(te, gr%mesh, hm, st%group%psib(ib, ik), ik, dt, &
            psib2 = zpsib_dt, deltat2 = M_TWO*dt)

          !use the dt propagation to calculate the density
          call density_calc_accumulate(dens_calc, ik, zpsib_dt)

          call batch_end(zpsib_dt)

        end do
      end do

   
    else

      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call exponential_apply_batch(te, gr%mesh, hm, st%group%psib(ib, ik), ik, dt)
          call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik))
        end do
      end do

    end if

    call density_calc_end(dens_calc)

    call profiling_out(prof)

    POP_SUB(propagation_ops_elec_fuse_density_exp_apply)

  end subroutine propagation_ops_elec_fuse_density_exp_apply

  ! ---------------------------------------------------------
  subroutine propagation_ops_elec_interpolate_get(mesh, hm, interp)
    type(mesh_t),                    intent(in)    :: mesh
    type(hamiltonian_elec_t),        intent(inout) :: hm
    type(potential_interpolation_t), intent(inout) :: interp

    PUSH_SUB(propagation_ops_elec_interpolate_get)

    if (family_is_mgga_with_exc(hm%xc)) then
      call potential_interpolation_get(interp, mesh%np, hm%d%nspin, 0, hm%vhxc, vtau = hm%vtau)
    else
      call potential_interpolation_get(interp, mesh%np, hm%d%nspin, 0, hm%vhxc)
    end if

    POP_SUB(propagation_ops_elec_interpolate_get)

  end subroutine propagation_ops_elec_interpolate_get

end module propagation_ops_elec_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

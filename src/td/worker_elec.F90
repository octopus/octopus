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

module worker_elec_oct_m
  use batch_oct_m
  use density_oct_m  
  use exponential_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use lda_u_oct_m
  use messages_oct_m
  use namespace_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use states_elec_oct_m
  use varinfo_oct_m
  use worker_abst_oct_m

  implicit none

  private
  public ::                            &
    worker_elec_t,                     &
    worker_elec_update_hamiltonian,    &
    worker_elec_exp_apply,             &
    worker_elec_fuse_density_exp_apply

  type, extends(worker_abst_t) :: worker_elec_t
    private

  contains

    procedure :: init => worker_elec_init
    procedure :: end => worker_elec_end
  end type worker_elec_t


contains

  subroutine worker_elec_init(wo)
    class(worker_elec_t),  intent(inout) :: wo

    PUSH_SUB(worker_elec_init)

    POP_SUB(worker_elec_init)
  end subroutine

  subroutine worker_elec_end(wo)
    class(worker_elec_t),  intent(inout) :: wo

    PUSH_SUB(worker_elec_end)

    POP_SUB(worker_elec_end)
  end subroutine

  subroutine worker_elec_update_hamiltonian(namespace, st, gr, hm, time)
    type(namespace_t),   intent(in)    :: namespace
    type(states_elec_t), intent(inout) :: st
    type(grid_t),        intent(inout) :: gr
    type(hamiltonian_t), intent(inout) :: hm
    FLOAT,               intent(in)    :: time

    type(profile_t), save :: prof

    PUSH_SUB(worker_elec_update_hamiltonian)

    call profiling_in(prof, 'ELEC_UPDATE_H')

    call hamiltonian_update(hm, gr%mesh, namespace, time = time)
    call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy)

    call profiling_out(prof)

    POP_SUB(worker_elec_update_hamiltonian)

  end subroutine worker_elec_update_hamiltonian

  subroutine worker_elec_exp_apply(te, st, gr, hm, psolver, dt)
    type(exponential_t), intent(inout) :: te 
    type(states_elec_t), intent(inout) :: st
    type(grid_t),        intent(inout) :: gr
    type(hamiltonian_t), intent(inout) :: hm
    type(poisson_t),     intent(in)    :: psolver
    FLOAT,               intent(in)    :: dt

    integer :: ik, ib
    type(profile_t), save :: prof

    PUSH_SUB(worker_elec_exp_apply)

    call profiling_in(prof, 'ELEC_EXP_APPLY')

    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call exponential_apply_batch(te, gr%der, hm, psolver, st%group%psib(ib, ik), ik, dt)
      end do
    end do

    call profiling_out(prof)

    POP_SUB(worker_elec_exp_apply)

  end subroutine worker_elec_exp_apply

  subroutine worker_elec_fuse_density_exp_apply(te, st, gr, hm, psolver, dt, dt2)
    type(exponential_t), intent(inout) :: te
    type(states_elec_t), intent(inout) :: st
    type(grid_t),        intent(inout) :: gr
    type(hamiltonian_t), intent(inout) :: hm
    type(poisson_t),     intent(in)    :: psolver
    FLOAT,               intent(in)    :: dt
    FLOAT,  optional,    intent(in)    :: dt2

    integer :: ik, ib
    type(batch_t) :: zpsib_dt
    type(density_calc_t) :: dens_calc
    type(profile_t), save :: prof

    PUSH_SUB(worker_elec_fuse_density_exp_apply)

    call profiling_in(prof, 'ELEC_FUSE_DENS_EXP_APPLY')

    call density_calc_init(dens_calc, st, gr, st%rho)

    if(present(dt2)) then
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end

          call batch_copy(st%group%psib(ib, ik), zpsib_dt)
          if(batch_is_packed(st%group%psib(ib, ik))) call batch_pack(zpsib_dt, copy = .false.)

          !propagate the state dt/2 and dt, simultaneously, with H(time - dt)
          call exponential_apply_batch(te, gr%der, hm, psolver, st%group%psib(ib, ik), ik, dt, &
            psib2 = zpsib_dt, deltat2 = M_TWO*dt)

          !use the dt propagation to calculate the density
          call density_calc_accumulate(dens_calc, ik, zpsib_dt)

          call batch_end(zpsib_dt)

        end do
      end do

   
    else

      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call exponential_apply_batch(te, gr%der, hm, psolver, st%group%psib(ib, ik), ik, dt)
          call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik))
        end do
      end do

    end if

    call density_calc_end(dens_calc)

    call profiling_out(prof)

    POP_SUB(worker_elec_fuse_density_exp_apply)

  end subroutine worker_elec_fuse_density_exp_apply



end module worker_elec_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

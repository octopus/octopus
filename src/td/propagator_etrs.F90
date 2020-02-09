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

module propagator_etrs_oct_m
  use accel_oct_m
  use batch_oct_m
  use density_oct_m
  use exponential_oct_m
  use grid_oct_m
  use geometry_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_elec_base_oct_m
  use ion_dynamics_oct_m
  use lalg_basic_oct_m
  use lda_u_oct_m
  use lda_u_io_oct_m
  use math_oct_m
  use messages_oct_m
  use mesh_function_oct_m
  use namespace_oct_m
  use parser_oct_m
  use potential_interpolation_oct_m
  use profiling_oct_m
  use propagator_base_oct_m
  use states_elec_dim_oct_m
  use states_elec_oct_m
  use types_oct_m
  use v_ks_oct_m
  use wfs_elec_oct_m
  use propagation_ops_elec_oct_m
  use xc_oct_m

  implicit none

  private

  public ::                       &
    td_etrs,                      &
    td_etrs_sc,                   &
    td_aetrs,                     &
    td_caetrs

contains

  ! ---------------------------------------------------------
  !> Propagator with enforced time-reversal symmetry
  subroutine td_etrs(ks, namespace, hm, gr, st, tr, time, dt, ionic_scale, ions, geo, move_ions)
    type(v_ks_t),             target, intent(inout) :: ks
    type(namespace_t),                intent(in)    :: namespace
    type(hamiltonian_elec_t), target, intent(inout) :: hm
    type(grid_t),             target, intent(inout) :: gr
    type(states_elec_t),      target, intent(inout) :: st
    type(propagator_t),       target, intent(inout) :: tr
    FLOAT,                            intent(in)    :: time
    FLOAT,                            intent(in)    :: dt
    FLOAT,                            intent(in)    :: ionic_scale
    type(ion_dynamics_t),             intent(inout) :: ions
    type(geometry_t),                 intent(inout) :: geo
    logical,                          intent(in)    :: move_ions

    FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)

    PUSH_SUB(td_etrs)

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then

      SAFE_ALLOCATE(vhxc_t1(1:gr%mesh%np, 1:st%d%nspin))
      SAFE_ALLOCATE(vhxc_t2(1:gr%mesh%np, 1:st%d%nspin))
      call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t1)

      call propagation_ops_elec_fuse_density_exp_apply(tr%te, namespace, st, gr, hm, CNST(0.5)*dt, dt)

      call v_ks_calc(ks, namespace, hm, st, geo, calc_current = .false.)

      call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t2)
      call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t1, hm%vhxc)
      call hamiltonian_elec_update(hm, gr%mesh, namespace, time = time - dt)

    else

      call propagation_ops_elec_exp_apply(tr%te, namespace, st, gr%mesh, hm, CNST(0.5)*dt)

    end if

    ! propagate dt/2 with H(t)

    ! first move the ions to time t
    call propagation_ops_elec_move_ions(tr%propagation_ops_elec, gr, hm, st, namespace, ions, geo, &
               time, ionic_scale*dt, move_ions = move_ions)

    call propagation_ops_elec_propagate_gauge_field(tr%propagation_ops_elec, namespace, hm, dt, time)

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t2, hm%vhxc)
    end if
    
    call propagation_ops_elec_update_hamiltonian(namespace, st, gr%mesh, hm, time)
    
    ! propagate dt/2 with H(time - dt)
    call propagation_ops_elec_fuse_density_exp_apply(tr%te, namespace, st, gr, hm, CNST(0.5)*dt)

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      SAFE_DEALLOCATE_A(vhxc_t1)
      SAFE_DEALLOCATE_A(vhxc_t2)
    end if

    POP_SUB(td_etrs)
  end subroutine td_etrs

  ! ---------------------------------------------------------
  !> Propagator with enforced time-reversal symmetry and self-consistency
  subroutine td_etrs_sc(ks, namespace, hm, gr, st, tr, time, dt, ionic_scale, ions, geo, move_ions, sctol, scsteps)
    type(v_ks_t),             target, intent(inout) :: ks
    type(namespace_t),                intent(in)    :: namespace
    type(hamiltonian_elec_t), target, intent(inout) :: hm
    type(grid_t),             target, intent(inout) :: gr
    type(states_elec_t),      target, intent(inout) :: st
    type(propagator_t),       target, intent(inout) :: tr
    FLOAT,                            intent(in)    :: time
    FLOAT,                            intent(in)    :: dt
    FLOAT,                            intent(in)    :: ionic_scale
    type(ion_dynamics_t),             intent(inout) :: ions
    type(geometry_t),                 intent(inout) :: geo
    logical,                          intent(in)    :: move_ions
    FLOAT,                            intent(in)    :: sctol
    integer,                optional, intent(out)   :: scsteps

    FLOAT :: diff
    FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
    integer :: ik, ib, iter, ip
    class(wfs_elec_t), allocatable :: psi2(:, :)
    ! these are hardcoded for the moment
    integer, parameter :: niter = 10

    PUSH_SUB(td_etrs_sc)

    ASSERT(hm%theory_level /= INDEPENDENT_PARTICLES)

    SAFE_ALLOCATE(vhxc_t1(1:gr%mesh%np, 1:st%d%nspin))
    SAFE_ALLOCATE(vhxc_t2(1:gr%mesh%np, 1:st%d%nspin))
    call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t1)

    call messages_new_line()
    call messages_write('        Self-consistency iteration:')
    call messages_info()

    !Propagate the states to t+dt/2 and compute the density at t+dt
    call propagation_ops_elec_fuse_density_exp_apply(tr%te, namespace, st, gr, hm, M_HALF*dt, dt)

    call v_ks_calc(ks, namespace, hm, st, geo, calc_current = .false.)

    call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t2)
    call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t1, hm%vhxc)

    call propagation_ops_elec_update_hamiltonian(namespace, st, gr%mesh, hm, time - dt)

    ! propagate dt/2 with H(t)

    ! first move the ions to time t
    call propagation_ops_elec_move_ions(tr%propagation_ops_elec, gr, hm, st, namespace, ions, geo, &
                time, ionic_scale*dt, move_ions = move_ions)

    call propagation_ops_elec_propagate_gauge_field(tr%propagation_ops_elec, namespace, hm, dt, time)

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t2, hm%vhxc)
    end if

    call propagation_ops_elec_update_hamiltonian(namespace, st, gr%mesh, hm, time)

    SAFE_ALLOCATE_TYPE_ARRAY(wfs_elec_t, psi2, (st%group%block_start:st%group%block_end, st%d%kpt%start:st%d%kpt%end))

    ! store the state at half iteration
    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call st%group%psib(ib, ik)%copy_to(psi2(ib, ik))
        if(st%group%psib(ib, ik)%is_packed()) call psi2(ib, ik)%do_pack(copy = .false.)
        call st%group%psib(ib, ik)%copy_data_to(gr%mesh%np, psi2(ib, ik))
      end do
    end do

    do iter = 1, niter

      call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t2)

      call propagation_ops_elec_fuse_density_exp_apply(tr%te, namespace, st, gr, hm, M_HALF*dt)

      call v_ks_calc(ks, namespace, hm, st, geo, time = time, calc_current = .false.)
      call lda_u_update_occ_matrices(hm%lda_u, namespace, gr%mesh, st, hm%hm_base, hm%energy )

      ! now check how much the potential changed
      do ip = 1, gr%mesh%np
        vhxc_t2(ip, 1) = sum(abs(vhxc_t2(ip, 1:st%d%nspin) - hm%vhxc(ip, 1:st%d%nspin)))
      end do
      diff = dmf_integrate(gr%mesh, vhxc_t2(:, 1))

      call messages_write('          step ')
      call messages_write(iter)
      call messages_write(', residue = ')
      call messages_write(abs(diff), fmt = '(1x,es9.2)')
      call messages_info()

      if(diff <= sctol) exit

      if(iter /= niter) then
        ! we are not converged, restore the states
        do ik = st%d%kpt%start, st%d%kpt%end
          do ib = st%group%block_start, st%group%block_end
            call psi2(ib, ik)%copy_data_to(gr%mesh%np, st%group%psib(ib, ik))
          end do
        end do
      end if

    end do

    if(hm%lda_u_level /= DFT_U_NONE) then 
      call lda_u_write_U(hm%lda_u, stdout) 
      call lda_u_write_V(hm%lda_u, stdout)
    end if

    ! print an empty line
    call messages_info()

    if(present(scsteps)) scsteps = iter

    SAFE_DEALLOCATE_A(vhxc_t1)
    SAFE_DEALLOCATE_A(vhxc_t2)

    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call psi2(ib, ik)%end()
      end do
    end do

    SAFE_DEALLOCATE_A(psi2)

    POP_SUB(td_etrs_sc)
  end subroutine td_etrs_sc

  ! ---------------------------------------------------------
  !> Propagator with approximate enforced time-reversal symmetry
  subroutine td_aetrs(namespace, hm, gr, st, tr, time, dt, ionic_scale, ions, geo, move_ions)
    type(namespace_t),                intent(in)    :: namespace
    type(hamiltonian_elec_t), target, intent(inout) :: hm
    type(grid_t),             target, intent(inout) :: gr
    type(states_elec_t),      target, intent(inout) :: st
    type(propagator_t),       target, intent(inout) :: tr
    FLOAT,                            intent(in)    :: time
    FLOAT,                            intent(in)    :: dt
    FLOAT,                            intent(in)    :: ionic_scale
    type(ion_dynamics_t),             intent(inout) :: ions
    type(geometry_t),                 intent(inout) :: geo
    logical,                          intent(in)    :: move_ions

    PUSH_SUB(td_aetrs)

    ! propagate half of the time step with H(time - dt)
    call propagation_ops_elec_exp_apply(tr%te, namespace, st, gr%mesh, hm, M_HALF*dt)

    !Get the potentials from the interpolator
    call propagation_ops_elec_interpolate_get(gr%mesh, hm, tr%vksold)

    ! move the ions to time t
    call propagation_ops_elec_move_ions(tr%propagation_ops_elec, gr, hm, st, namespace, ions, &
              geo, time, ionic_scale*dt, move_ions = move_ions)

    !Propagate gauge field
    call propagation_ops_elec_propagate_gauge_field(tr%propagation_ops_elec, namespace, hm, dt, time)

    !Update Hamiltonian
    call propagation_ops_elec_update_hamiltonian(namespace, st, gr%mesh, hm, time)

    !Do the time propagation for the second half of the time step
    call propagation_ops_elec_fuse_density_exp_apply(tr%te, namespace, st, gr, hm, M_HALF*dt)

    POP_SUB(td_aetrs)
  end subroutine td_aetrs

  ! ---------------------------------------------------------
  !> Propagator with approximate enforced time-reversal symmetry
  subroutine td_caetrs(ks, namespace, hm, gr, st, tr, time, dt, ionic_scale, ions, geo, move_ions)
    type(v_ks_t),             target, intent(inout) :: ks
    type(namespace_t),                intent(in)    :: namespace
    type(hamiltonian_elec_t), target, intent(inout) :: hm
    type(grid_t),             target, intent(inout) :: gr
    type(states_elec_t),      target, intent(inout) :: st
    type(propagator_t),       target, intent(inout) :: tr
    FLOAT,                            intent(in)    :: time
    FLOAT,                            intent(in)    :: dt
    FLOAT,                            intent(in)    :: ionic_scale
    type(ion_dynamics_t),             intent(inout) :: ions
    type(geometry_t),                 intent(inout) :: geo
    logical,                          intent(in)    :: move_ions

    integer :: ik, ispin, ip, ist, ib
    FLOAT :: vv
    CMPLX :: phase
    type(density_calc_t)  :: dens_calc
    type(profile_t), save :: phase_prof
    integer               :: pnp, iprange
    FLOAT, allocatable    :: vold(:, :), vtauold(:, :)
    type(accel_mem_t)    :: phase_buff

    PUSH_SUB(td_caetrs)

    SAFE_ALLOCATE(vold(1:gr%mesh%np, 1:st%d%nspin))
    if (family_is_mgga_with_exc(hm%xc)) then 
      SAFE_ALLOCATE(vtauold(1:gr%mesh%np, 1:st%d%nspin))
      call potential_interpolation_get(tr%vksold, gr%mesh%np, st%d%nspin, 2, vold, vtau = vtauold)
      call hamiltonian_elec_set_vhxc(hm, gr%mesh, vold, vtauold)
    else
      call potential_interpolation_get(tr%vksold, gr%mesh%np, st%d%nspin, 2, vold)
      call hamiltonian_elec_set_vhxc(hm, gr%mesh, vold)
    endif

    call propagation_ops_elec_update_hamiltonian(namespace, st, gr%mesh, hm, time - dt) 

    call v_ks_calc_start(ks, namespace, hm, st, geo, time = time - dt, calc_energy = .false., &
           calc_current = .false.)

    ! propagate half of the time step with H(time - dt)
    call propagation_ops_elec_exp_apply(tr%te, namespace, st, gr%mesh, hm, M_HALF*dt)

    call v_ks_calc_finish(ks, hm, namespace)

    if (family_is_mgga_with_exc(hm%xc)) then 
      call potential_interpolation_set(tr%vksold, gr%mesh%np, st%d%nspin, 1, hm%vhxc, vtau = hm%vtau)
      call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), &
         tr%vksold%v_old(:, :, 1:3), time, tr%vksold%v_old(:, :, 0))
      call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), &
         tr%vksold%vtau_old(:, :, 1:3), time, tr%vksold%vtau_old(:, :, 0))
      forall(ispin = 1:st%d%nspin, ip = 1:gr%mesh%np)
        vold(ip, ispin) =  CNST(0.5)*dt*(hm%vhxc(ip, ispin) - vold(ip, ispin))
        vtauold(ip, ispin) =  CNST(0.5)*dt*(hm%vtau(ip, ispin) - vtauold(ip, ispin))
      end forall      
    else
      call potential_interpolation_set(tr%vksold, gr%mesh%np, st%d%nspin, 1, hm%vhxc)
      call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), &
         tr%vksold%v_old(:, :, 1:3), time, tr%vksold%v_old(:, :, 0))

      forall(ispin = 1:st%d%nspin, ip = 1:gr%mesh%np)
        vold(ip, ispin) =  CNST(0.5)*dt*(hm%vhxc(ip, ispin) - vold(ip, ispin))
      end forall
    end if

    ! copy vold to a cl buffer
    if(accel_is_enabled() .and. hamiltonian_elec_apply_packed(hm)) then
      if(family_is_mgga_with_exc(hm%xc)) then
        call messages_not_implemented('CAETRS propagator with accel and MGGA with energy functionals', namespace=namespace)
      end if
      pnp = accel_padded_size(gr%mesh%np)
      call accel_create_buffer(phase_buff, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, pnp*st%d%nspin)
      ASSERT(ubound(vold, dim = 1) == gr%mesh%np)
      do ispin = 1, st%d%nspin
        call accel_write_buffer(phase_buff, gr%mesh%np, vold(:, ispin), offset = (ispin - 1)*pnp)
      end do
    end if

    !Get the potentials from the interpolator
    call propagation_ops_elec_interpolate_get(gr%mesh, hm, tr%vksold)

    ! move the ions to time t
    call propagation_ops_elec_move_ions(tr%propagation_ops_elec, gr, hm, st, namespace, ions, &
              geo, time, ionic_scale*dt, move_ions = move_ions)

    call propagation_ops_elec_propagate_gauge_field(tr%propagation_ops_elec, namespace, hm, dt, time)

    call propagation_ops_elec_update_hamiltonian(namespace, st, gr%mesh, hm, time)

    call density_calc_init(dens_calc, st, gr, st%rho)

    ! propagate the other half with H(t)
    do ik = st%d%kpt%start, st%d%kpt%end
      ispin = states_elec_dim_get_spin_index(st%d, ik)

      do ib = st%group%block_start, st%group%block_end
        if (hamiltonian_elec_apply_packed(hm)) then
          call st%group%psib(ib, ik)%do_pack()
          if (hamiltonian_elec_inh_term(hm)) call hm%inh_st%group%psib(ib, ik)%do_pack()
        end if

        call profiling_in(phase_prof, "CAETRS_PHASE")
        select case(st%group%psib(ib, ik)%status())
        case(BATCH_NOT_PACKED)
          do ip = 1, gr%mesh%np
            vv = vold(ip, ispin)
            phase = TOCMPLX(cos(vv), -sin(vv))
            forall(ist = 1:st%group%psib(ib, ik)%nst_linear)
              st%group%psib(ib, ik)%states_linear(ist)%zpsi(ip) = st%group%psib(ib, ik)%states_linear(ist)%zpsi(ip)*phase
            end forall
          end do
        case(BATCH_PACKED)
          do ip = 1, gr%mesh%np
            vv = vold(ip, ispin)
            phase = TOCMPLX(cos(vv), -sin(vv))
            forall(ist = 1:st%group%psib(ib, ik)%nst_linear)
              st%group%psib(ib, ik)%pack%zpsi(ist, ip) = st%group%psib(ib, ik)%pack%zpsi(ist, ip)*phase
            end forall
          end do
        case(BATCH_DEVICE_PACKED)
          call accel_set_kernel_arg(kernel_phase, 0, pnp*(ispin - 1))
          call accel_set_kernel_arg(kernel_phase, 1, phase_buff)
          call accel_set_kernel_arg(kernel_phase, 2, st%group%psib(ib, ik)%pack%buffer)
          call accel_set_kernel_arg(kernel_phase, 3, log2(st%group%psib(ib, ik)%pack%size(1)))

          iprange = accel_max_workgroup_size()/st%group%psib(ib, ik)%pack%size(1)

          call accel_kernel_run(kernel_phase, (/st%group%psib(ib, ik)%pack%size(1), pnp/), &
            (/st%group%psib(ib, ik)%pack%size(1), iprange/))
        end select
        call profiling_out(phase_prof)

        call hamiltonian_elec_base_set_phase_corr(hm%hm_base, gr%mesh, st%group%psib(ib, ik))
        if (hamiltonian_elec_inh_term(hm)) then
          call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, st%group%psib(ib, ik), CNST(0.5)*dt, &
            inh_psib = hm%inh_st%group%psib(ib, ik))
        else
          call exponential_apply_batch(tr%te, namespace, gr%mesh, hm, st%group%psib(ib, ik), CNST(0.5)*dt)
        end if
        call hamiltonian_elec_base_unset_phase_corr(hm%hm_base, gr%mesh, st%group%psib(ib, ik))

        call density_calc_accumulate(dens_calc, st%group%psib(ib, ik))

        if (hamiltonian_elec_apply_packed(hm)) then
          call st%group%psib(ib, ik)%do_unpack()
          if (hamiltonian_elec_inh_term(hm)) call hm%inh_st%group%psib(ib, ik)%do_unpack()
        end if
      end do
    end do

    call density_calc_end(dens_calc)

    if(accel_is_enabled() .and. hamiltonian_elec_apply_packed(hm)) then
      call accel_release_buffer(phase_buff)
    end if

    SAFE_DEALLOCATE_A(vold)
    SAFE_DEALLOCATE_A(vtauold)

    POP_SUB(td_caetrs)
  end subroutine td_caetrs

end module propagator_etrs_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

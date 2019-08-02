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
  use gauge_field_oct_m
  use grid_oct_m
  use geometry_oct_m
  use global_oct_m
  use hamiltonian_oct_m
  use ion_dynamics_oct_m
  use lalg_basic_oct_m
  use lda_u_oct_m
  use lda_u_io_oct_m
  use math_oct_m
  use messages_oct_m
  use mesh_function_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use potential_interpolation_oct_m
  use profiling_oct_m
  use propagator_base_oct_m
  use states_dim_oct_m
  use states_oct_m
  use types_oct_m
  use v_ks_oct_m

  implicit none

  private

  public ::                       &
    td_etrs,                      &
    td_etrs_sc,                   &
    td_aetrs

contains

  ! ---------------------------------------------------------
  !> Propagator with enforced time-reversal symmetry
  subroutine td_etrs(ks, namespace, hm, psolver, gr, st, tr, time, dt, ionic_scale, ions, geo, move_ions)
    type(v_ks_t), target,            intent(inout) :: ks
    type(namespace_t),               intent(in)    :: namespace
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(poisson_t),                 intent(in)    :: psolver
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    type(propagator_t),  target,     intent(inout) :: tr
    FLOAT,                           intent(in)    :: time
    FLOAT,                           intent(in)    :: dt
    FLOAT,                           intent(in)    :: ionic_scale
    type(ion_dynamics_t),            intent(inout) :: ions
    type(geometry_t),                intent(inout) :: geo
    logical,                         intent(in)    :: move_ions

    FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
    integer :: ik, ib
    type(batch_t) :: zpsib_dt
    type(density_calc_t) :: dens_calc

    PUSH_SUB(td_etrs)

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then

      SAFE_ALLOCATE(vhxc_t1(1:gr%mesh%np, 1:st%d%nspin))
      SAFE_ALLOCATE(vhxc_t2(1:gr%mesh%np, 1:st%d%nspin))
      call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t1)

      call density_calc_init(dens_calc, st, gr, st%rho)

      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end

          call batch_copy(st%group%psib(ib, ik), zpsib_dt)
          if(batch_is_packed(st%group%psib(ib, ik))) call batch_pack(zpsib_dt, copy = .false.)

          !propagate the state dt/2 and dt, simultaneously, with H(time - dt)
          call exponential_apply_batch(tr%te, gr%der, hm, psolver, st%group%psib(ib, ik), ik, CNST(0.5)*dt, &
            psib2 = zpsib_dt, deltat2 = dt)

          !use the dt propagation to calculate the density
          call density_calc_accumulate(dens_calc, ik, zpsib_dt)

          call batch_end(zpsib_dt)

        end do
      end do

      call density_calc_end(dens_calc)

      call v_ks_calc(ks, namespace, hm, st, geo, calc_current = .false.)

      call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t2)
      call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t1, hm%vhxc)
      call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = time - dt)

    else

      ! propagate dt/2 with H(time - dt)
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call exponential_apply_batch(tr%te, gr%der, hm, psolver, st%group%psib(ib, ik), ik, CNST(0.5)*dt)
        end do
      end do

    end if

    ! propagate dt/2 with H(t)

    ! first move the ions to time t
    if(move_ions .and. ion_dynamics_ions_move(ions)) then
      call ion_dynamics_propagate(ions, gr%sb, geo, time, ionic_scale*dt)
      call hamiltonian_epot_generate(hm, namespace,  gr, geo, st, psolver, time = time)
    end if

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_propagate(hm%ep%gfield, dt, time)
    end if

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t2, hm%vhxc)
    end if
    call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = time)
    !We update the occupation matrices
    call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy )

    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call exponential_apply_batch(tr%te, gr%der, hm, psolver, st%group%psib(ib, ik), ik, CNST(0.5)*dt)
      end do
    end do

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      SAFE_DEALLOCATE_A(vhxc_t1)
      SAFE_DEALLOCATE_A(vhxc_t2)
    end if

    call density_calc(st, gr, st%rho)

    POP_SUB(td_etrs)
  end subroutine td_etrs

  ! ---------------------------------------------------------
  !> Propagator with enforced time-reversal symmetry and self-consistency
  subroutine td_etrs_sc(ks, namespace, hm, psolver, gr, st, tr, time, dt, ionic_scale, ions, geo, move_ions, sctol, scsteps)
    type(v_ks_t), target,            intent(inout) :: ks
    type(namespace_t),               intent(in)    :: namespace
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(poisson_t),                 intent(in)    :: psolver
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    type(propagator_t),  target,     intent(inout) :: tr
    FLOAT,                           intent(in)    :: time
    FLOAT,                           intent(in)    :: dt
    FLOAT,                           intent(in)    :: ionic_scale
    type(ion_dynamics_t),            intent(inout) :: ions
    type(geometry_t),                intent(inout) :: geo
    logical,                         intent(in)    :: move_ions
    FLOAT,                           intent(in)    :: sctol
    integer,              optional,  intent(out)   :: scsteps

    FLOAT :: diff
    FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
    integer :: ik, ib, iter, ip
    type(batch_t) :: zpsib_dt
    type(density_calc_t) :: dens_calc
    type(batch_t), allocatable :: psi2(:, :)
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
    
    call density_calc_init(dens_calc, st, gr, st%rho)

    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end

        call batch_copy(st%group%psib(ib, ik), zpsib_dt)
        if(batch_is_packed(st%group%psib(ib, ik))) call batch_pack(zpsib_dt, copy = .false.)

        !propagate the state dt/2 and dt, simultaneously, with H(time - dt)
        call exponential_apply_batch(tr%te, gr%der, hm, psolver, st%group%psib(ib, ik), ik, CNST(0.5)*dt, &
          psib2 = zpsib_dt, deltat2 = dt)

        !use the dt propagation to calculate the density
        call density_calc_accumulate(dens_calc, ik, zpsib_dt)

        call batch_end(zpsib_dt)

      end do
    end do

    call density_calc_end(dens_calc)

    call v_ks_calc(ks, namespace, hm, st, geo, calc_current = .false.)

    call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t2)
    call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t1, hm%vhxc)
    call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = time - dt)
    call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy )

    ! propagate dt/2 with H(t)

    ! first move the ions to time t
    if(move_ions .and. ion_dynamics_ions_move(ions)) then
      call ion_dynamics_propagate(ions, gr%sb, geo, time, ionic_scale*dt)
      call hamiltonian_epot_generate(hm, namespace,  gr, geo, st, psolver, time = time)
    end if

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_propagate(hm%ep%gfield, dt, time)
    end if

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t2, hm%vhxc)
    end if

    call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = time)
    call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy )

    SAFE_ALLOCATE(psi2(st%group%block_start:st%group%block_end, st%d%kpt%start:st%d%kpt%end))

    ! store the state at half iteration
    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call batch_copy(st%group%psib(ib, ik), psi2(ib, ik))
        if(batch_is_packed(st%group%psib(ib, ik))) call batch_pack(psi2(ib, ik), copy = .false.)
        call batch_copy_data(gr%mesh%np, st%group%psib(ib, ik), psi2(ib, ik))
      end do
    end do

    do iter = 1, niter

      call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t2)

      call density_calc_init(dens_calc, st, gr, st%rho)
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call exponential_apply_batch(tr%te, gr%der, hm, psolver, st%group%psib(ib, ik), ik, CNST(0.5)*dt)
          call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik))
        end do
      end do
      call density_calc_end(dens_calc)

      call v_ks_calc(ks, namespace, hm, st, geo, time = time, calc_current = .false.)
      call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy )

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
            call batch_copy_data(gr%mesh%np, psi2(ib, ik), st%group%psib(ib, ik))
          end do
        end do
      end if

    end do

    if(hm%lda_u_level /= DFT_U_NONE) then 
      call lda_u_write_U(hm%lda_u, stdout) 
    end if

    ! print an empty line
    call messages_info()

    if(present(scsteps)) scsteps = iter

    SAFE_DEALLOCATE_A(vhxc_t1)
    SAFE_DEALLOCATE_A(vhxc_t2)

    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call batch_end(psi2(ib, ik))
      end do
    end do

    SAFE_DEALLOCATE_A(psi2)

    POP_SUB(td_etrs_sc)
  end subroutine td_etrs_sc

  ! ---------------------------------------------------------
  !> Propagator with approximate enforced time-reversal symmetry
  subroutine td_aetrs(ks, namespace, hm, psolver, gr, st, tr, time, dt, ionic_scale, ions, geo, move_ions)
    type(v_ks_t), target,            intent(inout) :: ks
    type(namespace_t),               intent(in)    :: namespace
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(poisson_t),                 intent(in)    :: psolver
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    type(propagator_t),  target,     intent(inout) :: tr
    FLOAT,                           intent(in)    :: time
    FLOAT,                           intent(in)    :: dt
    FLOAT,                           intent(in)    :: ionic_scale
    type(ion_dynamics_t),            intent(inout) :: ions
    type(geometry_t),                intent(inout) :: geo
    logical,                         intent(in)    :: move_ions

    integer :: ik, ispin, ip, ist, ib
    FLOAT :: vv
    CMPLX :: phase
    type(density_calc_t)  :: dens_calc
    type(profile_t), save :: phase_prof
    integer               :: pnp, iprange
    FLOAT, allocatable    :: vold(:, :), vtauold(:, :)
    type(accel_mem_t)    :: phase_buff

    PUSH_SUB(td_aetrs)

    if(tr%method == PROP_CAETRS) then
      SAFE_ALLOCATE(vold(1:gr%mesh%np, 1:st%d%nspin))
      if(hm%family_is_mgga_with_exc) then 
        call potential_interpolation_get(tr%vksold, gr%mesh%np, st%d%nspin, 2, vold, vtau = vtauold)
        call lalg_copy(gr%mesh%np, st%d%nspin, vold, hm%vhxc)
        call lalg_copy(gr%mesh%np, st%d%nspin, vtauold, hm%vtau)
      else
        call potential_interpolation_get(tr%vksold, gr%mesh%np, st%d%nspin, 2, vold)
        call lalg_copy(gr%mesh%np, st%d%nspin, vold, hm%vhxc)
      endif

      call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = time - dt)
      !We update the occupation matrices
      call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy )
      call v_ks_calc_start(ks, namespace, hm, st, geo, time = time - dt, calc_energy = .false., &
             calc_current = .false.)
    end if

    ! propagate half of the time step with H(time - dt)
    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call exponential_apply_batch(tr%te, gr%der, hm, psolver, st%group%psib(ib, ik), ik, CNST(0.5)*dt)
      end do
    end do

    if(tr%method == PROP_CAETRS) then
      call v_ks_calc_finish(ks, hm)

      if(hm%family_is_mgga_with_exc) then 
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
      if(accel_is_enabled() .and. hamiltonian_apply_packed(hm, gr%mesh)) then
        if(hm%family_is_mgga_with_exc) then
          call messages_not_implemented('CAETRS propagator with accel and MGGA with energy functionals')
        end if
        pnp = accel_padded_size(gr%mesh%np)
        call accel_create_buffer(phase_buff, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, pnp*st%d%nspin)
        ASSERT(ubound(vold, dim = 1) == gr%mesh%np)
        do ispin = 1, st%d%nspin
          call accel_write_buffer(phase_buff, gr%mesh%np, vold(:, ispin), offset = (ispin - 1)*pnp)
        end do
      end if

    end if

    if(hm%family_is_mgga_with_exc) then
      call potential_interpolation_get(tr%vksold, gr%mesh%np, st%d%nspin, 0, hm%vhxc, vtau = hm%vtau)
    else
      call potential_interpolation_get(tr%vksold, gr%mesh%np, st%d%nspin, 0, hm%vhxc)
    end if

    ! move the ions to time t
    if(move_ions .and. ion_dynamics_ions_move(ions)) then
      call ion_dynamics_propagate(ions, gr%sb, geo, time, ionic_scale*dt)
      call hamiltonian_epot_generate(hm, namespace, gr, geo, st, psolver, time = time)
    end if

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_propagate(hm%ep%gfield, dt, time)
    end if

    call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = time)
    !We update the occupation matrices
    call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy )

    call density_calc_init(dens_calc, st, gr, st%rho)

    ! propagate the other half with H(t)
    do ik = st%d%kpt%start, st%d%kpt%end
      ispin = states_dim_get_spin_index(st%d, ik)

      do ib = st%group%block_start, st%group%block_end

        if(tr%method == PROP_CAETRS) then
          call profiling_in(phase_prof, "CAETRS_PHASE")
          select case(batch_status(st%group%psib(ib, ik)))
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
        end if

        call exponential_apply_batch(tr%te, gr%der, hm, psolver, st%group%psib(ib, ik), ik, CNST(0.5)*dt)
        call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik))

      end do
    end do

    if(tr%method == PROP_CAETRS .and. accel_is_enabled() .and. hamiltonian_apply_packed(hm, gr%mesh)) then
      call accel_release_buffer(phase_buff)
    end if

    call density_calc_end(dens_calc)

    if(tr%method == PROP_CAETRS) then
      SAFE_DEALLOCATE_A(vold)
    end if

    POP_SUB(td_aetrs)
  end subroutine td_aetrs

end module propagator_etrs_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

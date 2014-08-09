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
    !> Propagator with enforced time-reversal symmetry
    subroutine td_etrs(ks, hm, gr, st, tr, time, dt, mu, ions, geo, gauge_force)
    type(v_ks_t), target,            intent(inout) :: ks
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    type(propagator_t),  target,     intent(inout) :: tr
    FLOAT,                           intent(in)    :: time
    FLOAT,                           intent(in)    :: dt
    FLOAT,                           intent(in)    :: mu
    type(ion_dynamics_t),            intent(inout) :: ions
    type(geometry_t),                intent(inout) :: geo
    type(gauge_force_t),  optional,  intent(inout) :: gauge_force

      FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
      integer :: ik, ib
      type(batch_t) :: zpsib_save
      type(density_calc_t) :: dens_calc

      PUSH_SUB(propagator_dt.td_etrs)

      if(hm%theory_level /= INDEPENDENT_PARTICLES) then

        SAFE_ALLOCATE(vhxc_t1(1:gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE(vhxc_t2(1:gr%mesh%np, 1:st%d%nspin))
        call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t1)

        call density_calc_init(dens_calc, st, gr, st%rho)

        do ik = st%d%kpt%start, st%d%kpt%end
          do ib = st%group%block_start, st%group%block_end

            !save the state
            call batch_copy(st%group%psib(ib, ik), zpsib_save, reference = .false.)
            if(batch_is_packed(st%group%psib(ib, ik))) call batch_pack(zpsib_save, copy = .false.)
            call batch_copy_data(gr%der%mesh%np, st%group%psib(ib, ik), zpsib_save)

            !propagate the state dt with H(time - dt)
            call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, dt/mu, time - dt)
            call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik))

            !restore the saved state
            call batch_copy_data(gr%der%mesh%np, zpsib_save, st%group%psib(ib, ik))

            call batch_end(zpsib_save)

          end do
        end do

        call density_calc_end(dens_calc)

        call v_ks_calc(ks, hm, st, geo)

        call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t2)
        call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t1, hm%vhxc)
        call hamiltonian_update(hm, gr%mesh, time = time - dt)
      end if

      ! propagate dt/2 with H(time - dt)
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, dt/(mu*M_TWO), time - dt)
        end do
      end do

      ! propagate dt/2 with H(t)

      ! first move the ions to time t
      if(ion_dynamics_ions_move(ions)) then
        call ion_dynamics_propagate(ions, gr%sb, geo, time, dt)
        call hamiltonian_epot_generate(hm, gr, geo, st, time = time)
      end if

      if(gauge_field_is_applied(hm%ep%gfield)) then
        call gauge_field_propagate(hm%ep%gfield, gauge_force, dt)
      end if

      if(hm%theory_level /= INDEPENDENT_PARTICLES) then
        call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t2, hm%vhxc)
      end if
      call hamiltonian_update(hm, gr%mesh, time = time)

      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, dt/(M_TWO*mu), time)
        end do
      end do

      if(hm%theory_level /= INDEPENDENT_PARTICLES) then
        SAFE_DEALLOCATE_A(vhxc_t1)
        SAFE_DEALLOCATE_A(vhxc_t2)
      end if

      if(.not. hm%cmplxscl%space) then
        call density_calc(st, gr, st%rho)
      else
        call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
      end if

      POP_SUB(propagator_dt.td_etrs)
    end subroutine td_etrs


    ! ---------------------------------------------------------
    !> Propagator with approximate enforced time-reversal symmetry
    subroutine td_aetrs(ks, hm, gr, st, tr, time, dt, mu, ions, geo, gauge_force)
    type(v_ks_t), target,            intent(inout) :: ks
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    type(propagator_t),  target,     intent(inout) :: tr
    FLOAT,                           intent(in)    :: time
    FLOAT,                           intent(in)    :: dt
    FLOAT,                           intent(in)    :: mu
    type(ion_dynamics_t),            intent(inout) :: ions
    type(geometry_t),                intent(inout) :: geo
    type(gauge_force_t),  optional,  intent(inout) :: gauge_force

      integer :: ik, ispin, ip, ist, ib
      FLOAT :: vv
      CMPLX :: phase
      type(density_calc_t)  :: dens_calc
      type(profile_t), save :: phase_prof
      integer               :: pnp, iprange
      FLOAT, allocatable    :: vold(:, :), imvold(:, :)
      type(opencl_mem_t)    :: phase_buff

      PUSH_SUB(propagator_dt.td_aetrs)

      if(tr%method == PROP_CAETRS) then
        SAFE_ALLOCATE(vold(1:gr%mesh%np, 1:st%d%nspin))
        if(hm%cmplxscl%space) SAFE_ALLOCATE(Imvold(1:gr%mesh%np, 1:st%d%nspin))
        if(hm%cmplxscl%space) then
          call vksinterp_get(tr%vksold, gr%mesh%np, st%d%nspin, 2, vold, imvold)
        else
          call vksinterp_get(tr%vksold, gr%mesh%np, st%d%nspin, 2, vold)
        end if
        call lalg_copy(gr%mesh%np, st%d%nspin, vold, hm%vhxc)
        if(hm%cmplxscl%space) call lalg_copy(gr%mesh%np, st%d%nspin, Imvold, hm%Imvhxc)
        call hamiltonian_update(hm, gr%mesh, time = time - dt)
        call v_ks_calc_start(ks, hm, st, geo, time = time - dt, calc_energy = .false.)
      end if

      ! propagate half of the time step with H(time - dt)
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, dt/(M_TWO*mu), time - dt)
        end do
      end do

      if(tr%method == PROP_CAETRS) then
        call v_ks_calc_finish(ks, hm)

        call vksinterp_set(tr%vksold, gr%mesh%np, st%d%nspin, 1, hm%vhxc)
        call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), &
          tr%vksold%v_old(:, :, 1:3), time, tr%vksold%v_old(:, :, 0))

        forall(ispin = 1:st%d%nspin, ip = 1:gr%mesh%np) 
          vold(ip, ispin) =  dt/(M_TWO*mu)*(hm%vhxc(ip, ispin) - vold(ip, ispin))
        end forall

        ! copy vold to a cl buffer
        if(opencl_is_enabled() .and. hamiltonian_apply_packed(hm, gr%mesh)) then
#ifdef HAVE_OPENCL
          pnp = opencl_padded_size(gr%mesh%np)
          call opencl_create_buffer(phase_buff, CL_MEM_READ_ONLY, TYPE_FLOAT, pnp*st%d%nspin)
          ASSERT(ubound(vold, dim = 1) == gr%mesh%np)
          do ispin = 1, st%d%nspin
            call opencl_write_buffer(phase_buff, gr%mesh%np, vold(:, ispin), offset = (ispin - 1)*pnp)
          end do
#endif
        end if

      end if

      call vksinterp_get(tr%vksold, gr%mesh%np, st%d%nspin, 0, hm%vhxc)

      ! move the ions to time t
      if(ion_dynamics_ions_move(ions)) then
        call ion_dynamics_propagate(ions, gr%sb, geo, time, dt)
        call hamiltonian_epot_generate(hm, gr, geo, st, time = time)
      end if

      if(gauge_field_is_applied(hm%ep%gfield)) then
        call gauge_field_propagate(hm%ep%gfield, gauge_force, dt)
      end if

      call hamiltonian_update(hm, gr%mesh, time = time)

      call density_calc_init(dens_calc, st, gr, st%rho)

      ! propagate the other half with H(t)
      do ik = st%d%kpt%start, st%d%kpt%end
        ispin = states_dim_get_spin_index(st%d, ik)

        do ib = st%group%block_start, st%group%block_end
          if(hamiltonian_apply_packed(hm, gr%mesh)) call batch_pack(st%group%psib(ib, ik))
          
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
            case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
              call opencl_set_kernel_arg(kernel_phase, 0, pnp*(ispin - 1))
              call opencl_set_kernel_arg(kernel_phase, 1, phase_buff)
              call opencl_set_kernel_arg(kernel_phase, 2, st%group%psib(ib, ik)%pack%buffer)
              call opencl_set_kernel_arg(kernel_phase, 3, log2(st%group%psib(ib, ik)%pack%size(1)))

              iprange = opencl_max_workgroup_size()/st%group%psib(ib, ik)%pack%size(1)

              call opencl_kernel_run(kernel_phase, (/st%group%psib(ib, ik)%pack%size(1), pnp/), &
                (/st%group%psib(ib, ik)%pack%size(1), iprange/))
#endif
            end select
            call profiling_out(phase_prof)
          end if

          call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, dt/(M_TWO*mu), time)
          call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik))

          if(hamiltonian_apply_packed(hm, gr%mesh)) call batch_unpack(st%group%psib(ib, ik))
        end do
      end do

#ifdef HAVE_OPENCL
      if(tr%method == PROP_CAETRS .and. opencl_is_enabled() .and. hamiltonian_apply_packed(hm, gr%mesh)) then
        call opencl_release_buffer(phase_buff)
      end if
#endif
      
      call density_calc_end(dens_calc)

      if(tr%method == PROP_CAETRS) then
        SAFE_DEALLOCATE_A(vold)
      end if

      POP_SUB(propagator_dt.td_aetrs)
    end subroutine td_aetrs



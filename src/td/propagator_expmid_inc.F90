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
  !> Exponential midpoint
  subroutine exponential_midpoint(ks, hm, gr, st, tr, time, dt, mu, ions, geo, gauge_force)
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

    integer :: ib, ik
    type(ion_state_t) :: ions_state
    FLOAT :: vecpot(1:MAX_DIM), vecpot_vel(1:MAX_DIM)
    CMPLX :: zt, zdt

    PUSH_SUB(propagator_dt.exponential_midpoint)
      
    if(hm%cmplxscl%time) then
      zt =  TOCMPLX(time, M_ZERO) *exp(M_zI * TOCMPLX(hm%cmplxscl%alphaR, M_ZERO))
      zdt = TOCMPLX(dt,   M_ZERO) *exp(M_zI * TOCMPLX(hm%cmplxscl%alphaR, M_ZERO))
        
      !FIXME: not adapted yet
      if(hm%theory_level /= INDEPENDENT_PARTICLES) then
        if(hm%cmplxscl%space) then
          call vksinterp_interpolate(tr%vksold, 3, &
            gr%mesh%np, st%d%nspin, time, dt, time - dt/M_TWO, hm%vhxc, hm%imvhxc)
        else
          call vksinterp_interpolate(tr%vksold, 3, &
            gr%mesh%np, st%d%nspin, time, dt, time - dt/M_TWO, hm%vhxc)
        end if
      end if
      !FIXME: not implemented yet
      !move the ions to time 'time - dt/2'
      if(ion_dynamics_ions_move(ions)) then
        call ion_dynamics_save_state(ions, geo, ions_state)
        call ion_dynamics_propagate(ions, gr%sb, geo, time - dt/M_TWO, M_HALF*dt)
        call hamiltonian_epot_generate(hm, gr, geo, st, time = time - dt/M_TWO)
      end if
        
      !FIXME: not implemented yet      
      if(gauge_field_is_applied(hm%ep%gfield)) then
        vecpot = gauge_field_get_vec_pot(hm%ep%gfield)
        vecpot_vel = gauge_field_get_vec_pot_vel(hm%ep%gfield)
        call gauge_field_propagate(hm%ep%gfield, gauge_force, M_HALF*dt)
      end if

      call hamiltonian_update(hm, gr%mesh, time = real(zt - zdt/M_z2, REAL_PRECISION), Imtime = aimag(zt - zdt/M_z2  ))

      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
           call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, &
             real(zdt/mu,REAL_PRECISION), real(zt - zdt/M_z2,REAL_PRECISION), &
             Imdeltat = aimag(zdt/mu), Imtime = aimag(zt -  zdt / M_z2 ) )

        end do
      end do
        
    else

      if(hm%theory_level /= INDEPENDENT_PARTICLES) then
        if(hm%cmplxscl%space) then
          call vksinterp_interpolate(tr%vksold, 3, &
            gr%mesh%np, st%d%nspin, time, dt, time - dt/M_TWO, hm%vhxc, hm%imvhxc)
        else
          call vksinterp_interpolate(tr%vksold, 3, &
            gr%mesh%np, st%d%nspin, time, dt, time - dt/M_TWO, hm%vhxc)
        end if
      end if
      !move the ions to time 'time - dt/2'
      if(ion_dynamics_ions_move(ions)) then
        call ion_dynamics_save_state(ions, geo, ions_state)
        call ion_dynamics_propagate(ions, gr%sb, geo, time - dt/M_TWO, M_HALF*dt)
        call hamiltonian_epot_generate(hm, gr, geo, st, time = time - dt/M_TWO)
      end if
     
      if(gauge_field_is_applied(hm%ep%gfield)) then
        vecpot = gauge_field_get_vec_pot(hm%ep%gfield)
        vecpot_vel = gauge_field_get_vec_pot_vel(hm%ep%gfield)
        call gauge_field_propagate(hm%ep%gfield, gauge_force, M_HALF*dt)
      end if
      call hamiltonian_update(hm, gr%mesh, time = time - M_HALF*dt)
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end

          call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, dt/mu, time - dt/M_TWO)

        end do
      end do
    end if

    if(hm%cmplxscl%space) then ! Propagate the left state
      ! (L(t+dt)| = (L|U(t-dt) = (L|e^{i H_\theta(t-dt/2) (dt)}
        
      if(hm%cmplxscl%time) then
        zt =  TOCMPLX(time,M_ZERO) *exp(M_zI * TOCMPLX(hm%cmplxscl%alphaL, M_ZERO))
        zdt = TOCMPLX(dt,  M_ZERO) *exp(M_zI * TOCMPLX(hm%cmplxscl%alphaL, M_ZERO))
        ! FIXME: check this interpolation!! 
        ! probably need some rethinking 
        if(hm%theory_level /= INDEPENDENT_PARTICLES) then
           call vksinterp_interpolate(tr%vksold, 3, gr%mesh%np, st%d%nspin, &
             time, dt, time+dt/M_TWO, hm%vhxc, hm%imvhxc)
        end if
      
        call hamiltonian_update(hm, gr%mesh, time = real(zt + zdt/M_z2, REAL_PRECISION), Imtime = aimag(zt + zdt/M_z2  ))
          
        do ik = st%d%kpt%start, st%d%kpt%end
          do ib = st%group%block_start, st%group%block_end
            call exponential_apply_batch(tr%te, gr%der, hm, st%psibL(ib, ik), ik,&
              real(-zdt/mu, REAL_PRECISION), real(zt + zdt/M_z2, REAL_PRECISION), &
              Imdeltat = aimag(-zdt/mu), Imtime = aimag(zt +  zdt / M_z2 ) )
           end do
        end do
       
      else
          
        ! FIXME: check this interpolation!! 
        ! probably need some rethinking 
        if(hm%theory_level /= INDEPENDENT_PARTICLES) then
          call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%vksold%v_old(:, :, 0:2), time + dt/M_TWO, hm%vhxc(:, :))
          if(hm%cmplxscl%space) &
            call interpolate( (/time, time - dt, time - M_TWO*dt/), &
              tr%vksold%Imv_old(:, :, 0:2), time + dt/M_TWO, hm%Imvhxc(:, :))
        end if

        call hamiltonian_update(hm, gr%mesh, time = time + M_HALF*dt)

        do ik = st%d%kpt%start, st%d%kpt%end
          do ib = st%group%block_start, st%group%block_end
            call exponential_apply_batch(tr%te, gr%der, hm, st%psibL(ib, ik), ik, -dt/mu, time + dt/M_TWO)
          end do
        end do
        
      end if
                
    end if

    !restore to time 'time - dt'
    if(ion_dynamics_ions_move(ions)) call ion_dynamics_restore_state(ions, geo, ions_state)

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_set_vec_pot(hm%ep%gfield, vecpot)
      call gauge_field_set_vec_pot_vel(hm%ep%gfield, vecpot_vel)
      call hamiltonian_update(hm, gr%mesh)
    end if

    if(.not. hm%cmplxscl%space) then
      call density_calc(st, gr, st%rho)
    else
      call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
    end if

    POP_SUB(propagator_dt.exponential_midpoint)
  end subroutine exponential_midpoint
  ! ---------------------------------------------------------


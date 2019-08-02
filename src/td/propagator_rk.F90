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

module propagator_rk_oct_m
  use batch_ops_oct_m
  use comm_oct_m
  use density_oct_m
  use forces_oct_m
  use gauge_field_oct_m
  use grid_oct_m
  use geometry_oct_m
  use global_oct_m
  use hamiltonian_oct_m
  use ion_dynamics_oct_m
  use lda_u_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use namespace_oct_m
  use oct_exchange_oct_m
  use opt_control_state_oct_m
  use parser_oct_m
  use pert_oct_m
  use poisson_oct_m
  use potential_interpolation_oct_m
  use profiling_oct_m
  use propagator_base_oct_m
  use species_oct_m
  use sparskit_oct_m
  use states_oct_m
  use v_ks_oct_m
  use xc_oct_m

  implicit none

  private

  public ::                    &
    td_explicit_runge_kutta4,  &
    td_runge_kutta2,           &
    td_runge_kutta4
  
  type(grid_t),            pointer, private :: grid_p
  type(hamiltonian_t),     pointer, private :: hm_p
  type(poisson_t),         pointer, private :: psolver_p
  type(states_t),          pointer, private :: st_p
  type(xc_t),              pointer, private :: xc_p
  type(propagator_t),      pointer, private :: tr_p
  integer,                 private :: dim_op
  FLOAT,                   private :: t_op, dt_op
  FLOAT, allocatable, private      :: vhxc1_op(:, :), vhxc2_op(:, :), vpsl1_op(:), vpsl2_op(:)
  logical :: move_ions_op
  
contains
  
  subroutine td_explicit_runge_kutta4(ks, namespace, hm, psolver, gr, st, time, dt, ions, geo, qcchi)
    type(v_ks_t), target,            intent(inout) :: ks
    type(namespace_t),               intent(in)    :: namespace
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(poisson_t),     target,     intent(in)    :: psolver
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    FLOAT,                           intent(in)    :: time
    FLOAT,                           intent(in)    :: dt
    type(ion_dynamics_t),            intent(inout) :: ions
    type(geometry_t),                intent(inout) :: geo
    type(opt_control_state_t), optional, target, intent(inout) :: qcchi

    type(states_t), pointer :: chi
    FLOAT, pointer :: q(:, :), p(:, :)

    integer :: np_part, np, kp1, kp2, st1, st2, nspin, ik, ist, iatom, ib
    CMPLX, allocatable :: zphi(:, :, :, :), zchi(:, :, :, :), dvpsi(:, :, :)
    type(states_t) :: hst, stphi, inh, hchi, stchi
    logical :: propagate_chi
    FLOAT, allocatable :: pos0(:, :), vel0(:, :), &
      posk(:, :), velk(:, :), &
      pos(:, :), vel(:, :), &
      posfinal(:, :), velfinal(:, :), &
      pos0t(:, :), vel0t(:, :), &
      poskt(:, :), velkt(:, :), &
      post(:, :), velt(:, :), &
      posfinalt(:, :), velfinalt(:, :), &
      coforce(:, :)

    PUSH_SUB(td_explicit_runge_kutta4)

    propagate_chi = present(qcchi)
    if(propagate_chi) then
      chi => opt_control_point_qs(qcchi)
      q => opt_control_point_q(qcchi)
      p => opt_control_point_p(qcchi)
    end if

    st1 = st%st_start
    st2 = st%st_end
    kp1 = st%d%kpt%start
    kp2 = st%d%kpt%end
    np_part = gr%mesh%np_part
    np = gr%mesh%np
    nspin = hm%d%nspin

    SAFE_ALLOCATE(zphi(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
    if(propagate_chi) then
      SAFE_ALLOCATE(zchi(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
    end if
    if(ion_dynamics_ions_move(ions)) then
      SAFE_ALLOCATE(pos(1:geo%space%dim, 1:geo%natoms))
      SAFE_ALLOCATE(vel(1:geo%space%dim, 1:geo%natoms))
      SAFE_ALLOCATE(pos0(1:geo%space%dim, 1:geo%natoms))
      SAFE_ALLOCATE(vel0(1:geo%space%dim, 1:geo%natoms))
      SAFE_ALLOCATE(posk(1:geo%space%dim, 1:geo%natoms))
      SAFE_ALLOCATE(velk(1:geo%space%dim, 1:geo%natoms))
      SAFE_ALLOCATE(posfinal(1:geo%space%dim, 1:geo%natoms))
      SAFE_ALLOCATE(velfinal(1:geo%space%dim, 1:geo%natoms))

      if(propagate_chi) then
        SAFE_ALLOCATE(post(1:geo%space%dim, 1:geo%natoms))
        SAFE_ALLOCATE(velt(1:geo%space%dim, 1:geo%natoms))
        SAFE_ALLOCATE(pos0t(1:geo%space%dim, 1:geo%natoms))
        SAFE_ALLOCATE(vel0t(1:geo%space%dim, 1:geo%natoms))
        SAFE_ALLOCATE(poskt(1:geo%space%dim, 1:geo%natoms))
        SAFE_ALLOCATE(velkt(1:geo%space%dim, 1:geo%natoms))
        SAFE_ALLOCATE(posfinalt(1:geo%space%dim, 1:geo%natoms))
        SAFE_ALLOCATE(velfinalt(1:geo%space%dim, 1:geo%natoms))
        SAFE_ALLOCATE(coforce(1:geo%natoms, 1:geo%space%dim))
      end if
    end if

    call states_copy(hst, st)
    call states_copy(stphi, st)
    call states_get_state(st, gr%mesh, zphi)

    if(propagate_chi) then
      call states_copy(hchi, chi)
      call states_copy(stchi, chi)
      call states_get_state(chi, gr%mesh, zchi)
    end if

    if(ion_dynamics_ions_move(ions)) then
      do iatom = 1, geo%natoms
        pos0(1:geo%space%dim, iatom) = geo%atom(iatom)%x(1:geo%space%dim)
        vel0(1:geo%space%dim, iatom) = geo%atom(iatom)%v(1:geo%space%dim)
      end do
      posfinal = pos0
      velfinal = vel0

      if(propagate_chi) then
        do iatom = 1, geo%natoms
          pos0t(1:geo%space%dim, iatom) = q(iatom, 1:geo%space%dim)
          vel0t(1:geo%space%dim, iatom) = p(iatom, 1:geo%space%dim) / species_mass(geo%atom(iatom)%species)
        end do
        posfinalt = pos0t
        velfinalt = vel0t
      end if
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Stage 1.
    !

    ! 
    call states_set_state(stphi, gr%mesh, zphi)
    if(propagate_chi) then
      call states_set_state(stchi, gr%mesh, zchi)    
    end if
    if(ion_dynamics_ions_move(ions)) then
      pos = pos0
      vel = vel0
      if(propagate_chi) then
        post = pos0t
        velt = vel0t
      end if
    end if

    call f_psi(time - dt)
    if(propagate_chi) call f_chi(time - dt)
    if(ion_dynamics_ions_move(ions)) call f_ions(time - dt)

    call update_state(M_ONE/CNST(6.0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Stage 2.
    
    call states_set_state(stphi, gr%mesh, zphi)
    if(propagate_chi) then
      call states_set_state(stchi, gr%mesh, zchi)    
    end if
    
    do ik = stphi%d%kpt%start, stphi%d%kpt%end
      do ib = stphi%group%block_start, stphi%group%block_end
        call batch_axpy(gr%mesh%np, -CNST(0.5)*M_zI*dt, hst%group%psib(ib, ik), stphi%group%psib(ib, ik))
        if(propagate_chi) then
          call batch_axpy(gr%mesh%np, -CNST(0.5)*M_zI*dt, hchi%group%psib(ib, ik), stchi%group%psib(ib, ik))
        end if
      end do
    end do
        
    if(ion_dynamics_ions_move(ions)) then
      pos = pos0 + M_HALF * posk
      vel = vel0 + M_HALF * velk
      if(propagate_chi) then
        post = pos0t + M_HALF * poskt
        velt = vel0t + M_HALF * velkt
      end if
    end if

    call f_psi(time-M_HALF*dt)
    if(propagate_chi) call f_chi(time-M_HALF*dt)
    if(ion_dynamics_ions_move(ions)) call f_ions(time-M_HALF*dt)
    call update_state(M_THIRD)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Stage 3.

    call states_set_state(stphi, gr%mesh, zphi)
    if(propagate_chi) then
      call states_set_state(stchi, gr%mesh, zchi)    
    end if
    
    do ik = stphi%d%kpt%start, stphi%d%kpt%end
      do ib = stphi%group%block_start, stphi%group%block_end
        call batch_axpy(gr%mesh%np, -CNST(0.5)*M_zI*dt, hst%group%psib(ib, ik), stphi%group%psib(ib, ik))
        if(propagate_chi) then
          call batch_axpy(gr%mesh%np, -CNST(0.5)*M_zI*dt, hchi%group%psib(ib, ik), stchi%group%psib(ib, ik))
        end if
      end do
    end do

    if(ion_dynamics_ions_move(ions)) then
      pos = pos0 + M_HALF * posk
      vel = vel0 + M_HALF * velk
      if(propagate_chi) then
        post = pos0t + M_HALF * poskt
        velt = vel0t + M_HALF * velkt
      end if
    end if

    call f_psi(time-M_HALF*dt)
    if(propagate_chi) call f_chi(time-M_HALF*dt)
    if(ion_dynamics_ions_move(ions)) call f_ions(time-M_HALF*dt)

    call update_state(M_THIRD)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Stage 4.

    call states_set_state(stphi, gr%mesh, zphi)
    if(propagate_chi) then
      call states_set_state(stchi, gr%mesh, zchi)    
    end if
    
    do ik = stphi%d%kpt%start, stphi%d%kpt%end
      do ib = stphi%group%block_start, stphi%group%block_end
        call batch_axpy(gr%mesh%np, -M_zI*dt, hst%group%psib(ib, ik), stphi%group%psib(ib, ik))
        if(propagate_chi) then
          call batch_axpy(gr%mesh%np, -M_zI*dt, hchi%group%psib(ib, ik), stchi%group%psib(ib, ik))
        end if
      end do
    end do

    if(ion_dynamics_ions_move(ions)) then
      pos = pos0 + posk
      vel = vel0 + velk
      if(propagate_chi) then
        post = pos0t + poskt
        velt = vel0t + velkt
      end if
    end if

    call f_psi(time)
    if(propagate_chi) call f_chi(time)
    if(ion_dynamics_ions_move(ions)) call f_ions(time)

    call update_state( M_ONE/CNST(6.0) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Collect the results.

    call density_calc(st, gr, st%rho)
    if(ion_dynamics_ions_move(ions)) then
      do iatom = 1, geo%natoms
        geo%atom(iatom)%x(1:geo%space%dim) = posfinal(:, iatom)
        geo%atom(iatom)%v(1:geo%space%dim) = velfinal(:, iatom)
      end do
      call hamiltonian_epot_generate(hm, namespace,  gr, geo, st, psolver, time)
      !call forces_calculate(gr, namespace, geo, hm, stphi, time, dt)
      geo%kinetic_energy = ion_dynamics_kinetic_energy(geo)

      if(propagate_chi) then
        do iatom = 1, geo%natoms
          q(iatom, 1:geo%space%dim) = posfinalt(1:geo%space%dim, iatom)
          p(iatom, 1:geo%space%dim) = species_mass(geo%atom(iatom)%species) * velfinalt(1:geo%space%dim, iatom)
        end do
      end if
    end if

    call states_end(hst)
    call states_end(stphi)
    SAFE_DEALLOCATE_A(zphi)
    if(propagate_chi) then
      call states_end(hchi)
      call states_end(stchi)
      SAFE_DEALLOCATE_A(zchi)
      nullify(chi)
      nullify(p)
      nullify(q)
    end if

    if(ion_dynamics_ions_move(ions)) then
      SAFE_DEALLOCATE_A(pos)
      SAFE_DEALLOCATE_A(vel)
      SAFE_DEALLOCATE_A(pos0)
      SAFE_DEALLOCATE_A(vel0)
      SAFE_DEALLOCATE_A(posk)
      SAFE_DEALLOCATE_A(velk)
      SAFE_DEALLOCATE_A(posfinal)
      SAFE_DEALLOCATE_A(velfinal)

      if(propagate_chi) then
        SAFE_DEALLOCATE_A(post)
        SAFE_DEALLOCATE_A(velt)
        SAFE_DEALLOCATE_A(pos0t)
        SAFE_DEALLOCATE_A(vel0t)
        SAFE_DEALLOCATE_A(poskt)
        SAFE_DEALLOCATE_A(velkt)
        SAFE_DEALLOCATE_A(posfinalt)
        SAFE_DEALLOCATE_A(velfinalt)
        SAFE_DEALLOCATE_A(coforce)
      end if
    end if
    POP_SUB(td_explicit_runge_kutta4)

  contains

    subroutine f_psi(tau)
      FLOAT, intent(in) :: tau

      if(ion_dynamics_ions_move(ions)) then
        do iatom = 1, geo%natoms
          geo%atom(iatom)%x(1:geo%space%dim) = pos(:, iatom)
          geo%atom(iatom)%v(1:geo%space%dim) = vel(:, iatom)
        end do
        call hamiltonian_epot_generate(hm, namespace,  gr, geo, stphi, psolver, time = tau)
      end if
      if(.not.oct_exchange_enabled(hm%oct_exchange)) then
        call density_calc(stphi, gr, stphi%rho)
        call v_ks_calc(ks, namespace, hm, stphi, geo, calc_current = gauge_field_is_applied(hm%ep%gfield), time = tau)
      else
        call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = tau)
      end if
      call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy)
      call zhamiltonian_apply_all(hm, ks%xc, gr%der, psolver, stphi, hst)
    end subroutine f_psi

    subroutine f_ions(tau)
      FLOAT, intent(in) :: tau

      call forces_calculate(gr, namespace, geo, hm, stphi, ks, t = tau, dt = dt)
      do iatom = 1, geo%natoms
        posk(:, iatom) = dt * vel(:, iatom)
        velk(:, iatom) = dt * geo%atom(iatom)%f(1:geo%space%dim) / species_mass(geo%atom(iatom)%species)
      end do
      if(propagate_chi) then
        call forces_costate_calculate(gr, namespace, geo, hm, stphi, stchi, coforce, transpose(post))
        do iatom = 1, geo%natoms
          poskt(:, iatom) = dt * velt(:, iatom)
          velkt(:, iatom) = dt * coforce(iatom, :) / species_mass(geo%atom(iatom)%species)
        end do
      end if
    end subroutine f_ions

    subroutine f_chi(tau)
      FLOAT, intent(in) :: tau

      if( hm%theory_level /= INDEPENDENT_PARTICLES) call oct_exchange_set(hm%oct_exchange, stphi, gr%mesh)
      call prepare_inh()
      call hamiltonian_adjoint(hm)
      call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = tau)
      call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy)
      call zhamiltonian_apply_all(hm, ks%xc, gr%der, psolver, stchi, hchi)
      call hamiltonian_not_adjoint(hm)


      call apply_inh()
      if( hm%theory_level /= INDEPENDENT_PARTICLES) call oct_exchange_remove(hm%oct_exchange)
      if(ion_dynamics_ions_move(ions)) call hamiltonian_remove_inh(hm)
    end subroutine f_chi

    subroutine apply_inh()
      integer :: ib

      if(hamiltonian_inh_term(hm)) then
        do ik = kp1, kp2
          do ib = 1, st%group%block_start, st%group%block_end
            call batch_axpy(np, M_ZI, hm%inh_st%group%psib(ib, ik), hchi%group%psib(ib, ik))
          end do
        end do
      end if
    end subroutine apply_inh

    subroutine prepare_inh()
      integer :: idir
      CMPLX, allocatable :: psi(:, :), inhpsi(:, :)
      type(pert_t) :: pert

      if(ion_dynamics_ions_move(ions)) then
        call states_copy(inh, st)

        SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:st%d%dim))
        SAFE_ALLOCATE(inhpsi(1:gr%mesh%np_part, 1:st%d%dim))
        SAFE_ALLOCATE(dvpsi(1:gr%mesh%np_part, 1:st%d%dim, 1:gr%sb%dim))      

        do ik = 1, st%d%nik
          do ist = 1, st%nst

            inhpsi = M_z0
            call states_get_state(stphi, gr%mesh, ist, ik, psi)

            do iatom = 1, geo%natoms
              do idir = 1, gr%sb%dim
                call pert_init(pert, namespace, PERTURBATION_IONIC, gr, geo)
                call pert_setup_atom(pert, iatom)
                call pert_setup_dir(pert, idir)
                call zpert_apply(pert, namespace, gr, geo, hm, psolver, ik, psi(:, :), dvpsi(:, :, idir))
                dvpsi(:, :, idir) = - dvpsi(:, :, idir)
                inhpsi(1:gr%mesh%np, 1:stphi%d%dim) = inhpsi(1:gr%mesh%np, 1:stphi%d%dim) &
                  + st%occ(ist, ik)*post(idir, iatom)*dvpsi(1:gr%mesh%np, 1:stphi%d%dim, idir)
                call pert_end(pert)
              end do
            end do

            call states_set_state(inh, gr%mesh, ist, ik, inhpsi)

          end do
        end do

        call hamiltonian_set_inh(hm, inh)
        call states_end(inh)

        SAFE_DEALLOCATE_A(psi)
        SAFE_DEALLOCATE_A(inhpsi)
        SAFE_DEALLOCATE_A(dvpsi)

      end if
    end subroutine prepare_inh

    subroutine update_state(epsilon)
      FLOAT, intent(in) :: epsilon

      do ik = stphi%d%kpt%start, stphi%d%kpt%end
        do ib = stphi%group%block_start, stphi%group%block_end
          call batch_axpy(gr%mesh%np, -M_zI*dt*epsilon, hst%group%psib(ib, ik), st%group%psib(ib, ik))
          if(propagate_chi) then
            call batch_axpy(gr%mesh%np, -M_zI*dt*epsilon, hchi%group%psib(ib, ik), chi%group%psib(ib, ik))
          end if
        end do
      end do
      
      if(ion_dynamics_ions_move(ions)) then
        posfinal = posfinal + posk * epsilon
        velfinal = velfinal + velk * epsilon
        if(propagate_chi) then
          posfinalt = posfinalt + poskt * epsilon
          velfinalt = velfinalt + velkt * epsilon
        end if
      end if
    end subroutine update_state

  end subroutine td_explicit_runge_kutta4


  subroutine td_runge_kutta2(ks, namespace, hm, psolver, gr, st, tr, time, dt, ions, geo)
    type(v_ks_t), target,            intent(inout) :: ks
    type(namespace_t),               intent(in)    :: namespace
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(poisson_t),     target,     intent(in)    :: psolver
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    type(propagator_t),  target,     intent(inout) :: tr
    FLOAT,                           intent(in)    :: time
    FLOAT,                           intent(in)    :: dt
    type(ion_dynamics_t),            intent(inout) :: ions
    type(geometry_t),                intent(inout) :: geo

    integer :: np_part, np, kp1, kp2, st1, st2, nspin, ik, ist, idim, j, ip
    integer :: i
    FLOAT :: dres
    CMPLX, allocatable :: zphi(:, :, :, :)
    CMPLX, allocatable :: zpsi(:), rhs(:)
    CMPLX, allocatable :: k2(:, :, :, :), oldk2(:, :, :, :), rhs1(:, :, :, :)
    CMPLX, allocatable :: inhpsi(:)
    type(ion_state_t) :: ions_state

    PUSH_SUB(td_runge_kutta2)

    st1 = st%st_start
    st2 = st%st_end
    kp1 = st%d%kpt%start
    kp2 = st%d%kpt%end
    np_part = gr%mesh%np_part
    np = gr%mesh%np
    nspin = hm%d%nspin
    move_ions_op = ion_dynamics_ions_move(ions)

    sp_np = np
    sp_dim = st%d%dim
    sp_st1 = st1
    sp_st2 = st2
    sp_kp1 = kp1
    sp_kp2 = kp2
    sp_parallel = st%parallel_in_states .or. st%d%kpt%parallel
    if(sp_parallel) sp_comm = st%st_kpt_mpi_grp%comm

    ! define pointer and variables for usage in td_rk2op, td_rk2opt routines
    grid_p    => gr
    hm_p      => hm
    psolver_p => psolver
    tr_p      => tr
    st_p      => st
    dt_op = dt
    t_op  = time - dt/M_TWO
    dim_op = st%d%dim
    xc_p      => ks%xc

    SAFE_ALLOCATE(k2(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
    SAFE_ALLOCATE(oldk2(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
    SAFE_ALLOCATE(zphi(1:gr%mesh%np_part, st%d%dim, st1:st2, kp1:kp2))
    SAFE_ALLOCATE(rhs1(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
    SAFE_ALLOCATE(rhs(1:tr%tdsk_size))
    SAFE_ALLOCATE(zpsi(1:tr%tdsk_size))
    SAFE_ALLOCATE(vhxc1_op(1:np, 1:nspin))
    SAFE_ALLOCATE(vpsl1_op(1:np))

    ! First, we get the state that we want to propagate. For the moment being, only one state.
    do ik = kp1, kp2
      do ist = st1, st2
        call states_get_state(st, gr%mesh, ist, ik, zphi(:, :, ist, ik))
      end do
    end do

    if (oct_exchange_enabled(hm%oct_exchange)) then
      call oct_exchange_prepare(hm%oct_exchange, gr%mesh, zphi, ks%xc, psolver)
    end if
    call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = time-dt)
    call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy)
    rhs1 = M_z0
    do ik = kp1, kp2
      do ist = st1, st2
        call zhamiltonian_apply(hm_p, grid_p%der, psolver, zphi(:, :, ist, ik), rhs1(:, :, ist, ik), ist, ik)
      end do
    end do
    do ik = kp1, kp2
      do ist = st1, st2
        if(oct_exchange_enabled(hm%oct_exchange)) then
          call zoct_exchange_operator(hm%oct_exchange, gr%der, rhs1(:, :, ist, ik), ist, ik)
        end if
      end do
    end do

    rhs1 = zphi - M_zI * M_HALF * dt * rhs1

    if(hamiltonian_inh_term(hm)) then
      SAFE_ALLOCATE(inhpsi(1:gr%mesh%np))
      do ik = kp1, kp2
        do ist = st1, st2
          do idim = 1, st%d%dim
            call states_get_state(hm%inh_st, gr%mesh, idim, ist, ik, inhpsi)
            forall(ip = 1:gr%mesh%np) rhs1(ip, idim, ist, ik) = rhs1(ip, idim, ist, ik) + dt * inhpsi(ip)
          end do
        end do
      end do
      SAFE_DEALLOCATE_A(inhpsi)
    end if

    k2 = zphi

    i = 1
    do
      oldk2 = k2

      ! Set the Hamiltonian at the final time of the propagation
      if(.not.oct_exchange_enabled(hm_p%oct_exchange)) then
        do ik = kp1, kp2
          do ist = st1, st2
            call states_set_state(st, gr%mesh, ist, ik, k2(:, :, ist, ik))
          end do
        end do
        call density_calc(st, gr, st%rho)
        call v_ks_calc(ks, namespace, hm, st, geo, calc_current = gauge_field_is_applied(hm%ep%gfield))
      end if
      if(ion_dynamics_ions_move(ions)) then
        call ion_dynamics_save_state(ions, geo, ions_state)
        call ion_dynamics_propagate(ions, gr%sb, geo, time, dt)
        call hamiltonian_epot_generate(hm, namespace,  gr, geo, st, psolver, time = time)
        vpsl1_op = hm%ep%vpsl
      end if
      call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = time)
      call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy)
      if(.not.oct_exchange_enabled(hm_p%oct_exchange)) then
        if (i==1) then
          call potential_interpolation_get(tr%vksold, gr%mesh%np, st%d%nspin, 0, vhxc1_op)
          i = i + 1
        else
          vhxc1_op = hm%vhxc
        end if
        t_op  = time
      else
        vhxc1_op = hm%vhxc
      end if

      if(ion_dynamics_ions_move(ions)) then
        call ion_dynamics_restore_state(ions, geo, ions_state)
      end if

      j = 1
      do ik = kp1, kp2
        do ist = st1, st2
          do idim = 1, st%d%dim
            rhs(j:j+np-1) = rhs1(1:np, idim, ist, ik)
            j = j + np
          end do
        end do
      end do

      ! Now we populate an initial guess for the equation.
      j = 1
      do ik = kp1, kp2
        do ist = st1, st2
          do idim = 1, st%d%dim
            zpsi(j:j+np-1) = k2(1:np, idim, ist, ik)
            j = j + np
          end do
        end do
      end do

      t_op  = time - dt
      call zsparskit_solver_run(tr%tdsk, td_rk2op, td_rk2opt, zpsi, rhs)
      
      k2 = M_z0
      j = 1
      do ik = kp1, kp2
        do ist = st1, st2
          do idim = 1, st%d%dim
            k2(1:np, idim, ist, ik) = zpsi(j:j+np-1)
            j = j + np
          end do
        end do
      end do

      dres = M_ZERO
      do ik = kp1, kp2
        do ist = st1, st2
          do idim = 1, st%d%dim
            dres = dres + zmf_nrm2(gr%mesh, k2(:, idim, ist, ik) - oldk2(:, idim, ist, ik), reduce = .false.)
          end do
        end do
      end do

      call comm_allreduce(st%dom_st_kpt_mpi_grp%comm, dres)
      
      if(dres < tr%scf_threshold) exit
    end do

    zphi = k2
    do ik = kp1, kp2
      do ist = st1, st2
        call states_set_state(st, gr%mesh, ist, ik, zphi(:, :, ist, ik))
      end do
    end do

    call density_calc(st, gr, st%rho)

    SAFE_DEALLOCATE_A(k2)
    SAFE_DEALLOCATE_A(oldk2)
    SAFE_DEALLOCATE_A(zphi)
    SAFE_DEALLOCATE_A(rhs1)
    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(rhs)
    SAFE_DEALLOCATE_A(vhxc1_op)
    SAFE_DEALLOCATE_A(vpsl1_op)

    POP_SUB(td_runge_kutta2)
  end subroutine td_runge_kutta2

  !----------------------------------------------------------------------------

  subroutine td_runge_kutta4(ks, namespace, hm, psolver, gr, st, tr, time, dt, ions, geo)
    type(v_ks_t), target,            intent(inout) :: ks
    type(namespace_t),               intent(in)    :: namespace    
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(poisson_t),     target,     intent(in)    :: psolver
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    type(propagator_t),  target,     intent(inout) :: tr
    FLOAT,                           intent(in)    :: time
    FLOAT,                           intent(in)    :: dt
    type(ion_dynamics_t),            intent(inout) :: ions
    type(geometry_t),                intent(inout) :: geo

    integer :: np_part, np, idim, ip, ist, ik, j, kp1, kp2, st1, st2, nspin
    FLOAT :: dres
    CMPLX, allocatable :: inhpsi(:)
    CMPLX, allocatable :: zpsi(:), rhs(:)
    CMPLX, allocatable :: zphi(:, :, :, :)
    type(ion_state_t) :: ions_state

    FLOAT :: a(2, 2), c(2), b(2)
    CMPLX, allocatable :: k1(:, :, :, :), k2(:, :, :, :), oldk1(:, :, :, :), &
      oldk2(:, :, :, :), yn1(:, :, :, :), yn2(:, :, :, :), &
      rhs1(:, :, :, :), rhs2(:, :, :, :)

    PUSH_SUB(td_runge_kutta4)
    
    c(1) = M_HALF - sqrt(M_THREE)/CNST(6.0)
    c(2) = M_HALF + sqrt(M_THREE)/CNST(6.0)
    b(1) = M_HALF
    b(2) = M_HALF
    a(1, 1) = M_FOURTH
    a(1, 2) = M_FOURTH - sqrt(M_THREE)/CNST(6.0)
    a(2, 1) = M_FOURTH + sqrt(M_THREE)/CNST(6.0)
    a(2, 2) = M_FOURTH

    st1 = st%st_start
    st2 = st%st_end
    kp1 = st%d%kpt%start
    kp2 = st%d%kpt%end
    np_part = gr%mesh%np_part
    np = gr%mesh%np
    nspin = hm%d%nspin
    move_ions_op = ion_dynamics_ions_move(ions)

    sp_np = np
    sp_dim = st%d%dim
    sp_st1 = st1
    sp_st2 = st2
    sp_kp1 = kp1
    sp_kp2 = kp2
    sp_parallel = st%parallel_in_states .or. st%d%kpt%parallel
    if(sp_parallel) sp_comm = st%st_kpt_mpi_grp%comm

    ! define pointer and variables for usage in td_rk4op, td_rk4opt routines
    grid_p    => gr
    hm_p      => hm
    psolver_p => psolver
    tr_p      => tr
    st_p      => st
    dt_op = dt
    t_op  = time - dt/M_TWO
    dim_op = st%d%dim

    SAFE_ALLOCATE(vhxc1_op(1:np, 1:nspin))
    SAFE_ALLOCATE(vhxc2_op(1:np, 1:nspin))
    SAFE_ALLOCATE(vpsl1_op(1:np))
    SAFE_ALLOCATE(vpsl2_op(1:np))
    SAFE_ALLOCATE(k1(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
    SAFE_ALLOCATE(k2(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
    SAFE_ALLOCATE(oldk1(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
    SAFE_ALLOCATE(oldk2(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
    SAFE_ALLOCATE(yn1(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
    SAFE_ALLOCATE(yn2(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
    SAFE_ALLOCATE(rhs1(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
    SAFE_ALLOCATE(rhs2(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
    SAFE_ALLOCATE(rhs(1:tr%tdsk_size))
    SAFE_ALLOCATE(zpsi(1:tr%tdsk_size))
    SAFE_ALLOCATE(zphi(1:gr%mesh%np_part, st%d%dim, st1:st2, kp1:kp2))

    ! First, we get the state that we want to propagate. For the moment being, only one state.
    do ik = kp1, kp2
      do ist = st1, st2
        call states_get_state(st, gr%mesh, ist, ik, zphi(:, :, ist, ik))
      end do
    end do
    k1 = M_z0
    k2 = M_z0

    do

      oldk1 = k1
      oldk2 = k2

      yn1 = zphi + a(1, 1) * k1 + a(1, 2) * k2
      yn2 = zphi + a(2, 1) * k1 + a(2, 2) * k2

      ! Set the Hamiltonian at time-dt + c(1) * dt
      do ik = kp1, kp2
        do ist = st1, st2
          call states_set_state(st, gr%mesh, ist, ik, yn1(:, :, ist, ik))
        end do
      end do
      call density_calc(st, gr, st%rho)
      call v_ks_calc(ks, namespace, hm, st, geo, calc_current = gauge_field_is_applied(hm%ep%gfield))
      if(ion_dynamics_ions_move(ions)) then
        call ion_dynamics_save_state(ions, geo, ions_state)
        call ion_dynamics_propagate(ions, gr%sb, geo, time - dt + c(1)*dt, c(1)*dt)
        call hamiltonian_epot_generate(hm, namespace,  gr, geo, st, psolver, time = time - dt + c(1)*dt)
        vpsl1_op = hm%ep%vpsl
      end if
      call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = time - dt + c(1)*dt)
      call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy)
      vhxc1_op = hm%vhxc
      t_op  = time - dt + c(1) * dt
      rhs1 = M_z0
      do ik = kp1, kp2
        do ist = st1, st2
          call zhamiltonian_apply(hm_p, grid_p%der, psolver_p, zphi(:, :, ist, ik), rhs1(:, :, ist, ik), ist, ik)
          if(hamiltonian_inh_term(hm)) then
            SAFE_ALLOCATE(inhpsi(1:gr%mesh%np))
            do idim = 1, st%d%dim
              call states_get_state(hm%inh_st, gr%mesh, idim, ist, ik, inhpsi)
              forall(ip = 1:gr%mesh%np) rhs1(ip, idim, ist, ik) = rhs1(ip, idim, ist, ik) + M_zI * inhpsi(ip)
            end do
            SAFE_DEALLOCATE_A(inhpsi)
          end if
        end do
      end do
      rhs1 = - M_zI * dt * rhs1
      if(ion_dynamics_ions_move(ions)) then
        call ion_dynamics_restore_state(ions, geo, ions_state)
      end if

      ! Set the Hamiltonian at time-dt + c(2) * dt
      do ik = kp1, kp2
        do ist = st1, st2
          call states_set_state(st, gr%mesh, ist, ik, yn2(:, :, ist, ik))
        end do
      end do
      call density_calc(st, gr, st%rho)
      call v_ks_calc(ks, namespace, hm, st, geo, calc_current = gauge_field_is_applied(hm%ep%gfield))
      if(ion_dynamics_ions_move(ions)) then
        call ion_dynamics_save_state(ions, geo, ions_state)
        call ion_dynamics_propagate(ions, gr%sb, geo, time - dt + c(2)*dt, c(2)*dt)
        call hamiltonian_epot_generate(hm, namespace, gr, geo, st, psolver, time = time - dt + c(2)*dt)
        vpsl2_op = hm%ep%vpsl
      end if
      call hamiltonian_update(hm, gr%mesh, gr%der%boundaries, time = time - dt + c(2)*dt)
      call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy)
      vhxc2_op = hm%vhxc
      t_op  = time - dt + c(2) * dt
      rhs2 = M_z0
      do ik = kp1, kp2
        do ist = st1, st2
          call zhamiltonian_apply(hm_p, grid_p%der, psolver_p, zphi(:, :, ist, ik), rhs2(:, :, ist, ik), ist, ik)
          if(hamiltonian_inh_term(hm)) then
            SAFE_ALLOCATE(inhpsi(1:gr%mesh%np))
            do idim = 1, st%d%dim
              call states_get_state(hm%inh_st, gr%mesh, idim, ist, ik, inhpsi)
              forall(ip = 1:gr%mesh%np) rhs2(ip, idim, ist, ik) = rhs2(ip, idim, ist, ik) + M_zI * inhpsi(ip)
            end do
            SAFE_DEALLOCATE_A(inhpsi)
          end if
        end do
      end do
      rhs2 = -M_zI * dt * rhs2
      if(ion_dynamics_ions_move(ions)) then
        call ion_dynamics_restore_state(ions, geo, ions_state)
      end if

      j = 1
      do ik = kp1, kp2
        do ist = st1, st2
          do idim = 1, st%d%dim
            rhs(j:j+np-1) = rhs1(1:np, idim, ist, ik)
            j = j + np
          end do
        end do
      end do
      do ik = kp1, kp2
        do ist = st1, st2
          do idim = 1, st%d%dim
            rhs(j:j+np-1) = rhs2(1:np, idim, ist, ik)
            j = j + np
          end do
        end do
      end do

      ! Now we populate an initial guess for the equation.
      j = 1
      do ik = kp1, kp2
        do ist = st1, st2
          do idim = 1, st%d%dim
            zpsi(j:j+np-1) = k1(1:np, idim, ist, ik)
            j = j + np
          end do
        end do
      end do
      do ik = kp1, kp2
        do ist = st1, st2
          do idim = 1, st%d%dim
            zpsi(j:j+np-1) = k2(1:np, idim, ist, ik)
            j = j + np
          end do
        end do
      end do

      t_op  = time - dt
      call zsparskit_solver_run(tr%tdsk, td_rk4op, td_rk4opt, zpsi, rhs)
      
      k1 = M_z0
      k2 = M_z0
      j = 1
      do ik = kp1, kp2
        do ist = st1, st2
          do idim = 1, st%d%dim
            k1(1:np, idim, ist, ik) = zpsi(j:j+np-1)
            j = j + np
          end do
        end do
      end do
      do ik = kp1, kp2
        do ist = st1, st2
          do idim = 1, st%d%dim
            k2(1:np, idim, ist, ik) = zpsi(j:j+np-1)
            j = j + np
          end do
        end do
      end do

      dres = M_ZERO
      do ik = kp1, kp2
        do ist = st1, st2
          do idim = 1, st%d%dim
            dres = dres + zmf_nrm2(gr%mesh, k1(:, idim, ist, ik) - oldk1(:, idim, ist, ik))
            dres = dres + zmf_nrm2(gr%mesh, k2(:, idim, ist, ik) - oldk2(:, idim, ist, ik))
          end do
        end do
      end do
      if(sp_parallel) call comm_allreduce(sp_comm, dres)
      !write(*, *) 'Residual = ', dres

      if(dres < tr%scf_threshold) exit
    end do


    zphi = zphi + b(1) * k1 + b(2) * k2
    do ik = kp1, kp2
      do ist = st1, st2
        call states_set_state(st, gr%mesh, ist, ik, zphi(:, :, ist, ik))
      end do
    end do

    call density_calc(st, gr, st%rho)

    SAFE_DEALLOCATE_A(rhs1)
    SAFE_DEALLOCATE_A(rhs2)
    SAFE_DEALLOCATE_A(k1)
    SAFE_DEALLOCATE_A(k2)
    SAFE_DEALLOCATE_A(oldk1)
    SAFE_DEALLOCATE_A(oldk2)
    SAFE_DEALLOCATE_A(yn1)
    SAFE_DEALLOCATE_A(yn2)
    SAFE_DEALLOCATE_A(vhxc1_op)
    SAFE_DEALLOCATE_A(vhxc2_op)
    SAFE_DEALLOCATE_A(vpsl1_op)
    SAFE_DEALLOCATE_A(vpsl2_op)
    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(rhs)

    POP_SUB(td_runge_kutta4)
  end subroutine td_runge_kutta4

  ! ---------------------------------------------------------
  !> operators for Crank-Nicolson scheme
  subroutine td_rk4op(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:)
    FLOAT, intent(in)  :: xim(:)
    FLOAT, intent(out) :: yre(:)
    FLOAT, intent(out) :: yim(:)

    integer :: idim, j, ik, ist, kp1, kp2, st1, st2, dim, k, jj
    CMPLX, allocatable :: zpsi(:, :)
    CMPLX, allocatable :: opzpsi(:, :)
    FLOAT :: a(2, 2), c(2)
    integer :: np_part, np

    PUSH_SUB(td_rk4op)

    np_part = grid_p%mesh%np_part
    np = grid_p%mesh%np
    st1 = st_p%st_start
    st2 = st_p%st_end
    kp1 = st_p%d%kpt%start
    kp2 = st_p%d%kpt%end
    dim = st_p%d%dim

    SAFE_ALLOCATE(zpsi(1:np_part, 1:dim))
    SAFE_ALLOCATE(opzpsi(1:np_part, 1:dim))

    a(1, 1) = M_FOURTH
    a(1, 2) = M_FOURTH - sqrt(M_THREE)/CNST(6.0)
    a(2, 1) = M_FOURTH + sqrt(M_THREE)/CNST(6.0)
    a(2, 2) = M_FOURTH

    c(1) = M_HALF - sqrt(M_THREE)/CNST(6.0)
    c(2) = M_HALF + sqrt(M_THREE)/CNST(6.0)

    zpsi = M_z0

    hm_p%vhxc = vhxc1_op
    if(move_ions_op) hm_p%ep%vpsl = vpsl1_op
    call hamiltonian_update(hm_p, grid_p%mesh, grid_p%der%boundaries, time = t_op + c(1)*dt_op)
    call lda_u_update_occ_matrices(hm_p%lda_u, grid_p%mesh, st_p, hm_p%hm_base, hm_p%energy)
    j = 1
    k = np * (kp2 - kp1 + 1) * (st2 - st1 + 1) * dim + 1
    do ik = kp1, kp2
      do ist = st1, st2
        jj = j
        do idim = 1, dim
          zpsi(1:np, idim) = a(1, 1) * cmplx(xre(j:j+np-1), xim(j:j+np-1), REAL_PRECISION) + &
            a(1, 2) * cmplx(xre(k:k+np-1), xim(k:k+np-1), REAL_PRECISION)
          j = j + np
          k = k + np
        end do

        call zhamiltonian_apply(hm_p, grid_p%der, psolver_p, zpsi, opzpsi, ist, ik)

        do idim = 1, dim
          yre(jj:jj+np-1) = xre(jj:jj+np-1) + real(M_zI * dt_op * opzpsi(1:np, idim))
          yim(jj:jj+np-1) = xim(jj:jj+np-1) + aimag(M_zI * dt_op * opzpsi(1:np, idim))
          jj = jj + np
        end do
      end do
    end do

    hm_p%vhxc = vhxc2_op
    if(move_ions_op) hm_p%ep%vpsl = vpsl2_op
    call hamiltonian_update(hm_p, grid_p%mesh, grid_p%der%boundaries, time = t_op + c(2)*dt_op)
    call lda_u_update_occ_matrices(hm_p%lda_u, grid_p%mesh, st_p, hm_p%hm_base, hm_p%energy)
    j = 1
    k = np * (kp2 - kp1 + 1) * (st2 - st1 + 1) * dim + 1
    do ik = kp1, kp2
      do ist = st1, st2
        jj = k
        do idim = 1, dim
          zpsi(1:np, idim) = a(2, 1) * cmplx(xre(j:j+np-1), xim(j:j+np-1), REAL_PRECISION) + &
            a(2, 2) * cmplx(xre(k:k+np-1), xim(k:k+np-1), REAL_PRECISION)
          j = j + np
          k = k + np
        end do

        call zhamiltonian_apply(hm_p, grid_p%der, psolver_p, zpsi, opzpsi, ist, ik)

        do idim = 1, dim
          yre(jj:jj+np-1) = xre(jj:jj+np-1) + real(M_zI * dt_op * opzpsi(1:np, idim))
          yim(jj:jj+np-1) = xim(jj:jj+np-1) + aimag(M_zI * dt_op * opzpsi(1:np, idim))
          jj = jj + np
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(opzpsi)
    POP_SUB(td_rk4op)
  end subroutine td_rk4op
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Transpose of H (called e.g. by bi-conjugate gradient solver)
  subroutine td_rk4opt(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:)
    FLOAT, intent(in)  :: xim(:)
    FLOAT, intent(out) :: yre(:)
    FLOAT, intent(out) :: yim(:)

    integer :: idim, j, ik, ist, kp1, kp2, st1, st2, dim, k, jj
    CMPLX, allocatable :: zpsi(:, :)
    CMPLX, allocatable :: opzpsi(:, :)
    FLOAT :: a(2, 2), c(2)
    integer :: np_part, np

    PUSH_SUB(td_rk4opt)

    np_part = grid_p%mesh%np_part
    np = grid_p%mesh%np
    st1 = st_p%st_start
    st2 = st_p%st_end
    kp1 = st_p%d%kpt%start
    kp2 = st_p%d%kpt%end
    dim = st_p%d%dim

    SAFE_ALLOCATE(zpsi(1:np_part, 1:dim))
    SAFE_ALLOCATE(opzpsi(1:np_part, 1:dim))

    a(1, 1) = M_FOURTH
    a(1, 2) = M_FOURTH - sqrt(M_THREE)/CNST(6.0)
    a(2, 1) = M_FOURTH + sqrt(M_THREE)/CNST(6.0)
    a(2, 2) = M_FOURTH

    c(1) = M_HALF - sqrt(M_THREE)/CNST(6.0)
    c(2) = M_HALF + sqrt(M_THREE)/CNST(6.0)

    zpsi = M_z0

    hm_p%vhxc = vhxc1_op
    if(move_ions_op) hm_p%ep%vpsl = vpsl1_op
    call hamiltonian_update(hm_p, grid_p%mesh, grid_p%der%boundaries, time = t_op + c(1)*dt_op)
    call lda_u_update_occ_matrices(hm_p%lda_u, grid_p%mesh, st_p, hm_p%hm_base, hm_p%energy)
    j = 1
    k = np * (kp2 - kp1 + 1) * (st2 - st1 + 1) * dim + 1
    do ik = kp1, kp2
      do ist = st1, st2
        jj = j
        do idim = 1, dim
          zpsi(1:np, idim) = a(1, 1) * cmplx(xre(j:j+np-1), -xim(j:j+np-1), REAL_PRECISION) + &
            a(1, 2) * cmplx(xre(k:k+np-1), -xim(k:k+np-1), REAL_PRECISION)
          j = j + np
          k = k + np
        end do

        call zhamiltonian_apply(hm_p, grid_p%der, psolver_p, zpsi, opzpsi, ist, ik)

        do idim = 1, dim
          yre(jj:jj+np-1) = xre(jj:jj+np-1) + real(M_zI * dt_op * opzpsi(1:np, idim))
          yim(jj:jj+np-1) = xim(jj:jj+np-1) - aimag(M_zI * dt_op * opzpsi(1:np, idim))
          jj = jj + np
        end do
      end do
    end do

    hm_p%vhxc = vhxc2_op
    if(move_ions_op) hm_p%ep%vpsl = vpsl2_op
    call hamiltonian_update(hm_p, grid_p%mesh, grid_p%der%boundaries, time = t_op + c(2)*dt_op)
    call lda_u_update_occ_matrices(hm_p%lda_u, grid_p%mesh, st_p, hm_p%hm_base, hm_p%energy)
    j = 1
    k = np * (kp2 - kp1 + 1) * (st2 - st1 + 1) * dim + 1
    do ik = kp1, kp2
      do ist = st1, st2
        jj = k
        do idim = 1, dim
          zpsi(1:np, idim) = a(2, 1) * cmplx(xre(j:j+np-1), -xim(j:j+np-1), REAL_PRECISION) + &
            a(2, 2) * cmplx(xre(k:k+np-1), -xim(k:k+np-1), REAL_PRECISION)
          j = j + np
          k = k + np
        end do

        call zhamiltonian_apply(hm_p, grid_p%der, psolver_p, zpsi, opzpsi, ist, ik)

        do idim = 1, dim
          yre(jj:jj+np-1) = xre(jj:jj+np-1) + real(M_zI * dt_op * opzpsi(1:np, idim))
          yim(jj:jj+np-1) = xim(jj:jj+np-1) - aimag(M_zI * dt_op * opzpsi(1:np, idim))
          jj = jj + np
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(opzpsi)
    POP_SUB(td_rk4opt)
  end subroutine td_rk4opt
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> operator for the RK2 propagator
  subroutine td_rk2op(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:)
    FLOAT, intent(in)  :: xim(:)
    FLOAT, intent(out) :: yre(:)
    FLOAT, intent(out) :: yim(:)

    integer :: np_part, np, st1, st2, kp1, kp2, dim, idim, ik, ist, jj, j
    CMPLX, allocatable :: zpsi(:, :)
    CMPLX, allocatable :: opzpsi(:, :)
    CMPLX, allocatable :: zpsi_(:, :, :, :)

    PUSH_SUB(td_rk2op)

    np_part = grid_p%mesh%np_part
    np = grid_p%mesh%np
    st1 = st_p%st_start
    st2 = st_p%st_end
    kp1 = st_p%d%kpt%start
    kp2 = st_p%d%kpt%end
    dim = st_p%d%dim

    SAFE_ALLOCATE(zpsi(1:np_part, 1:dim))
    SAFE_ALLOCATE(opzpsi(1:np_part, 1:dim))
    SAFE_ALLOCATE(zpsi_(1:np_part, 1:dim, st1:st2, kp1:kp2))

    zpsi = M_z0
    opzpsi = M_z0

    hm_p%vhxc = vhxc1_op
    if(move_ions_op) hm_p%ep%vpsl = vpsl1_op
    call hamiltonian_update(hm_p, grid_p%mesh, grid_p%der%boundaries, time = t_op + dt_op)
   call lda_u_update_occ_matrices(hm_p%lda_u, grid_p%mesh, st_p, hm_p%hm_base, hm_p%energy)

    if(oct_exchange_enabled(hm_p%oct_exchange)) then
      zpsi_ = M_z0
      j = 1
      do ik = kp1, kp2
        do ist = st1, st2
          jj = j
          do idim = 1, dim
            zpsi_(1:np, idim, ist, ik) = cmplx(xre(j:j+np-1), xim(j:j+np-1), REAL_PRECISION)
            j = j + np
          end do
        end do
      end do
      call oct_exchange_prepare(hm_p%oct_exchange, grid_p%mesh, zpsi_, xc_p, psolver_p)
    end if

    j = 1
    do ik = kp1, kp2
      do ist = st1, st2
        jj = j
        do idim = 1, dim
          zpsi(1:np, idim) = cmplx(xre(j:j+np-1), xim(j:j+np-1), REAL_PRECISION)
          j = j + np
        end do
        call zhamiltonian_apply(hm_p, grid_p%der, psolver_p, zpsi, opzpsi, ist, ik)
        do idim = 1, dim
          yre(jj:jj+np-1) = xre(jj:jj+np-1) + real(M_zI * dt_op * M_HALF * opzpsi(1:np, idim))
          yim(jj:jj+np-1) = xim(jj:jj+np-1) + aimag(M_zI * dt_op * M_HALF * opzpsi(1:np, idim))
          jj = jj + np
        end do
      end do
    end do

    if(oct_exchange_enabled(hm_p%oct_exchange)) then
      j = 1
      do ik = kp1, kp2
        do ist = st1, st2
          jj = j
          do idim = 1, dim
            zpsi(1:np, idim) = cmplx(xre(j:j+np-1), xim(j:j+np-1), REAL_PRECISION)
            j = j + np
          end do
          opzpsi = M_z0
          call zoct_exchange_operator(hm_p%oct_exchange, grid_p%der, opzpsi, ist, ik)

          do idim = 1, dim
            yre(jj:jj+np-1) = yre(jj:jj+np-1) + real(M_zI * dt_op * M_HALF * opzpsi(1:np, idim))
            yim(jj:jj+np-1) = yim(jj:jj+np-1) + aimag(M_zI * dt_op * M_HALF * opzpsi(1:np, idim))
            jj = jj + np
          end do
        end do
      end do
    end if

    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(opzpsi)
    
    POP_SUB(td_rk2op)
  end subroutine td_rk2op
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> operator for the RK2 propagator
  subroutine td_rk2opt(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:)
    FLOAT, intent(in)  :: xim(:)
    FLOAT, intent(out) :: yre(:)
    FLOAT, intent(out) :: yim(:)

    integer :: np_part, np, st1, st2, kp1, kp2, dim, idim, ik, ist, jj, j
    CMPLX, allocatable :: zpsi(:, :)
    CMPLX, allocatable :: opzpsi(:, :)
    CMPLX, allocatable :: zpsi_(:, :, :, :)

    PUSH_SUB(td_rk2opt)

    np_part = grid_p%mesh%np_part
    np = grid_p%mesh%np
    st1 = st_p%st_start
    st2 = st_p%st_end
    kp1 = st_p%d%kpt%start
    kp2 = st_p%d%kpt%end
    dim = st_p%d%dim

    SAFE_ALLOCATE(zpsi(1:np_part, 1:dim))
    SAFE_ALLOCATE(opzpsi(1:np_part, 1:dim))
    SAFE_ALLOCATE(zpsi_(1:np_part, 1:dim, st1:st2, kp1:kp2))

    zpsi = M_z0
    opzpsi = M_z0

    hm_p%vhxc = vhxc1_op
    if(move_ions_op) hm_p%ep%vpsl = vpsl1_op
    call hamiltonian_update(hm_p, grid_p%mesh, grid_p%der%boundaries, time = t_op + dt_op)
    call lda_u_update_occ_matrices(hm_p%lda_u, grid_p%mesh, st_p, hm_p%hm_base, hm_p%energy)

    if(oct_exchange_enabled(hm_p%oct_exchange)) then
      zpsi_ = M_z0
      j = 1
      do ik = kp1, kp2
        do ist = st1, st2
          jj = j
          do idim = 1, dim
            zpsi_(1:np, idim, ist, ik) = cmplx(xre(j:j+np-1), -xim(j:j+np-1), REAL_PRECISION)
            j = j + np
          end do
        end do
      end do
      call oct_exchange_prepare(hm_p%oct_exchange, grid_p%mesh, zpsi_, xc_p, psolver_p)
    end if

    j = 1
    do ik = kp1, kp2
      do ist = st1, st2
        jj = j
        do idim = 1, dim
          zpsi(1:np, idim) = cmplx(xre(j:j+np-1), -xim(j:j+np-1), REAL_PRECISION)
          j = j + np
        end do
        call zhamiltonian_apply(hm_p, grid_p%der, psolver_p, zpsi, opzpsi, ist, ik)

        do idim = 1, dim
          yre(jj:jj+np-1) = xre(jj:jj+np-1) + real(M_zI * dt_op * M_HALF * opzpsi(1:np, idim))
          yim(jj:jj+np-1) = xim(jj:jj+np-1) - aimag(M_zI * dt_op * M_HALF * opzpsi(1:np, idim))
          jj = jj + np
        end do
      end do
    end do

    if(oct_exchange_enabled(hm_p%oct_exchange)) then
      j = 1
      do ik = kp1, kp2
        do ist = st1, st2
          jj = j
          do idim = 1, dim
            zpsi(1:np, idim) = cmplx(xre(j:j+np-1), xim(j:j+np-1), REAL_PRECISION)
            j = j + np
          end do
          opzpsi = M_z0
          call zoct_exchange_operator(hm_p%oct_exchange, grid_p%der, opzpsi, ist, ik)

          do idim = 1, dim
            yre(jj:jj+np-1) = yre(jj:jj+np-1) + real(M_zI * dt_op * M_HALF * opzpsi(1:np, idim))
            yim(jj:jj+np-1) = yim(jj:jj+np-1) - aimag(M_zI * dt_op * M_HALF * opzpsi(1:np, idim))
            jj = jj + np
          end do
        end do
      end do
    end if

    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(opzpsi)
    PUSH_SUB(td_rk2opt)
  end subroutine td_rk2opt
  ! ---------------------------------------------------------

end module propagator_rk_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

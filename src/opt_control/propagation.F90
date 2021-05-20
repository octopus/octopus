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

module propagation_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use boundary_op_oct_m
  use controlfunction_oct_m
  use density_oct_m
  use electrons_oct_m
  use energy_calc_oct_m
  use epot_oct_m
  use excited_states_oct_m
  use forces_oct_m
  use gauge_field_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use ion_dynamics_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use lasers_oct_m
  use loct_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use oct_exchange_oct_m
  use opt_control_state_oct_m
  use parser_oct_m
  use pert_oct_m
  use propagator_elec_oct_m
  use propagator_base_oct_m
  use profiling_oct_m
  use restart_oct_m
  use space_oct_m
  use species_oct_m
  use states_elec_oct_m
  use states_elec_restart_oct_m
  use target_oct_m
  use td_oct_m
  use td_write_oct_m
  use v_ks_oct_m

  implicit none

  private

  public :: propagation_mod_init, &
            propagate_forward,    &
            propagate_backward,   &
            fwd_step,             &
            bwd_step,             &
            bwd_step_2,           &
            oct_prop_t,           &
            oct_prop_init,        &
            oct_prop_check,       &
            oct_prop_end


  type oct_prop_t
    private
    integer :: number_checkpoints
    integer, allocatable :: iter(:)
    integer :: niter
    character(len=100) :: dirname
    type(restart_t) :: restart_load
    type(restart_t) :: restart_dump
  end type oct_prop_t


  !> Module variables 
  integer :: niter_
  integer :: number_checkpoints_
  FLOAT   :: eta_
  FLOAT   :: delta_
  logical :: zbr98_
  logical :: gradients_

contains

  !> This subroutine must be called before any QOCT propagations are
  !! done. It simply stores in the module some data that is needed for
  !! the propagations, and which should stay invariant during the whole
  !! run.
  !! There is no need for any propagation_mod_close.
  subroutine propagation_mod_init(niter, eta, delta, number_checkpoints, zbr98, gradients)
    integer, intent(in) :: niter
    FLOAT,   intent(in) :: eta
    FLOAT,   intent(in) :: delta
    integer, intent(in) :: number_checkpoints
    logical, intent(in) :: zbr98
    logical, intent(in) :: gradients

    ASSERT(.not. (zbr98 .and. gradients))

    PUSH_SUB(propagation_mod_init)

    niter_              = niter
    eta_                = eta
    delta_              = delta
    number_checkpoints_ = number_checkpoints
    zbr98_              = zbr98
    gradients_          = gradients

    POP_SUB(propagation_mod_init)
  end subroutine propagation_mod_init
  ! ---------------------------------------------------------


  !> ---------------------------------------------------------
  !! Performs a full propagation of state psi, with the laser
  !! field specified in par. If write_iter is present and is
  !! set to .true., writes down through the td_write module.
  !! ---------------------------------------------------------
  subroutine propagate_forward(sys, td, par, tg, qcpsi, prop, write_iter)
    type(electrons_t),          intent(inout)  :: sys
    type(td_t),                 intent(inout)  :: td
    type(controlfunction_t),    intent(in)     :: par
    type(target_t),             intent(inout)  :: tg
    type(opt_control_state_t),  intent(inout)  :: qcpsi
    type(oct_prop_t), optional, intent(inout)  :: prop
    logical, optional,          intent(in)     :: write_iter

    integer :: ii, istep, ierr
    logical :: write_iter_ = .false.
    type(td_write_t)           :: write_handler
    FLOAT, allocatable :: x_initial(:,:)
    logical :: vel_target_ = .false.
    type(states_elec_t), pointer :: psi

    FLOAT :: init_time, final_time

    PUSH_SUB(propagate_forward)

    message(1) = "Info: Forward propagation."
    call messages_info(1)

    call controlfunction_to_h(par, sys%hm%ext_lasers)

    write_iter_ = .false.
    if(present(write_iter)) write_iter_ = write_iter

    psi => opt_control_point_qs(qcpsi)
    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call psi%pack()
    call opt_control_get_classical(sys%ions, qcpsi)

    if(write_iter_) then
      call td_write_init(write_handler, sys%namespace, sys%space, sys%outp, sys%gr, sys%st, sys%hm, sys%ions, sys%ks, &
        ion_dynamics_ions_move(td%ions_dyn), gauge_field_is_applied(sys%hm%ep%gfield), sys%hm%ep%kick, td%iter, td%max_iter, &
        td%dt, sys%mc)
      call td_write_data(write_handler)
    end if

    call hamiltonian_elec_not_adjoint(sys%hm)

    ! setup the Hamiltonian
    call density_calc(psi, sys%gr, psi%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, psi, sys%ions, time = M_ZERO)
    call propagator_elec_run_zero_iter(sys%hm, sys%gr, td%tr)
    if(ion_dynamics_ions_move(td%ions_dyn)) then
      call hamiltonian_elec_epot_generate(sys%hm, sys%namespace,  sys%space, sys%gr, sys%ions, psi, time = M_ZERO)
      call forces_calculate(sys%gr, sys%namespace, sys%ions, sys%hm, psi, sys%ks, t = M_ZERO, dt = td%dt)
    end if


    if(target_type(tg) == oct_tg_hhgnew) &
      call target_init_propagation(tg)
    
    if(target_type(tg) == oct_tg_velocity .or. target_type(tg) == oct_tg_hhgnew) then
      SAFE_ALLOCATE_SOURCE_A(x_initial, sys%ions%pos)
      vel_target_ = .true.
      sys%ions%vel = M_ZERO
      sys%ions%tot_force = M_ZERO
    end if

    if(.not.target_move_ions(tg)) call epot_precalc_local_potential(sys%hm%ep, sys%namespace, sys%gr, sys%ions)

    call target_tdcalc(tg, sys%namespace, sys%space, sys%hm, sys%gr, sys%ions, psi, 0, td%max_iter)

    if (present(prop)) then
      call oct_prop_dump_states(prop, sys%space, 0, psi, sys%gr%mesh, sys%kpoints, ierr)
      if (ierr /= 0) then
        message(1) = "Unable to write OCT states restart."
        call messages_warning(1)
      end if
    end if

    init_time = loct_clock()
    if (mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, td%max_iter)

    ii = 1
    do istep = 1, td%max_iter
      ! time-iterate wavefunctions

      call propagator_elec_dt(sys%ks, sys%namespace, sys%space, sys%hm, sys%gr, psi, td%tr, istep*td%dt, td%dt, td%mu, istep, &
        td%ions_dyn, sys%ions, sys%outp)

      if(present(prop)) then
        call oct_prop_dump_states(prop, sys%space, istep, psi, sys%gr%mesh, sys%kpoints, ierr)
        if (ierr /= 0) then
          message(1) = "Unable to write OCT states restart."
          call messages_warning(1)
        end if
      end if

      ! update
      call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, psi, sys%ions, time = istep*td%dt)
      call energy_calc_total(sys%namespace, sys%space, sys%hm, sys%gr, psi)

      if(sys%hm%bc%abtype == MASK_ABSORBING) call zvmask(sys%gr%mesh, sys%hm, psi)

      ! if td_target
      call target_tdcalc(tg, sys%namespace, sys%space, sys%hm, sys%gr, sys%ions, psi, istep, td%max_iter)

      ! only write in final run
      if(write_iter_) then
        call td_write_iter(write_handler, sys%namespace, sys%space, sys%outp, sys%gr, psi, sys%hm,  sys%ions, sys%hm%ep%kick, &
          td%dt, istep)
        ii = ii + 1
        if (any(ii == sys%outp%output_interval + 1) .or. istep == td%max_iter) then ! output
          if (istep == td%max_iter) sys%outp%output_interval = ii - 1
          ii = istep
          call td_write_data(write_handler) 
        end if
      end if

      if ((mod(istep, 100) == 0) .and. mpi_grp_is_root(mpi_world)) call loct_progress_bar(istep, td%max_iter)
    end do
    if(mpi_grp_is_root(mpi_world)) write(stdout, '(1x)')

    final_time = loct_clock()
    write(message(1),'(a,f12.2,a)') 'Propagation time: ', final_time - init_time, ' seconds.'
    call messages_info(1)

    if(vel_target_) then
      sys%ions%pos = x_initial
       SAFE_DEALLOCATE_A(x_initial)
    end if

    call opt_control_set_classical(sys%ions, qcpsi)

    if(write_iter_) call td_write_end(write_handler)
    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call psi%unpack()
    nullify(psi)
    POP_SUB(propagate_forward)
  end subroutine propagate_forward
  ! ---------------------------------------------------------


  !> ---------------------------------------------------------
  !! Performs a full backward propagation of state psi, with the
  !! external fields specified in Hamiltonian h.
  !! ---------------------------------------------------------
  subroutine propagate_backward(sys, td, qcpsi, prop)
    type(electrons_t),         intent(inout) :: sys
    type(td_t),                intent(inout) :: td
    type(opt_control_state_t), intent(inout) :: qcpsi
    type(oct_prop_t),          intent(inout) :: prop

    integer :: istep, ierr
    type(states_elec_t), pointer :: psi

    PUSH_SUB(propagate_backward)
    
    message(1) = "Info: Backward propagation."
    call messages_info(1)

    psi => opt_control_point_qs(qcpsi)
    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call psi%pack()

    call hamiltonian_elec_adjoint(sys%hm)

    ! setup the Hamiltonian
    call density_calc(psi, sys%gr, psi%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, psi, sys%ions)
    call propagator_elec_run_zero_iter(sys%hm, sys%gr, td%tr)

    call oct_prop_dump_states(prop, sys%space, td%max_iter, psi, sys%gr%mesh, sys%kpoints, ierr)
    if (ierr /= 0) then
      message(1) = "Unable to write OCT states restart."
      call messages_warning(1)
    end if
    
    if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, td%max_iter)

    do istep = td%max_iter, 1, -1
      call propagator_elec_dt(sys%ks, sys%namespace, sys%space, sys%hm, sys%gr, psi, td%tr, &
        (istep - 1)*td%dt, -td%dt, td%mu, istep-1, td%ions_dyn, sys%ions, sys%outp)

      call oct_prop_dump_states(prop, sys%space, istep - 1, psi, sys%gr%mesh, sys%kpoints, ierr)
      if (ierr /= 0) then
        message(1) = "Unable to write OCT states restart."
        call messages_warning(1)
      end if

      call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, psi, sys%ions)
      if(mod(istep, 100) == 0 .and. mpi_grp_is_root(mpi_world)) call loct_progress_bar(td%max_iter - istep + 1, td%max_iter)
    end do

    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call psi%unpack()
    nullify(psi)
    POP_SUB(propagate_backward)
  end subroutine propagate_backward
  ! ---------------------------------------------------------


  !> ---------------------------------------------------------
  !! Performs a forward propagation on the state psi and on the
  !! Lagrange-multiplier state chi. It also updates the control
  !! function par,  according to the following scheme:
  !! 
  !! |chi> --> U[par_chi](T, 0)|chi>
  !! par = par[|psi>, |chi>]
  !! |psi> --> U[par](T, 0)|psi>
  !!
  !! Note that the control functions "par" are updated on the
  !! fly, so that the propagation of psi is performed with the
  !! "new" control functions.
  !! --------------------------------------------------------
  subroutine fwd_step(sys, td, tg, par, par_chi, qcpsi, prop_chi, prop_psi)
    type(electrons_t),         intent(inout) :: sys
    type(td_t),                intent(inout) :: td
    type(target_t),            intent(inout) :: tg
    type(controlfunction_t),   intent(inout) :: par
    type(controlfunction_t),   intent(in)    :: par_chi
    type(opt_control_state_t), intent(inout) :: qcpsi
    type(oct_prop_t),          intent(inout) :: prop_chi
    type(oct_prop_t),          intent(inout) :: prop_psi

    integer :: i, ierr
    logical :: aux_fwd_propagation
    type(states_elec_t) :: psi2
    type(opt_control_state_t) :: qcchi
    type(controlfunction_t) :: par_prev
    type(propagator_base_t) :: tr_chi
    type(propagator_base_t) :: tr_psi2
    type(states_elec_t), pointer :: psi, chi

    PUSH_SUB(fwd_step)

    message(1) = "Info: Forward propagation."
    call messages_info(1)

    call controlfunction_to_realtime(par)

    call opt_control_state_null(qcchi)
    call opt_control_state_copy(qcchi, qcpsi)

    psi => opt_control_point_qs(qcpsi)
    chi => opt_control_point_qs(qcchi)
    call propagator_elec_copy(tr_chi, td%tr)
    ! The propagation of chi should not be self-consistent, because the Kohn-Sham
    ! potential used is the one created by psi. Note, however, that it is likely that
    ! the first two iterations are done self-consistently nonetheless.
    call propagator_elec_remove_scf_prop(tr_chi)

    aux_fwd_propagation = (target_mode(tg) == oct_targetmode_td .or. &
                           (sys%hm%theory_level /= INDEPENDENT_PARTICLES .and. &
                            .not.sys%ks%frozen_hxc))
    if(aux_fwd_propagation) then
      call states_elec_copy(psi2, psi)
      call controlfunction_copy(par_prev, par)
      if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call psi2%pack()
    end if

    ! setup forward propagation
    call density_calc(psi, sys%gr, psi%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, psi, sys%ions)
    call propagator_elec_run_zero_iter(sys%hm, sys%gr, td%tr)
    call propagator_elec_run_zero_iter(sys%hm, sys%gr, tr_chi)
    if(aux_fwd_propagation) then
      call propagator_elec_copy(tr_psi2, td%tr)
      call propagator_elec_run_zero_iter(sys%hm, sys%gr, tr_psi2)
    end if

    call oct_prop_dump_states(prop_psi, sys%space, 0, psi, sys%gr%mesh, sys%kpoints, ierr)
    if (ierr /= 0) then
      message(1) = "Unable to write OCT states restart."
      call messages_warning(1)
    end if
    call oct_prop_load_states(prop_chi, sys%namespace, sys%space, chi, sys%gr%mesh, sys%kpoints, 0, ierr)
    if (ierr /= 0) then
      message(1) = "Unable to read OCT states restart."
      call messages_fatal(1)
    end if

    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call psi%pack()
    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call chi%pack()

    do i = 1, td%max_iter
      call update_field(i, par, sys%gr, sys%hm, sys%ions, qcpsi, qcchi, par_chi, dir = 'f')
      call update_hamiltonian_elec_chi(i, sys%namespace, sys%gr, sys%ks, sys%hm, td, tg, par_chi, sys%ions, psi2)
      call hamiltonian_elec_update(sys%hm, sys%gr%mesh, sys%namespace, sys%space, time = (i - 1)*td%dt)
      call propagator_elec_dt(sys%ks, sys%namespace, sys%space, sys%hm, sys%gr, chi, tr_chi, i*td%dt, td%dt, td%mu, i, &
        td%ions_dyn, sys%ions, sys%outp)
      if(aux_fwd_propagation) then
        call update_hamiltonian_elec_psi(i, sys%namespace, sys%space, sys%gr, sys%ks, sys%hm, td, tg, par_prev, psi2, sys%ions)
        call propagator_elec_dt(sys%ks, sys%namespace, sys%space, sys%hm, sys%gr, psi2, tr_psi2, i*td%dt, td%dt, td%mu, i, &
          td%ions_dyn, sys%ions, sys%outp)
      end if
      call update_hamiltonian_elec_psi(i, sys%namespace, sys%space, sys%gr, sys%ks, sys%hm, td, tg, par, psi, sys%ions)
      call hamiltonian_elec_update(sys%hm, sys%gr%mesh, sys%namespace, sys%space, time = (i - 1)*td%dt)
      call propagator_elec_dt(sys%ks, sys%namespace, sys%space, sys%hm, sys%gr, psi, td%tr, i*td%dt, td%dt, td%mu, i, &
        td%ions_dyn, sys%ions, sys%outp)
      call target_tdcalc(tg, sys%namespace, sys%space, sys%hm, sys%gr, sys%ions, psi, i, td%max_iter) 

      call oct_prop_dump_states(prop_psi, sys%space, i, psi, sys%gr%mesh, sys%kpoints, ierr)
      if (ierr /= 0) then
        message(1) = "Unable to write OCT states restart."
        call messages_warning(1)
      end if
      call oct_prop_check(prop_chi, sys%namespace, sys%space, chi, sys%gr%mesh, sys%kpoints, i)
    end do
    call update_field(td%max_iter+1, par, sys%gr, sys%hm, sys%ions, qcpsi, qcchi, par_chi, dir = 'f')

    call density_calc(psi, sys%gr, psi%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, psi, sys%ions)

    if (target_mode(tg) == oct_targetmode_td .or. &
        (sys%hm%theory_level /= INDEPENDENT_PARTICLES .and. (.not.sys%ks%frozen_hxc))) then
      call states_elec_end(psi2)
      call controlfunction_end(par_prev)
    end if

    call controlfunction_to_basis(par)
    if(aux_fwd_propagation) call propagator_elec_end(tr_psi2)
    call states_elec_end(chi)
    call propagator_elec_end(tr_chi)
    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call psi%unpack()
    nullify(psi)
    nullify(chi)
    POP_SUB(fwd_step)
  end subroutine fwd_step
  ! ---------------------------------------------------------


  !> --------------------------------------------------------
  !! Performs a backward propagation on the state psi and on the
  !! Lagrange-multiplier state chi, according to the following
  !! scheme:
  !!
  !! |psi> --> U[par](0, T)|psi>
  !! par_chi = par_chi[|psi>, |chi>]
  !! |chi> --> U[par_chi](0, T)|chi>
  !! --------------------------------------------------------
  subroutine bwd_step(sys, td, tg, par, par_chi, qcchi, prop_chi, prop_psi) 
    type(electrons_t),         intent(inout) :: sys
    type(td_t),                intent(inout) :: td
    type(target_t),            intent(inout) :: tg
    type(controlfunction_t),   intent(in)    :: par
    type(controlfunction_t),   intent(inout) :: par_chi
    type(opt_control_state_t), intent(inout) :: qcchi
    type(oct_prop_t),          intent(inout) :: prop_chi
    type(oct_prop_t),          intent(inout) :: prop_psi

    integer :: i, ierr
    type(propagator_base_t) :: tr_chi
    type(opt_control_state_t) :: qcpsi
    type(states_elec_t), pointer :: chi, psi

    PUSH_SUB(bwd_step)

    message(1) = "Info: Backward propagation."
    call messages_info(1)

    call controlfunction_to_realtime(par_chi)

    chi => opt_control_point_qs(qcchi)
    psi => opt_control_point_qs(qcpsi)

    call propagator_elec_copy(tr_chi, td%tr)
    ! The propagation of chi should not be self-consistent, because the Kohn-Sham
    ! potential used is the one created by psi. Note, however, that it is likely that
    ! the first two iterations are done self-consistently nonetheless.
    call propagator_elec_remove_scf_prop(tr_chi)

    call states_elec_copy(psi, chi)
    call oct_prop_load_states(prop_psi, sys%namespace, sys%space, psi, sys%gr%mesh, sys%kpoints, td%max_iter, ierr)
    if (ierr /= 0) then
      message(1) = "Unable to read OCT states restart."
      call messages_fatal(1)
    end if
    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call psi%pack()
    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call chi%pack()

    call density_calc(psi, sys%gr, psi%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, psi, sys%ions)
    call hamiltonian_elec_update(sys%hm, sys%gr%mesh, sys%namespace, sys%space)
    call propagator_elec_run_zero_iter(sys%hm, sys%gr, td%tr)
    call propagator_elec_run_zero_iter(sys%hm, sys%gr, tr_chi)

    td%dt = -td%dt
    call oct_prop_dump_states(prop_chi, sys%space, td%max_iter, chi, sys%gr%mesh, sys%kpoints, ierr)
    if (ierr /= 0) then
      message(1) = "Unable to write OCT states restart."
      call messages_warning(1)
    end if

    do i = td%max_iter, 1, -1
      call oct_prop_check(prop_psi, sys%namespace, sys%space, psi, sys%gr%mesh, sys%kpoints, i)
      call update_field(i, par_chi, sys%gr, sys%hm, sys%ions, qcpsi, qcchi, par, dir = 'b')
      call update_hamiltonian_elec_chi(i-1, sys%namespace, sys%gr, sys%ks, sys%hm, td, tg, par_chi, sys%ions, psi)
      call hamiltonian_elec_update(sys%hm, sys%gr%mesh, sys%namespace, sys%space, time = abs(i*td%dt))
      call propagator_elec_dt(sys%ks, sys%namespace, sys%space, sys%hm, sys%gr, chi, tr_chi, abs((i-1)*td%dt), td%dt, td%mu, &
        i-1, td%ions_dyn, sys%ions, sys%outp)
      call oct_prop_dump_states(prop_chi, sys%space, i-1, chi, sys%gr%mesh, sys%kpoints, ierr)
      if (ierr /= 0) then
        message(1) = "Unable to write OCT states restart."
        call messages_warning(1)
      end if
      call update_hamiltonian_elec_psi(i-1, sys%namespace, sys%space, sys%gr, sys%ks, sys%hm, td, tg, par, psi, sys%ions)
      call hamiltonian_elec_update(sys%hm, sys%gr%mesh, sys%namespace, sys%space, time = abs(i*td%dt))
      call propagator_elec_dt(sys%ks, sys%namespace, sys%space, sys%hm, sys%gr, psi, td%tr, abs((i-1)*td%dt), td%dt, td%mu, &
        i-1, td%ions_dyn, sys%ions, sys%outp)
    end do
    td%dt = -td%dt
    call update_field(0, par_chi, sys%gr, sys%hm, sys%ions, qcpsi, qcchi, par, dir = 'b')

    call density_calc(psi, sys%gr, psi%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, psi, sys%ions)
    call hamiltonian_elec_update(sys%hm, sys%gr%mesh, sys%namespace, sys%space)

    call controlfunction_to_basis(par_chi)
    call states_elec_end(psi)
    call propagator_elec_end(tr_chi)
    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call chi%unpack()
    nullify(chi)
    nullify(psi)
    POP_SUB(bwd_step)
  end subroutine bwd_step
  ! ---------------------------------------------------------


  !> --------------------------------------------------------
  !! Performs a backward propagation on the state psi and on the
  !! Lagrange-multiplier state chi, according to the following
  !! scheme:
  !!
  !! |psi> --> U[par](0, T)|psi>
  !! |chi> --> U[par](0, T)|chi>
  !! 
  !! It also calculates during the propagation, a new "output" field:
  !!
  !! par_chi = par_chi[|psi>, |chi>]
  !! --------------------------------------------------------
  subroutine bwd_step_2(sys, td, tg, par, par_chi, qcchi, prop_chi, prop_psi) 
    type(electrons_t),                 intent(inout) :: sys
    type(td_t),                        intent(inout) :: td
    type(target_t),                    intent(inout) :: tg
    type(controlfunction_t),           intent(in)    :: par
    type(controlfunction_t),           intent(inout) :: par_chi
    type(opt_control_state_t),         intent(inout) :: qcchi
    type(oct_prop_t),                  intent(inout) :: prop_chi
    type(oct_prop_t),                  intent(inout) :: prop_psi

    integer :: i, ierr, ik, ib
    logical :: freeze
    type(propagator_base_t) :: tr_chi
    type(opt_control_state_t) :: qcpsi
    type(states_elec_t) :: st_ref
    type(states_elec_t), pointer :: chi, psi
    FLOAT, pointer :: q(:, :), p(:, :)
    FLOAT, allocatable :: qtildehalf(:, :), qinitial(:, :)
    FLOAT, allocatable :: vhxc(:, :)
    FLOAT, allocatable :: fold(:, :), fnew(:, :)
    type(ion_state_t) :: ions_state_initial, ions_state_final

    FLOAT :: init_time, final_time

    PUSH_SUB(bwd_step_2)

    chi => opt_control_point_qs(qcchi)
    q => opt_control_point_q(qcchi)
    p => opt_control_point_p(qcchi)
    SAFE_ALLOCATE(qtildehalf(1:sys%ions%natoms, 1:sys%ions%space%dim))
    SAFE_ALLOCATE(qinitial(1:sys%ions%space%dim, 1:sys%ions%natoms))

    call propagator_elec_copy(tr_chi, td%tr)
    ! The propagation of chi should not be self-consistent, because the Kohn-Sham
    ! potential used is the one created by psi. Note, however, that it is likely that
    ! the first two iterations are done self-consistently nonetheless.
    call propagator_elec_remove_scf_prop(tr_chi)

    call opt_control_state_null(qcpsi)
    call opt_control_state_copy(qcpsi, qcchi)
    psi => opt_control_point_qs(qcpsi)
    call oct_prop_load_states(prop_psi, sys%namespace, sys%space, psi, sys%gr%mesh, sys%kpoints, td%max_iter, ierr)
    if (ierr /= 0) then
      message(1) = "Unable to read OCT states restart."
      call messages_fatal(1)
    end if

    SAFE_ALLOCATE(vhxc(1:sys%gr%mesh%np, 1:sys%hm%d%nspin))

    call density_calc(psi, sys%gr, psi%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, psi, sys%ions)
    call hamiltonian_elec_update(sys%hm, sys%gr%mesh, sys%namespace, sys%space)
    call propagator_elec_run_zero_iter(sys%hm, sys%gr, td%tr)
    call propagator_elec_run_zero_iter(sys%hm, sys%gr, tr_chi)
    td%dt = -td%dt
    call oct_prop_dump_states(prop_chi, sys%space, td%max_iter, chi, sys%gr%mesh, sys%kpoints, ierr)
    if (ierr /= 0) then
      message(1) = "Unable to write OCT states restart."
      call messages_warning(1)
    end if

    call states_elec_copy(st_ref, psi)
    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call psi%pack()
    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call st_ref%pack()

    if (ion_dynamics_ions_move(td%ions_dyn)) then
      call forces_calculate(sys%gr, sys%namespace, sys%ions, sys%hm, psi, sys%ks, t = td%max_iter*abs(td%dt), dt = td%dt)
    end if

    message(1) = "Info: Backward propagation."
    call messages_info(1)
    if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, td%max_iter)

    init_time = loct_clock()

    do i = td%max_iter, 1, -1

      call oct_prop_check(prop_psi, sys%namespace, sys%space, psi, sys%gr%mesh, sys%kpoints, i)
      call update_field(i, par_chi, sys%gr, sys%hm, sys%ions, qcpsi, qcchi, par, dir = 'b')

      select case(td%tr%method)

      case(PROP_EXPLICIT_RUNGE_KUTTA4)

        call update_hamiltonian_elec_psi(i-1, sys%namespace, sys%space, sys%gr, sys%ks, sys%hm, td, tg, par, psi, sys%ions)
        call propagator_elec_dt(sys%ks, sys%namespace, sys%space, sys%hm, sys%gr, psi, td%tr, abs((i-1)*td%dt), td%dt, td%mu, &
          i-1, td%ions_dyn, sys%ions, sys%outp, qcchi = qcchi)

      case default

        if(ion_dynamics_ions_move(td%ions_dyn)) then
          qtildehalf = q

          call ion_dynamics_save_state(td%ions_dyn, sys%ions, ions_state_initial)
          call ion_dynamics_propagate(td%ions_dyn, sys%ions, abs((i-1)*td%dt), M_HALF * td%dt, sys%namespace)
          call sys%ions%get_positions(qinitial)
          call ion_dynamics_restore_state(td%ions_dyn, sys%ions, ions_state_initial)

          SAFE_ALLOCATE(fold(1:sys%ions%natoms, 1:sys%gr%sb%dim))
          SAFE_ALLOCATE(fnew(1:sys%ions%natoms, 1:sys%gr%sb%dim))
          call forces_costate_calculate(sys%gr, sys%namespace, sys%ions, sys%hm, psi, chi, fold, q)

          call ion_dynamics_verlet_step1(sys%ions, qtildehalf, p, fold, M_HALF * td%dt)
          call ion_dynamics_verlet_step1(sys%ions, q, p, fold, td%dt)
        end if

        ! Here propagate psi one full step, and then simply interpolate to get the state
        ! at half the time interval. Perhaps one could gain some accuracy by performing two
        ! successive propagations of half time step.
        call update_hamiltonian_elec_psi(i-1, sys%namespace, sys%space, sys%gr, sys%ks, sys%hm, td, tg, par, psi, sys%ions)

        do ik = psi%d%kpt%start, psi%d%kpt%end
          do ib = psi%group%block_start, psi%group%block_end
            call psi%group%psib(ib, ik)%copy_data_to(sys%gr%mesh%np, st_ref%group%psib(ib, ik))
          end do
        end do

        vhxc(:, :) = sys%hm%vhxc(:, :)
        call propagator_elec_dt(sys%ks, sys%namespace, sys%space, sys%hm, sys%gr, psi, td%tr, abs((i-1)*td%dt), td%dt, td%mu, &
          i-1, td%ions_dyn, sys%ions, sys%outp)

        if(ion_dynamics_ions_move(td%ions_dyn)) then
          call ion_dynamics_save_state(td%ions_dyn, sys%ions, ions_state_final)
          call sys%ions%set_positions(qinitial)
          call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, sys%space, sys%gr, sys%ions, psi, time = abs((i-1)*td%dt))
        end if

        do ik = psi%d%kpt%start, psi%d%kpt%end
          do ib = psi%group%block_start, psi%group%block_end
            call batch_scal(sys%gr%mesh%np, TOCMPLX(M_HALF, M_ZERO), &
              st_ref%group%psib(ib, ik))
            call batch_axpy(sys%gr%mesh%np, TOCMPLX(M_HALF, M_ZERO), &
              psi%group%psib(ib, ik), st_ref%group%psib(ib, ik))
          end do
        end do

        sys%hm%vhxc(:, :) = M_HALF * (sys%hm%vhxc(:, :) + vhxc(:, :))
        call update_hamiltonian_elec_chi(i-1, sys%namespace, sys%gr, sys%ks, sys%hm, td, tg, par, sys%ions, st_ref, &
          qtildehalf)
        freeze = ion_dynamics_freeze(td%ions_dyn)
        call propagator_elec_dt(sys%ks, sys%namespace, sys%space, sys%hm, sys%gr, chi, tr_chi, abs((i-1)*td%dt), td%dt, td%mu, &
          i-1, td%ions_dyn, sys%ions, sys%outp)
        if(freeze) call ion_dynamics_unfreeze(td%ions_dyn)

        if(ion_dynamics_ions_move(td%ions_dyn)) then
          call ion_dynamics_restore_state(td%ions_dyn, sys%ions, ions_state_final)
          call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, sys%space, sys%gr, sys%ions, psi, time = abs((i-1)*td%dt))
          call forces_calculate(sys%gr, sys%namespace, sys%ions, sys%hm, psi, sys%ks, t = abs((i-1)*td%dt), dt = td%dt)
          call forces_costate_calculate(sys%gr, sys%namespace, sys%ions, sys%hm, psi, chi, fnew, q)
          call ion_dynamics_verlet_step2(sys%ions, p, fold, fnew, td%dt)
          SAFE_DEALLOCATE_A(fold)
          SAFE_DEALLOCATE_A(fnew)
        end if

        sys%hm%vhxc(:, :) = vhxc(:, :)

      end select

      call oct_prop_dump_states(prop_chi, sys%space, i-1, chi, sys%gr%mesh, sys%kpoints, ierr)
      if (ierr /= 0) then
        message(1) = "Unable to write OCT states restart."
        call messages_warning(1)
      end if

      if ((mod(i, 100) == 0).and. mpi_grp_is_root(mpi_world)) call loct_progress_bar(td%max_iter-i, td%max_iter)
    end do
    if(mpi_grp_is_root(mpi_world)) then
      call loct_progress_bar(td%max_iter, td%max_iter)
      write(stdout, '(1x)')
    end if

    final_time = loct_clock()
    write(message(1),'(a,f12.2,a)') 'Propagation time: ', final_time - init_time, ' seconds.'
    call messages_info(1)

    call states_elec_end(st_ref)

    td%dt = -td%dt
    call update_hamiltonian_elec_psi(0, sys%namespace, sys%space, sys%gr, sys%ks, sys%hm, td, tg, par, psi, sys%ions)
    call update_field(0, par_chi, sys%gr, sys%hm, sys%ions, qcpsi, qcchi, par, dir = 'b')

    call density_calc(psi, sys%gr, psi%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%space, sys%hm, psi, sys%ions)
    call hamiltonian_elec_update(sys%hm, sys%gr%mesh, sys%namespace, sys%space)

    call propagator_elec_end(tr_chi)

    SAFE_DEALLOCATE_A(vhxc)
    call states_elec_end(psi)

    nullify(chi)
    nullify(psi)
    nullify(q)
    nullify(p)
    SAFE_DEALLOCATE_A(qtildehalf)
    SAFE_DEALLOCATE_A(qinitial)
    POP_SUB(bwd_step_2)
  end subroutine bwd_step_2
  ! ----------------------------------------------------------


  ! ----------------------------------------------------------
  !
  ! ----------------------------------------------------------
  subroutine update_hamiltonian_elec_chi(iter, namespace, gr, ks, hm, td, tg, par_chi, ions, st, qtildehalf)
    integer,                  intent(in)    :: iter
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),             intent(inout) :: gr
    type(v_ks_t),             intent(inout) :: ks
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(td_t),               intent(inout) :: td
    type(target_t),           intent(inout) :: tg
    type(controlfunction_t),  intent(in)    :: par_chi
    type(ions_t),             intent(in)    :: ions
    type(states_elec_t),      intent(inout) :: st
    FLOAT,          optional, intent(in)    :: qtildehalf(:, :)

    type(states_elec_t) :: inh
    type(pert_t) :: pert
    integer :: j, iatom, idim
    CMPLX, allocatable :: dvpsi(:, :, :), zpsi(:, :), inhzpsi(:, :)
    integer :: ist, ik, ib

    PUSH_SUB(update_hamiltonian_elec_chi)

    if(target_mode(tg) == oct_targetmode_td) then
      call states_elec_copy(inh, st)
      call target_inh(st, gr, hm%kpoints, tg, abs(td%dt)*iter, inh, iter)
      call hamiltonian_elec_set_inh(hm, inh)
      call states_elec_end(inh)
    end if

    if(ion_dynamics_ions_move(td%ions_dyn)) then
      call states_elec_copy(inh, st)
      SAFE_ALLOCATE(dvpsi(1:gr%mesh%np_part, 1:st%d%dim, 1:gr%sb%dim))
      do ik = inh%d%kpt%start, inh%d%kpt%end
        do ib = inh%group%block_start, inh%group%block_end
          call batch_set_zero(inh%group%psib(ib, ik))
        end do
      end do
      SAFE_ALLOCATE(zpsi(1:gr%mesh%np_part, 1:st%d%dim))
      SAFE_ALLOCATE(inhzpsi(1:gr%mesh%np_part, 1:st%d%dim))
      do ist = 1, st%nst
        do ik = 1, st%d%nik

          call states_elec_get_state(st, gr%mesh, ist, ik, zpsi)
          call states_elec_get_state(inh, gr%mesh, ist, ik, inhzpsi)

          do iatom = 1, ions%natoms
            do idim = 1, gr%sb%dim
              call pert_init(pert, namespace, PERTURBATION_IONIC, gr, ions)
              call pert_setup_atom(pert, iatom)
              call pert_setup_dir(pert, idim)
              call zpert_apply(pert, namespace, gr, ions, hm, ik, zpsi(:, :), dvpsi(:, :, idim))
              dvpsi(:, :, idim) = - dvpsi(:, :, idim)
              inhzpsi(:,  :)  = &
                inhzpsi(:, :) + st%occ(ist, ik)*qtildehalf(iatom, idim)*dvpsi(:, :, idim)
              call pert_end(pert)
            end do
          end do
          call states_elec_set_state(inh, gr%mesh, ist, ik, inhzpsi)
        end do
      end do

      SAFE_DEALLOCATE_A(zpsi)
      SAFE_DEALLOCATE_A(inhzpsi)
      SAFE_DEALLOCATE_A(dvpsi)
      call hamiltonian_elec_set_inh(hm, inh)
      call states_elec_end(inh)
    end if

    if (hm%theory_level /= INDEPENDENT_PARTICLES .and. (.not. ks%frozen_hxc)) then
      call density_calc(st, gr, st%rho)
      call oct_exchange_set(hm%oct_exchange, st, gr%mesh)
    end if

    call hamiltonian_elec_adjoint(hm)

    do j = iter - 2, iter + 2
      if(j >= 0 .and. j<=td%max_iter) then
        call controlfunction_to_h_val(par_chi, hm%ext_lasers, j+1)
      end if
    end do

    POP_SUB(update_hamiltonian_elec_chi)
  end subroutine update_hamiltonian_elec_chi
  ! ---------------------------------------------------------


  ! ----------------------------------------------------------
  !
  ! ----------------------------------------------------------
  subroutine update_hamiltonian_elec_psi(iter, namespace, space, gr, ks, hm, td, tg, par, st, ions)
    integer,                  intent(in)    :: iter
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(grid_t),             intent(inout) :: gr
    type(v_ks_t),             intent(inout) :: ks
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(td_t),               intent(inout) :: td
    type(target_t),           intent(inout) :: tg
    type(controlfunction_t),  intent(in)    :: par
    type(states_elec_t),      intent(inout) :: st
    type(ions_t),             intent(in)    :: ions

    integer :: j

    PUSH_SUB(update_hamiltonian_elec_psi)

    if (target_mode(tg) == oct_targetmode_td) then
      call hamiltonian_elec_remove_inh(hm)
    end if

    if(ion_dynamics_ions_move(td%ions_dyn)) then
      call hamiltonian_elec_remove_inh(hm)
    end if

    if(hm%theory_level /= INDEPENDENT_PARTICLES .and. (.not. ks%frozen_hxc)) then
      call oct_exchange_remove(hm%oct_exchange)
    end if

    call hamiltonian_elec_not_adjoint(hm)

    do j = iter - 2, iter + 2
      if(j >= 0 .and. j<=td%max_iter) then
        call controlfunction_to_h_val(par, hm%ext_lasers, j+1)
      end if
    end do
    if (hm%theory_level /= INDEPENDENT_PARTICLES .and. (.not. ks%frozen_hxc)) then
      call density_calc(st, gr, st%rho)
      call v_ks_calc(ks, namespace, space, hm, st, ions)
      call hamiltonian_elec_update(hm, gr%mesh, namespace, space)
    end if

    POP_SUB(update_hamiltonian_elec_psi)
  end subroutine update_hamiltonian_elec_psi
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine calculate_g(gr, hm, psi, chi, dl, dq)
    type(grid_t),                   intent(inout) :: gr
    type(hamiltonian_elec_t),       intent(in)    :: hm
    type(states_elec_t),            intent(inout) :: psi
    type(states_elec_t),            intent(inout) :: chi
    CMPLX,                          intent(inout) :: dl(:), dq(:)

    CMPLX, allocatable :: zpsi(:, :), zoppsi(:, :)
    integer :: no_parameters, j, ik, p

    PUSH_SUB(calculate_g)

    no_parameters = hm%ext_lasers%no_lasers

    SAFE_ALLOCATE(zpsi(1:gr%mesh%np_part, 1:chi%d%dim))
    SAFE_ALLOCATE(zoppsi(1:gr%mesh%np_part, 1:chi%d%dim))
    
    do j = 1, no_parameters

      dl(j) = M_z0
      do ik = 1, psi%d%nik
        do p = 1, psi%nst

          call states_elec_get_state(psi, gr%mesh, p, ik, zpsi)
          
          zoppsi = M_z0
          if (allocated(hm%ep%a_static)) then
            call vlaser_operator_linear(hm%ext_lasers%lasers(j), gr%der, hm%d, zpsi, &
              zoppsi, ik, hm%ep%gyromagnetic_ratio, hm%ep%a_static)
          else
            call vlaser_operator_linear(hm%ext_lasers%lasers(j), gr%der, hm%d, zpsi, &
              zoppsi, ik, hm%ep%gyromagnetic_ratio)
          end if

          call states_elec_get_state(chi, gr%mesh, p, ik, zpsi)
          dl(j) = dl(j) + zmf_dotp(gr%mesh, psi%d%dim, zpsi, zoppsi)
        end do
      end do

      ! The quadratic part should only be computed if necessary.
      if (laser_kind(hm%ext_lasers%lasers(j)) == E_FIELD_MAGNETIC) then

        dq(j) = M_z0
        do ik = 1, psi%d%nik
          do p = 1, psi%nst
            zoppsi = M_z0

            call states_elec_get_state(psi, gr%mesh, p, ik, zpsi)
            call vlaser_operator_quadratic(hm%ext_lasers%lasers(j), gr%der, zpsi, zoppsi)

            call states_elec_get_state(chi, gr%mesh, p, ik, zpsi)
            dq(j) = dq(j) + zmf_dotp(gr%mesh, psi%d%dim, zpsi, zoppsi)
            
          end do
        end do

      else
        dq(j) = M_z0
      end if
    end do

    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(zoppsi)
    
    POP_SUB(calculate_g)
  end subroutine calculate_g
  ! ---------------------------------------------------------




  !> Calculates the value of the control functions at iteration
  !! iter, from the state psi and the Lagrange-multiplier chi.
  !!
  !! If dir = 'f', the field must be updated for a forward
  !! propagation. In that case, the propagation step that is
  !! going to be done moves from (iter-1)*|dt| to iter*|dt|.
  !!
  !! If dir = 'b', the field must be updated for a backward
  !! propagation. In taht case, the propagation step that is
  !! going to be done moves from iter*|dt| to (iter-1)*|dt|.
  !!
  !! cp = (1-eta)*cpp - (eta/alpha) * <chi|V|Psi>
  subroutine update_field(iter, cp, gr, hm, ions, qcpsi, qcchi, cpp, dir)
    integer,                   intent(in)    :: iter
    type(controlfunction_t),   intent(inout) :: cp
    type(grid_t),              intent(inout) :: gr
    type(hamiltonian_elec_t),       intent(in)    :: hm
    type(ions_t),              intent(in)    :: ions
    type(opt_control_state_t), intent(inout) :: qcpsi
    type(opt_control_state_t), intent(inout) :: qcchi
    type(controlfunction_t),   intent(in)    :: cpp
    character(len=1),          intent(in)    :: dir

    CMPLX :: d1, pol(MAX_DIM)
    CMPLX, allocatable  :: dl(:), dq(:), zpsi(:, :), zchi(:, :)
    FLOAT, allocatable :: d(:)
    integer :: j, no_parameters, iatom
    type(states_elec_t), pointer :: psi, chi
    FLOAT, pointer :: q(:, :)

    PUSH_SUB(update_field)

    psi => opt_control_point_qs(qcpsi)
    chi => opt_control_point_qs(qcchi)
    q => opt_control_point_q(qcchi)

    no_parameters = controlfunction_number(cp)
    
    SAFE_ALLOCATE(dl(1:no_parameters))
    SAFE_ALLOCATE(dq(1:no_parameters))
    SAFE_ALLOCATE( d(1:no_parameters))

    call calculate_g(gr, hm, psi, chi, dl, dq)
    d1 = M_z1
    if(zbr98_) then
      SAFE_ALLOCATE(zpsi(1:gr%mesh%np, 1:psi%d%dim))
      SAFE_ALLOCATE(zchi(1:gr%mesh%np, 1:chi%d%dim))

      call states_elec_get_state(psi, gr%mesh, 1, 1, zpsi)
      call states_elec_get_state(chi, gr%mesh, 1, 1, zchi)

      d1 = zmf_dotp(gr%mesh, psi%d%dim, zpsi, zchi)
      do j = 1, no_parameters
        d(j) = aimag(d1*dl(j)) / controlfunction_alpha(cp, j)
      end do

      SAFE_DEALLOCATE_A(zpsi)
      SAFE_DEALLOCATE_A(zchi)

    elseif(gradients_) then
      do j = 1, no_parameters
        d(j) = M_TWO * aimag(dl(j))
      end do
    else
      do j = 1, no_parameters
        d(j) = aimag(dl(j)) / controlfunction_alpha(cp, j)
      end do
    end if

    ! This is for the classical target.
    if(dir == 'b') then
      pol = laser_polarization(hm%ext_lasers%lasers(1))
      do iatom = 1, ions%natoms
        d(1) = d(1) - species_zval(ions%atom(iatom)%species) * &
          TOFLOAT(sum(pol(1:gr%sb%dim)*q(iatom, 1:gr%sb%dim)))
      end do
    end if


    if(dir == 'f') then
      call controlfunction_update(cp, cpp, dir, iter, delta_, d, dq)
    else
      call controlfunction_update(cp, cpp, dir, iter, eta_, d, dq)
    end if

    nullify(q)
    nullify(psi)
    nullify(chi)
    SAFE_DEALLOCATE_A(d)
    SAFE_DEALLOCATE_A(dl)
    SAFE_DEALLOCATE_A(dq)
    POP_SUB(update_field)
  end subroutine update_field
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_prop_init(prop, namespace, dirname, mesh, mc)
    type(oct_prop_t),  intent(inout) :: prop
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: dirname
    type(mesh_t),      intent(in)    :: mesh
    type(multicomm_t), intent(in)    :: mc
    
    integer :: j, ierr

    PUSH_SUB(oct_prop_init)

    prop%dirname = dirname
    prop%niter = niter_
    prop%number_checkpoints = number_checkpoints_

    ! The OCT_DIR//trim(dirname) will be used to write and read information during the calculation,
    ! so they need to use the same path.
    call restart_init(prop%restart_dump, namespace, RESTART_OCT, RESTART_TYPE_DUMP, mc, ierr, mesh=mesh)
    call restart_init(prop%restart_load, namespace, RESTART_OCT, RESTART_TYPE_LOAD, mc, ierr, mesh=mesh)

    SAFE_ALLOCATE(prop%iter(1:prop%number_checkpoints+2))
    prop%iter(1) = 0
    do j = 1, prop%number_checkpoints
      prop%iter(j+1) = nint(TOFLOAT(niter_)/(prop%number_checkpoints+1) * j)
    end do
    prop%iter(prop%number_checkpoints+2) = niter_

    POP_SUB(oct_prop_init)
  end subroutine oct_prop_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_prop_end(prop)
    type(oct_prop_t), intent(inout) :: prop

    PUSH_SUB(oct_prop_end)

    call restart_end(prop%restart_load)
    call restart_end(prop%restart_dump)

    SAFE_DEALLOCATE_A(prop%iter)
    ! This routine should maybe delete the files?

    POP_SUB(oct_prop_end)
  end subroutine oct_prop_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_prop_check(prop, namespace, space, psi, mesh, kpoints, iter)
    type(oct_prop_t),    intent(inout) :: prop
    type(namespace_t),   intent(in)    :: namespace
    type(space_t),       intent(in)    :: space
    type(states_elec_t), intent(inout) :: psi
    type(mesh_t),        intent(in)    :: mesh
    type(kpoints_t),     intent(in)    :: kpoints
    integer,             intent(in)    :: iter

    type(states_elec_t) :: stored_st
    character(len=80) :: dirname
    integer :: j, ierr
    CMPLX :: overlap, prev_overlap
    FLOAT, parameter :: WARNING_THRESHOLD = CNST(1.0e-2)

    PUSH_SUB(oct_prop_check)

    do j = 1, prop%number_checkpoints + 2
     if(prop%iter(j)  ==  iter) then
       call states_elec_copy(stored_st, psi)
       write(dirname,'(a, i4.4)') trim(prop%dirname), j
       call restart_open_dir(prop%restart_load, dirname, ierr)
       if (ierr == 0) then
         call states_elec_load(prop%restart_load, namespace, space, stored_st, mesh, kpoints, ierr, verbose=.false.)
       end if
       if(ierr /= 0) then
         message(1) = "Unable to read wavefunctions from '"//trim(dirname)//"'."
         call messages_fatal(1)
       end if
       call restart_close_dir(prop%restart_load)
       prev_overlap = zstates_elec_mpdotp(namespace, mesh, stored_st, stored_st)
       overlap = zstates_elec_mpdotp(namespace, mesh, stored_st, psi)
       if (abs(overlap - prev_overlap) > WARNING_THRESHOLD) then
          write(message(1), '(a,es13.4)') &
            "Forward-backward propagation produced an error of", abs(overlap-prev_overlap)
          write(message(2), '(a,i8)') "Iter = ", iter
          call messages_warning(2)
       end if
       ! Restore state only if the number of checkpoints is larger than zero.
       if(prop%number_checkpoints > 0) then
         call states_elec_end(psi)
         call states_elec_copy(psi, stored_st)
       end if
       call states_elec_end(stored_st)
     end if
    end do
    POP_SUB(oct_prop_check)
  end subroutine oct_prop_check
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_prop_dump_states(prop, space, iter, psi, mesh, kpoints, ierr)
    type(oct_prop_t),    intent(inout) :: prop
    type(space_t),       intent(in)    :: space
    integer,             intent(in)    :: iter
    type(states_elec_t), intent(inout) :: psi
    type(mesh_t),        intent(inout) :: mesh
    type(kpoints_t),     intent(in)    :: kpoints
    integer,             intent(out)   :: ierr

    integer :: j, err
    character(len=80) :: dirname

    PUSH_SUB(oct_prop_dump_states)

    ierr = 0

    if (restart_skip(prop%restart_dump)) then
      POP_SUB(oct_prop_dump_states)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing OCT propagation states restart."
      call messages_info(1)
    end if

    do j = 1, prop%number_checkpoints + 2
      if(prop%iter(j)  ==  iter) then
        write(dirname,'(a,i4.4)') trim(prop%dirname), j
        call restart_open_dir(prop%restart_dump, dirname, err)
        if (err == 0) then
          call states_elec_dump(prop%restart_dump, space, psi, mesh, kpoints, err, iter, verbose = .false.)
        end if
        if(err /= 0) then
          message(1) = "Unable to write wavefunctions to '"//trim(dirname)//"'."
          call messages_warning(1)
          ierr = ierr + 2**j
        end if
        call restart_close_dir(prop%restart_dump)
      end if
    end do

    if (debug%info) then
      message(1) = "Debug: Writing OCT propagation states restart done."
      call messages_info(1)
    end if

    POP_SUB(oct_prop_dump_states)
  end subroutine oct_prop_dump_states
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_prop_load_states(prop, namespace, space, psi, mesh, kpoints, iter, ierr)
    type(oct_prop_t),    intent(inout) :: prop
    type(namespace_t),   intent(in)    :: namespace
    type(space_t),       intent(in)    :: space
    type(states_elec_t), intent(inout) :: psi
    type(mesh_t),        intent(in)    :: mesh
    type(kpoints_t),     intent(in)    :: kpoints
    integer,             intent(in)    :: iter
    integer,             intent(out)   :: ierr

    integer :: j, err
    character(len=80) :: dirname

    PUSH_SUB(oct_prop_load_states)

    ierr = 0

    if (restart_skip(prop%restart_load)) then
      ierr = -1
      POP_SUB(oct_prop_load_states)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading OCT propagation states restart."
      call messages_info(1)
    end if

    do j = 1, prop%number_checkpoints + 2
      if (prop%iter(j)  ==  iter) then
        write(dirname,'(a, i4.4)') trim(prop%dirname), j
        call restart_open_dir(prop%restart_load, dirname, err)
        if (err == 0) then
          call states_elec_load(prop%restart_load, namespace, space, psi, mesh, kpoints, err, verbose=.false.)
        end if
        if(err /= 0) then
          message(1) = "Unable to read wavefunctions from '"//trim(dirname)//"'."
          call messages_warning(1)
          ierr = ierr + 2**j
        end if
        call restart_close_dir(prop%restart_load)
      end if
    end do

    if (debug%info) then
      message(1) = "Debug: Reading OCT propagation states restart done."
      call messages_info(1)
    end if

    POP_SUB(oct_prop_load_states)
  end subroutine oct_prop_load_states
  ! ---------------------------------------------------------


end module propagation_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

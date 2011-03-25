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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: opt_control.F90 2873 2007-04-29 22:05:29Z acastro $

#include "global.h"

module opt_control_propagation_m
  use controlfunction_m
  use datasets_m
  use density_m
  use energy_m
  use epot_m
  use excited_states_m
  use forces_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use ion_dynamics_m 
  use lasers_m
  use loct_m
  use mesh_function_m
  use messages_m
  use opt_control_target_m
  use profiling_m
  use restart_m
  use species_m
  use states_m
  use system_m
  use td_m
  use propagator_m
  use td_write_m
  use v_ks_m
  use varinfo_m

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
            oct_prop_end,         &
            oct_prop_output


  type oct_prop_t
    private
    integer :: number_checkpoints
    integer, pointer :: iter(:)
    integer :: niter
    character(len=100) :: dirname
  end type oct_prop_t


  ! Module variables 
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

    ASSERT(.not. (zbr98 .and. gradients) )

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


  ! ---------------------------------------------------------
  ! Performs a full propagation of state psi, with the laser
  ! field specified in par. If write_iter is present and is
  ! set to .true., writes down through the td_write module.
  ! ---------------------------------------------------------
  subroutine propagate_forward(sys, hm, td, par, target, psi, prop, write_iter)
    type(system_t),             intent(inout)  :: sys
    type(hamiltonian_t),        intent(inout)  :: hm
    type(td_t),                 intent(inout)  :: td
    type(controlfunction_t),    intent(in)     :: par
    type(target_t),             intent(inout)  :: target
    type(states_t),             intent(inout)  :: psi
    type(oct_prop_t), optional, intent(in)     :: prop
    logical, optional,          intent(in)     :: write_iter

    integer :: ii, i
    logical :: write_iter_ = .false.
    type(grid_t),  pointer :: gr
    type(td_write_t)           :: write_handler
    FLOAT, allocatable :: x_initial(:,:)
    logical :: vel_target_ = .false. , move_ions_ = .false.
    integer :: iatom

    PUSH_SUB(propagate_forward)

    message(1) = "Info: Forward propagation."
    call messages_info(1)

    call controlfunction_to_h(par, hm%ep)

    write_iter_ = .false.
    if(present(write_iter)) write_iter_ = write_iter

    gr => sys%gr

    if(write_iter_) then
      call td_write_init(write_handler, gr, sys%st, hm, sys%geo, &
        ion_dynamics_ions_move(td%ions), gauge_field_is_applied(hm%ep%gfield), td%kick, td%iter, td%max_iter, td%dt)
      call td_write_data(write_handler, gr, psi, hm, sys%outp, sys%geo, 0)
    end if

    call hamiltonian_not_adjoint(hm)

    ! setup the Hamiltonian
    call density_calc(psi, gr, psi%rho)
    call v_ks_calc(sys%ks, hm, psi, time = M_ZERO)
    call propagator_run_zero_iter(hm, td%tr)

    if(target_type(target) .eq. oct_tg_velocity) then
       SAFE_ALLOCATE(x_initial(1:sys%geo%natoms,1:MAX_DIM))
       vel_target_ = .true.
       do iatom=1, sys%geo%natoms
          sys%geo%atom(iatom)%f(1:MAX_DIM) = M_ZERO
          sys%geo%atom(iatom)%v(1:MAX_DIM) = M_ZERO
          x_initial(iatom,1:MAX_DIM) = sys%geo%atom(iatom)%x(1:MAX_DIM)
       end do
       if(target_move_ions(target)) then
         move_ions_ = .true.
       else
         call epot_precalc_local_potential(hm%ep, sys%gr, sys%geo, time = M_ZERO)
       end if
    end if

    call target_tdcalc(target, gr, psi, 0)

    if(present(prop)) call oct_prop_output(prop, 0, psi, gr)
    ii = 1
    do i = 1, td%max_iter
      ! time-iterate wavefunctions
      call propagator_dt(sys%ks, hm, gr, psi, td%tr, i*td%dt, td%dt, td%mu, td%max_iter, i)

      if(present(prop)) call oct_prop_output(prop, i, psi, gr)

      ! update
      call density_calc(psi, gr, psi%rho)
      call v_ks_calc(sys%ks, hm, psi, time = i*td%dt)
      call total_energy(hm, sys%gr, psi, -1)

      ! if td_target
      call target_tdcalc(target, gr, psi, i)

      ! calculate velocity and new position of each atom
      if(vel_target_) then
         call forces_calculate(gr, sys%geo, hm%ep, psi, i*td%dt)
         do iatom=1, sys%geo%natoms
           if(i.ne.td%max_iter) then
             sys%geo%atom(iatom)%v(1:MAX_DIM) = sys%geo%atom(iatom)%v(1:MAX_DIM) + &
               sys%geo%atom(iatom)%f(1:MAX_DIM)*td%dt/species_weight(sys%geo%atom(iatom)%spec)
           else
             sys%geo%atom(iatom)%v(1:MAX_DIM) = sys%geo%atom(iatom)%v(1:MAX_DIM) + &
               M_HALF * sys%geo%atom(iatom)%f(1:MAX_DIM)*td%dt/species_weight(sys%geo%atom(iatom)%spec)
           end if
           if(move_ions_) then
             sys%geo%atom(iatom)%x(1:MAX_DIM) = sys%geo%atom(iatom)%x(1:MAX_DIM) + &
               sys%geo%atom(iatom)%v(1:MAX_DIM)*td%dt
           end if
         end do
         call hamiltonian_epot_generate(hm, gr, sys%geo, psi, i*td%dt)
      end if

      ! only write in final run
      if(write_iter_) then
        call td_write_iter(write_handler, gr, psi, hm, sys%geo, td%kick, td%dt, i)
        ii = ii + 1 
        if(ii==sys%outp%iter+1 .or. i == td%max_iter) then ! output
          if(i == td%max_iter) sys%outp%iter = ii - 1
          ii = i
          call td_write_data(write_handler, gr, psi, hm, sys%outp, sys%geo, i) 
        end if
      end if
    end do

    if(vel_target_) then
       do iatom=1, sys%geo%natoms
          sys%geo%atom(iatom)%x(1:MAX_DIM) = x_initial(iatom,1:MAX_DIM)
       end do
       SAFE_DEALLOCATE_A(x_initial)
    end if

    if(write_iter_) call td_write_end(write_handler)
    POP_SUB(propagate_forward)
  end subroutine propagate_forward
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Performs a full backward propagation of state psi, with the
  ! external fields specified in Hamiltonian h.
  ! ---------------------------------------------------------
  subroutine propagate_backward(sys, hm, td, psi, prop)
    type(system_t),          intent(inout) :: sys
    type(hamiltonian_t),     intent(inout) :: hm
    type(td_t),              intent(inout) :: td
    type(states_t),          intent(inout) :: psi
    type(oct_prop_t),        intent(in)    :: prop

    integer :: i
    type(grid_t),  pointer :: gr

    PUSH_SUB(propagate_backward)
    
    message(1) = "Info: Backward propagation."
    call messages_info(1)

    gr => sys%gr

    call hamiltonian_adjoint(hm)

    ! setup the Hamiltonian
    call density_calc(psi, gr, psi%rho)
    call v_ks_calc(sys%ks, hm, psi)
    call propagator_run_zero_iter(hm, td%tr)

    call oct_prop_output(prop, td%max_iter, psi, gr)
    do i = td%max_iter, 1, -1
      call propagator_dt(sys%ks, hm, gr, psi, td%tr, (i-1)*td%dt, -td%dt, td%mu, td%max_iter, i)
      call oct_prop_output(prop, i-1, psi, gr)
      call density_calc(psi, gr, psi%rho)
      call v_ks_calc(sys%ks, hm, psi)
    end do

    POP_SUB(propagate_backward)
  end subroutine propagate_backward
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Performs a forward propagation on the state psi and on the
  ! Lagrange-multiplier state chi. It also updates the control
  ! function par,  according to the following scheme:
  ! 
  ! |chi> --> U[par_chi](T, 0)|chi>
  ! par = par[|psi>, |chi>]
  ! |psi> --> U[par](T, 0)|psi>
  !
  ! Note that the control functions "par" are updated on the
  ! fly, so that the propagation of psi is performed with the
  ! "new" control functions.
  ! --------------------------------------------------------
  subroutine fwd_step(sys, td, hm, target, par, par_chi, psi, prop_chi, prop_psi)
    type(system_t), intent(inout)                 :: sys
    type(td_t), intent(inout)                     :: td
    type(hamiltonian_t), intent(inout)            :: hm
    type(target_t), intent(inout)                 :: target
    type(controlfunction_t), intent(inout)        :: par
    type(controlfunction_t), intent(in)           :: par_chi
    type(states_t), intent(inout)                 :: psi
    type(oct_prop_t), intent(in)                  :: prop_chi
    type(oct_prop_t), intent(in)                  :: prop_psi

    integer :: i
    logical :: aux_fwd_propagation
    type(states_t) :: psi2
    type(states_t) :: chi
    type(controlfunction_t) :: par_prev
    type(grid_t), pointer :: gr
    type(propagator_t) :: tr_chi
    type(propagator_t) :: tr_psi2

    PUSH_SUB(fwd_step)

    message(1) = "Info: Forward propagation."
    call messages_info(1)

    gr => sys%gr
    call propagator_copy(tr_chi, td%tr)
    ! The propagation of chi should not be self-consistent, because the Kohn-Sham
    ! potential used is the one created by psi. Note, however, that it is likely that
    ! the first two iterations are done self-consistently nonetheless.
    call propagator_remove_scf_prop(tr_chi)

    aux_fwd_propagation = ( target_mode(target) == oct_targetmode_td .or. &
                           (hm%theory_level.ne.INDEPENDENT_PARTICLES .and. &
                            .not.sys%ks%frozen_hxc ) )
    if(aux_fwd_propagation) then
      call states_copy(psi2, psi)
      call controlfunction_copy(par_prev, par)
    end if

    
    ! setup forward propagation
    call density_calc(psi, gr, psi%rho)
    call v_ks_calc(sys%ks, hm, psi)
    call propagator_run_zero_iter(hm, td%tr)
    call propagator_run_zero_iter(hm, tr_chi)
    if(aux_fwd_propagation) then
      call propagator_copy(tr_psi2, td%tr)
      call propagator_run_zero_iter(hm, tr_psi2)
    end if

    call oct_prop_output(prop_psi, 0, psi, gr)
    call states_copy(chi, psi)
    call oct_prop_read_state(prop_chi, chi, gr, sys%geo, 0)

    do i = 1, td%max_iter
      call update_field(i, par, gr, hm, psi, chi, par_chi, dir = 'f')
      call update_hamiltonian_chi(i, gr, sys%ks, hm, td, target, par_chi, psi2)
      call hamiltonian_update(hm, gr%mesh, time = (i - 1)*td%dt)
      call propagator_dt(sys%ks, hm, gr, chi, tr_chi, i*td%dt, td%dt, td%mu, td%max_iter, i)
      if(aux_fwd_propagation) then
        call update_hamiltonian_psi(i, gr, sys%ks, hm, td, target, par_prev, psi2)
        call propagator_dt(sys%ks, hm, gr, psi2, tr_psi2, i*td%dt, td%dt, td%mu, td%max_iter, i)
      end if
      call update_hamiltonian_psi(i, gr, sys%ks, hm, td, target, par, psi)
      call hamiltonian_update(hm, gr%mesh, time = (i - 1)*td%dt)
      call propagator_dt(sys%ks, hm, gr, psi, td%tr, i*td%dt, td%dt, td%mu, td%max_iter, i)
      call target_tdcalc(target, gr, psi, i) 
      call oct_prop_output(prop_psi, i, psi, gr)
      call oct_prop_check(prop_chi, chi, gr, sys%geo, i)
    end do
    call update_field(td%max_iter+1, par, gr, hm, psi, chi, par_chi, dir = 'f')

    call density_calc(psi, gr, psi%rho)
    call v_ks_calc(sys%ks, hm, psi)

    if( target_mode(target) == oct_targetmode_td .or. &
        (hm%theory_level.ne.INDEPENDENT_PARTICLES .and. (.not.sys%ks%frozen_hxc) ) ) then
      call states_end(psi2)
      call controlfunction_end(par_prev)
    end if

    if(aux_fwd_propagation) call propagator_end(tr_psi2)
    call states_end(chi)
    call propagator_end(tr_chi)
    POP_SUB(fwd_step)
  end subroutine fwd_step
  ! ---------------------------------------------------------


  ! --------------------------------------------------------
  ! Performs a backward propagation on the state psi and on the
  ! Lagrange-multiplier state chi, according to the following
  ! scheme:
  !
  ! |psi> --> U[par](0, T)|psi>
  ! par_chi = par_chi[|psi>, |chi>]
  ! |chi> --> U[par_chi](0, T)|chi>
  ! --------------------------------------------------------
  subroutine bwd_step(sys, td, hm, target, par, par_chi, chi, prop_chi, prop_psi) 
    type(system_t), intent(inout)                 :: sys
    type(td_t), intent(inout)                     :: td
    type(hamiltonian_t), intent(inout)            :: hm
    type(target_t), intent(inout)                 :: target
    type(controlfunction_t), intent(in)           :: par
    type(controlfunction_t), intent(inout)        :: par_chi
    type(states_t), intent(inout)                 :: chi
    type(oct_prop_t), intent(in)                  :: prop_chi
    type(oct_prop_t), intent(in)                  :: prop_psi

    integer :: i
    type(grid_t), pointer :: gr
    type(propagator_t) :: tr_chi
    type(states_t) :: psi

    PUSH_SUB(bwd_step)

    message(1) = "Info: Backward propagation."
    call messages_info(1)

    gr => sys%gr

    call propagator_copy(tr_chi, td%tr)
    ! The propagation of chi should not be self-consistent, because the Kohn-Sham
    ! potential used is the one created by psi. Note, however, that it is likely that
    ! the first two iterations are done self-consistently nonetheless.
    call propagator_remove_scf_prop(tr_chi)

    call states_copy(psi, chi)
    call oct_prop_read_state(prop_psi, psi, gr, sys%geo, td%max_iter)

    call density_calc(psi, gr, psi%rho)
    call v_ks_calc(sys%ks, hm, psi)
    call hamiltonian_update(hm, gr%mesh)
    call propagator_run_zero_iter(hm, td%tr)
    call propagator_run_zero_iter(hm, tr_chi)

    td%dt = -td%dt
    call oct_prop_output(prop_chi, td%max_iter, chi, gr)
    do i = td%max_iter, 1, -1
      call oct_prop_check(prop_psi, psi, gr, sys%geo, i)
      call update_field(i, par_chi, gr, hm, psi, chi, par, dir = 'b')
      call update_hamiltonian_chi(i-1, gr, sys%ks, hm, td, target, par_chi, psi)
      call hamiltonian_update(hm, gr%mesh, abs(i*td%dt))
      call propagator_dt(sys%ks, hm, gr, chi, tr_chi, abs((i-1)*td%dt), td%dt, td%mu, td%max_iter, i)
      call oct_prop_output(prop_chi, i-1, chi, gr)
      call update_hamiltonian_psi(i-1, gr, sys%ks, hm, td, target, par, psi)
      call hamiltonian_update(hm, gr%mesh, abs(i*td%dt))
      call propagator_dt(sys%ks, hm, gr, psi, td%tr, abs((i-1)*td%dt), td%dt, td%mu, td%max_iter, i)
    end do
    td%dt = -td%dt
    call update_field(0, par_chi, gr, hm, psi, chi, par, dir = 'b')

    call density_calc(psi, gr, psi%rho)
    call v_ks_calc(sys%ks, hm, psi)
    call hamiltonian_update(hm, gr%mesh)

    call states_end(psi)
    call propagator_end(tr_chi)
    POP_SUB(bwd_step)
  end subroutine bwd_step
  ! ---------------------------------------------------------


  ! --------------------------------------------------------
  ! Performs a backward propagation on the state psi and on the
  ! Lagrange-multiplier state chi, according to the following
  ! scheme:
  !
  ! |psi> --> U[par](0, T)|psi>
  ! |chi> --> U[par](0, T)|chi>
  ! 
  ! It also calculates during the propagation, a new "output" field:
  !
  ! par_chi = par_chi[|psi>, |chi>]
  ! --------------------------------------------------------
  subroutine bwd_step_2(sys, td, hm, target, par, par_chi, chi, prop_chi, prop_psi) 
    type(system_t), intent(inout)                 :: sys
    type(td_t), intent(inout)                     :: td
    type(hamiltonian_t), intent(inout)            :: hm
    type(target_t), intent(inout)                 :: target
    type(controlfunction_t), intent(in)           :: par
    type(controlfunction_t), intent(inout)        :: par_chi
    type(states_t), intent(inout)                 :: chi
    type(oct_prop_t), intent(in)                  :: prop_chi
    type(oct_prop_t), intent(in)                  :: prop_psi

    integer :: i
    type(grid_t), pointer :: gr
    type(propagator_t) :: tr_chi
    type(states_t) :: psi
    type(states_t) :: psi_aux

    PUSH_SUB(bwd_step_2)

    message(1) = "Info: Backward propagation."
    call messages_info(1)

    gr => sys%gr

    call propagator_copy(tr_chi, td%tr)
    ! The propagation of chi should not be self-consistent, because the Kohn-Sham
    ! potential used is the one created by psi. Note, however, that it is likely that
    ! the first two iterations are done self-consistently nonetheless.
    call propagator_remove_scf_prop(tr_chi)

    call states_copy(psi, chi)
    call oct_prop_read_state(prop_psi, psi, gr, sys%geo, td%max_iter)

    call density_calc(psi, gr, psi%rho)
    call v_ks_calc(sys%ks, hm, psi)
    call hamiltonian_update(hm, gr%mesh)
    call propagator_run_zero_iter(hm, td%tr)
    call propagator_run_zero_iter(hm, tr_chi)
    td%dt = -td%dt
    call oct_prop_output(prop_chi, td%max_iter, chi, gr)
    do i = td%max_iter, 1, -1
      call oct_prop_check(prop_psi, psi, gr, sys%geo, i)
      call update_field(i, par_chi, gr, hm, psi, chi, par, dir = 'b')

      if(target_mode(target) == oct_targetmode_td) then
        call states_copy(psi_aux, psi)
        call update_hamiltonian_psi(i-1, gr, sys%ks, hm, td, target, par, psi_aux)
        call propagator_dt(sys%ks, hm, gr, psi_aux, td%tr, abs((i-1)*td%dt), td%dt*M_HALF, td%mu, td%max_iter, i)
        call update_hamiltonian_chi(i-1, gr, sys%ks, hm, td, target, par, psi, st_aux = psi_aux)
        call propagator_dt(sys%ks, hm, gr, chi, tr_chi, abs((i-1)*td%dt), td%dt, td%mu, td%max_iter, i)
        call update_hamiltonian_psi(i-1, gr, sys%ks, hm, td, target, par, psi)
        call propagator_dt(sys%ks, hm, gr, psi, td%tr, abs((i-1)*td%dt), td%dt, td%mu, td%max_iter, i)
        call oct_prop_output(prop_chi, i-1, chi, gr)
        call states_end(psi_aux)

      else

        call update_hamiltonian_chi(i-1, gr, sys%ks, hm, td, target, par, psi)
        call propagator_dt(sys%ks, hm, gr, chi, tr_chi, abs((i-1)*td%dt), td%dt, td%mu, td%max_iter, i)
        call oct_prop_output(prop_chi, i-1, chi, gr)
        call update_hamiltonian_psi(i-1, gr, sys%ks, hm, td, target, par, psi)
        call propagator_dt(sys%ks, hm, gr, psi, td%tr, abs((i-1)*td%dt), td%dt, td%mu, td%max_iter, i)

      end if

    end do
    td%dt = -td%dt
    call update_field(0, par_chi, gr, hm, psi, chi, par, dir = 'b')

    call density_calc(psi, gr, psi%rho)
    call v_ks_calc(sys%ks, hm, psi)
    call hamiltonian_update(hm, gr%mesh)

    call propagator_end(tr_chi)

    if(target_mode(target) == oct_targetmode_td) call states_end(psi_aux)
    call states_end(psi)

    POP_SUB(bwd_step_2)
  end subroutine bwd_step_2
  ! ----------------------------------------------------------


  ! ----------------------------------------------------------
  !
  ! ----------------------------------------------------------
  subroutine update_hamiltonian_chi(iter, gr, ks, hm, td, target, par_chi, st, st_aux)
    integer, intent(in)                        :: iter
    type(grid_t), intent(inout)                :: gr
    type(v_ks_t), intent(inout)                :: ks
    type(hamiltonian_t), intent(inout)         :: hm
    type(td_t), intent(inout)                  :: td
    type(target_t), intent(inout)              :: target
    type(controlfunction_t), intent(in)        :: par_chi
    type(states_t), intent(inout)              :: st
    type(states_t), optional, intent(inout)    :: st_aux

    type(states_t)                             :: inh
    integer :: j

    PUSH_SUB(update_hamiltonian_chi)

    if(target_mode(target) == oct_targetmode_td) then
      if(present(st_aux)) then
        call states_copy(inh, st_aux)
        call target_inh(st_aux, gr, target, td%dt*iter, inh)
      else
        call states_copy(inh, st)
        call target_inh(st, gr, target, td%dt*iter, inh)
      end if
      call hamiltonian_set_inh(hm, inh)
      call states_end(inh)
    end if

    if( hm%theory_level.ne.INDEPENDENT_PARTICLES .and. (.not.ks%frozen_hxc) ) then
      call hamiltonian_set_oct_exchange(hm, st, gr, ks%xc)
    end if

    call hamiltonian_adjoint(hm)

    do j = iter - 2, iter + 2
      if(j >= 0 .and. j<=td%max_iter) then
        call controlfunction_to_h_val(par_chi, hm%ep, j+1)
      end if
    end do
    if(hm%theory_level.ne.INDEPENDENT_PARTICLES .and. (.not.ks%frozen_hxc) ) then
      call density_calc(st, gr, st%rho)
      call v_ks_calc(ks, hm, st)
      call hamiltonian_update(hm, gr%mesh)
    end if

    POP_SUB(update_hamiltonian_chi)
  end subroutine update_hamiltonian_chi
  ! ---------------------------------------------------------


  ! ----------------------------------------------------------
  !
  ! ----------------------------------------------------------
  subroutine update_hamiltonian_psi(iter, gr, ks, hm, td, target, par, st)
    integer, intent(in)                        :: iter
    type(grid_t), intent(inout)                :: gr
    type(v_ks_t), intent(inout)                :: ks
    type(hamiltonian_t), intent(inout)         :: hm
    type(td_t), intent(inout)                  :: td
    type(target_t), intent(inout)              :: target
    type(controlfunction_t), intent(in)        :: par
    type(states_t), intent(inout)              :: st

    integer :: j

    PUSH_SUB(update_hamiltonian_psi)

    if(target_mode(target) == oct_targetmode_td) then
      call hamiltonian_remove_inh(hm)
    end if

    if(hm%theory_level.ne.INDEPENDENT_PARTICLES .and. (.not.ks%frozen_hxc) ) then
      call hamiltonian_remove_oct_exchange(hm)
    end if

    call hamiltonian_not_adjoint(hm)

    do j = iter - 2, iter + 2
      if(j >= 0 .and. j<=td%max_iter) then
        call controlfunction_to_h_val(par, hm%ep, j+1)
      end if
    end do
    if(hm%theory_level.ne.INDEPENDENT_PARTICLES .and. (.not.ks%frozen_hxc) ) then
      call density_calc(st, gr, st%rho)
      call v_ks_calc(ks, hm, st)
      call hamiltonian_update(hm, gr%mesh)
    end if

    POP_SUB(update_hamiltonian_psi)
  end subroutine update_hamiltonian_psi
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine calculate_g(gr, hm, psi, chi, dl, dq)
    type(grid_t),                   intent(inout) :: gr
    type(hamiltonian_t),            intent(in)    :: hm
    type(states_t),                 intent(inout) :: psi
    type(states_t),                 intent(inout) :: chi
    CMPLX,                          intent(inout) :: dl(:), dq(:)

    type(states_t) :: oppsi
    integer :: no_parameters, j, ik, p

    PUSH_SUB(calculate_g)

    no_parameters = hm%ep%no_lasers

    do j = 1, no_parameters
      call states_copy(oppsi, psi)
      dl(j) = M_z0
      do ik = 1, psi%d%nik
        do p = 1, psi%nst
          oppsi%zpsi(:, :, p, ik) = M_z0
          if(associated(hm%ep%a_static)) then
            call zvlaser_operator_linear(hm%ep%lasers(j), gr%der, hm%d, psi%zpsi(:, :, p, ik), &
              oppsi%zpsi(:, :, p, ik), ik, hm%ep%gyromagnetic_ratio, hm%ep%a_static)
          else
            call zvlaser_operator_linear(hm%ep%lasers(j), gr%der, hm%d, psi%zpsi(:, :, p, ik), &
              oppsi%zpsi(:, :, p, ik), ik, hm%ep%gyromagnetic_ratio)
          end if
          dl(j) = dl(j) + psi%occ(p, ik) * zmf_dotp(gr%mesh, psi%d%dim, chi%zpsi(:, :, p, ik), &
            oppsi%zpsi(:, :, p, ik))
        end do
      end do
      call states_end(oppsi)

      ! The quadratic part should only be computed if necessary.
      if(laser_kind(hm%ep%lasers(j)).eq.E_FIELD_MAGNETIC ) then
        call states_copy(oppsi, psi)
        dq(j) = M_z0
        do ik = 1, psi%d%nik
          do p = 1, psi%nst
            oppsi%zpsi(:, :, p, ik) = M_z0
            call zvlaser_operator_quadratic(hm%ep%lasers(j), gr%der, hm%d, &
              psi%zpsi(:, :, p, ik), oppsi%zpsi(:, :, p, ik))
            dq(j) = dq(j) + psi%occ(p, ik)*zmf_dotp(gr%mesh, psi%d%dim, &
              chi%zpsi(:, :, p, ik), oppsi%zpsi(:, :, p, ik))
          end do
        end do
        call states_end(oppsi)
      else
        dq(j) = M_z0
      end if
    end do

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
  subroutine update_field(iter, cp, gr, hm, psi, chi, cpp, dir)
    integer, intent(in)        :: iter
    type(controlfunction_t), intent(inout) :: cp
    type(grid_t), intent(inout)   :: gr
    type(hamiltonian_t), intent(in) :: hm
    type(states_t), intent(inout) :: psi
    type(states_t), intent(inout) :: chi
    type(controlfunction_t), intent(in) :: cpp
    character(len=1),intent(in) :: dir

    CMPLX :: d1
    CMPLX, allocatable  :: dl(:), dq(:)
    FLOAT, allocatable :: d(:)
    integer :: j, no_parameters

    PUSH_SUB(update_field)

    no_parameters = controlfunction_number(cp)
    
    SAFE_ALLOCATE(dl(1:no_parameters))
    SAFE_ALLOCATE(dq(1:no_parameters))
    SAFE_ALLOCATE( d(1:no_parameters))

    call calculate_g(gr, hm, psi, chi, dl, dq)
    d1 = M_z1
    if(zbr98_) then
      d1 = zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, 1, 1), chi%zpsi(:, :, 1, 1))
      forall(j = 1:no_parameters) d(j) = aimag(d1*dl(j)) / controlfunction_alpha(cp, j) 
    elseif(gradients_) then
      forall(j = 1:no_parameters) d(j) = aimag(dl(j))
    else
      forall(j = 1:no_parameters) d(j) = aimag(dl(j)) / controlfunction_alpha(cp, j) 
    end if

    if(dir == 'f') then
      call controlfunction_update(cp, cpp, dir, iter, delta_, d, dq)
    else
      call controlfunction_update(cp, cpp, dir, iter, eta_, d, dq)
    end if

    SAFE_DEALLOCATE_A(d)
    SAFE_DEALLOCATE_A(dl)
    SAFE_DEALLOCATE_A(dq)
    POP_SUB(update_field)
  end subroutine update_field
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_prop_init(prop, dirname)
    type(oct_prop_t), intent(inout) :: prop
    character(len=*), intent(in)    :: dirname

    integer :: j

    PUSH_SUB(oct_prop_init)

    prop%dirname = OCT_DIR//trim(dirname)
    call io_mkdir(trim(prop%dirname))
    prop%niter = niter_
    prop%number_checkpoints = number_checkpoints_

    SAFE_ALLOCATE(prop%iter(1:prop%number_checkpoints+2))
    prop%iter(1) = 0
    do j = 1, prop%number_checkpoints
      prop%iter(j+1) = nint( real(niter_)/(prop%number_checkpoints+1) * j)
    end do
    prop%iter(prop%number_checkpoints+2) = niter_

    POP_SUB(oct_prop_init)
  end subroutine oct_prop_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_prop_end(prop)
    type(oct_prop_t), intent(inout) :: prop

    PUSH_SUB(oct_prop_end)

    SAFE_DEALLOCATE_P(prop%iter)
    ! This routine should maybe delete the files?

    POP_SUB(oct_prop_end)
  end subroutine oct_prop_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_prop_check(prop, psi, gr, geo, iter)
    type(oct_prop_t),  intent(in)    :: prop
    type(states_t),    intent(inout) :: psi
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(in)    :: geo
    integer,           intent(in)    :: iter

    type(states_t) :: stored_st
    character(len=100) :: filename
    integer :: j, ierr
    CMPLX :: overlap, prev_overlap
    FLOAT, parameter :: WARNING_THRESHOLD = CNST(1.0e-2)

    PUSH_SUB(oct_prop_check)

    do j = 1, prop%number_checkpoints + 2
     if(prop%iter(j) .eq. iter) then
       call states_copy(stored_st, psi)
       write(filename,'(a,i4.4)') trim(prop%dirname)//'/', j
       call restart_read(trim(filename), stored_st, gr, geo, ierr)
       prev_overlap = zstates_mpdotp(gr%mesh, stored_st, stored_st)
       overlap = zstates_mpdotp(gr%mesh, stored_st, psi)
       if( abs(overlap - prev_overlap) > WARNING_THRESHOLD ) then
          write(message(1), '(a,es13.4)') &
            "WARNING: forward-backward propagation produced an error of", abs(overlap-prev_overlap)
          write(message(2), '(a,i8)') "Iter = ", iter
          call messages_warning(2)
       end if
       ! Restore state only if the number of checkpoints is larger than zero.
       if(prop%number_checkpoints > 0) then
         call states_end(psi)
         call states_copy(psi, stored_st)
       end if
       call states_end(stored_st)
     end if
    end do
    POP_SUB(oct_prop_check)
  end subroutine oct_prop_check
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_prop_read_state(prop, psi, gr, geo, iter)
    type(oct_prop_t),  intent(in)    :: prop
    type(states_t),    intent(inout) :: psi
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(in)    :: geo
    integer,           intent(in)    :: iter

    character(len=100) :: filename
    integer :: j, ierr

    PUSH_SUB(oct_prop_read_state)

    do j = 1, prop%number_checkpoints + 2
     if(prop%iter(j) .eq. iter) then
       write(filename,'(a,i4.4)') trim(prop%dirname)//'/', j
       call restart_read(trim(filename), psi, gr, geo, ierr)
     end if
    end do

    POP_SUB(oct_prop_read_state)
  end subroutine oct_prop_read_state
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_prop_output(prop, iter, psi, gr)
    type(oct_prop_t), intent(in) :: prop
    integer, intent(in)           :: iter
    type(states_t), intent(inout) :: psi
    type(grid_t), intent(inout)   :: gr

    integer :: j, ierr
    character(len=100) :: filename

    PUSH_SUB(oct_prop_output)

    do j = 1, prop%number_checkpoints + 2
      if(prop%iter(j) .eq. iter) then
        write(filename,'(a,i4.4)') trim(prop%dirname)//'/', j
        call restart_write(io_workpath(filename), psi, gr, ierr, iter)
      end if
    end do

    POP_SUB(oct_prop_output)
  end subroutine oct_prop_output
  ! ---------------------------------------------------------


end module opt_control_propagation_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

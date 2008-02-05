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
  use datasets_m
  use varinfo_m
  use global_m
  use loct_m
  use io_m
  use messages_m
  use units_m
  use grid_m
  use mesh_function_m
  use states_m
  use excited_states_m
  use hamiltonian_m
  use geometry_m
  use system_m
  use timedep_m
  use restart_m
  use td_rti_m
  use td_write_m
  use opt_control_constants_m
  use opt_control_target_m
  use opt_control_parameters_m
  use opt_control_output_m
  use lasers_m
  use v_ks_m
  use tdf_m

  implicit none

  private
  public :: propagate_forward,  &
            propagate_backward, &
            fwd_step,           &
            bwd_step,           &
            bwd_step_2,         &
            oct_prop_t,         &
            oct_prop_init,      &
            oct_prop_check,     &
            oct_prop_end,       &
            oct_prop_output


  type oct_prop_t
    integer :: number_checkpoints
    integer, pointer :: iter(:)
    integer :: niter
    character(len=100) :: dirname
  end type oct_prop_t
  

  contains

  ! ---------------------------------------------------------
  ! Performs a full propagation of state psi, with the laser
  ! field specified in hamiltonian h. If write_iter is present
  ! and is set to .true., writes down through the td_write
  ! module.
  ! ---------------------------------------------------------
  subroutine propagate_forward(sys, h, td, target, psi, prop, write_iter)
    type(system_t),             intent(inout) :: sys
    type(hamiltonian_t),        intent(inout) :: h
    type(td_t),                 intent(inout) :: td
    type(target_t),             intent(inout) :: target
    type(states_t),             intent(inout) :: psi
    type(oct_prop_t), optional, intent(in)    :: prop
    logical, optional,          intent(in)    :: write_iter

    integer :: ii, i
    logical :: write_iter_ = .false.
    type(grid_t),  pointer :: gr
    type(td_write_t)           :: write_handler

    call push_sub('propagation.propagate_forward')

    message(1) = "Info: Forward propagation."
    call write_info(1)

    write_iter_ = .false.
    if(present(write_iter)) write_iter_ = write_iter

    gr => sys%gr

    if(write_iter_) then
      call td_write_init(write_handler, gr, sys%st, sys%geo, &
        (td%move_ions>0), h%ep%with_gauge_field, td%iter, td%dt)
      call td_write_data(write_handler, gr, psi, h, sys%outp, sys%geo, 0)
    end if

    call hamiltonian_not_adjoint(h)

    ! setup the hamiltonian
    call states_calc_dens(psi, NP_PART, psi%rho)
    call v_ks_calc(gr, sys%ks, h, psi)
    call td_rti_run_zero_iter(h, td%tr)

    if(present(prop)) call oct_prop_output(prop, 0, psi, gr)
    ii = 1
    do i = 1, td%max_iter
      ! time iterate wavefunctions
      call td_rti_dt(sys%ks, h, gr, psi, td%tr, i*td%dt, td%dt, td%max_iter)

      if(present(prop)) call oct_prop_output(prop, i, psi, gr)

      ! if td_target
      if(target%mode==oct_targetmode_td) call calc_tdfitness(target, gr, psi, i)
      
      ! update
      call states_calc_dens(psi, NP_PART, psi%rho)
      call v_ks_calc(gr, sys%ks, h, psi)
      call hamiltonian_energy(h, sys%gr, sys%geo, psi, -1)

      ! only write in final run
      if(write_iter_) then
        call td_write_iter(write_handler, gr, psi, h, sys%geo, td%kick, td%dt, i)
        ii = ii + 1 
        if(ii==sys%outp%iter+1 .or. i == td%max_iter) then ! output
          if(i == td%max_iter) sys%outp%iter = ii - 1
          ii = i
          call td_write_data(write_handler, gr, psi, h, sys%outp, sys%geo, i) 
        end if
      end if
    end do

    if(write_iter_) call td_write_end(write_handler)
    call pop_sub()
  end subroutine propagate_forward
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Performs a full bacward propagation of state psi, with the
  ! external fields specified in hamiltonian h.
  ! ---------------------------------------------------------
  subroutine propagate_backward(sys, h, td, psi, prop)
    type(system_t),          intent(inout) :: sys
    type(hamiltonian_t),     intent(inout) :: h
    type(td_t),              intent(inout) :: td
    type(states_t),          intent(inout) :: psi
    type(oct_prop_t),        intent(in)    :: prop

    integer :: i
    type(grid_t),  pointer :: gr

    call push_sub('propagation.propagate_backward')
    
    message(1) = "Info: Backward propagation."
    call write_info(1)

    gr => sys%gr

    call hamiltonian_adjoint(h)

    ! setup the hamiltonian
    call states_calc_dens(psi, NP_PART, psi%rho)
    call v_ks_calc(gr, sys%ks, h, psi)
    call td_rti_run_zero_iter(h, td%tr)

    call oct_prop_output(prop, td%max_iter, psi, gr)
    do i = td%max_iter, 1, -1
      call td_rti_dt(sys%ks, h, gr, psi, td%tr, (i-1)*td%dt, -td%dt, td%max_iter)
      call oct_prop_output(prop, i-1, psi, gr)
      call states_calc_dens(psi, NP_PART, psi%rho)
      call v_ks_calc(gr, sys%ks, h, psi)
    end do

    call pop_sub()
  end subroutine propagate_backward
  ! ---------------------------------------------------------


  ! /*---------------------------------------------------------
  ! Performs a forward propagation on the state psi and on the
  ! lagrange-multiplier state chi. It also updates the control
  ! parameter par,  according to the following scheme:
  ! 
  ! |chi> --> U[par_chi](T, 0)|chi>
  ! par = par[|psi>, |chi>]
  ! |psi> --> U[par](T, 0)|psi>
  !
  ! Note that the control parameters "par" are updated "on the
  ! fly", so that the propagation of psi is performed with the
  ! "new" parameters.
  ! */--------------------------------------------------------
  subroutine fwd_step(oct, sys, td, h, target, par, par_chi, psi, prop_chi, prop_psi)
    type(oct_t), intent(in)                       :: oct
    type(system_t), intent(inout)                 :: sys
    type(td_t), intent(inout)                     :: td
    type(hamiltonian_t), intent(inout)            :: h
    type(target_t), intent(inout)                 :: target
    type(oct_control_parameters_t), intent(inout) :: par
    type(oct_control_parameters_t), intent(in)    :: par_chi
    type(states_t), intent(inout)                 :: psi
    type(oct_prop_t), intent(in)                  :: prop_chi
    type(oct_prop_t), intent(in)                  :: prop_psi

    integer :: i
    logical :: aux_fwd_propagation
    type(states_t) :: psi2
    type(states_t) :: chi
    type(oct_control_parameters_t) :: par_prev
    type(grid_t), pointer :: gr
    type(td_rti_t) :: tr_chi
    type(td_rti_t) :: tr_psi2

    call push_sub('propagation.fwd_step')

    message(1) = "Info: Forward propagation."
    call write_info(1)

    gr => sys%gr
    call td_rti_copy(tr_chi, td%tr)

    aux_fwd_propagation = (target%mode == oct_targetmode_td .or. (h%theory_level.ne.INDEPENDENT_PARTICLES))
    if(aux_fwd_propagation) then
      call states_copy(psi2, psi)
      call parameters_copy(par_prev, par)
    end if

    
    ! setup forward propagation
    call states_densities_init(psi, gr, sys%geo)
    call states_calc_dens(psi, NP_PART, psi%rho)
    call v_ks_calc(gr, sys%ks, h, psi)
    call td_rti_run_zero_iter(h, td%tr)
    call td_rti_run_zero_iter(h, tr_chi)
    if(aux_fwd_propagation) then
      call td_rti_copy(tr_psi2, td%tr)
      call td_rti_run_zero_iter(h, tr_psi2)
    end if

    call oct_prop_output(prop_psi, 0, psi, gr)
    call states_copy(chi, psi)
    call oct_prop_read_state(prop_chi, chi, gr, sys%geo, 0)

    do i = 1, td%max_iter
      if(aux_fwd_propagation) then
        call update_hamiltonian_psi(i, gr, sys%ks, h, td, target, par_prev, psi2)
        call td_rti_dt(sys%ks, h, gr, psi2, tr_psi2, i*td%dt, td%dt, td%max_iter)
      end if
      call update_field(oct, i, par, gr, td, h, psi, chi, par_chi, dir = 'f')
      call update_hamiltonian_chi(i, gr, sys%ks, h, td, target, par_chi, psi2)
      call td_rti_dt(sys%ks, h, gr, chi, tr_chi, i*td%dt, td%dt, td%max_iter)
      call update_hamiltonian_psi(i, gr, sys%ks, h, td, target, par, psi)
      call td_rti_dt(sys%ks, h, gr, psi, td%tr, i*td%dt, td%dt, td%max_iter)
      if(target%mode == oct_targetmode_td) call calc_tdfitness(target, gr, psi, i) 
      call oct_prop_output(prop_psi, i, psi, gr)
      call oct_prop_check(prop_chi, chi, gr, sys%geo, i)
    end do

    call states_calc_dens(psi, NP_PART, psi%rho)
    call v_ks_calc(gr, sys%ks, h, psi)

    if(target%mode == oct_targetmode_td .or. (h%theory_level.ne.INDEPENDENT_PARTICLES) ) then
      call states_end(psi2)
      call parameters_end(par_prev)
    end if

    if(aux_fwd_propagation) call td_rti_end(tr_psi2)
    call states_end(chi)
    call td_rti_end(tr_chi)
    call pop_sub()
  end subroutine fwd_step
  ! ---------------------------------------------------------


  ! --------------------------------------------------------
  ! Performs a backward propagation on the state psi and on the
  ! lagrange-multiplier state chi, according to the following
  ! scheme:
  !
  ! |psi> --> U[par](0, T)|psi>
  ! par_chi = par_chi[|psi>, |chi>]
  ! |chi> --> U[par_chi](0, T)|chi>
  ! --------------------------------------------------------
  subroutine bwd_step(oct, sys, td, h, target, par, par_chi, chi, prop_chi, prop_psi) 
    type(oct_t), intent(in)                       :: oct
    type(system_t), intent(inout)                 :: sys
    type(td_t), intent(inout)                     :: td
    type(hamiltonian_t), intent(inout)            :: h
    type(target_t), intent(inout)                 :: target
    type(oct_control_parameters_t), intent(in)    :: par
    type(oct_control_parameters_t), intent(inout) :: par_chi
    type(states_t), intent(inout)                 :: chi
    type(oct_prop_t), intent(in)                  :: prop_chi
    type(oct_prop_t), intent(in)                  :: prop_psi

    integer :: i
    type(grid_t), pointer :: gr
    type(td_rti_t) :: tr_chi
    type(states_t) :: psi

    call push_sub('propagation.bwd_step')

    message(1) = "Info: Backward propagation."
    call write_info(1)

    gr => sys%gr

    call td_rti_copy(tr_chi, td%tr)

    call states_copy(psi, chi)
    call oct_prop_read_state(prop_psi, psi, gr, sys%geo, td%max_iter)

    call states_calc_dens(psi, NP_PART, psi%rho)
    call v_ks_calc(gr, sys%ks, h, psi)
    call td_rti_run_zero_iter(h, td%tr)
    call td_rti_run_zero_iter(h, tr_chi)

    td%dt = -td%dt
    call oct_prop_output(prop_chi, td%max_iter, chi, gr)
    do i = td%max_iter, 1, -1
      call oct_prop_check(prop_psi, psi, gr, sys%geo, i)
      call update_field(oct, i, par_chi, gr, td, h, psi, chi, par, dir = 'b')
      call update_hamiltonian_chi(i-1, gr, sys%ks, h, td, target, par_chi, psi)
      call td_rti_dt(sys%ks, h, gr, chi, tr_chi, abs((i-1)*td%dt), td%dt, td%max_iter)
      call oct_prop_output(prop_chi, i-1, chi, gr)
      call update_hamiltonian_psi(i-1, gr, sys%ks, h, td, target, par, psi)
      call td_rti_dt(sys%ks, h, gr, psi, td%tr, abs((i-1)*td%dt), td%dt, td%max_iter)
    end do
    td%dt = -td%dt

    call states_calc_dens(psi, NP_PART, psi%rho)
    call v_ks_calc(gr, sys%ks, h, psi)

    call states_end(psi)
    call td_rti_end(tr_chi)
    call pop_sub()
  end subroutine bwd_step
  ! ---------------------------------------------------------


  ! --------------------------------------------------------
  ! Performs a backward propagation on the state psi and on the
  ! lagrange-multiplier state chi, according to the following
  ! scheme:
  !
  ! |psi> --> U[par](0, T)|psi>
  ! |chi> --> U[par](0, T)|chi>
  ! 
  ! It also calculates during the propagation, a new "output" field:
  !
  ! par_chi = par_chi[|psi>, |chi>]
  ! --------------------------------------------------------
  subroutine bwd_step_2(oct, sys, td, h, target, par, par_chi, chi, prop_chi, prop_psi) 
    type(oct_t), intent(in)                       :: oct
    type(system_t), intent(inout)                 :: sys
    type(td_t), intent(inout)                     :: td
    type(hamiltonian_t), intent(inout)            :: h
    type(target_t), intent(inout)                 :: target
    type(oct_control_parameters_t), intent(in)    :: par
    type(oct_control_parameters_t), intent(inout) :: par_chi
    type(states_t), intent(inout)                 :: chi
    type(oct_prop_t), intent(in)                  :: prop_chi
    type(oct_prop_t), intent(in)                  :: prop_psi

    integer :: i
    type(grid_t), pointer :: gr
    type(td_rti_t) :: tr_chi
    type(states_t) :: psi

    call push_sub('propagation.bwd_step_2')

    message(1) = "Info: Backward propagation."
    call write_info(1)

    gr => sys%gr

    call td_rti_copy(tr_chi, td%tr)

    call states_copy(psi, chi)
    call oct_prop_read_state(prop_psi, psi, gr, sys%geo, td%max_iter)

    call states_calc_dens(psi, NP_PART, psi%rho)
    call v_ks_calc(gr, sys%ks, h, psi)
    call td_rti_run_zero_iter(h, td%tr)
    call td_rti_run_zero_iter(h, tr_chi)

    td%dt = -td%dt
    call oct_prop_output(prop_chi, td%max_iter, chi, gr)
    do i = td%max_iter, 1, -1
      call oct_prop_check(prop_psi, psi, gr, sys%geo, i)
      call update_field(oct, i, par_chi, gr, td, h, psi, chi, par, dir = 'b')
      call update_hamiltonian_chi(i-1, gr, sys%ks, h, td, target, par, psi)
      call td_rti_dt(sys%ks, h, gr, chi, tr_chi, abs((i-1)*td%dt), td%dt, td%max_iter)
      call oct_prop_output(prop_chi, i-1, chi, gr)
      call update_hamiltonian_psi(i-1, gr, sys%ks, h, td, target, par, psi)
      call td_rti_dt(sys%ks, h, gr, psi, td%tr, abs((i-1)*td%dt), td%dt, td%max_iter)
    end do
    td%dt = -td%dt

    call states_calc_dens(psi, NP_PART, psi%rho)
    call v_ks_calc(gr, sys%ks, h, psi)

    call td_rti_end(tr_chi)

    call pop_sub()
  end subroutine bwd_step_2


  ! ----------------------------------------------------------
  !
  ! ----------------------------------------------------------
  subroutine update_hamiltonian_chi(iter, gr, ks, h, td, target, par_chi, st)
    integer, intent(in)                        :: iter
    type(grid_t), intent(inout)                :: gr
    type(v_ks_t), intent(inout)                :: ks
    type(hamiltonian_t), intent(inout)         :: h
    type(td_t), intent(inout)                  :: td
    type(target_t), intent(inout)              :: target
    type(oct_control_parameters_t), intent(in) :: par_chi
    type(states_t), intent(inout)              :: st
    type(states_t)                             :: inh

    integer :: j
    call push_sub('propagation.update_hamiltonian_chi')

    if(target%mode == oct_targetmode_td) then
      call calc_inh(st, gr, target, td%dt*iter, inh)
      call hamiltonian_set_inh(h, inh)
    end if

    if(h%theory_level.ne.INDEPENDENT_PARTICLES) then
      call hamiltonian_set_oct_exchange(h, st)
    end if

    call hamiltonian_adjoint(h)

    do j = iter - 2, iter + 2
      if(j >= 0 .and. j<=td%max_iter) then
        call parameters_to_h_val(par_chi, h%ep, j+1)
      end if
    end do
    if(h%theory_level.ne.INDEPENDENT_PARTICLES) then
      call states_calc_dens(st, NP_PART, st%rho)
      call v_ks_calc(gr, ks, h, st)
    end if

    call pop_sub()
  end subroutine update_hamiltonian_chi
  ! ---------------------------------------------------------


  ! ----------------------------------------------------------
  !
  ! ----------------------------------------------------------
  subroutine update_hamiltonian_psi(iter, gr, ks, h, td, target, par, st)
    integer, intent(in)                        :: iter
    type(grid_t), intent(inout)                :: gr
    type(v_ks_t), intent(inout)                :: ks
    type(hamiltonian_t), intent(inout)         :: h
    type(td_t), intent(inout)                  :: td
    type(target_t), intent(inout)              :: target
    type(oct_control_parameters_t), intent(in) :: par
    type(states_t), intent(inout)              :: st

    integer :: j
    call push_sub('propagation.update_hamiltonian_psi')

    if(target%mode == oct_targetmode_td) then
      call hamiltonian_remove_inh(h)
    end if

    if(h%theory_level.ne.INDEPENDENT_PARTICLES) then
      call hamiltonian_remove_oct_exchange(h)
    end if

    call hamiltonian_not_adjoint(h)

    do j = iter - 2, iter + 2
      if(j >= 0 .and. j<=td%max_iter) then
        call parameters_to_h_val(par, h%ep, j+1)
      end if
    end do
    if(h%theory_level.ne.INDEPENDENT_PARTICLES) then
      call states_calc_dens(st, NP_PART, st%rho)
      call v_ks_calc(gr, ks, h, st)
    end if

    call pop_sub()
  end subroutine update_hamiltonian_psi
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Calculates the value of the control parameters at iteration
  ! iter, from the state psi and the Lagrange-multiplier chi.
  !
  ! If dir = 'f', the field must be updated for a forward
  ! propagation. In that case, the propagation step that is
  ! going to be done moves from (iter-1)*|dt| to iter*|dt|.
  !
  ! If dir = 'b', the field must be updated for a backward
  ! propagation. In taht case, the propagation step that is
  ! going to be done moves from iter*|dt| to (iter-1)*|dt|.
  !
  ! cp = (1-eta)*cpp - (eta/alpha) * <chi|V|Psi>
  ! ---------------------------------------------------------
  subroutine update_field(oct, iter, cp, gr, td, h, psi, chi, cpp, dir)
    type(oct_t), intent(in) :: oct
    integer, intent(in)        :: iter
    type(oct_control_parameters_t), intent(inout) :: cp
    type(grid_t), intent(inout)   :: gr
    type(td_t), intent(in)     :: td
    type(hamiltonian_t), intent(in) :: h
    type(states_t), intent(inout) :: psi
    type(states_t), intent(inout) :: chi
    type(oct_control_parameters_t), intent(in) :: cpp
    character(len=1),intent(in) :: dir

    type(states_t) :: oppsi
    
    CMPLX :: d1
    CMPLX, allocatable  :: dl(:), dq(:)
    integer :: ik, p, j
    FLOAT :: value

    call push_sub('propagation.update_field')
    
    ALLOCATE(dl(cp%no_parameters), cp%no_parameters)
    ALLOCATE(dq(cp%no_parameters), cp%no_parameters)

    do j = 1, cp%no_parameters
      oppsi = psi
      dl(j) = M_z0
      do ik = 1, psi%d%nik
        do p = 1, psi%nst
          oppsi%zpsi(:, :, p, ik) = M_z0
          call zvlaser_operator_linear(gr, h, psi%zpsi(:, :, p, ik), &
                                       oppsi%zpsi(:, :, p, ik), ik, laser_number = j)
          dl(j) = dl(j) + psi%occ(p, ik) * zstates_dotp(gr%m, psi%d%dim, chi%zpsi(:, :, p, ik), &
            oppsi%zpsi(:, :, p, ik))
        end do
      end do
      call states_end(oppsi)

      ! The quadratic part should only be computed if necessary.
      if(h%ep%lasers(j)%field.eq.E_FIELD_MAGNETIC ) then
        oppsi = psi
        dq(j) = M_z0
        do ik = 1, psi%d%nik
          do p = 1, psi%nst
            oppsi%zpsi(:, :, p, ik) = M_z0
            call zvlaser_operator_quadratic(gr, h, psi%zpsi(:, :, p, ik), &
                                         oppsi%zpsi(:, :, p, ik), ik, laser_number = j)
            dq(j) = dq(j) + psi%occ(p, ik) * &
                    zstates_dotp(gr%m, psi%d%dim, chi%zpsi(:, :, p, ik), &
              oppsi%zpsi(:, :, p, ik))
          end do
        end do
        call states_end(oppsi)
      else
        dq(j) = M_z0
      end if
    end do

    d1 = M_z1
    if(oct%algorithm_type .eq. oct_algorithm_zbr98) d1 = zstates_mpdotp(gr%m, psi, chi)

    select case(dir)
      case('f')
        do j = 1, cp%no_parameters
          value = (M_ONE / cp%alpha(j)) * aimag(d1*dl(j)) / &
           ( tdf(cp%td_penalty(j), iter) - M_TWO*aimag(dq(j)) )
          value = (M_ONE - oct%delta)*tdf(cpp%f(j), iter) + oct%delta * value
          call tdf_set_numerical(cp%f(j), iter, value)
          if(iter+1 <= td%max_iter + 1)  call tdf_set_numerical(cp%f(j), iter+1, value)
          if(iter+2 <= td%max_iter + 1)  call tdf_set_numerical(cp%f(j), iter+2, value)
        end do

      case('b')
        do j = 1, cp%no_parameters
          value = (M_ONE / cp%alpha(j)) * aimag(d1*dl(j)) / &
           ( tdf(cp%td_penalty(j), iter+1) - M_TWO*aimag(dq(j)) ) 
          value = (M_ONE - oct%eta)*tdf(cpp%f(j), iter+1) + oct%eta * value
          call tdf_set_numerical(cp%f(j), iter+1, value)
          if(iter > 0) call tdf_set_numerical(cp%f(j), iter, value)
          if(iter-1 > 0) call tdf_set_numerical(cp%f(j), iter-1, value)
        end do
    end select

    deallocate(dl, dq)
    call pop_sub()
  end subroutine update_field
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_prop_init(prop, dirname, niter, number_checkpoints)
    type(oct_prop_t), intent(inout) :: prop
    character(len=*), intent(in)    :: dirname
    integer,          intent(in)    :: niter
    integer,          intent(in)    :: number_checkpoints

    integer :: j

    prop%dirname = 'opt-control/'//trim(dirname)
    call io_mkdir(trim(prop%dirname))
    prop%niter = niter
    prop%number_checkpoints = number_checkpoints

    ALLOCATE(prop%iter(prop%number_checkpoints+2), prop%number_checkpoints+2)
    prop%iter(1) = 0
    do j = 1, prop%number_checkpoints
      prop%iter(j+1) = nint( real(niter)/(prop%number_checkpoints+1) * j)
    end do
    prop%iter(prop%number_checkpoints+2) = niter

  end subroutine oct_prop_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_prop_end(prop)
    type(oct_prop_t), intent(inout) :: prop

    deallocate(prop%iter)
    ! This routine should maybe delete the files?

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

    do j = 1, prop%number_checkpoints + 2
     if(prop%iter(j) .eq. iter) then
       call states_copy(stored_st, psi)
       write(filename,'(a,i4.4)') trim(prop%dirname)//'/', j
       call restart_read(trim(filename), stored_st, gr, geo, ierr)
       prev_overlap = zstates_mpdotp(gr%m, stored_st, stored_st)
       overlap = zstates_mpdotp(gr%m, stored_st, psi)
       if( abs(overlap - prev_overlap) > WARNING_THRESHOLD ) then
          write(message(1), '(a,es13.4)') &
            "WARNING: forward-backward propagation produced an error of", abs(overlap-prev_overlap)
          write(message(2), '(a,i8)') "Iter = ", iter
          call write_warning(2)
       end if
       ! Restore state only if the number of checkpoints is larger than zero.
       if(prop%number_checkpoints > 0) then
         call states_end(psi)
         call states_copy(psi, stored_st)
       end if
       call states_end(stored_st)
     end if
    end do

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

    do j = 1, prop%number_checkpoints + 2
     if(prop%iter(j) .eq. iter) then
       write(filename,'(a,i4.4)') trim(prop%dirname)//'/', j
       call restart_read(trim(filename), psi, gr, geo, ierr)
     end if
    end do

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

    do j = 1, prop%number_checkpoints + 2
      if(prop%iter(j) .eq. iter) then
        write(filename,'(a,i4.4)') trim(prop%dirname)//'/', j
        call restart_write(trim(filename), psi, gr, ierr, iter)
      end if
    end do

  end subroutine oct_prop_output
  ! ---------------------------------------------------------


end module opt_control_propagation_m

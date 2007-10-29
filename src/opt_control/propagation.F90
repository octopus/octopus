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
  use lib_oct_m
  use messages_m
  use units_m
  use grid_m
  use mesh_function_m
  use states_m
  use excited_states_m
  use hamiltonian_m
  use system_m
  use timedep_m
  use td_rti_m
  use td_write_m
  use opt_control_constants_m
  use opt_control_target_m
  use opt_control_parameters_m
  use lasers_m
  use v_ks_m
  use tdf_m

  implicit none

  private
  public :: propagate_forward,  &
            propagate_backward, &
            fwd_step,           &
            bwd_step

  contains

  ! ---------------------------------------------------------
  ! Performs a full propagation of state psi_n, with the laser
  ! field specified in hamiltonian h. If write_iter is present
  ! and is set to .true., writes down through the td_write
  ! module.
  ! ---------------------------------------------------------
  subroutine propagate_forward(sys, h, td, target, psi_n, write_iter)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h
    type(td_t),          intent(inout) :: td
    type(target_t),      intent(inout) :: target
    type(states_t),      intent(inout) :: psi_n
    logical, optional, intent(in)      :: write_iter

    integer :: ii, i
    logical :: write_iter_ = .false.
    type(grid_t),  pointer :: gr
    type(td_write_t)           :: write_handler

    call push_sub('propagation.propagate_forward')

    write_iter_ = .false.
    if(present(write_iter)) write_iter_ = write_iter

    gr => sys%gr

    if(write_iter_) then
      call td_write_init(write_handler, gr, sys%st, sys%geo, &
        (td%move_ions>0), h%ep%with_gauge_field, td%iter, td%dt)
    end if

    ! setup the hamiltonian
    call states_calc_dens(psi_n, NP_PART, psi_n%rho)
    call v_ks_calc(gr, sys%ks, h, psi_n)
    call td_rti_run_zero_iter(h, td%tr)

    ii = 1
    do i = 1, td%max_iter
      ! time iterate wavefunctions
      call td_rti_dt(sys%ks, h, gr, psi_n, td%tr, i*td%dt, td%dt, td%max_iter)

      ! if td_target
      if(target%targetmode==oct_targetmode_td) call calc_tdfitness(target, gr, psi_n, i)
      
      ! update
      call states_calc_dens(psi_n, NP_PART, psi_n%rho)
      call v_ks_calc(gr, sys%ks, h, psi_n)
      call hamiltonian_energy(h, sys%gr, sys%geo, psi_n, -1)

      ! only write in final run
      if(write_iter_) then
        ii = ii + 1 
        if(ii==sys%outp%iter+1 .or. i == td%max_iter) then ! output 
          if(i == td%max_iter) sys%outp%iter = ii - 1 
          ii = 1 
          call td_write_iter(write_handler, gr, psi_n, h, sys%geo, td%kick, td%dt, i)
          call td_write_data(write_handler, gr, psi_n, h, sys%outp, sys%geo, td%dt, i) 
        end if
      end if
    end do

    if(write_iter_) call td_write_end(write_handler)
    call pop_sub()
  end subroutine propagate_forward


  ! ---------------------------------------------------------
  ! Performs a full bacward propagation of state psi_n, with the
  ! external fields specified in hamiltonian h.
  ! ---------------------------------------------------------
  subroutine propagate_backward(sys, h, td, psi_n)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h
    type(td_t),          intent(inout) :: td
    type(states_t), intent(inout)      :: psi_n

    integer :: i
    type(grid_t),  pointer :: gr

    call push_sub('propagation.propagate_backward')
    
    message(1) = "Info: Backward propagating Chi"
    call write_info(1)

    gr => sys%gr

    ! setup the hamiltonian
    call states_calc_dens(psi_n, NP_PART, psi_n%rho)
    call v_ks_calc(gr, sys%ks, h, psi_n)
    call td_rti_run_zero_iter(h, td%tr)

    do i = td%max_iter, 1, -1
      call td_rti_dt(sys%ks, h, gr, psi_n, td%tr, (i-1)*td%dt, -td%dt, td%max_iter)
      call states_calc_dens(psi_n, NP_PART, psi_n%rho)
      call v_ks_calc(gr, sys%ks, h, psi_n)
    end do
    
    call pop_sub()
  end subroutine propagate_backward


  ! /*---------------------------------------------------------
  ! Performs a forward propagation on the state psi and on the
  ! lagrange-multiplier state chi. It also updates the control
  ! parameter par,  according to the following scheme:
  ! 
  ! |chi> --> U[par_tmp](T, 0)|chi>
  ! par = par[|psi>, |chi>]
  ! |psi> --> U[par](T, 0)|psi>
  !
  ! Note that the control parameters "par" are updated "on the
  ! fly", so that the propagation of psi is performed with the
  ! "new" parameters.
  ! */--------------------------------------------------------
  subroutine fwd_step(oct, sys, td, h, target, par, par_tmp, chi, psi)
    type(oct_t), intent(in) :: oct
    type(system_t), intent(inout)                 :: sys
    type(td_t), intent(inout)                     :: td
    type(hamiltonian_t), intent(inout)            :: h
    type(target_t), intent(inout)          :: target
    type(oct_control_parameters_t), intent(inout) :: par
    type(oct_control_parameters_t), intent(in)    :: par_tmp
    type(states_t), intent(inout)                 :: chi
    type(states_t), intent(inout)                 :: psi

    integer :: i
    type(states_t) :: psi2
    type(grid_t), pointer :: gr

    call push_sub('propagation.fwd_step')

    gr => sys%gr
    psi2 = psi
    
    ! setup forward propagation
    call states_densities_init(psi, gr, sys%geo)
    call states_calc_dens(psi, NP_PART, psi%rho)
    call v_ks_calc(gr, sys%ks, h, psi)
    call td_rti_run_zero_iter(h, td%tr)

    message(1) = "Info: Propagating forward"
    call write_info(1)

    do i = 1, td%max_iter

      if(target%targetmode==oct_targetmode_td) then
        call calc_inh(psi2, gr, target, td%dt*i, chi)
        call parameters_to_h(par, h%ep)
        call states_calc_dens(psi2, NP_PART, psi2%rho)
        call v_ks_calc(gr, sys%ks, h, psi2)
        call td_rti_dt(sys%ks, h, gr, psi2, td%tr, abs(i*td%dt), abs(td%dt), td%max_iter)
      end if

      call update_field(oct, i, par, gr, td, h, psi, chi, par_tmp, dir = 'f')
      call prop_iter(i, gr, sys%ks, h, td, par_tmp, chi)
      call prop_iter(i, gr, sys%ks, h, td, par, psi)

      if(target%targetmode==oct_targetmode_td) call calc_tdfitness(target, gr, psi, i) 

    end do

    call states_calc_dens(psi, NP_PART, psi%rho)
    call v_ks_calc(gr, sys%ks, h, psi)

    call states_end(psi2)
    nullify(gr)
    call pop_sub()
  end subroutine fwd_step


  ! --------------------------------------------------------
  ! Performs a backward propagation on the state psi and on the
  ! lagrange-multiplier state chi, according to the following
  ! scheme:
  !
  ! |psi> --> U[par](0, T)|psi>
  ! par_tmp = par_tmp[|psi>, |chi>]
  ! |chi> --> U[par_tmp](0, T)|chi>
  ! --------------------------------------------------------
  subroutine bwd_step(oct, sys, td, h, target, par, par_tmp, chi, psi) 
    type(oct_t), intent(in) :: oct
    type(system_t), intent(inout) :: sys
    type(td_t), intent(inout)                     :: td
    type(hamiltonian_t), intent(inout)            :: h
    type(target_t), intent(inout)          :: target
    type(oct_control_parameters_t), intent(in)    :: par
    type(oct_control_parameters_t), intent(inout) :: par_tmp
    type(states_t), intent(inout)                 :: chi
    type(states_t), intent(inout)                 :: psi

    integer :: i
    type(grid_t), pointer :: gr

    call push_sub('propagation.bwd_step')

    message(1) = "Info: Propagating backward"
    call write_info(1)

    gr => sys%gr

    call states_calc_dens(psi, NP_PART, psi%rho)
    call v_ks_calc(gr, sys%ks, h, psi)
    call td_rti_run_zero_iter(h, td%tr)

    td%dt = -td%dt
    do i = td%max_iter, 1, -1
      if(target%targetmode==oct_targetmode_td) &
        call calc_inh(psi, gr, target, (i-1)*td%dt, chi)
      call update_field(oct, i, par_tmp, gr, td, h, psi, chi, par, dir = 'b')
      call prop_iter(i-1, gr, sys%ks, h, td, par_tmp, chi)
      call prop_iter(i-1, gr, sys%ks, h, td, par, psi)
    end do
    td%dt = -td%dt

    call states_calc_dens(psi, NP_PART, psi%rho)
    call v_ks_calc(gr, sys%ks, h, psi)

    call pop_sub()
  end subroutine bwd_step


  ! ----------------------------------------------------------
  ! Performs one propagation step:
  ! o If td%dt > 0, from iter*td%dt-td%dt to iter*td%dt
  ! o If td%dt < 0, from iter*|td%dt|+|td%dt| to iter*|td%dt|
  ! ----------------------------------------------------------
  subroutine prop_iter(iter, gr, ks, h, td, par, st)
    integer, intent(in)                        :: iter
    type(grid_t), intent(inout)                :: gr
    type(v_ks_t), intent(inout)                :: ks
    type(hamiltonian_t), intent(inout)         :: h
    type(td_t), intent(inout)                  :: td
    type(oct_control_parameters_t), intent(in) :: par
    type(states_t), intent(inout)              :: st

    integer :: j
    call push_sub('propagation.prop_iter')

    do j = iter - 2, iter + 2
      if(j >= 0 .and. j<=td%max_iter) then
        call parameters_to_h_val(par, h%ep, j+1)
      end if
    end do
    if(.not.h%ip_app) then
      call states_calc_dens(st, NP_PART, st%rho)
      call v_ks_calc(gr, ks, h, st)
    end if
    call td_rti_dt(ks, h, gr, st, td%tr, abs(iter*td%dt), td%dt, td%max_iter)

    call pop_sub()
  end subroutine prop_iter


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

end module opt_control_propagation_m

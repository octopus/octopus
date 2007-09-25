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

    integer :: ierr, ii, i
    logical :: write_iter_ = .false.
    FLOAT, allocatable :: dens(:,:)
    type(grid_t),  pointer :: gr
    type(td_write_t)           :: write_handler

    call push_sub('propagation.propagate_forward')

    message(1) = "Info: Forward propagating Psi"
    call write_info(1)

    write_iter_ = .false.
    if(present(write_iter)) write_iter_ = write_iter

    gr => sys%gr

    if(write_iter_) then
      call td_write_init(write_handler, gr, sys%st, sys%geo, &
        (td%move_ions>0), h%ep%with_gauge_field, td%iter, td%dt)
    end if

    ALLOCATE(dens(NP_PART, psi_n%d%nspin), NP_PART*psi_n%d%nspin) 

    ! setup the hamiltonian
    call states_calc_dens(psi_n, NP_PART, dens)
    psi_n%rho = dens
    call v_ks_calc(gr, sys%ks, h, psi_n)
    call td_rti_run_zero_iter(h, td%tr)

    ii = 1

    do i = 1, td%max_iter
      ! time iterate wavefunctions
      call td_rti_dt(sys%ks, h, gr, psi_n, td%tr, abs(i*td%dt), abs(td%dt), td%max_iter)

      ! if td_target
      if(target%targetmode==oct_targetmode_td) call calc_tdfitness(target, gr, psi_n, i)
      
      ! update
      call states_calc_dens(psi_n, NP_PART, dens)
      psi_n%rho = dens
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
    deallocate(dens)
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
    FLOAT, allocatable :: dens(:,:)
    type(grid_t),  pointer :: gr

    call push_sub('propagation.propagate_backward')
    
    message(1) = "Info: Backward propagating Chi"
    call write_info(1)

    gr => sys%gr

    ALLOCATE(dens(NP_PART, psi_n%d%nspin), NP_PART*psi_n%d%nspin) 

    ! setup the hamiltonian
    call states_calc_dens(psi_n, NP_PART, dens)
    psi_n%rho = dens
    call v_ks_calc(gr, sys%ks, h, psi_n)!, calc_eigenval=.true.)
    call td_rti_run_zero_iter(h, td%tr)

    td%dt = -td%dt
    do i = td%max_iter-1, 0, -1
      ! time iterate wavefunctions
      call td_rti_dt(sys%ks, h, gr, psi_n, td%tr, abs(i*td%dt), td%dt, td%max_iter)
      ! update
      call states_calc_dens(psi_n, NP_PART, dens)
      psi_n%rho = dens
      call v_ks_calc(gr, sys%ks, h, psi_n)!, calc_eigenval=.true.)
      
    end do
    td%dt = -td%dt
    
    deallocate(dens)
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
  subroutine fwd_step(method, sys, td, h, target, par, par_tmp, chi, psi)
    integer, intent(in)                           :: method
    type(system_t), intent(inout)                 :: sys
    type(td_t), intent(inout)                     :: td
    type(hamiltonian_t), intent(inout)            :: h
    type(target_t), intent(inout)          :: target
    type(oct_control_parameters_t), intent(inout) :: par
    type(oct_control_parameters_t), intent(in)    :: par_tmp
    type(states_t), intent(inout)                 :: chi
    type(states_t), intent(inout)                 :: psi

    integer :: i, nspin
    FLOAT, allocatable :: dens_tmp(:,:)
    type(states_t) :: psi2
    type(grid_t), pointer :: gr

    call push_sub('propagation.fwd_step')

    gr => sys%gr
    nspin = psi%d%nspin
    ALLOCATE(dens_tmp(NP_PART, nspin), NP_PART*nspin) 

    psi2 = psi
    
    ! setup forward propagation
    call states_densities_init(psi, gr, sys%geo)
    call states_calc_dens(psi, NP_PART, dens_tmp)
    psi%rho = dens_tmp
    call v_ks_calc(gr, sys%ks, h, psi)!, calc_eigenval=.true.)
    call td_rti_run_zero_iter(h, td%tr)

    message(1) = "Info: Propagating forward"
    call write_info(1)

    do i = 1, td%max_iter

      if(target%targetmode==oct_targetmode_td) then
        call calc_inh(psi2, gr, target, i, td%max_iter, td%dt, chi)
        call parameters_to_h(par, h%ep)
        call states_calc_dens(psi2, NP_PART, dens_tmp)
        psi2%rho = dens_tmp
        call v_ks_calc(gr, sys%ks, h, psi2)!, calc_eigenval=.true.)
        call td_rti_dt(sys%ks, h, gr, psi2, td%tr, abs(i*td%dt), abs(td%dt), td%max_iter)
      end if

      call update_field(method, i-1, par, gr, td, h, psi, chi)
      call prop_iter_fwd(i, gr, sys%ks, h, td, par_tmp, chi)
      call prop_iter_fwd(i, gr, sys%ks, h, td, par, psi)

      if(target%targetmode==oct_targetmode_td) call calc_tdfitness(target, gr, psi, i) 

    end do

    call states_end(psi2)

    nullify(gr)
    deallocate(dens_tmp) 
     
    call pop_sub()
  end subroutine fwd_step


  ! ----------------------------------------------------------
  ! Performs one step of forward propagation.
  ! ----------------------------------------------------------
  subroutine prop_iter_fwd(iter, gr, ks, h, td, par, st)
    integer, intent(in)                        :: iter
    type(grid_t), intent(inout)                :: gr
    type(v_ks_t), intent(inout)                :: ks
    type(hamiltonian_t), intent(inout)         :: h
    type(td_t), intent(inout)                  :: td
    type(oct_control_parameters_t), intent(in) :: par
    type(states_t), intent(inout)              :: st

    FLOAT, allocatable :: dens_tmp(:,:)
    integer :: nspin
    
    call push_sub('propagation.prop_iter_fwd_')

    nspin = st%d%nspin
    ALLOCATE(dens_tmp(NP_PART, nspin), NP_PART*nspin) 

    call parameters_to_h(par, h%ep)
    call states_calc_dens(st, NP_PART, dens_tmp)
    st%rho = dens_tmp
    call v_ks_calc(gr, ks, h, st)
    call td_rti_dt(ks, h, gr, st, td%tr, abs(iter*td%dt), abs(td%dt), td%max_iter)

    deallocatE(dens_tmp)
    call pop_sub()
  end subroutine prop_iter_fwd


  ! --------------------------------------------------------
  ! Performs a forward propagation on the state psi and on the
  ! lagrange-multiplier state chi, according to the following
  ! scheme:
  !
  ! |psi> --> U[par](0, T)|psi>
  ! par_tmp = par_tmp[|psi>, |chi>]
  ! |chi> --> U[par_tmp](0, T)|chi>
  ! --------------------------------------------------------
  subroutine bwd_step(method, sys, td, h, target, par, par_tmp, chi, psi) 
    integer, intent(in) :: method
    type(system_t), intent(inout) :: sys
    type(td_t), intent(inout)                     :: td
    type(hamiltonian_t), intent(inout)            :: h
    type(target_t), intent(inout)          :: target
    type(oct_control_parameters_t), intent(in)    :: par
    type(oct_control_parameters_t), intent(inout) :: par_tmp
    type(states_t), intent(inout)                 :: chi
    type(states_t), intent(inout)                 :: psi

    integer :: i, nspin
    FLOAT, allocatable :: dens_tmp(:,:)
    type(grid_t), pointer :: gr

    call push_sub('propagation.bwd_step')

    message(1) = "Info: Propagating backward"
    call write_info(1)

    gr => sys%gr
    nspin = psi%d%nspin

    ALLOCATE(dens_tmp(NP_PART, nspin), NP_PART*nspin) 

    ! setup backward propagation
    call states_calc_dens(chi, NP_PART, dens_tmp)
    chi%rho = dens_tmp
    call v_ks_calc(gr, sys%ks, h, chi)
    call td_rti_run_zero_iter(h, td%tr)

    td%dt = -td%dt
    do i = td%max_iter-1, 0, -1
      if(target%targetmode==oct_targetmode_td) &
        call calc_inh(psi, gr, target, i, td%max_iter, td%dt, chi)
      call update_field(method, i+1, par_tmp, gr, td, h, psi, chi)
      call prop_iter_bwd(i, gr, sys%ks, h, td, par_tmp, chi)
      call prop_iter_bwd(i, gr, sys%ks, h, td, par, psi)
    end do
    td%dt = -td%dt

    deallocate(dens_tmp)
    call pop_sub()
  end subroutine bwd_step


  ! ----------------------------------------------------------
  ! Performs one step of forward propagation.
  ! ----------------------------------------------------------
  subroutine prop_iter_bwd(iter, gr, ks, h, td, par, st)
    integer, intent(in) :: iter
    type(grid_t), intent(inout) :: gr
    type(v_ks_t), intent(inout) :: ks
    type(hamiltonian_t), intent(inout) :: h
    type(td_t), intent(inout) :: td
    type(oct_control_parameters_t), intent(in) :: par
    type(states_t), intent(inout) :: st

    FLOAT, allocatable :: dens_tmp(:,:)

    integer :: nspin
    
    call push_sub('propagation.prop_iter_bwd')

    nspin = st%d%nspin

    ALLOCATE(dens_tmp(NP_PART, nspin), NP_PART*nspin) 

    call parameters_to_h(par, h%ep)
    call states_calc_dens(st, NP_PART, dens_tmp)
    st%rho = dens_tmp
    call v_ks_calc(gr, ks, h, st)
    call td_rti_dt(ks, h, gr, st, td%tr, abs(iter*td%dt), td%dt, td%max_iter)

    deallocate(dens_tmp)
    call pop_sub()
  end subroutine prop_iter_bwd


  ! ---------------------------------------------------------
  ! Calculates the value of the control parameters at iteration
  ! iter, from the state psi and the Lagrange-multiplier chi.
  ! ---------------------------------------------------------
  subroutine update_field(algorithm_type, iter, cp, gr, td, h, psi, chi)
    integer, intent(in) :: algorithm_type
    integer, intent(in)        :: iter
    type(oct_control_parameters_t), intent(inout) :: cp
    type(grid_t), intent(inout)   :: gr
    type(td_t), intent(in)     :: td
    type(hamiltonian_t), intent(in) :: h
    type(states_t), intent(inout) :: psi
    type(states_t), intent(inout) :: chi

    type(states_t) :: oppsi
    
    CMPLX :: d1
    CMPLX, allocatable  :: d2(:), dq(:)
    integer :: ik, p, dim, i, j
    FLOAT :: value

    call push_sub('propagation.update_field')
    
    ALLOCATE(d2(cp%no_parameters), cp%no_parameters)
    ALLOCATE(dq(cp%no_parameters), cp%no_parameters)

    do j = 1, cp%no_parameters
      oppsi = psi
      do ik = 1, psi%d%nik
        do p = 1, psi%nst
          oppsi%zpsi(:, :, p, ik) = M_z0
          call zvlaser_operator_linear(gr, h, psi%zpsi(:, :, p, ik), &
                                       oppsi%zpsi(:, :, p, ik), ik, laser_number = j)
        end do
      end do
      d2(j) = zstates_mpmatrixelement(gr%m, chi, psi, oppsi)
      call states_end(oppsi)

      ! The quadratic part should only be computed if necessary.
      if(h%ep%lasers(j)%field.eq.E_FIELD_MAGNETIC ) then
        oppsi = psi
        do ik = 1, psi%d%nik
          do p = 1, psi%nst
            oppsi%zpsi(:, :, p, ik) = M_z0
            call zvlaser_operator_quadratic(gr, h, psi%zpsi(:, :, p, ik), &
                                         oppsi%zpsi(:, :, p, ik), ik, laser_number = j)
          end do
        end do
        dq(j) = zstates_mpmatrixelement(gr%m, chi, psi, oppsi)
        call states_end(oppsi)
      else
        dq(j) = M_z0
      end if
    end do

    d1 = M_z1
    if(algorithm_type .eq. oct_algorithm_zbr98) d1 = zstates_mpdotp(gr%m, psi, chi)

    do j = 1, cp%no_parameters
      value = aimag(d1*d2(j)) / ( tdf(cp%td_penalty(j), iter+1) - M_TWO*aimag(dq(j)) ) 
      call tdf_set_numerical(cp%f(j), iter+1, value)
      i = int(sign(M_ONE, td%dt))
      if(iter==0.or.iter==td%max_iter) then
        call tdf_set_numerical(cp%f(j), iter+1, value)
      else
        value = M_HALF*( M_FOUR*tdf(cp%f(j), iter+1) - M_TWO*tdf(cp%f(j), iter-i+1))
        call tdf_set_numerical(cp%f(j), iter+i+1, value)
      end if
    end do

    deallocate(d2, dq)
    call pop_sub()
  end subroutine update_field

end module opt_control_propagation_m

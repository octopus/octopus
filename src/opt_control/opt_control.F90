!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or
!! modify
!! it under the terms of the GNU General Public License as published
!! by
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
!! $Id$

#include "global.h"

module opt_control_m
  use lalg_basic_m
  use lalg_adv_m
  use loct_math_m
  use excited_states_m
  use datasets_m
  use varinfo_m
  use global_m
  use filter_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lasers_m
  use loct_parser_m
  use loct_m
  use messages_m
  use simul_box_m
  use mesh_m
  use mesh_function_m
  use functions_m
  use h_sys_output_m
  use geometry_m
  use states_m
  use string_m
  use system_m
  use td_rti_m
  use td_write_m
  use timedep_m
  use units_m
  use v_ks_m
  use external_pot_m
  use restart_m
  use tdf_m
  use opt_control_constants_m
  use opt_control_propagation_m
  use opt_control_parameters_m
  use opt_control_iter_m
  use opt_control_output_m
  use opt_control_target_m

  implicit none

  private
  public :: opt_control_run

  ! For the algorithm_direct scheme:
  type(oct_control_parameters_t) :: par_
  type(system_t), pointer :: sys_
  type(hamiltonian_t), pointer :: h_
  type(td_t), pointer :: td_
  type(target_t), pointer :: target_
  type(states_t), pointer :: psi_
  type(oct_iterator_t), pointer :: iterator_
  type(oct_t), pointer :: oct_

contains

  ! ---------------------------------------------------------
  subroutine direct_opt_calc(n, x, f)
    integer, intent(in)  :: n
    real(8), intent(in)  :: x(n)
    real(8), intent(out) :: f

    integer :: j
    FLOAT :: j1, fluence
    type(states_t) :: psi

    call parameters_x_to_par(par_, x)
    call parameters_to_realtime(par_)
    call parameters_to_h(par_, h_%ep)
    call states_copy(psi, psi_)
    call propagate_forward(sys_, h_, td_, target_, psi)
    call parameters_set_rep(par_)

    j1 = j1_functional(sys_%gr, psi, target_)
    fluence = parameters_fluence(par_)
    f = - j1 !+ par_%alpha(1) * (fluence - par_%targetfluence)**2

    call states_end(psi)
  end subroutine direct_opt_calc
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine direct_opt_write_info(geom_iter, n, energy, maxdx, x)
    integer, intent(in) :: geom_iter, n
    real(8), intent(in) :: energy, maxdx
    real(8), intent(in) :: x(n)

    FLOAT :: fluence, j, j1, j2

    call parameters_x_to_par(par_, x)

    iterator_%ctr_iter = geom_iter
    call iteration_manager_direct(energy, par_, iterator_, maxdx)
    if(oct_%dump_intermediate) call iterator_write(iterator_, psi_, par_, sys_%gr, sys_%outp)

  end subroutine direct_opt_write_info
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine opt_control_run(sys, h)
    type(system_t), target,      intent(inout) :: sys
    type(hamiltonian_t), target, intent(inout) :: h

    type(oct_t), target            :: oct
    type(oct_iterator_t), target   :: iterator
    type(td_t), target             :: td
    type(states_t)                 :: psi
    type(states_t), target         :: initial_st
    type(target_t), target         :: target
    type(filter_t)                 :: filter
    type(oct_control_parameters_t) :: par, par_new, par_prev
    logical                        :: stop_loop
    FLOAT                          :: j1
    type(oct_prop_t)               :: prop_chi, prop_psi;

    call push_sub('opt_control.opt_control_run')

    call io_mkdir('opt-control')

    call td_init(sys, h, td)
    if(h%theory_level .ne. INDEPENDENT_PARTICLES ) then
      call td_rti_set_scf_prop(td%tr)
    end if

    call states_allocate_wfns(sys%st, sys%gr%m, M_CMPLX)

    call oct_read_inp(oct)

    call parameters_set_initial(par, h%ep, sys%gr%m, td%dt, td%max_iter, &
                                oct%mode_fixed_fluence, oct%mode_basis_set)

    call oct_iterator_init(iterator, par)

    if(oct%use_mixing) call parameters_mixing_init(par)

    call filter_init(td%max_iter, td%dt, filter)
    call filter_write(filter)

    call initial_state_init(sys%gr, sys%geo, sys%st, initial_st)
    call target_init(sys%gr, sys%geo, sys%st, td, target)

    call check_faulty_runmodes(oct, sys, h, target, td%tr)

    call h_sys_output_states(initial_st, sys%gr, 'opt-control/initial', sys%outp)
    call target_output(target, sys%gr, 'opt-control/target', sys%outp)

    ! psi is the "working state".
    call states_copy(psi, initial_st)

    ! mode switcher
    select case(oct%algorithm)
      case(oct_algorithm_zbr98)
        message(1) = "Info: Starting OCT iteration using scheme: ZBR98"
        call write_info(1)
        call scheme_zbr98
      case(oct_algorithm_wg05)
        message(1) = "Info: Starting OCT iteration using scheme: WG05"
        call write_info(1)
        call scheme_wg05
      case(oct_algorithm_zr98)
        message(1) = "Info: Starting OCT iteration using scheme: ZR98"
        call write_info(1)
        call scheme_mt03
      case(oct_algorithm_mt03)
        message(1) = "Info: Starting OCT iteration using scheme: MT03"
        call write_info(1)
        call scheme_mt03
      case(oct_algorithm_krotov)
        message(1) = "Info: Starting OCT iteration using scheme: KROTOV"
        call write_info(1)
        call scheme_mt03
      case(oct_algorithm_str_iter)
        message(1) = "Info: Starting OCT iterations using scheme: STRAIGHT ITERATION"
        call write_info(1)
        call scheme_straight_iteration
      case(oct_algorithm_direct)
        message(1) = "Info: Starting OCT iterations using scheme: DIRECT OPTIMIZATION"
        call write_info(1)
        call scheme_direct
    case default
      call input_error('OCTScheme')
    end select

    ! Some informative output.
    call oct_output(iterator, sys%gr, sys%outp, psi)

    ! do final test run: propagate initial state with optimal field
    call oct_finalcheck(oct, initial_st, target, iterator%best_par, sys, h, td)

    ! clean up
    call parameters_end(par)
    call oct_iterator_end(iterator)
    if(oct%use_mixing) call parameters_mixing_end()
    call filter_end(filter)
    call td_end(td)
    call states_end(psi)
    call states_end(initial_st)
    call target_end(target)
   
    call pop_sub()

  contains


    ! ---------------------------------------------------------
    subroutine scheme_straight_iteration
      call push_sub('opt_control.scheme_mt03')

      call oct_prop_init(prop_chi, "chi", td%max_iter, oct%number_checkpoints)
      call oct_prop_init(prop_psi, "psi", td%max_iter, oct%number_checkpoints)

      call parameters_set_rep(par)
      call parameters_copy(par_new, par)
      ctr_loop: do
        call parameters_copy(par_prev, par)
        call f_striter(oct, sys, h, td, filter, psi, initial_st, target, par, prop_psi, prop_chi, j1)
        if(oct%dump_intermediate) call iterator_write(iterator, psi, par_prev, sys%gr, sys%outp)
        stop_loop = iteration_manager(j1, par_prev, par, iterator)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if(oct%use_mixing) then
          call parameters_mixing(iterator%ctr_iter, par_prev, par, par_new)
          call parameters_copy(par, par_new)
          if(oct%mode_fixed_fluence) call parameters_set_fluence(par)
        end if
      end do ctr_loop

      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call parameters_end(par_new)
      call parameters_end(par_prev)
      call pop_sub()
    end subroutine scheme_straight_iteration
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_mt03
      call push_sub('opt_control.scheme_mt03')

      call oct_prop_init(prop_chi, "chi", td%max_iter, oct%number_checkpoints)
      call oct_prop_init(prop_psi, "psi", td%max_iter, oct%number_checkpoints)

      call parameters_copy(par_new, par)
      ctr_loop: do
        call parameters_copy(par_prev, par)
        call f_iter(oct, iterator, sys, h, td, psi, initial_st, target, par, prop_psi, prop_chi, j1)
        if(oct%dump_intermediate) call iterator_write(iterator, psi, par, sys%gr, sys%outp)
        stop_loop = iteration_manager(j1, par, par_prev, iterator)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if( oct%use_mixing .and. (iterator%ctr_iter > 1) ) then
          ! We do not mix if it is the first iteration, since in that case f_iter only propagates
          ! with the input field, and does not generate any output field.
          call parameters_mixing(iterator%ctr_iter-1, par_prev, par, par_new)
          call parameters_copy(par, par_new)
        end if
      end do ctr_loop

      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call parameters_end(par_new)
      call parameters_end(par_prev)
      call pop_sub()
    end subroutine scheme_mt03
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_wg05
      call push_sub('opt_control.scheme_wg05')

      call oct_prop_init(prop_chi, "chi", td%max_iter, oct%number_checkpoints)
      call oct_prop_init(prop_psi, "psi", td%max_iter, oct%number_checkpoints)

      if (oct%mode_fixed_fluence) then
        par%alpha(1) = (M_ONE/sqrt(par%targetfluence)) * sqrt ( parameters_fluence(par) )
      end if

      call parameters_copy(par_new, par)      
      ctr_loop: do
        call parameters_copy(par_prev, par)
        call f_wg05(oct, iterator, sys, h, td, filter, psi, initial_st, target, par, prop_psi, prop_chi, j1)
        if(oct%dump_intermediate) call iterator_write(iterator, psi, par, sys%gr, sys%outp)
        stop_loop = iteration_manager(j1, par, par_prev, iterator)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if(oct%use_mixing) then
          call parameters_mixing(iterator%ctr_iter, par_prev, par, par_new)
          call parameters_copy(par, par_new)
        end if
      end do ctr_loop

      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call parameters_end(par_new)
      call parameters_end(par_prev)
      call pop_sub()
    end subroutine scheme_wg05


    ! ---------------------------------------------------------
    subroutine scheme_zbr98
      call push_sub('opt_control.scheme_zbr98')

      call oct_prop_init(prop_chi, "chi", td%max_iter, oct%number_checkpoints)
      call oct_prop_init(prop_psi, "psi", td%max_iter, oct%number_checkpoints)

      call parameters_copy(par_prev, par)
      call propagate_forward(sys, h, td, target, psi, prop_psi)
      j1 = j1_functional(sys%gr, psi, target)
      if(oct%dump_intermediate) call iterator_write(iterator, psi, par, sys%gr, sys%outp)
      stop_loop = iteration_manager(j1, par, par_prev, iterator)
      call parameters_end(par_prev)

      call parameters_copy(par_new, par)
      ctr_loop: do
        call parameters_copy(par_prev, par)
        call f_zbr98(oct, sys, h, td, psi, initial_st, target, prop_psi, prop_chi, par)
        j1 = j1_functional(sys%gr, psi, target)
        if(oct%dump_intermediate) call iterator_write(iterator, psi, par, sys%gr, sys%outp)
        stop_loop = iteration_manager(j1, par, par_prev, iterator)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if(oct%use_mixing) then
          call parameters_mixing(iterator%ctr_iter-1, par_prev, par, par_new)
          call parameters_copy(par, par_new)
        end if
      end do ctr_loop

      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call parameters_end(par_new)
      call parameters_end(par_prev)
      call pop_sub()
    end subroutine scheme_zbr98
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_direct
      integer :: ierr
      FLOAT :: minvalue, step
      FLOAT, allocatable :: x(:)
      call push_sub('opt_control.scheme_direct')

      call parameters_set_rep(par)

      ALLOCATE(x(par%nfreqs-1), par%nfreqs-1)

      ! Set the module pointers, so that the calc_point and write_iter_info routines
      ! can use them.
      call parameters_copy(par_, par)
      sys_      => sys
      h_        => h
      td_       => td
      target_   => target
      psi_      => initial_st
      iterator_ => iterator
      oct_      => oct

      ! Do a zero iteration, with the input field.
      ! (This could be removed in a final version, since the minimization algorithm itself
      ! has to repeat this run. Now it stays for clarity).
      call parameters_to_realtime(par)
      call states_copy(psi, initial_st)
      call parameters_to_h(par, h%ep)
      call propagate_forward(sys, h, td, target, psi)
      call parameters_set_rep(par)
      if(oct%dump_intermediate) call iterator_write(iterator, psi, par, sys%gr, sys%outp)

      j1 = j1_functional(sys%gr, psi, target)
      call iteration_manager_direct(j1, par, iterator)

      call parameters_par_to_x(par, x)
      step = oct%direct_step * M_PI
      ierr = loct_minimize_direct(MINMETHOD_NMSIMPLEX, par%nfreqs-1, x(1), step,&
               iterator%eps, iterator%ctr_iter_max, &
               direct_opt_calc, direct_opt_write_info, minvalue)

      if(ierr.ne.0) then
        if(ierr <= 1024) then
          message(1) = "Error occurred during the GSL minimization procedure:"
          call loct_strerror(ierr, message(2))
          call write_fatal(2)
        else
          message(1) = "The OCT direct optimization did not meet the convergence criterion."
          call write_info(1)
        end if
      end if

      deallocate(x)
      call pop_sub()
    end subroutine scheme_direct
    ! ---------------------------------------------------------

  end subroutine opt_control_run


  ! ---------------------------------------------------------
  subroutine f_zbr98(oct, sys, h, td, psi, initial_st, target, prop_psi, prop_chi, par)
    type(oct_t), intent(in)                       :: oct
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: h
    type(td_t), intent(inout)                     :: td
    type(states_t), intent(inout)                 :: psi
    type(states_t), intent(in)                    :: initial_st
    type(target_t), intent(inout)                 :: target
    type(oct_prop_t), intent(inout)               :: prop_psi, prop_chi
    type(oct_control_parameters_t), intent(inout) :: par

    type(states_t) :: chi
    type(oct_control_parameters_t) :: par_chi

    call push_sub('opt_control.f_zbr98')

    call parameters_copy(par_chi, par)

    if(oct%use_mixing) then
      call states_end(psi)
      call states_copy(psi, initial_st)
      call parameters_to_h(par, h%ep)
      call propagate_forward(sys, h, td, target, psi, prop_psi)
    end if

    call states_copy(chi, target%st)
    call bwd_step(oct, sys, td, h, target, par, par_chi, chi, prop_chi, prop_psi)

    call states_end(psi)
    call states_copy(psi, initial_st)
    call fwd_step(oct, sys, td, h, target, par, par_chi, psi, prop_chi, prop_psi)

    call states_end(chi)
    call parameters_end(par_chi)
    call pop_sub()
  end subroutine f_zbr98


  ! ---------------------------------------------------------
  subroutine f_wg05(oct, iterator, sys, h, td, filter, psi, initial_st, target, par, prop_psi, prop_chi, j1)
    type(oct_t), intent(in)                       :: oct
    type(oct_iterator_t), intent(in)              :: iterator
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: h
    type(td_t), intent(inout)                     :: td
    type(filter_t), intent(inout)                 :: filter
    type(states_t), intent(inout)                 :: psi
    type(states_t), intent(in)                    :: initial_st
    type(target_t), intent(inout)                 :: target
    type(oct_control_parameters_t), intent(inout) :: par
    type(oct_prop_t), intent(inout)               :: prop_psi, prop_chi
    FLOAT, intent(out)                            :: j1

    integer :: j
    FLOAT :: fluence, old_penalty, new_penalty
    type(states_t) :: chi
    type(oct_control_parameters_t) :: parp

    call push_sub('opt_control.f_wg05')

    if( iterator%ctr_iter .eq. 0) then
      call states_end(psi)
      call states_copy(psi, initial_st)
      call parameters_to_h(par, h%ep)
      call propagate_forward(sys, h, td, target, psi, prop_psi)
      j1 = j1_functional(sys%gr, psi, target)
      call pop_sub()
      return
    end if

    call parameters_copy(parp, par)

    call states_copy(chi, psi)
    call calc_chi(oct, sys%gr, target, psi, chi)
    call bwd_step(oct, sys, td, h, target, par, parp, chi, prop_chi, prop_psi)

    do j = 1, parp%no_parameters
      call filter_apply(parp%f(j), filter)
    end do

    ! recalc field
    if (oct%mode_fixed_fluence) then
      fluence = parameters_fluence(parp) 
      old_penalty = par%alpha(1)
      new_penalty = sqrt( fluence * old_penalty**2 / par%targetfluence )
      par%alpha(:) = new_penalty
      parp%alpha(:) = new_penalty
      call parameters_set_fluence(parp)
    end if

    call parameters_end(par)
    call parameters_copy(par, parp)
    call parameters_apply_envelope(par)

    call states_end(psi)
    call states_copy(psi, initial_st)
    call parameters_to_h(par, h%ep)
    call propagate_forward(sys, h, td, target, psi, prop_psi)

    j1 = j1_functional(sys%gr, psi, target)

    call states_end(chi)
    call parameters_end(parp)
    call pop_sub()
  end subroutine f_wg05
  ! ---------------------------------------------------------


   ! ---------------------------------------------------------
  subroutine f_striter(oct, sys, h, td, filter, psi, initial_st, target, par, prop_psi, prop_chi, j1)
    type(oct_t), intent(in)                       :: oct
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: h
    type(td_t), intent(inout)                     :: td
    type(filter_t), intent(inout)                 :: filter
    type(states_t), intent(inout)                 :: psi
    type(states_t), intent(in)                    :: initial_st
    type(target_t), intent(inout)                 :: target
    type(oct_control_parameters_t), intent(inout) :: par
    type(oct_prop_t), intent(inout)               :: prop_psi, prop_chi
    FLOAT, intent(out)                            :: j1

    integer :: j
    type(states_t) :: chi
    type(oct_control_parameters_t) :: par_chi

    call parameters_to_realtime(par)

    call parameters_copy(par_chi, par)

    ! First, a forward propagation with the input field.
    call states_end(psi)
    call states_copy(psi, initial_st)
    call parameters_to_h(par, h%ep)
    call propagate_forward(sys, h, td, target, psi, prop_psi)

    ! Check the performance.
    j1 = j1_functional(sys%gr, psi, target)

    ! Set the boundary condition for the backward propagation.
    call states_copy(chi, psi)
    call calc_chi(oct, sys%gr, target, psi, chi)

    ! Backward propagation, while at the same time finding the output field, 
    ! which is placed at par_chi
    call bwd_step_2(oct, sys, td, h, target, par, par_chi, chi, prop_chi, prop_psi)

    call parameters_set_rep(par_chi)

    do j = 1, par_chi%no_parameters
      call filter_apply(par_chi%f(j), filter)
    end do

    ! Fix the fluence, in case it is needed.
    if(oct%mode_fixed_fluence) call parameters_set_fluence(par_chi)

    ! Copy par_chi to par
    call parameters_end(par)
    call parameters_copy(par, par_chi)
    call pop_sub()
  end subroutine f_striter
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine f_iter(oct, iterator, sys, h, td, psi, initial_st, target, par, prop_psi, prop_chi, j1)
    type(oct_t), intent(in)                       :: oct
    type(oct_iterator_t), intent(in)              :: iterator
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: h
    type(td_t), intent(inout)                     :: td
    type(states_t), intent(inout)                 :: psi
    type(states_t), intent(in)                    :: initial_st
    type(target_t), intent(inout)                 :: target
    type(oct_control_parameters_t), intent(inout) :: par
    type(oct_prop_t), intent(inout)               :: prop_psi, prop_chi
    FLOAT, intent(out)                            :: j1

    type(states_t) :: chi
    type(oct_control_parameters_t) :: par_chi

    call push_sub('opt_control.f_zbr98')

    if( iterator%ctr_iter .eq. 0) then
      call states_end(psi)
      call states_copy(psi, initial_st)
      call parameters_to_h(par, h%ep)
      call propagate_forward(sys, h, td, target, psi, prop_psi)
      j1 = j1_functional(sys%gr, psi, target)
      call pop_sub()
      return
    end if

    call parameters_copy(par_chi, par)

    if(oct%use_mixing .and. (iterator%ctr_iter > 1) ) then
      ! No need to do this auxiliary propagation if no mixing has been done previously.
      call states_end(psi)
      call states_copy(psi, initial_st)
      call parameters_to_h(par, h%ep)
      call propagate_forward(sys, h, td, target, psi, prop_psi)
    end if
    
    call states_copy(chi, psi)
    call calc_chi(oct, sys%gr, target, psi, chi)
    call bwd_step(oct, sys, td, h, target, par, par_chi, chi, prop_chi, prop_psi)

    call states_end(psi)
    call states_copy(psi, initial_st)
    call fwd_step(oct, sys, td, h, target, par, par_chi, psi, prop_chi, prop_psi)

    j1 = j1_functional(sys%gr, psi, target)

    call states_end(chi)
    call parameters_end(par_chi)
    call pop_sub()
  end subroutine f_iter
  ! ---------------------------------------------------------


#include "read.F90"
#include "aux.F90"
#include "defstates.F90"
#include "finalcheck.F90"

end module opt_control_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

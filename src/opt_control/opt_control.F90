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
  use h_sys_output_m
  use geometry_m
  use states_m
  use states_dim_m
  use states_calc_m
  use string_m
  use system_m
  use td_m
  use td_rti_m
  use td_write_m
  use units_m
  use v_ks_m
  use external_pot_m
  use restart_m
  use opt_control_global_m
  use opt_control_propagation_m
  use opt_control_parameters_m
  use opt_control_iter_m
  use opt_control_target_m
  use opt_control_initst_m
#if defined(HAVE_NEWUOA)
  use newuoa_m
#endif
  use profiling_m

  implicit none

  private
  public :: opt_control_run

  ! Module variables
  type(filter_t), save       :: filter
  type(oct_t), save          :: oct
  type(oct_iterator_t), save :: iterator
  type(target_t), save       :: target
  type(states_t), save       :: initial_st
  

  ! For the algorithm_direct scheme:
  type(oct_control_parameters_t), save :: par_
  type(system_t), pointer :: sys_
  type(hamiltonian_t), pointer :: hm_
  type(td_t), pointer :: td_

contains


  ! ---------------------------------------------------------
  subroutine direct_opt_calc(n, x, f)
    integer, intent(in)  :: n
    REAL_DOUBLE, intent(in)  :: x(n)
    REAL_DOUBLE, intent(out) :: f

    FLOAT :: j1, delta
    type(states_t) :: psi
    type(oct_control_parameters_t) :: par_new

    call parameters_set_theta(par_, x)
    call parameters_theta_to_basis(par_)

    if(oct%delta == M_ZERO) then
      ! We only need the value of the target functional.
      call states_copy(psi, initial_st)
      call propagate_forward(sys_, hm_, td_, par_, target, psi)
      f = - j1_functional(target, sys_%gr, psi) - parameters_j2(par_)
      if(oct%dump_intermediate) call iterator_write(iterator, par_)
      call iteration_manager_direct(real(-f, REAL_PRECISION), par_, iterator)
      call states_end(psi)
    else
      call parameters_copy(par_new, par_)
      call f_striter(sys_, hm_, td_, par_new, j1)
      delta = parameters_diff(par_, par_new)
      f = - oct%eta * j1 + oct%delta * delta
      if(oct%dump_intermediate) call iterator_write(iterator, par_)
      call iteration_manager_direct(real(-f, REAL_PRECISION), par_, iterator, delta)
      call parameters_end(par_new)
    end if


  end subroutine direct_opt_calc
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine direct_opt_write_info(iter, n, val, maxdx, x)
    integer, intent(in) :: iter, n
    REAL_DOUBLE, intent(in) :: val, maxdx
    REAL_DOUBLE, intent(in) :: x(n)

    FLOAT :: fluence, j1, j2, j

    call parameters_set_theta(par_, x)
    call parameters_theta_to_basis(par_)

    j = - val
    fluence = parameters_fluence(par_)
    j2 = parameters_j2(par_)
    j1 = j - j2

    write(message(1), '(a,i5)') 'Direct optimization iteration #', iter
    call messages_print_stress(stdout, trim(message(1)))

    write(message(1), '(6x,a,f12.5)')    " => J1       = ", j1
    write(message(2), '(6x,a,f12.5)')    " => J        = ", j
    write(message(3), '(6x,a,f12.5)')    " => J2       = ", j2
    write(message(4), '(6x,a,f12.5)')    " => Fluence  = ", fluence
    write(message(5), '(6x,a,f12.5)')    " => Delta    = ", maxdx
    call write_info(5)
    call messages_print_stress(stdout)

  end subroutine direct_opt_write_info
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine opt_control_run(sys, hm)
    type(system_t), target,      intent(inout) :: sys
    type(hamiltonian_t), target, intent(inout) :: hm

    type(td_t), target             :: td
    type(oct_control_parameters_t) :: par, par_new, par_prev
    logical                        :: stop_loop
    FLOAT                          :: j1
    type(oct_prop_t)               :: prop_chi, prop_psi;

    call push_sub('opt_control.opt_control_run')

    call io_mkdir('opt-control')


    ! Initialize the time propagator.
    call td_init(td, sys, hm)
    if(hm%theory_level .ne. INDEPENDENT_PARTICLES ) then
      call td_rti_set_scf_prop(td%tr)
    end if


    ! Read general information about how the OCT run will be made, from inp file.
    call oct_read_inp(oct)


    ! Read info about, and prepare, the control functions (a.k.a. parameters).
    call parameters_mod_init(hm%ep, td%dt, td%max_iter, oct%mode_fixed_fluence, &
                             (oct%ctr_function_rep == oct_ctr_function_parametrized) )
    call parameters_init(par, td%dt, td%max_iter)
    call parameters_set(par, hm%ep)
      ! This prints the initial control parameters, exactly as described in the inp file,
      ! that is, without applying any envelope or filter.
    call parameters_write('opt-control/initial_laser_inp', par)
    call parameters_prepare_initial(par)
    call parameters_to_h(par, hm%ep)
    call messages_print_stress(stdout, "TD ext. fields after processing")
    call laser_write_info(hm%ep%lasers, stdout)
    call messages_print_stress(stdout)
    call parameters_write('opt-control/initial_laser', par)


    ! Startup of the iterator data type (takes care of counting iterations, stopping, etc).
    call oct_iterator_init(iterator, par)


    ! Initialization of the propagation_m module.
    call propagation_mod_init(td%max_iter, oct%eta, oct%delta, oct%number_checkpoints, &
      (oct%algorithm == oct_algorithm_zbr98))


    ! If mixing is required, the mixing machinery has to be initialized -- inside the parameters module.
    if(oct%use_mixing) call parameters_mixing_init(par)


    ! If filters are to be used, they also have to be initialized.
    call filter_init(td%max_iter, td%dt, filter)
    call filter_write(filter)


    ! Figure out how is the starting wave function(s), and the target.
    call initial_state_init(sys%gr, sys%geo, sys%st, initial_st)
    call target_init(sys%gr, sys%geo, sys%st, td, parameters_w0(par), target)


    ! Sanity checks.
    call check_faulty_runmodes(sys, hm, td%tr)


    ! Informative output.
    call h_sys_output_states(initial_st, sys%gr, sys%geo, 'opt-control/initial', sys%outp)
    call target_output(target, sys%gr, 'opt-control/target', sys%geo, sys%outp)


    ! mode switcher; here is where the real run is made.
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
        message(1) = "Info: Starting OCT iterations using scheme: DIRECT OPTIMIZATION (NELDER-MEAD)"
        call write_info(1)
        call scheme_direct
      case(oct_algorithm_newuoa)
        message(1) = "Info: Starting OCT iterations using scheme: DIRECT OPTIMIZATION (NEWUOA)"
        call write_info(1)
        call scheme_newuoa
    case default
      call input_error('OCTScheme')
    end select

    ! do final test run: propagate initial state with optimal field
    call oct_finalcheck(sys, hm, td)

    ! clean up
    call parameters_end(par)
    call oct_iterator_end(iterator)
    if(oct%use_mixing) call parameters_mixing_end()
    call filter_end(filter)
    call td_end(td)
    call states_end(initial_st)
    call target_end(target)
    call parameters_mod_close()
   
    call pop_sub()

  contains


    ! ---------------------------------------------------------
    subroutine scheme_straight_iteration
      call push_sub('opt_control.scheme_mt03')

      call parameters_set_rep(par)
      call parameters_copy(par_new, par)
      ctr_loop: do
        call parameters_copy(par_prev, par)
        call f_striter(sys, hm, td, par, j1)
        if(oct%dump_intermediate) call iterator_write(iterator, par_prev)
        stop_loop = iteration_manager(j1, par_prev, par, iterator)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if(oct%use_mixing) then
          call parameters_mixing(oct_iterator_current(iterator), par_prev, par, par_new)
          call parameters_copy(par, par_new)
          if(oct%mode_fixed_fluence) call parameters_set_fluence(par)
        end if
      end do ctr_loop

      call parameters_end(par_new)
      call parameters_end(par_prev)
      call pop_sub()
    end subroutine scheme_straight_iteration
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_mt03
      type(states_t) :: psi
      call push_sub('opt_control.scheme_mt03')

      call states_copy(psi, initial_st)
      call oct_prop_init(prop_chi, "chi")
      call oct_prop_init(prop_psi, "psi")

      call parameters_copy(par_new, par)
      ctr_loop: do
        call parameters_copy(par_prev, par)
        call f_iter(sys, hm, td, psi, par, prop_psi, prop_chi, j1)
        if(oct%dump_intermediate) call iterator_write(iterator, par)
        stop_loop = iteration_manager(j1, par, par_prev, iterator)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if( oct%use_mixing .and. (oct_iterator_current(iterator) > 1) ) then
          ! We do not mix if it is the first iteration, since in that case f_iter only propagates
          ! with the input field, and does not generate any output field.
          call parameters_mixing(oct_iterator_current(iterator) - 1, par_prev, par, par_new)
          call parameters_copy(par, par_new)
        end if
      end do ctr_loop

      call states_end(psi)
      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call parameters_end(par_new)
      call parameters_end(par_prev)
      call pop_sub()
    end subroutine scheme_mt03
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_wg05
      type(states_t) :: psi
      call push_sub('opt_control.scheme_wg05')

      call states_copy(psi, initial_st)
      call oct_prop_init(prop_chi, "chi")
      call oct_prop_init(prop_psi, "psi")

      if (oct%mode_fixed_fluence) then
        call parameters_set_alpha(par, sqrt( parameters_fluence(par) / parameters_targetfluence()))
      end if

      call parameters_copy(par_new, par)      
      ctr_loop: do
        call parameters_copy(par_prev, par)
        call f_wg05(sys, hm, td, psi, par, prop_psi, prop_chi, j1)
        if(oct%dump_intermediate) call iterator_write(iterator, par)
        stop_loop = iteration_manager(j1, par, par_prev, iterator)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if(oct%use_mixing) then
          call parameters_mixing(oct_iterator_current(iterator), par_prev, par, par_new)
          call parameters_copy(par, par_new)
        end if
      end do ctr_loop

      call states_end(psi)
      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call parameters_end(par_new)
      call parameters_end(par_prev)
      call pop_sub()
    end subroutine scheme_wg05


    ! ---------------------------------------------------------
    subroutine scheme_zbr98
      type(states_t) :: psi
      call push_sub('opt_control.scheme_zbr98')

      call states_copy(psi, initial_st)
      call oct_prop_init(prop_chi, "chi")
      call oct_prop_init(prop_psi, "psi")

      call parameters_copy(par_prev, par)
      call propagate_forward(sys, hm, td, par, target, psi, prop_psi)
      j1 = j1_functional(target, sys%gr, psi)
      if(oct%dump_intermediate) call iterator_write(iterator, par)
      stop_loop = iteration_manager(j1, par, par_prev, iterator)
      call parameters_end(par_prev)

      call parameters_copy(par_new, par)
      ctr_loop: do
        call parameters_copy(par_prev, par)
        call f_zbr98(sys, hm, td, psi, prop_psi, prop_chi, par)
        j1 = j1_functional(target, sys%gr, psi)
        if(oct%dump_intermediate) call iterator_write(iterator, par)
        stop_loop = iteration_manager(j1, par, par_prev, iterator)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if(oct%use_mixing) then
          call parameters_mixing(oct_iterator_current(iterator) - 1, par_prev, par, par_new)
          call parameters_copy(par, par_new)
        end if
      end do ctr_loop

      call states_end(psi)
      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call parameters_end(par_new)
      call parameters_end(par_prev)
      call pop_sub()
    end subroutine scheme_zbr98
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_direct
      integer :: ierr, maxiter
      REAL_DOUBLE :: minvalue, step
      REAL_DOUBLE, allocatable :: x(:)
      FLOAT :: f
      integer :: dim
      type(states_t) :: psi

      call push_sub('opt_control.scheme_direct')

      call parameters_set_rep(par)

      dim = parameters_dof(par)

      ALLOCATE(x(dim), dim)

      call states_copy(psi, initial_st)
      call propagate_forward(sys, hm, td, par, target, psi)
      f = - j1_functional(target, sys%gr, psi) - parameters_j2(par)
      if(oct%dump_intermediate) call iterator_write(iterator, par)
      call iteration_manager_direct(-f, par, iterator)
      call states_end(psi)

      if(oct%random_initial_guess) call parameters_randomize(par)

      ! Set the module pointers, so that the direct_opt_calc and direct_opt_write_info routines
      ! can use them.
      call parameters_copy(par_, par)
      sys_      => sys
      hm_       => hm
      td_       => td

      call parameters_basis_to_theta(par)
      call parameters_get_theta(par, x)

      step = oct%direct_step * M_PI
      maxiter = oct_iterator_maxiter(iterator)

      ierr = loct_minimize_direct(MINMETHOD_NMSIMPLEX, dim, x(1), step,&
               real(oct_iterator_tolerance(iterator), 8), maxiter, &
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


    ! ---------------------------------------------------------
    subroutine scheme_newuoa
#if defined(HAVE_NEWUOA)
      integer :: iprint, npt, maxfun, sizeofw, dim
      REAL_DOUBLE :: rhobeg, rhoend
      REAL_DOUBLE, allocatable :: x(:), w(:), xl(:), xu(:)
      FLOAT :: f
      type(states_t) :: psi
      call push_sub('opt_control.scheme_direct')

      call parameters_set_rep(par)

      dim = parameters_dof(par)
      ALLOCATE(x(dim), dim)
      ALLOCATE(xl(dim), dim)
      ALLOCATE(xu(dim), dim)
      call parameters_bounds(par, xl, xu)

      call states_copy(psi, initial_st)
      call propagate_forward(sys, hm, td, par, target, psi)
      f = - j1_functional(target, sys%gr, psi) - parameters_j2(par)
      if(oct%dump_intermediate) call iterator_write(iterator, par)
      call iteration_manager_direct(-f, par, iterator)      
      call states_end(psi)

      if(oct%random_initial_guess) call parameters_randomize(par)

      ! Set the module pointers, so that the calc_point and write_iter_info routines
      ! can use them.
      call parameters_copy(par_, par)
      sys_      => sys
      hm_        => hm
      td_       => td

      call parameters_basis_to_theta(par)
      call parameters_get_theta(par, x)

      iprint = 2
      npt = 2*dim + 1
      rhoend = oct_iterator_tolerance(iterator)
      rhobeg = oct%direct_step * M_PI
      maxfun = oct_iterator_maxiter(iterator)
      sizeofw = (npt + 13)*(npt + dim) + 3 * dim*(dim + 3)/2 
      ALLOCATE(w(sizeofw), sizeofw)
      w = M_ZERO
      call newuoa(dim, npt, x, rhobeg, rhoend, iprint, maxfun, w, direct_opt_calc)

      deallocate(x, xl, xu, w)
      call pop_sub()
#endif
    end subroutine scheme_newuoa
    ! ---------------------------------------------------------

  end subroutine opt_control_run
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine f_zbr98(sys, hm, td, psi, prop_psi, prop_chi, par)
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: hm
    type(td_t), intent(inout)                     :: td
    type(states_t), intent(inout)                 :: psi
    type(oct_prop_t), intent(inout)               :: prop_psi, prop_chi
    type(oct_control_parameters_t), intent(inout) :: par

    type(states_t) :: chi
    type(oct_control_parameters_t) :: par_chi

    call push_sub('opt_control.f_zbr98')

    call parameters_copy(par_chi, par)

    if(oct%use_mixing) then
      call states_end(psi)
      call states_copy(psi, initial_st)
      call propagate_forward(sys, hm, td, par, target, psi, prop_psi)
    end if

    call target_get_state(target, chi)
    call bwd_step(sys, td, hm, target, par, par_chi, chi, prop_chi, prop_psi)

    call states_end(psi)
    call states_copy(psi, initial_st)
    call fwd_step(sys, td, hm, target, par, par_chi, psi, prop_chi, prop_psi)

    call states_end(chi)
    call parameters_end(par_chi)
    call pop_sub()
  end subroutine f_zbr98


  ! ---------------------------------------------------------
  subroutine f_wg05(sys, hm, td, psi, par, prop_psi, prop_chi, j1)
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: hm
    type(td_t), intent(inout)                     :: td
    type(states_t), intent(inout)                 :: psi
    type(oct_control_parameters_t), intent(inout) :: par
    type(oct_prop_t), intent(inout)               :: prop_psi, prop_chi
    FLOAT, intent(out)                            :: j1

    FLOAT :: new_penalty
    type(states_t) :: chi
    type(oct_control_parameters_t) :: parp

    call push_sub('opt_control.f_wg05')

    if( oct_iterator_current(iterator) .eq. 0) then
      call states_end(psi)
      call states_copy(psi, initial_st)
      call propagate_forward(sys, hm, td, par, target, psi, prop_psi)
      j1 = j1_functional(target, sys%gr, psi)
      call pop_sub(); return
    end if

    call parameters_copy(parp, par)

    call states_copy(chi, psi)
    call calc_chi(target, sys%gr, psi, chi)
    call bwd_step(sys, td, hm, target, par, parp, chi, prop_chi, prop_psi)

    call parameters_filter(parp, filter)

    ! recalc field
    if (oct%mode_fixed_fluence) then
      new_penalty = parameters_alpha(par, 1) * sqrt( parameters_fluence(parp) / parameters_targetfluence() )
      call parameters_set_alpha(parp, new_penalty)
      call parameters_set_fluence(parp)
    end if

    call parameters_end(par)
    call parameters_copy(par, parp)
    call parameters_apply_envelope(par)

    call states_end(psi)
    call states_copy(psi, initial_st)
    call propagate_forward(sys, hm, td, par, target, psi, prop_psi)

    j1 = j1_functional(target, sys%gr, psi)

    call states_end(chi)
    call parameters_end(parp)
    call pop_sub()
  end subroutine f_wg05
  ! ---------------------------------------------------------


   ! ---------------------------------------------------------
  subroutine f_striter(sys, hm, td, par, j1)
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: hm
    type(td_t), intent(inout)                     :: td
    type(oct_control_parameters_t), intent(inout) :: par
    FLOAT, intent(out)                            :: j1

    type(states_t) :: chi
    type(states_t) :: psi
    type(oct_control_parameters_t) :: par_chi
    type(oct_prop_t)               :: prop_chi, prop_psi;

    call oct_prop_init(prop_chi, "chi")
    call oct_prop_init(prop_psi, "psi")

    call parameters_to_realtime(par)

    call parameters_copy(par_chi, par)

    ! First, a forward propagation with the input field.
    call states_copy(psi, initial_st)
    call propagate_forward(sys, hm, td, par, target, psi, prop_psi)

    ! Check the performance.
    j1 = j1_functional(target, sys%gr, psi)

    ! Set the boundary condition for the backward propagation.
    call states_copy(chi, psi)
    call calc_chi(target, sys%gr, psi, chi)

    ! Backward propagation, while at the same time finding the output field, 
    ! which is placed at par_chi
    call bwd_step_2(sys, td, hm, target, par, par_chi, chi, prop_chi, prop_psi)

    call parameters_set_rep(par_chi)

    ! Fix the fluence, in case it is needed.
    if(oct%mode_fixed_fluence) call parameters_set_fluence(par_chi)

    ! Copy par_chi to par
    call parameters_end(par)
    call parameters_copy(par, par_chi)

    call parameters_end(par_chi)
    call oct_prop_end(prop_chi)
    call oct_prop_end(prop_psi)
    call pop_sub()
  end subroutine f_striter
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine f_iter(sys, hm, td, psi, par, prop_psi, prop_chi, j1)
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: hm
    type(td_t), intent(inout)                     :: td
    type(states_t), intent(inout)                 :: psi
    type(oct_control_parameters_t), intent(inout) :: par
    type(oct_prop_t), intent(inout)               :: prop_psi, prop_chi
    FLOAT, intent(out)                            :: j1

    type(states_t) :: chi
    type(oct_control_parameters_t) :: par_chi

    call push_sub('opt_control.f_zbr98')

    if( oct_iterator_current(iterator) .eq. 0) then
      call states_end(psi)
      call states_copy(psi, initial_st)
      call propagate_forward(sys, hm, td, par, target, psi, prop_psi)
      j1 = j1_functional(target, sys%gr, psi)
      call pop_sub(); return
    end if

    call parameters_copy(par_chi, par)

    if(oct%use_mixing .and. (oct_iterator_current(iterator) > 1) ) then
      ! No need to do this auxiliary propagation if no mixing has been done previously.
      call states_end(psi)
      call states_copy(psi, initial_st)
      call propagate_forward(sys, hm, td, par, target, psi, prop_psi)
    end if
    
    call states_copy(chi, psi)
    call calc_chi(target, sys%gr, psi, chi)
    call bwd_step(sys, td, hm, target, par, par_chi, chi, prop_chi, prop_psi)

    call states_end(psi)
    call states_copy(psi, initial_st)
    call fwd_step(sys, td, hm, target, par, par_chi, psi, prop_chi, prop_psi)

    j1 = j1_functional(target, sys%gr, psi)

    call states_end(chi)
    call parameters_end(par_chi)
    call pop_sub()
  end subroutine f_iter
  ! ---------------------------------------------------------


#include "check_input.F90"
#include "finalcheck.F90"

end module opt_control_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

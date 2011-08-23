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

!> This module contains the main procedure ("opt_control_run") that is 
!! used when optimal control runs are requested.
module opt_control_m
  use messages_m
  use controlfunction_m
  use excited_states_m
  use exponential_m
  use filter_m
  use geometry_m
  use global_m
  use grid_m
  use output_m
  use hamiltonian_m
  use io_m
  use lasers_m
  use loct_m
  use loct_math_m
  use mesh_m
#if defined(HAVE_NEWUOA)
  use newuoa_m
#endif
  use opt_control_global_m
  use opt_control_state_m
  use opt_control_propagation_m
  use opt_control_iter_m
  use opt_control_target_m
  use opt_control_initst_m
  use profiling_m
  use propagator_m
  use restart_m
  use simul_box_m
  use states_m
  use states_dim_m
  use system_m
  use td_m
  
  implicit none

  private
  public :: opt_control_run,                  &
            opt_control_cg_calc,              &
            opt_control_cg_write_info,     &
            opt_control_direct_calc,          &
            opt_control_direct_message_info,  &
            opt_control_function_forward


  ! Module variables
  type(filter_t), save       :: filter
  type(oct_t), save          :: oct
  type(oct_iterator_t), save :: iterator
  type(target_t), save       :: target
  type(states_t), save       :: initial_st
  

  ! For the direct, newuoa, and cg schemes:
  type(controlfunction_t), save :: par_
  type(system_t), pointer :: sys_
  type(hamiltonian_t), pointer :: hm_
  type(td_t), pointer :: td_
  FLOAT, allocatable :: x_(:)
  integer :: index_

contains


  !> This is the main procedure for all types of optimal control runs.
  !! It is called from the "run" procedure in the "run_m" module.
  subroutine opt_control_run(sys, hm)
    type(system_t), target,      intent(inout) :: sys
    type(hamiltonian_t), target, intent(inout) :: hm

    type(td_t), target             :: td
    type(controlfunction_t)        :: par, par_new, par_prev
    logical                        :: stop_loop
    FLOAT                          :: j1
    type(oct_prop_t)               :: prop_chi, prop_psi

    PUSH_SUB(opt_control_run)

    ! Creates a directory where the optimal control stuff will be written. The name of the directory
    ! is stored in the preprocessor macro OCT_DIR, which should be defined in src/include/global.h
    call io_mkdir(OCT_DIR)

    ! Initializes the time propagator. Then, it forces the propagation to be self consistent, in case
    ! the theory level is not "independent particles".
    call td_init(td, sys, hm)
    if(hm%theory_level .ne. INDEPENDENT_PARTICLES ) call propagator_set_scf_prop(td%tr)

    ! Read general information about how the OCT run will be made, from inp file. "oct_read_inp" is
    ! in the opt_control_global_m module (like the definition of the oct_t data type)
    call oct_read_inp(oct)

    ! Read info about, and prepare, the control functions
    call controlfunction_mod_init(hm%ep, td%dt, td%max_iter, oct%mode_fixed_fluence, &
                             (oct%ctr_function_rep == oct_ctr_function_parametrized) )
    call controlfunction_init(par, td%dt, td%max_iter)
    call controlfunction_set(par, hm%ep)
      ! This prints the initial control parameters, exactly as described in the inp file,
      ! that is, without applying any envelope or filter.
    call controlfunction_write(OCT_DIR//'initial_laser_inp', par)
    call controlfunction_prepare_initial(par)
    call controlfunction_to_h(par, hm%ep)
    call messages_print_stress(stdout, "TD ext. fields after processing")
    call laser_write_info(hm%ep%lasers, stdout)
    call messages_print_stress(stdout)
    call controlfunction_write(OCT_DIR//'initial_laser', par)


    ! Startup of the iterator data type (takes care of counting iterations, stopping, etc).
    call oct_iterator_init(iterator, par)


    ! Initialization of the propagation_m module.
    call propagation_mod_init(td%max_iter, oct%eta, oct%delta, oct%number_checkpoints, &
      (oct%algorithm == oct_algorithm_zbr98), (oct%algorithm == oct_algorithm_cg) )


    ! If mixing is required, the mixing machinery has to be initialized -- inside the controlfunction_m module.
    if(oct%use_mixing) call controlfunction_mixing_init(par)


    ! If filters are to be used, they also have to be initialized.
    call filter_init(td%max_iter, td%dt, filter)
    call filter_write(filter)


    ! Figure out the starting wavefunction(s), and the target.
    call initial_state_init(sys, hm, initial_st)
    call target_init(sys%gr, sys%geo, initial_st, td, controlfunction_w0(par), target, oct)

    ! Sanity checks.
    call check_faulty_runmodes(sys, hm, td%tr)


    ! Informative output.
    call output_states(initial_st, sys%gr, sys%geo, OCT_DIR//'initial', sys%outp)
    call target_output(target, sys%gr, OCT_DIR//'target', sys%geo, sys%outp)


    ! mode switcher; here is where the real run is made.
    select case(oct%algorithm)
      case(oct_algorithm_zbr98)
        message(1) = "Info: Starting OCT iteration using scheme: ZBR98"
        call messages_info(1)
        call scheme_zbr98
      case(oct_algorithm_wg05)
        message(1) = "Info: Starting OCT iteration using scheme: WG05"
        call messages_info(1)
        call scheme_wg05
      case(oct_algorithm_zr98)
        message(1) = "Info: Starting OCT iteration using scheme: ZR98"
        call messages_info(1)
        call scheme_mt03
      case(oct_algorithm_mt03)
        message(1) = "Info: Starting OCT iteration using scheme: MT03"
        call messages_info(1)
        call scheme_mt03
      case(oct_algorithm_krotov)
        message(1) = "Info: Starting OCT iteration using scheme: KROTOV"
        call messages_info(1)
        call scheme_mt03
      case(oct_algorithm_str_iter)
        message(1) = "Info: Starting OCT iterations using scheme: STRAIGHT ITERATION"
        call messages_info(1)
        call scheme_straight_iteration
      case(oct_algorithm_cg)
        message(1) = "Info: Starting OCT iterations using scheme: CONJUGATE GRADIENTS"
        call messages_info(1)
        call scheme_cg
      case(oct_algorithm_direct)
        message(1) = "Info: Starting OCT iterations using scheme: DIRECT OPTIMIZATION (NELDER-MEAD)"
        call messages_info(1)
        call scheme_direct
      case(oct_algorithm_newuoa)
        message(1) = "Info: Starting OCT iterations using scheme: DIRECT OPTIMIZATION (NEWUOA)"
        call messages_info(1)
        call scheme_newuoa
    case default
      call input_error('OCTScheme')
    end select

    ! do final test run: propagate initial state with optimal field
    call oct_finalcheck(sys, hm, td)

    ! clean up
    call controlfunction_end(par)
    call oct_iterator_end(iterator)
    if(oct%use_mixing) call controlfunction_mixing_end()
    call filter_end(filter)
    call td_end(td)
    call states_end(initial_st)
    call target_end(target, oct)
    call controlfunction_mod_close()
   
    POP_SUB(opt_control_run)

  contains


    ! ---------------------------------------------------------
    subroutine scheme_straight_iteration
      PUSH_SUB(opt_control_run.scheme_straight_iteration)

      call controlfunction_set_rep(par)
      call controlfunction_copy(par_new, par)
      ctr_loop: do
        call controlfunction_copy(par_prev, par)
        call f_striter(sys, hm, td, par, j1)
        stop_loop = iteration_manager(j1, par_prev, par, iterator)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if(oct%use_mixing) then
          call controlfunction_mixing(oct_iterator_current(iterator), par_prev, par, par_new)
          call controlfunction_copy(par, par_new)
          if(oct%mode_fixed_fluence) call controlfunction_set_fluence(par)
        end if
      end do ctr_loop

      call controlfunction_end(par_new)
      call controlfunction_end(par_prev)
      POP_SUB(opt_control_run.scheme_straight_iteration)
    end subroutine scheme_straight_iteration
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_mt03
      type(states_t) :: psi
      PUSH_SUB(opt_control_run.scheme_mt03)

      call states_copy(psi, initial_st)
      call oct_prop_init(prop_chi, "chi")
      call oct_prop_init(prop_psi, "psi")

      call controlfunction_copy(par_new, par)
      ctr_loop: do
        call controlfunction_copy(par_prev, par)
        call f_iter(sys, hm, td, psi, par, prop_psi, prop_chi, j1)
        stop_loop = iteration_manager(j1, par, par_prev, iterator)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if( oct%use_mixing .and. (oct_iterator_current(iterator) > 1) ) then
          ! We do not mix if it is the first iteration, since in that case f_iter only propagates
          ! with the input field, and does not generate any output field.
          call controlfunction_mixing(oct_iterator_current(iterator) - 1, par_prev, par, par_new)
          call controlfunction_copy(par, par_new)
        end if
      end do ctr_loop

      call states_end(psi)
      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call controlfunction_end(par_new)
      call controlfunction_end(par_prev)
      POP_SUB(opt_control_run.scheme_mt03)
    end subroutine scheme_mt03
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_wg05
      type(states_t) :: psi
      PUSH_SUB(opt_control_run.scheme_wg05)

      call states_copy(psi, initial_st)
      call oct_prop_init(prop_chi, "chi")
      call oct_prop_init(prop_psi, "psi")

      if (oct%mode_fixed_fluence) then
        call controlfunction_set_alpha(par, sqrt( controlfunction_fluence(par) / controlfunction_targetfluence()))
      end if

      call controlfunction_copy(par_new, par)      
      ctr_loop: do
        call controlfunction_copy(par_prev, par)
        call f_wg05(sys, hm, td, psi, par, prop_psi, prop_chi, j1)
        stop_loop = iteration_manager(j1, par, par_prev, iterator)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if(oct%use_mixing) then
          call controlfunction_mixing(oct_iterator_current(iterator), par_prev, par, par_new)
          call controlfunction_copy(par, par_new)
        end if
      end do ctr_loop

      call states_end(psi)
      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call controlfunction_end(par_new)
      call controlfunction_end(par_prev)
      POP_SUB(opt_control_run.scheme_wg05)
    end subroutine scheme_wg05
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_zbr98
      type(states_t) :: psi
      PUSH_SUB(opt_control_run.scheme_zbr98)

      call states_copy(psi, initial_st)
      call oct_prop_init(prop_chi, "chi")
      call oct_prop_init(prop_psi, "psi")

      call controlfunction_copy(par_prev, par)
      call propagate_forward(sys, hm, td, par, target, psi, prop_psi)
      j1 = target_j1(target, sys%gr, psi)
      stop_loop = iteration_manager(j1, par, par_prev, iterator)
      call controlfunction_end(par_prev)
      if(clean_stop() .or. stop_loop) then
        call states_end(psi)
        call oct_prop_end(prop_chi)
        call oct_prop_end(prop_psi)
        POP_SUB(opt_control_run.scheme_zbr98)
        return        
      end if

      call controlfunction_copy(par_new, par)
      ctr_loop: do
        call controlfunction_copy(par_prev, par)
        call f_zbr98(sys, hm, td, psi, prop_psi, prop_chi, par)
        j1 = target_j1(target, sys%gr, psi)
        stop_loop = iteration_manager(j1, par, par_prev, iterator)
        if(clean_stop() .or. stop_loop) exit ctr_loop
        if(oct%use_mixing) then
          call controlfunction_mixing(oct_iterator_current(iterator) - 1, par_prev, par, par_new)
          call controlfunction_copy(par, par_new)
        end if
      end do ctr_loop

      call states_end(psi)
      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call controlfunction_end(par_new)
      call controlfunction_end(par_prev)
      POP_SUB(opt_control_run.scheme_zbr98)
    end subroutine scheme_zbr98
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_cg
      integer :: dof, ierr, maxiter
      REAL_DOUBLE   :: step, minvalue
      FLOAT, allocatable :: theta(:)
      REAL_DOUBLE, allocatable :: x(:)
      FLOAT   :: f
      type(states_t) :: psi
      PUSH_SUB(opt_control_run.scheme_cg)

      call controlfunction_set_rep(par)

      call states_copy(psi, initial_st)
      call propagate_forward(sys, hm, td, par, target, psi)
      f = - target_j1(target, sys%gr, psi, sys%geo) - controlfunction_j2(par)
      call iteration_manager_direct(-f, par, iterator, sys)
      call states_end(psi)
      if(oct_iterator_maxiter(iterator).eq.0) then
        ! Nothing to do.
        POP_SUB(opt_control_run.scheme_cg)
        return
      end if

      ! Set the module pointers, so that the direct_opt_calc and direct_opt_write_info routines
      ! can use them.
      call controlfunction_copy(par_, par)
      sys_      => sys
      hm_       => hm
      td_       => td

      dof = controlfunction_dof(par)
      SAFE_ALLOCATE(x(1:dof))
      SAFE_ALLOCATE(theta(1:dof))
      call controlfunction_get_theta(par, theta)
      x = theta

      step = oct%direct_step * M_PI
      maxiter = oct_iterator_maxiter(iterator) - 1

      ierr = loct_minimize(MINMETHOD_BFGS2, dof, x(1), step, &
           real(oct_iterator_tolerance(iterator), 8), real(oct_iterator_tolerance(iterator), 8), &
           maxiter, opt_control_cg_calc, opt_control_cg_write_info, minvalue)

      if(ierr.ne.0) then
        if(ierr <= 1024) then
          message(1) = "Error occurred during the GSL minimization procedure:"
          call loct_strerror(ierr, message(2))
          call messages_fatal(2, only_root_writes = .true.)
        else
          message(1) = "The CG optimization did not meet the convergence criterion."
          call messages_info(1)
        end if
      end if

      SAFE_DEALLOCATE_A(x)
      SAFE_DEALLOCATE_A(theta)
      POP_SUB(opt_control_run.scheme_cg)
    end subroutine scheme_cg
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_direct
      integer :: ierr, maxiter
      REAL_DOUBLE :: minvalue, step
      FLOAT, allocatable :: theta(:)
      REAL_DOUBLE, allocatable :: x(:)
      FLOAT :: f
      integer :: dim
      type(states_t) :: psi

      PUSH_SUB(opt_control_run.scheme_direct)

      call controlfunction_set_rep(par)
      dim = controlfunction_dof(par)

      call states_copy(psi, initial_st)
      call propagate_forward(sys, hm, td, par, target, psi)
      f = - target_j1(target, sys%gr, psi, sys%geo) - controlfunction_j2(par)
      call iteration_manager_direct(-f, par, iterator, sys)
      call states_end(psi)
      if(oct_iterator_maxiter(iterator).eq.0) then
        ! Nothing to do.
        POP_SUB(opt_control_run.scheme_cg)
        return
      end if

      SAFE_ALLOCATE(x(1:dim))
      SAFE_ALLOCATE(theta(1:dim))

      if(oct%random_initial_guess) call controlfunction_randomize(par)

      ! Set the module pointers, so that the direct_opt_calc and direct_opt_write_info routines
      ! can use them.
      call controlfunction_copy(par_, par)
      sys_      => sys
      hm_       => hm
      td_       => td

      call controlfunction_basis_to_theta(par)
      ! theta may be in single precision, whereas x is always double precision.
      call controlfunction_get_theta(par, theta)
      x = theta

      step = oct%direct_step * M_PI
      maxiter = oct_iterator_maxiter(iterator)

      ierr = loct_minimize_direct(MINMETHOD_NMSIMPLEX, dim, x(1), step, &
               real(oct_iterator_tolerance(iterator), 8), maxiter, &
               opt_control_direct_calc, opt_control_direct_message_info, minvalue)

      if(ierr.ne.0) then
        if(ierr <= 1024) then
          message(1) = "Error occurred during the GSL minimization procedure:"
          call loct_strerror(ierr, message(2))
          call messages_fatal(2)
        else
          message(1) = "The OCT direct optimization did not meet the convergence criterion."
          call messages_info(1)
        end if
      end if

      SAFE_DEALLOCATE_A(x)
      SAFE_DEALLOCATE_A(theta)
      POP_SUB(opt_control_run.scheme_direct)
    end subroutine scheme_direct
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_newuoa
#if defined(HAVE_NEWUOA)
      integer :: iprint, npt, maxfun, sizeofw, dim
      REAL_DOUBLE :: rhobeg, rhoend
      FLOAT, allocatable :: xl(:), xu(:)
      REAL_DOUBLE, allocatable :: x(:), w(:)
      FLOAT, allocatable :: theta(:)
      FLOAT :: f
      type(states_t) :: psi
      PUSH_SUB(opt_control_run.scheme_newuoa)

      call controlfunction_set_rep(par)

      call states_copy(psi, initial_st)
      call propagate_forward(sys, hm, td, par, target, psi)
      f = - target_j1(target, sys%gr, psi, sys%geo) - controlfunction_j2(par)
      call iteration_manager_direct(-f, par, iterator, sys)      
      call states_end(psi)
      if(oct_iterator_maxiter(iterator).eq.0) then
        ! Nothing to do.
        POP_SUB(opt_control_run.scheme_cg)
        return
      end if

      dim = controlfunction_dof(par)
      SAFE_ALLOCATE( x(1:dim))
      SAFE_ALLOCATE(theta(1:dim))
      SAFE_ALLOCATE(xl(1:dim))
      SAFE_ALLOCATE(xu(1:dim))
      call controlfunction_bounds(par, xl, xu)

      if(oct%random_initial_guess) call controlfunction_randomize(par)

      ! Set the module pointers, so that the calc_point and write_iter_info routines
      ! can use them.
      call controlfunction_copy(par_, par)
      sys_      => sys
      hm_        => hm
      td_       => td

      call controlfunction_basis_to_theta(par)
      call controlfunction_get_theta(par, theta)
      x = theta

      iprint = 2
      npt = 2*dim + 1
      rhoend = oct_iterator_tolerance(iterator)
      rhobeg = oct%direct_step * M_PI
      maxfun = oct_iterator_maxiter(iterator)
      sizeofw = (npt + 13)*(npt + dim) + 3 * dim*(dim + 3)/2 
      SAFE_ALLOCATE(w(1:sizeofw))
      w = M_ZERO
      call newuoa(dim, npt, x, rhobeg, rhoend, iprint, maxfun, w, opt_control_direct_calc)

      SAFE_DEALLOCATE_A(x)
      SAFE_DEALLOCATE_A(theta)
      SAFE_DEALLOCATE_A(xl)
      SAFE_DEALLOCATE_A(xu)
      SAFE_DEALLOCATE_A(w)
      POP_SUB(opt_control_run.scheme_newuoa)
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
    type(controlfunction_t), intent(inout)        :: par

    type(states_t) :: chi
    type(controlfunction_t) :: par_chi

    PUSH_SUB(f_zbr98)

    call controlfunction_copy(par_chi, par)

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
    call controlfunction_end(par_chi)
    POP_SUB(f_zbr98)
  end subroutine f_zbr98


  ! ---------------------------------------------------------
  subroutine f_wg05(sys, hm, td, psi, par, prop_psi, prop_chi, j1)
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: hm
    type(td_t), intent(inout)                     :: td
    type(states_t), intent(inout)                 :: psi
    type(controlfunction_t), intent(inout)        :: par
    type(oct_prop_t), intent(inout)               :: prop_psi, prop_chi
    FLOAT, intent(out)                            :: j1

    FLOAT :: new_penalty
    type(states_t) :: chi
    type(controlfunction_t) :: parp

    PUSH_SUB(f_wg05)

    if( oct_iterator_current(iterator) .eq. 0) then
      call states_end(psi)
      call states_copy(psi, initial_st)
      call propagate_forward(sys, hm, td, par, target, psi, prop_psi)
      j1 = target_j1(target, sys%gr, psi)
      POP_SUB(f_wg05)
      return
    end if

    call controlfunction_copy(parp, par)

    call states_copy(chi, psi)
    call target_chi(target, sys%gr, psi, chi, sys%geo)
    call bwd_step(sys, td, hm, target, par, parp, chi, prop_chi, prop_psi)

    call controlfunction_filter(parp, filter)

    ! recalc field
    if (oct%mode_fixed_fluence) then
      new_penalty = controlfunction_alpha(par, 1) * sqrt( controlfunction_fluence(parp) / controlfunction_targetfluence() )
      call controlfunction_set_alpha(parp, new_penalty)
      call controlfunction_set_fluence(parp)
    end if

    call controlfunction_end(par)
    call controlfunction_copy(par, parp)
    call controlfunction_apply_envelope(par)

    call states_end(psi)
    call states_copy(psi, initial_st)
    call propagate_forward(sys, hm, td, par, target, psi, prop_psi)

    j1 = target_j1(target, sys%gr, psi)

    call states_end(chi)
    call controlfunction_end(parp)
    POP_SUB(f_wg05)
  end subroutine f_wg05
  ! ---------------------------------------------------------


   ! ---------------------------------------------------------
  subroutine f_striter(sys, hm, td, par, j1)
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: hm
    type(td_t), intent(inout)                     :: td
    type(controlfunction_t), intent(inout)        :: par
    FLOAT, intent(out)                            :: j1

    type(states_t) :: chi
    type(states_t) :: psi
    type(controlfunction_t) :: par_chi
    type(oct_prop_t)        :: prop_chi, prop_psi;

    PUSH_SUB(f_striter)

    call oct_prop_init(prop_chi, "chi")
    call oct_prop_init(prop_psi, "psi")

    call controlfunction_to_realtime(par)

    call controlfunction_copy(par_chi, par)

    ! First, a forward propagation with the input field.
    call states_copy(psi, initial_st)
    call propagate_forward(sys, hm, td, par, target, psi, prop_psi)

    ! Check the performance.
    j1 = target_j1(target, sys%gr, psi, sys%geo)

    ! Set the boundary condition for the backward propagation.
    call states_copy(chi, psi)
    call target_chi(target, sys%gr, psi, chi, sys%geo)

    ! Backward propagation, while at the same time finding the output field, 
    ! which is placed at par_chi
    call bwd_step_2(sys, td, hm, target, par, par_chi, chi, prop_chi, prop_psi)

    call controlfunction_set_rep(par_chi)

    ! Fix the fluence, in case it is needed.
    if(oct%mode_fixed_fluence) call controlfunction_set_fluence(par_chi)

    ! Copy par_chi to par
    call controlfunction_end(par)
    call controlfunction_copy(par, par_chi)

    call controlfunction_end(par_chi)
    call oct_prop_end(prop_chi)
    call oct_prop_end(prop_psi)
    POP_SUB(f_striter)
  end subroutine f_striter
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine f_iter(sys, hm, td, psi, par, prop_psi, prop_chi, j1)
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: hm
    type(td_t), intent(inout)                     :: td
    type(states_t), intent(inout)                 :: psi
    type(controlfunction_t), intent(inout)        :: par
    type(oct_prop_t), intent(inout)               :: prop_psi, prop_chi
    FLOAT, intent(out)                            :: j1

    type(states_t) :: chi
    type(controlfunction_t) :: par_chi

    PUSH_SUB(f_iter)

    if( oct_iterator_current(iterator) .eq. 0) then
      call states_end(psi)
      call states_copy(psi, initial_st)
      call propagate_forward(sys, hm, td, par, target, psi, prop_psi)
      j1 = target_j1(target, sys%gr, psi)
      POP_SUB(f_iter)
      return
    end if

    call controlfunction_copy(par_chi, par)

    if(oct%use_mixing .and. (oct_iterator_current(iterator) > 1) ) then
      ! No need to do this auxiliary propagation if no mixing has been done previously.
      call states_end(psi)
      call states_copy(psi, initial_st)
      call propagate_forward(sys, hm, td, par, target, psi, prop_psi)
    end if
    
    call states_copy(chi, psi)
    call target_chi(target, sys%gr, psi, chi, sys%geo)
    call bwd_step(sys, td, hm, target, par, par_chi, chi, prop_chi, prop_psi)

    call states_end(psi)
    call states_copy(psi, initial_st)
    call fwd_step(sys, td, hm, target, par, par_chi, psi, prop_chi, prop_psi)

    j1 = target_j1(target, sys%gr, psi)

    call states_end(chi)
    call controlfunction_end(par_chi)
    POP_SUB(f_iter)
  end subroutine f_iter
  ! ---------------------------------------------------------

#include "opt_control_c.F90"
#include "check_input_inc.F90"
#include "finalcheck_inc.F90"


end module opt_control_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

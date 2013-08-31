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

#include "global.h"

!> This module contains the main procedure ("opt_control_run") that is 
!! used when optimal control runs are requested.
module opt_control_m
  use messages_m
  use controlfunction_m
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
  use minimizer_m
  use opt_control_global_m
  use opt_control_state_m
  use propagation_m
  use opt_control_iter_m
  use target_m
  use initst_m
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
            opt_control_cg_write_info,        &
            opt_control_direct_calc,          &
            opt_control_direct_message_info,  &
            opt_control_function_forward


  !> Module variables
  type(filter_t), save       :: filter
  type(oct_t), save          :: oct
  type(oct_iterator_t), save :: iterator
  type(target_t), save       :: oct_target
  type(opt_control_state_t), save :: initial_st
  

  !> For the direct, newuoa, and cg schemes:
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
    type(states_t)                 :: psi

    PUSH_SUB(opt_control_run)

    ! Creates a directory where the optimal control stuff will be written. The name of the directory
    ! is stored in the preprocessor macro OCT_DIR, which should be defined in src/include/global.h
    call io_mkdir(OCT_DIR)

    ! Initializes the time propagator. Then, it forces the propagation to be self consistent, in case
    ! the theory level is not "independent particles".
    call td_init(td, sys, hm)
    if(hm%theory_level /= INDEPENDENT_PARTICLES ) call propagator_set_scf_prop(td%tr, threshold = CNST(1.0e-14))

    ! Read general information about how the OCT run will be made, from inp file. "oct_read_inp" is
    ! in the opt_control_global_m module (like the definition of the oct_t data type)
    call oct_read_inp(oct)

    ! Read info about, and prepare, the control functions
    call controlfunction_mod_init(hm%ep, td%dt, td%max_iter, oct%mode_fixed_fluence)
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


    ! If filters are to be used, they also have to be initialized.
    call filter_init(td%max_iter, td%dt, filter)
    call filter_write(filter)


    ! Figure out the starting wavefunction(s), and the target.
    call initial_state_init(sys, hm, initial_st)
    call target_init(sys%gr, sys%geo, initial_st, td, controlfunction_w0(par), oct_target, oct, hm%ep)

    ! Sanity checks.
    call check_faulty_runmodes(sys, hm, td%tr)


    ! Informative output.
    call opt_control_get_qs(psi, initial_st)
    call output_states(psi, sys%gr, sys%geo, OCT_DIR//'initial', sys%outp)
    call target_output(oct_target, sys%gr, OCT_DIR//'target', sys%geo, sys%outp)
    call states_end(psi)


    ! mode switcher; here is where the real run is made.
    select case(oct%algorithm)
      case(oct_algorithm_zbr98)
        message(1) = "Info: Starting OCT iteration using scheme: ZBR98"
        call messages_info(1)
        call scheme_zbr98()
      case(oct_algorithm_wg05)
        message(1) = "Info: Starting OCT iteration using scheme: WG05"
        call messages_info(1)
        call scheme_wg05()
      case(oct_algorithm_zr98)
        message(1) = "Info: Starting OCT iteration using scheme: ZR98"
        call messages_info(1)
        call scheme_mt03()
      case(oct_algorithm_mt03)
        message(1) = "Info: Starting OCT iteration using scheme: MT03"
        call messages_info(1)
        call scheme_mt03()
      case(oct_algorithm_krotov)
        message(1) = "Info: Starting OCT iteration using scheme: KROTOV"
        call messages_info(1)
        call scheme_mt03()
      case(oct_algorithm_str_iter)
        message(1) = "Info: Starting OCT iterations using scheme: STRAIGHT ITERATION"
        call messages_info(1)
        call scheme_straight_iteration()
      case(oct_algorithm_cg)
        message(1) = "Info: Starting OCT iterations using scheme: CONJUGATE GRADIENTS"
        call messages_info(1)
        call scheme_cg()
      case(oct_algorithm_direct)
        message(1) = "Info: Starting OCT iterations using scheme: DIRECT OPTIMIZATION (NELDER-MEAD)"
        call messages_info(1)
        call scheme_direct()
      case(oct_algorithm_newuoa)
        message(1) = "Info: Starting OCT iterations using scheme: DIRECT OPTIMIZATION (NEWUOA)"
        call messages_info(1)
        call scheme_newuoa()
    case default
      call input_error('OCTScheme')
    end select

    ! do final test run: propagate initial state with optimal field
    call oct_finalcheck(sys, hm, td)

    ! clean up
    call controlfunction_end(par)
    call oct_iterator_end(iterator)
    call filter_end(filter)
    call td_end(td)
    call opt_control_state_end(initial_st)
    call target_end(oct_target, oct)
    call controlfunction_mod_close()
   
    POP_SUB(opt_control_run)

  contains


    ! ---------------------------------------------------------
    subroutine scheme_straight_iteration()
      PUSH_SUB(opt_control_run.scheme_straight_iteration)

      call controlfunction_set_rep(par)
      call controlfunction_copy(par_new, par)
      ctr_loop: do
        call controlfunction_copy(par_prev, par)
        call f_striter(sys, hm, td, par, j1)
        stop_loop = iteration_manager(j1, par_prev, par, iterator)
        if(clean_stop(sys%mc%master_comm) .or. stop_loop) exit ctr_loop
      end do ctr_loop

      call controlfunction_end(par_new)
      call controlfunction_end(par_prev)
      POP_SUB(opt_control_run.scheme_straight_iteration)
    end subroutine scheme_straight_iteration
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_mt03()
      type(opt_control_state_t) :: psi
      PUSH_SUB(opt_control_run.scheme_mt03)

      call opt_control_state_copy(psi, initial_st)
      call oct_prop_init(prop_chi, "chi")
      call oct_prop_init(prop_psi, "psi")

      call controlfunction_copy(par_new, par)
      ctr_loop: do
        call controlfunction_copy(par_prev, par)
        call f_iter(sys, hm, td, psi, par, prop_psi, prop_chi, j1)
        stop_loop = iteration_manager(j1, par, par_prev, iterator)
        if(clean_stop(sys%mc%master_comm) .or. stop_loop) exit ctr_loop
      end do ctr_loop

      call opt_control_state_end(psi)
      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call controlfunction_end(par_new)
      call controlfunction_end(par_prev)
      POP_SUB(opt_control_run.scheme_mt03)
    end subroutine scheme_mt03
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_wg05()
      type(opt_control_state_t) :: psi
      PUSH_SUB(opt_control_run.scheme_wg05)

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
        if(clean_stop(sys%mc%master_comm) .or. stop_loop) exit ctr_loop
      end do ctr_loop

      call opt_control_state_end(psi)
      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call controlfunction_end(par_new)
      call controlfunction_end(par_prev)
      POP_SUB(opt_control_run.scheme_wg05)
    end subroutine scheme_wg05
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_zbr98()
      type(opt_control_state_t) :: qcpsi
      PUSH_SUB(opt_control_run.scheme_zbr98)

      call opt_control_state_copy(qcpsi, initial_st)
      call oct_prop_init(prop_chi, "chi")
      call oct_prop_init(prop_psi, "psi")

      call controlfunction_copy(par_prev, par)
      call propagate_forward(sys, hm, td, par, oct_target, qcpsi, prop_psi)
      j1 = target_j1(oct_target, sys%gr, qcpsi)
      stop_loop = iteration_manager(j1, par, par_prev, iterator)
      if(clean_stop(sys%mc%master_comm) .or. stop_loop) then
        call opt_control_state_end(qcpsi)
        call oct_prop_end(prop_chi)
        call oct_prop_end(prop_psi)
        POP_SUB(opt_control_run.scheme_zbr98)
        return        
      end if

      call controlfunction_copy(par_new, par)
      ctr_loop: do
        call controlfunction_copy(par_prev, par)
        call f_zbr98(sys, hm, td, qcpsi, prop_psi, prop_chi, par)
        j1 = target_j1(oct_target, sys%gr, qcpsi)
        stop_loop = iteration_manager(j1, par, par_prev, iterator)
        if(clean_stop(sys%mc%master_comm) .or. stop_loop) exit ctr_loop
      end do ctr_loop

      call opt_control_state_end(qcpsi)
      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call controlfunction_end(par_new)
      call controlfunction_end(par_prev)
      POP_SUB(opt_control_run.scheme_zbr98)
    end subroutine scheme_zbr98
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_cg()
      integer :: dof, ierr, maxiter
      REAL_DOUBLE   :: step, minvalue
      FLOAT, allocatable :: theta(:)
      REAL_DOUBLE, allocatable :: x(:)
      FLOAT   :: f
      type(opt_control_state_t) :: qcpsi
      PUSH_SUB(opt_control_run.scheme_cg)

      call controlfunction_set_rep(par)

      call opt_control_state_copy(qcpsi, initial_st)
      call propagate_forward(sys, hm, td, par, oct_target, qcpsi)
      f = - target_j1(oct_target, sys%gr, qcpsi, sys%geo) - controlfunction_j2(par)
      call opt_control_state_end(qcpsi)
      call iteration_manager_direct(-f, par, iterator, sys)
      if(oct_iterator_maxiter(iterator) == 0) then
        ! Nothing to do.
        POP_SUB(opt_control_run.scheme_cg)
        return
      end if

      if(oct%random_initial_guess) call controlfunction_randomize(par)

      ! Set the module pointers, so that the opt_control_cg_calc and opt_control_cg_write_info routines
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

      call minimize_multidim(MINMETHOD_BFGS2, dof, x, step, real(0.1, 8), &
        real(oct_iterator_tolerance(iterator), 8), real(oct_iterator_tolerance(iterator), 8), &
        maxiter, opt_control_cg_calc, opt_control_cg_write_info, minvalue, ierr)

      if(ierr /= 0) then
        if(ierr <= 1024) then
          message(1) = "Error occurred during the GSL minimization procedure:"
          call loct_strerror(ierr, message(2))
          call messages_fatal(2)
        else
          message(1) = "The CG optimization did not meet the convergence criterion."
          call messages_info(1)
        end if
      end if

      call controlfunction_end(par_)
      SAFE_DEALLOCATE_A(x)
      SAFE_DEALLOCATE_A(theta)
      POP_SUB(opt_control_run.scheme_cg)
    end subroutine scheme_cg
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_direct()
      integer :: ierr, maxiter
      REAL_DOUBLE :: minvalue, step
      FLOAT, allocatable :: theta(:)
      REAL_DOUBLE, allocatable :: x(:)
      FLOAT :: f
      integer :: dim
      type(opt_control_state_t) :: qcpsi

      PUSH_SUB(opt_control_run.scheme_direct)

      call controlfunction_set_rep(par)
      dim = controlfunction_dof(par)

      call opt_control_state_copy(qcpsi, initial_st)
      call propagate_forward(sys, hm, td, par, oct_target, qcpsi)
      f = - target_j1(oct_target, sys%gr, qcpsi, sys%geo) - controlfunction_j2(par)
      call opt_control_state_end(qcpsi)
      call iteration_manager_direct(-f, par, iterator, sys)
      if(oct_iterator_maxiter(iterator) == 0) then
        ! Nothing to do.
        POP_SUB(opt_control_run.scheme_direct)
        return
      end if

      SAFE_ALLOCATE(x(1:dim))
      SAFE_ALLOCATE(theta(1:dim))

      if(oct%random_initial_guess) call controlfunction_randomize(par)

      ! Set the module pointers, so that the opt_control_direct_calc and opt_control_direct_message_info routines
      ! can use them.
      call controlfunction_copy(par_, par)
      sys_      => sys
      hm_       => hm
      td_       => td

      ! theta may be in single precision, whereas x is always double precision.
      call controlfunction_get_theta(par, theta)
      x = theta

      step = oct%direct_step * M_PI
      maxiter = oct_iterator_maxiter(iterator)

      call minimize_multidim_nograd(MINMETHOD_NMSIMPLEX, dim, x, step, &
        real(oct_iterator_tolerance(iterator), 8), maxiter, &
        opt_control_direct_calc, opt_control_direct_message_info, minvalue, ierr)

      if(ierr /= 0) then
        if(ierr <= 1024) then
          message(1) = "Error occurred during the GSL minimization procedure:"
          call loct_strerror(ierr, message(2))
          call messages_fatal(2)
        else
          message(1) = "The OCT direct optimization did not meet the convergence criterion."
          call messages_info(1)
        end if
      end if

      call controlfunction_end(par_)
      SAFE_DEALLOCATE_A(x)
      SAFE_DEALLOCATE_A(theta)
      POP_SUB(opt_control_run.scheme_direct)
    end subroutine scheme_direct
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_newuoa()
#if defined(HAVE_NEWUOA)
      integer :: dim, ierr, maxfun
      REAL_DOUBLE :: minvalue, step
      FLOAT, allocatable :: xl(:), xu(:)
      REAL_DOUBLE, allocatable :: x(:), w(:)
      FLOAT, allocatable :: theta(:)
      FLOAT :: f
      type(opt_control_state_t) :: qcpsi
      PUSH_SUB(opt_control_run.scheme_newuoa)

      call controlfunction_set_rep(par)

      call opt_control_state_copy(qcpsi, initial_st)
      call propagate_forward(sys, hm, td, par, oct_target,qcpsi)
      f = - target_j1(oct_target, sys%gr, qcpsi, sys%geo) - controlfunction_j2(par)
      call opt_control_state_end(qcpsi)
      call iteration_manager_direct(-f, par, iterator, sys)      
      if(oct_iterator_maxiter(iterator) == 0) then
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
      hm_       => hm
      td_       => td

      call controlfunction_get_theta(par, theta)
      x = theta

      maxfun = oct_iterator_maxiter(iterator)
      step = oct%direct_step * M_PI
      call minimize_multidim_nograd(MINMETHOD_NEWUOA, dim, x, step, &
        real(oct_iterator_tolerance(iterator), 8), maxfun, &
        opt_control_direct_calc, opt_control_direct_message_info, minvalue, ierr)

      call controlfunction_end(par_)
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
  subroutine f_zbr98(sys, hm, td, qcpsi, prop_psi, prop_chi, par)
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: hm
    type(td_t), intent(inout)                     :: td
    type(opt_control_state_t), intent(inout)      :: qcpsi
    type(oct_prop_t), intent(inout)               :: prop_psi, prop_chi
    type(controlfunction_t), intent(inout)        :: par

    type(states_t) :: chi
    type(opt_control_state_t) :: qcchi
    type(controlfunction_t) :: par_chi

    PUSH_SUB(f_zbr98)

    call controlfunction_copy(par_chi, par)

    call target_get_state(oct_target, chi)
    call opt_control_state_init(qcchi, chi, sys%geo)
    call bwd_step(sys, td, hm, oct_target, par, par_chi, qcchi, prop_chi, prop_psi)
    call opt_control_state_copy(qcpsi, initial_st)
    call fwd_step(sys, td, hm, oct_target, par, par_chi, qcpsi, prop_chi, prop_psi)

    call states_end(chi)
    call opt_control_state_end(qcchi)
    call controlfunction_end(par_chi)
    POP_SUB(f_zbr98)
  end subroutine f_zbr98


  ! ---------------------------------------------------------
  subroutine f_wg05(sys, hm, td, qcpsi, par, prop_psi, prop_chi, j1)
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: hm
    type(td_t), intent(inout)                     :: td
    type(opt_control_state_t), intent(inout)      :: qcpsi
    type(controlfunction_t), intent(inout)        :: par
    type(oct_prop_t), intent(inout)               :: prop_psi, prop_chi
    FLOAT, intent(out)                            :: j1

    FLOAT :: new_penalty
    type(opt_control_state_t) :: qcchi
    type(controlfunction_t) :: parp

    PUSH_SUB(f_wg05)

    if( oct_iterator_current(iterator)  ==  0) then
      call opt_control_state_copy(qcpsi, initial_st)
      call propagate_forward(sys, hm, td, par, oct_target, qcpsi, prop_psi)
      j1 = target_j1(oct_target, sys%gr, qcpsi)
      POP_SUB(f_wg05)
      return
    end if

    call controlfunction_copy(parp, par)

    call opt_control_state_copy(qcchi, qcpsi)
    call target_chi(oct_target, sys%gr, qcpsi, qcchi, sys%geo)
    call bwd_step(sys, td, hm, oct_target, par, parp, qcchi, prop_chi, prop_psi)

    call controlfunction_filter(parp, filter)

    ! recalc field
    if (oct%mode_fixed_fluence) then
      new_penalty = controlfunction_alpha(par, 1) * sqrt( controlfunction_fluence(parp) / controlfunction_targetfluence() )
      call controlfunction_set_alpha(parp, new_penalty)
      call controlfunction_set_fluence(parp)
    end if

    call controlfunction_copy(par, parp)
    call controlfunction_apply_envelope(par)

    call opt_control_state_copy(qcpsi, initial_st)
    call propagate_forward(sys, hm, td, par, oct_target, qcpsi, prop_psi)

    j1 = target_j1(oct_target, sys%gr, qcpsi)

    call opt_control_state_end(qcchi)
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

    type(opt_control_state_t) :: qcpsi
    type(opt_control_state_t) :: qcchi
    type(controlfunction_t) :: par_chi
    type(oct_prop_t)        :: prop_chi, prop_psi;

    PUSH_SUB(f_striter)

    call oct_prop_init(prop_chi, "chi")
    call oct_prop_init(prop_psi, "psi")

    call controlfunction_to_realtime(par)

    call controlfunction_copy(par_chi, par)

    ! First, a forward propagation with the input field.
    call opt_control_state_copy(qcpsi, initial_st)
    call propagate_forward(sys, hm, td, par, oct_target, qcpsi, prop_psi)

    ! Check the performance.
    j1 = target_j1(oct_target, sys%gr, qcpsi, sys%geo)

    ! Set the boundary condition for the backward propagation.
    call opt_control_state_copy(qcchi, qcpsi)
    call target_chi(oct_target, sys%gr, qcpsi, qcchi, sys%geo)

    ! Backward propagation, while at the same time finding the output field, 
    ! which is placed at par_chi
    call bwd_step_2(sys, td, hm, oct_target, par, par_chi, qcchi, prop_chi, prop_psi)
    !if(oct%mode_fixed_fluence) call controlfunction_set_fluence(par_chi)

    ! Copy par_chi to par
    call controlfunction_copy(par, par_chi)

    call opt_control_state_end(qcpsi)
    call opt_control_state_end(qcchi)
    call controlfunction_end(par_chi)
    call oct_prop_end(prop_chi)
    call oct_prop_end(prop_psi)

    POP_SUB(f_striter)
  end subroutine f_striter
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine f_iter(sys, hm, td, qcpsi, par, prop_psi, prop_chi, j1)
    type(system_t), intent(inout)                 :: sys
    type(hamiltonian_t), intent(inout)            :: hm
    type(td_t), intent(inout)                     :: td
    type(opt_control_state_t), intent(inout)      :: qcpsi
    type(controlfunction_t), intent(inout)        :: par
    type(oct_prop_t), intent(inout)               :: prop_psi, prop_chi
    FLOAT, intent(out)                            :: j1

    type(opt_control_state_t) :: qcchi
    type(controlfunction_t) :: par_chi

    PUSH_SUB(f_iter)

    if( oct_iterator_current(iterator)  ==  0) then
      call opt_control_state_copy(qcpsi, initial_st)
      call propagate_forward(sys, hm, td, par, oct_target, qcpsi, prop_psi)
      j1 = target_j1(oct_target, sys%gr, qcpsi)
      POP_SUB(f_iter)
      return
    end if

    call controlfunction_copy(par_chi, par)

    call opt_control_state_copy(qcchi, qcpsi)
    call target_chi(oct_target, sys%gr, qcpsi, qcchi, sys%geo)
    call bwd_step(sys, td, hm, oct_target, par, par_chi, qcchi, prop_chi, prop_psi)

    call opt_control_state_copy(qcpsi, initial_st)
    call fwd_step(sys, td, hm, oct_target, par, par_chi, qcpsi, prop_chi, prop_psi)

    j1 = target_j1(oct_target, sys%gr, qcpsi)

    call opt_control_state_end(qcchi)
    call controlfunction_end(par_chi)
    POP_SUB(f_iter)
  end subroutine f_iter
  ! ---------------------------------------------------------

#include "opt_control_c_inc.F90"
#include "check_input_inc.F90"
#include "finalcheck_inc.F90"


end module opt_control_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

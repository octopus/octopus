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

!> This module contains the main procedure ("opt_control_run") that is 
!! used when optimal control runs are requested.
module opt_control_oct_m
  use boundary_op_oct_m
  use controlfunction_oct_m
  use exponential_oct_m
  use filter_oct_m
  use global_oct_m
  use grid_oct_m
  use initst_oct_m
  use ions_oct_m
  use iso_c_binding
  use output_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use lasers_oct_m
  use loct_oct_m
  use math_oct_m
  use messages_oct_m
  use minimizer_oct_m
  use multisystem_basic_oct_m
  use namespace_oct_m
  use opt_control_global_oct_m
  use opt_control_iter_oct_m
  use opt_control_state_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagation_oct_m
  use propagator_elec_oct_m
  use propagator_base_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use electrons_oct_m
  use target_oct_m
  use td_oct_m

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
  

  !> For the direct, nlopt, and cg schemes:
  type(controlfunction_t), save :: par_
  type(electrons_t), pointer :: sys_
  type(hamiltonian_elec_t), pointer :: hm_
  type(td_t), pointer :: td_
  FLOAT, allocatable :: x_(:)
  integer :: index_

contains

  ! ---------------------------------------------------------
  subroutine opt_control_run(system)
    class(*), intent(inout) :: system

    PUSH_SUB(opt_control_run)

    select type (system)
    class is (multisystem_basic_t)
      message(1) = "CalculationMode = opt_control not implemented for multi-system calculations"
      call messages_fatal(1)
    type is (electrons_t)
      call opt_control_run_legacy(system)
    end select

    POP_SUB(opt_control_run)
  end subroutine opt_control_run

  !> This is the main procedure for all types of optimal control runs.
  !! It is called from the "run" procedure in the "run_m" module.
  subroutine opt_control_run_legacy(sys)
    type(electrons_t), target,      intent(inout) :: sys

    type(td_t), target             :: td
    type(controlfunction_t)        :: par, par_new, par_prev
    logical                        :: stop_loop
    FLOAT                          :: j1
    type(oct_prop_t)               :: prop_chi, prop_psi
    type(states_elec_t)            :: psi

    PUSH_SUB(opt_control_run_legacy)

    if (sys%hm%pcm%run_pcm) then
      call messages_not_implemented("PCM for CalculationMode /= gs or td")
    end if

    if (sys%kpoints%use_symmetries) then
      call messages_experimental("KPoints symmetries with CalculationMode = opt_control")
    end if

    ! Creates a directory where the optimal control stuff will be written. The name of the directory
    ! is stored in the preprocessor macro OCT_DIR, which should be defined in src/include/global.h
    call io_mkdir(OCT_DIR, sys%namespace)

    ! Initializes the time propagator. Then, it forces the propagation to be self consistent, in case
    ! the theory level is not "independent particles".
    call td_init(td, sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm, sys%outp)
    if(sys%hm%theory_level /= INDEPENDENT_PARTICLES ) call propagator_elec_set_scf_prop(td%tr, threshold = CNST(1.0e-14))

    ! Read general information about how the OCT run will be made, from inp file. "oct_read_inp" is
    ! in the opt_control_global_oct_m module (like the definition of the oct_t data type)
    call oct_read_inp(oct, sys%namespace)

    ! Read info about, and prepare, the control functions
    call controlfunction_mod_init(sys%hm%ext_lasers, sys%namespace, td%dt, td%max_iter, oct%mode_fixed_fluence)
    call controlfunction_init(par, td%dt, td%max_iter)
    call controlfunction_set(par, sys%hm%ext_lasers)
      ! This prints the initial control parameters, exactly as described in the inp file,
      ! that is, without applying any envelope or filter.
    call controlfunction_write(OCT_DIR//'initial_laser_inp', par, sys%namespace)
    call controlfunction_prepare_initial(par)
    call controlfunction_to_h(par, sys%hm%ext_lasers)
    call messages_print_stress(stdout, "TD ext. fields after processing")
    call laser_write_info(sys%hm%ext_lasers%lasers, stdout)
    call messages_print_stress(stdout)
    call controlfunction_write(OCT_DIR//'initial_laser', par, sys%namespace)


    ! Startup of the iterator data type (takes care of counting iterations, stopping, etc).
    call oct_iterator_init(iterator, sys%namespace, par)


    ! Initialization of the propagation_oct_m module.
    call propagation_mod_init(td%max_iter, oct%eta, oct%delta, oct%number_checkpoints, &
      (oct%algorithm == OPTION__OCTSCHEME__OCT_ZBR98), &
      (oct%algorithm == OPTION__OCTSCHEME__OCT_CG) .or. &
      (oct%algorithm == OPTION__OCTSCHEME__OCT_BFGS) .or. &
      (oct%algorithm == OPTION__OCTSCHEME__OCT_NLOPT_LBFGS) )


    ! If filters are to be used, they also have to be initialized.
    call filter_init(td%max_iter, sys%namespace, td%dt, filter)
    call filter_write(filter, sys%namespace)


    ! Figure out the starting wavefunction(s), and the target.
    call initial_state_init(sys, initial_st)
    call target_init(sys%gr, sys%kpoints, sys%namespace, sys%space, sys%ions, initial_st, td, &
               controlfunction_w0(par), oct_target, oct, sys%hm%ep, sys%mc)

    ! Sanity checks.
    call check_faulty_runmodes(sys, td%tr)


    ! Informative output.
    call opt_control_get_qs(psi, initial_st)
    call output_states(sys%outp, sys%namespace, sys%space, OCT_DIR//'initial', psi, sys%gr, sys%ions, sys%hm, -1)
    call target_output(oct_target, sys%namespace, sys%space, sys%gr, OCT_DIR//'target', sys%ions, sys%hm, sys%outp)
    call states_elec_end(psi)


    ! mode switcher; here is where the real run is made.
    select case(oct%algorithm)
      case(OPTION__OCTSCHEME__OCT_ZBR98)
        message(1) = "Info: Starting OCT iteration using scheme: ZBR98"
        call messages_info(1)
        call scheme_zbr98()
      case(OPTION__OCTSCHEME__OCT_WG05)
        message(1) = "Info: Starting OCT iteration using scheme: WG05"
        call messages_info(1)
        call scheme_wg05()
      case(OPTION__OCTSCHEME__OCT_ZR98)
        message(1) = "Info: Starting OCT iteration using scheme: ZR98"
        call messages_info(1)
        call scheme_mt03()
      case(OPTION__OCTSCHEME__OCT_MT03)
        message(1) = "Info: Starting OCT iteration using scheme: MT03"
        call messages_info(1)
        call scheme_mt03()
      case(OPTION__OCTSCHEME__OCT_KROTOV)
        message(1) = "Info: Starting OCT iteration using scheme: KROTOV"
        call messages_info(1)
        call scheme_mt03()
      case(OPTION__OCTSCHEME__OCT_STRAIGHT_ITERATION)
        message(1) = "Info: Starting OCT iterations using scheme: STRAIGHT ITERATION"
        call messages_info(1)
        call scheme_straight_iteration()
      case(OPTION__OCTSCHEME__OCT_CG)
        message(1) = "Info: Starting OCT iterations using scheme: CONJUGATE GRADIENTS"
        call messages_info(1)
        call scheme_cg()
      case(OPTION__OCTSCHEME__OCT_BFGS)
        message(1) = "Info: Starting OCT iterations using scheme: BFGS"
        call messages_info(1)
        call scheme_cg()
      case(OPTION__OCTSCHEME__OCT_DIRECT)
        message(1) = "Info: Starting OCT iterations using scheme: DIRECT OPTIMIZATION (NELDER-MEAD)"
        call messages_info(1)
        call scheme_direct()
      case(OPTION__OCTSCHEME__OCT_NLOPT_BOBYQA)
        message(1) = "Info: Starting OCT iterations using scheme: DIRECT OPTIMIZATION (NLOPT - BOBYQA)"
        call messages_info(1)
        call scheme_nlopt()
      case(OPTION__OCTSCHEME__OCT_NLOPT_LBFGS)
        message(1) = "Info: Starting OCT iterations using scheme: DIRECT OPTIMIZATION (NLOPT - LBFGS)"
        call messages_info(1)
        call scheme_nlopt()
    case default
      call messages_input_error(sys%namespace, 'OCTScheme')
    end select

    ! do final test run: propagate initial state with optimal field
    call oct_finalcheck(sys, td)

    ! clean up
    call controlfunction_end(par)
    call oct_iterator_end(iterator, sys%namespace)
    call filter_end(filter)
    call td_end(td)
    call opt_control_state_end(initial_st)
    call target_end(oct_target, oct)
    call controlfunction_mod_close()
   
    POP_SUB(opt_control_run_legacy)

  contains


    ! ---------------------------------------------------------
    subroutine scheme_straight_iteration()
      PUSH_SUB(opt_control_run_legacy.scheme_straight_iteration)

      call controlfunction_set_rep(par)
      call controlfunction_copy(par_new, par)
      ctr_loop: do
        call controlfunction_copy(par_prev, par)
        call f_striter(sys, td, par, j1)
        stop_loop = iteration_manager(sys%namespace, j1, par_prev, par, iterator)
        if(clean_stop(sys%mc%master_comm) .or. stop_loop) exit ctr_loop
      end do ctr_loop

      call controlfunction_end(par_new)
      call controlfunction_end(par_prev)
      POP_SUB(opt_control_run_legacy.scheme_straight_iteration)
    end subroutine scheme_straight_iteration
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_mt03()
      type(opt_control_state_t) :: psi
      PUSH_SUB(opt_control_run_legacy.scheme_mt03)

      call opt_control_state_null(psi)
      call opt_control_state_copy(psi, initial_st)
      call oct_prop_init(prop_chi, sys%namespace, "chi", sys%gr%mesh, sys%mc)
      call oct_prop_init(prop_psi, sys%namespace, "psi", sys%gr%mesh, sys%mc)

      call controlfunction_copy(par_new, par)
      ctr_loop: do
        call controlfunction_copy(par_prev, par)
        call f_iter(sys, td, psi, par, prop_psi, prop_chi, j1)
        stop_loop = iteration_manager(sys%namespace, j1, par, par_prev, iterator)
        if(clean_stop(sys%mc%master_comm) .or. stop_loop) exit ctr_loop
      end do ctr_loop

      call opt_control_state_end(psi)
      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call controlfunction_end(par_new)
      call controlfunction_end(par_prev)
      POP_SUB(opt_control_run_legacy.scheme_mt03)
    end subroutine scheme_mt03
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_wg05()
      type(opt_control_state_t) :: psi
      PUSH_SUB(opt_control_run_legacy.scheme_wg05)

      call oct_prop_init(prop_chi, sys%namespace, "chi", sys%gr%mesh, sys%mc)
      call oct_prop_init(prop_psi, sys%namespace, "psi", sys%gr%mesh, sys%mc)

      if (oct%mode_fixed_fluence) then
        call controlfunction_set_alpha(par, sqrt( controlfunction_fluence(par) / controlfunction_targetfluence()))
      end if

      call opt_control_state_null(psi)

      call controlfunction_copy(par_new, par)
      ctr_loop: do
        call controlfunction_copy(par_prev, par)
        call f_wg05(sys, td, psi, par, prop_psi, prop_chi, j1)
        stop_loop = iteration_manager(sys%namespace, j1, par, par_prev, iterator)
        if(clean_stop(sys%mc%master_comm) .or. stop_loop) exit ctr_loop
      end do ctr_loop

      call opt_control_state_end(psi)
      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call controlfunction_end(par_new)
      call controlfunction_end(par_prev)
      POP_SUB(opt_control_run_legacy.scheme_wg05)
    end subroutine scheme_wg05
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_zbr98()
      type(opt_control_state_t) :: qcpsi
      PUSH_SUB(opt_control_run_legacy.scheme_zbr98)

      call opt_control_state_null(qcpsi)
      call opt_control_state_copy(qcpsi, initial_st)
      call oct_prop_init(prop_chi, sys%namespace, "chi", sys%gr%mesh, sys%mc)
      call oct_prop_init(prop_psi, sys%namespace, "psi", sys%gr%mesh, sys%mc)

      call controlfunction_copy(par_prev, par)
      call propagate_forward(sys, td, par, oct_target, qcpsi, prop_psi)
      j1 = target_j1(oct_target, sys%namespace, sys%gr, sys%kpoints, qcpsi)
      stop_loop = iteration_manager(sys%namespace, j1, par, par_prev, iterator)
      if(clean_stop(sys%mc%master_comm) .or. stop_loop) then
        call opt_control_state_end(qcpsi)
        call oct_prop_end(prop_chi)
        call oct_prop_end(prop_psi)
        POP_SUB(opt_control_run_legacy.scheme_zbr98)
        return        
      end if

      call controlfunction_copy(par_new, par)
      ctr_loop: do
        call controlfunction_copy(par_prev, par)
        call f_zbr98(sys, td, qcpsi, prop_psi, prop_chi, par)
        j1 = target_j1(oct_target, sys%namespace, sys%gr, sys%kpoints, qcpsi)
        stop_loop = iteration_manager(sys%namespace, j1, par, par_prev, iterator)
        if(clean_stop(sys%mc%master_comm) .or. stop_loop) exit ctr_loop
      end do ctr_loop

      call opt_control_state_end(qcpsi)
      call oct_prop_end(prop_chi)
      call oct_prop_end(prop_psi)
      call controlfunction_end(par_new)
      call controlfunction_end(par_prev)
      POP_SUB(opt_control_run_legacy.scheme_zbr98)
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
      PUSH_SUB(opt_control_run_legacy.scheme_cg)

      call controlfunction_set_rep(par)

      call opt_control_state_null(qcpsi)
      call opt_control_state_copy(qcpsi, initial_st)
      call propagate_forward(sys, td, par, oct_target, qcpsi)
      f = - target_j1(oct_target, sys%namespace, sys%gr, sys%kpoints, qcpsi, sys%ions) - controlfunction_j2(par)
      call opt_control_state_end(qcpsi)
      call iteration_manager_direct(-f, par, iterator, sys)
      if(oct_iterator_maxiter(iterator) == 0) then
        ! Nothing to do.
        POP_SUB(opt_control_run_legacy.scheme_cg)
        return
      end if

      if(oct%random_initial_guess) call controlfunction_randomize(par)

      ! Set the module pointers, so that the opt_control_cg_calc and opt_control_cg_write_info routines
      ! can use them.
      call controlfunction_copy(par_, par)
      sys_      => sys
      hm_       => sys%hm
      td_       => td

      dof = controlfunction_dof(par)
      SAFE_ALLOCATE(x(1:dof))
      SAFE_ALLOCATE(theta(1:dof))
      call controlfunction_get_theta(par, theta)
      x = theta

      step = oct%direct_step * M_PI
      maxiter = oct_iterator_maxiter(iterator) - 1

      select case(oct%algorithm)
      case(OPTION__OCTSCHEME__OCT_BFGS)
        call minimize_multidim(MINMETHOD_BFGS2, dof, x, step, CNST(0.1), &
          TOFLOAT(oct_iterator_tolerance(iterator)), TOFLOAT(oct_iterator_tolerance(iterator)), &
          maxiter, opt_control_cg_calc, opt_control_cg_write_info, minvalue, ierr)
      case(OPTION__OCTSCHEME__OCT_CG)
        call minimize_multidim(MINMETHOD_FR_CG, dof, x, step, CNST(0.1), &
          TOFLOAT(oct_iterator_tolerance(iterator)), TOFLOAT(oct_iterator_tolerance(iterator)), &
          maxiter, opt_control_cg_calc, opt_control_cg_write_info, minvalue, ierr)
      end select

      if(ierr /= 0) then
        if(ierr <= 1024) then
          message(1) = "Error occurred during the GSL minimization procedure:"
          call loct_strerror(ierr, message(2))
          call messages_fatal(2)
        else
          message(1) = "The optimization did not meet the convergence criterion."
          call messages_info(1)
        end if
      end if

      call controlfunction_end(par_)
      SAFE_DEALLOCATE_A(x)
      SAFE_DEALLOCATE_A(theta)
      POP_SUB(opt_control_run_legacy.scheme_cg)
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

      PUSH_SUB(opt_control_run_legacy.scheme_direct)

      call controlfunction_set_rep(par)
      dim = controlfunction_dof(par)

      call opt_control_state_null(qcpsi)
      call opt_control_state_copy(qcpsi, initial_st)
      call propagate_forward(sys, td, par, oct_target, qcpsi)
      f = - target_j1(oct_target, sys%namespace, sys%gr, sys%kpoints, qcpsi, sys%ions) - controlfunction_j2(par)
      call opt_control_state_end(qcpsi)
      call iteration_manager_direct(-f, par, iterator, sys)
      if(oct_iterator_maxiter(iterator) == 0) then
        ! Nothing to do.
        POP_SUB(opt_control_run_legacy.scheme_direct)
        return
      end if

      SAFE_ALLOCATE(x(1:dim))
      SAFE_ALLOCATE(theta(1:dim))

      if(oct%random_initial_guess) call controlfunction_randomize(par)

      ! Set the module pointers, so that the opt_control_direct_calc and opt_control_direct_message_info routines
      ! can use them.
      call controlfunction_copy(par_, par)
      sys_      => sys
      hm_       => sys%hm
      td_       => td

      ! theta may be in single precision, whereas x is always double precision.
      call controlfunction_get_theta(par, theta)
      x = theta

      step = oct%direct_step * M_PI
      maxiter = oct_iterator_maxiter(iterator)

      call minimize_multidim_nograd(MINMETHOD_NMSIMPLEX, dim, x, step, &
        TOFLOAT(oct_iterator_tolerance(iterator)), maxiter, &
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
      POP_SUB(opt_control_run_legacy.scheme_direct)
    end subroutine scheme_direct
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine scheme_nlopt()
#if defined(HAVE_NLOPT)
      integer :: method, dim, maxiter, ierr
      FLOAT, allocatable :: x(:), xl(:), xu(:)
      FLOAT :: step, toldr, minimum, f
      type(opt_control_state_t) :: qcpsi
      PUSH_SUB(opt_control_run_legacy.scheme_nlopt)

      call controlfunction_set_rep(par)

      call opt_control_state_null(qcpsi)
      call opt_control_state_copy(qcpsi, initial_st)
      call propagate_forward(sys, td, par, oct_target, qcpsi)
      f = - target_j1(oct_target, sys%namespace, sys%gr, sys%kpoints, qcpsi, sys%ions) - controlfunction_j2(par)
      call opt_control_state_end(qcpsi)
      call iteration_manager_direct(-f, par, iterator, sys)      
      if(oct_iterator_maxiter(iterator) == 0) then
        ! Nothing to do.
        POP_SUB(opt_control_run_legacy.scheme_cg)
        return
      end if

      dim = controlfunction_dof(par)
      SAFE_ALLOCATE(x(dim))
      SAFE_ALLOCATE(xl(1:dim))
      SAFE_ALLOCATE(xu(1:dim))
      call controlfunction_bounds(par, xl, xu)

      if(oct%random_initial_guess) call controlfunction_randomize(par)

      ! Set the module pointers, so that the calc_point and write_iter_info routines
      ! can use them.
      call controlfunction_copy(par_, par)
      sys_      => sys
      hm_       => sys%hm
      td_       => td

      call controlfunction_get_theta(par, x)

      maxiter = oct_iterator_maxiter(iterator)
      step = oct%direct_step
      select case(oct%algorithm)
      case(OPTION__OCTSCHEME__OCT_NLOPT_BOBYQA)
        method = MINMETHOD_NLOPT_BOBYQA
      case(OPTION__OCTSCHEME__OCT_NLOPT_LBFGS)
        method = MINMETHOD_NLOPT_LBFGS
      end select
      toldr = oct_iterator_tolerance(iterator)

      call minimize_multidim_nlopt(ierr, method, dim, x, step, toldr, maxiter, opt_control_nlopt_func, minimum, &
        xl, xu)
      if(ierr < 1 .or. ierr > 4) then
         message(1) = "The nlopt minimization procedure did not find convergence, or found an error"
         write(message(2),'(a,i5)') "Error code =", ierr
         call messages_info(2)
      end if

      call controlfunction_end(par_)
      SAFE_DEALLOCATE_A(xl)
      SAFE_DEALLOCATE_A(xu)
      SAFE_DEALLOCATE_A(x)
      POP_SUB(opt_control_run_legacy.scheme_nlopt)
#endif
    end subroutine scheme_nlopt
    ! ---------------------------------------------------------

  end subroutine opt_control_run_legacy
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine f_zbr98(sys, td, qcpsi, prop_psi, prop_chi, par)
    type(electrons_t),         intent(inout) :: sys
    type(td_t),                intent(inout) :: td
    type(opt_control_state_t), intent(inout) :: qcpsi
    type(oct_prop_t),          intent(inout) :: prop_psi
    type(oct_prop_t),          intent(inout) :: prop_chi
    type(controlfunction_t),   intent(inout) :: par

    type(states_elec_t) :: chi
    type(opt_control_state_t) :: qcchi
    type(controlfunction_t) :: par_chi

    PUSH_SUB(f_zbr98)

    call controlfunction_copy(par_chi, par)

    call target_get_state(oct_target, chi)
    call opt_control_state_init(qcchi, chi, sys%ions)
    call bwd_step(sys, td, oct_target, par, par_chi, qcchi, prop_chi, prop_psi)
    call opt_control_state_copy(qcpsi, initial_st)
    call fwd_step(sys, td, oct_target, par, par_chi, qcpsi, prop_chi, prop_psi)

    call states_elec_end(chi)
    call opt_control_state_end(qcchi)
    call controlfunction_end(par_chi)
    POP_SUB(f_zbr98)
  end subroutine f_zbr98


  ! ---------------------------------------------------------
  subroutine f_wg05(sys, td, qcpsi, par, prop_psi, prop_chi, j1)
    type(electrons_t),         intent(inout) :: sys
    type(td_t),                intent(inout) :: td
    type(opt_control_state_t), intent(inout) :: qcpsi
    type(controlfunction_t),   intent(inout) :: par
    type(oct_prop_t),          intent(inout) :: prop_psi
    type(oct_prop_t),          intent(inout) :: prop_chi
    FLOAT,                     intent(out)   :: j1

    FLOAT :: new_penalty
    type(opt_control_state_t) :: qcchi
    type(controlfunction_t) :: parp

    PUSH_SUB(f_wg05)

    if( oct_iterator_current(iterator)  ==  0) then
      call opt_control_state_copy(qcpsi, initial_st)
      call propagate_forward(sys, td, par, oct_target, qcpsi, prop_psi)
      j1 = target_j1(oct_target, sys%namespace, sys%gr, sys%kpoints, qcpsi)
      POP_SUB(f_wg05)
      return
    end if

    call controlfunction_copy(parp, par)

    call opt_control_state_null(qcchi)
    call opt_control_state_copy(qcchi, qcpsi)
    call target_chi(oct_target, sys%namespace, sys%gr, sys%kpoints, qcpsi, qcchi, sys%ions)
    call bwd_step(sys, td, oct_target, par, parp, qcchi, prop_chi, prop_psi)

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
    call propagate_forward(sys, td, par, oct_target, qcpsi, prop_psi)

    j1 = target_j1(oct_target, sys%namespace, sys%gr, sys%kpoints, qcpsi)

    call opt_control_state_end(qcchi)
    call controlfunction_end(parp)
    POP_SUB(f_wg05)
  end subroutine f_wg05
  ! ---------------------------------------------------------


   ! ---------------------------------------------------------
  subroutine f_striter(sys, td, par, j1)
    type(electrons_t),       intent(inout) :: sys
    type(td_t),              intent(inout) :: td
    type(controlfunction_t), intent(inout) :: par
    FLOAT,                   intent(out)   :: j1

    type(opt_control_state_t) :: qcpsi
    type(opt_control_state_t) :: qcchi
    type(controlfunction_t) :: par_chi
    type(oct_prop_t)        :: prop_chi, prop_psi;

    PUSH_SUB(f_striter)

    call oct_prop_init(prop_chi, sys%namespace, "chi", sys%gr%mesh, sys%mc)
    call oct_prop_init(prop_psi, sys%namespace, "psi", sys%gr%mesh, sys%mc)

    call controlfunction_to_realtime(par)

    call controlfunction_copy(par_chi, par)

    ! First, a forward propagation with the input field.
    call opt_control_state_null(qcpsi)
    call opt_control_state_copy(qcpsi, initial_st)
    call propagate_forward(sys, td, par, oct_target, qcpsi, prop_psi)

    ! Check the performance.
    j1 = target_j1(oct_target, sys%namespace, sys%gr, sys%kpoints, qcpsi, sys%ions)

    ! Set the boundary condition for the backward propagation.
    call opt_control_state_null(qcchi)
    call opt_control_state_copy(qcchi, qcpsi)
    call target_chi(oct_target, sys%namespace, sys%gr, sys%kpoints, qcpsi, qcchi, sys%ions)

    ! Backward propagation, while at the same time finding the output field, 
    ! which is placed at par_chi
    call bwd_step_2(sys, td, oct_target, par, par_chi, qcchi, prop_chi, prop_psi)
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
  subroutine f_iter(sys, td, qcpsi, par, prop_psi, prop_chi, j1)
    type(electrons_t),         intent(inout) :: sys
    type(td_t),                intent(inout) :: td
    type(opt_control_state_t), intent(inout) :: qcpsi
    type(controlfunction_t),   intent(inout) :: par
    type(oct_prop_t),          intent(inout) :: prop_psi
    type(oct_prop_t),          intent(inout) :: prop_chi
    FLOAT,                     intent(out)   :: j1

    type(opt_control_state_t) :: qcchi
    type(controlfunction_t) :: par_chi

    PUSH_SUB(f_iter)

    if( oct_iterator_current(iterator)  ==  0) then
      call opt_control_state_copy(qcpsi, initial_st)
      call propagate_forward(sys, td, par, oct_target, qcpsi, prop_psi)
      j1 = target_j1(oct_target, sys%namespace, sys%gr, sys%kpoints, qcpsi)
      POP_SUB(f_iter)
      return
    end if

    call controlfunction_copy(par_chi, par)

    call opt_control_state_null(qcchi)
    call opt_control_state_copy(qcchi, qcpsi)
    call target_chi(oct_target, sys%namespace, sys%gr, sys%kpoints, qcpsi, qcchi, sys%ions)
    call bwd_step(sys, td, oct_target, par, par_chi, qcchi, prop_chi, prop_psi)

    call opt_control_state_copy(qcpsi, initial_st)
    call fwd_step(sys, td, oct_target, par, par_chi, qcpsi, prop_chi, prop_psi)

    j1 = target_j1(oct_target, sys%namespace, sys%gr, sys%kpoints, qcpsi)

    call opt_control_state_end(qcchi)
    call controlfunction_end(par_chi)
    POP_SUB(f_iter)
  end subroutine f_iter
  ! ---------------------------------------------------------

#include "opt_control_c_inc.F90"
#include "check_input_inc.F90"
#include "finalcheck_inc.F90"


end module opt_control_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

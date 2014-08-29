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

module propagator_m
  use batch_m
  use batch_ops_m
  use blas_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use cmplxscl_m
  use comm_m
  use cube_function_m
  use datasets_m
  use density_m
  use energy_calc_m
  use exponential_m
  use forces_m
  use gauge_field_m
  use grid_m
  use geometry_m
  use global_m
  use hamiltonian_m
  use ion_dynamics_m
  use lalg_basic_m
  use lasers_m
  use loct_pointer_m
  use parser_m
  use math_m
  use mesh_function_m
  use messages_m
  use multicomm_m
  use ob_mem_m
  use ob_propagator_m
  use ob_terms_m
  use opencl_m
  use opt_control_state_m
  use output_m
  use profiling_m
  use restart_m
  use scf_m
  use species_m
  use states_dim_m
  use solvers_m
  use sparskit_m
  use states_m
  use types_m
  use v_ks_m
  use varinfo_m
  use xc_m

  implicit none

  private
  public ::                         &
    propagator_t,                   &
    propagator_init,                &
    propagator_end,                 &
    propagator_copy,                &
    propagator_run_zero_iter,       &
    propagator_dt,                  &
    propagator_set_scf_prop,        &
    propagator_remove_scf_prop,     &
    propagator_ions_are_propagated, &
    propagator_dt_bo,               &
    vksinterp_t,                    &
    vksinterp_dump,                 &
    vksinterp_load


  integer, public, parameter ::        &
    PROP_ETRS                    = 2,  &
    PROP_AETRS                   = 3,  &
    PROP_EXPONENTIAL_MIDPOINT    = 4,  &
    PROP_CRANK_NICOLSON          = 5,  &
    PROP_CRANK_NICOLSON_SPARSKIT = 6,  &
    PROP_MAGNUS                  = 7,  &
    PROP_CRANK_NICOLSON_SRC_MEM  = 8,  &
    PROP_QOCT_TDDFT_PROPAGATOR   = 10, &
    PROP_CAETRS                  = 12, &
    PROP_RUNGE_KUTTA4            = 13, &
    PROP_RUNGE_KUTTA2            = 14, & 
    PROP_EXPLICIT_RUNGE_KUTTA4   = 15

  type vksinterp_t
    private
    FLOAT, pointer      :: v_old(:, :, :) => null()
    FLOAT, pointer      :: Imv_old(:, :, :) => null()
  end type vksinterp_t

  type propagator_t
    integer             :: method           !< Which evolution method to use.
    type(exponential_t) :: te               !< How to apply the propagator \f$ e^{-i H \Delta t} \f$.
    !> Storage of the KS potential of previous iterations.
    type(vksinterp_t) :: vksold
    !> Auxiliary function to store the Magnus potentials.
    FLOAT, pointer      :: vmagnus(:, :, :) => null() 
    type(ob_terms_t)    :: ob               !< For open boundaries: leads, memory
    integer             :: scf_propagation_steps 
    logical             :: first
#ifdef HAVE_SPARSKIT
    type(sparskit_solver_t), pointer :: tdsk
    integer             :: tdsk_size
#endif
    FLOAT              :: scf_threshold
  end type propagator_t

  type(grid_t),            pointer, private :: grid_p
  type(hamiltonian_t),     pointer, private :: hm_p
  type(propagator_t),      pointer, private :: tr_p
  type(states_t),          pointer, private :: st_p
  type(xc_t),              pointer, private :: xc_p
  integer,                 private :: ik_op, ist_op, idim_op, dim_op, nst_op
  type(states_t),          private :: st_op
  FLOAT,                   private :: t_op, dt_op
  FLOAT, allocatable, private      :: vhxc1_op(:, :), vhxc2_op(:, :), vpsl1_op(:), vpsl2_op(:)
  logical :: move_ions_op

contains

#include "propagator_vksold_inc.F90"

  ! ---------------------------------------------------------
  subroutine propagator_nullify(this)
    type(propagator_t), intent(out) :: this

    PUSH_SUB(propagator_nullify)

    !this%method
    !call exponential_nullify(this%te)
    call vksinterp_nullify(this%vksold)
    this%vmagnus=>null() 
    !call ob_terms_nullify(this%ob)
    !this%scf_propagation_steps 
    !this%first

#ifdef HAVE_SPARSKIT
    nullify(this%tdsk)
#endif

    POP_SUB(propagator_nullify)
  end subroutine propagator_nullify
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine propagator_copy(tro, tri)
    type(propagator_t), intent(inout) :: tro
    type(propagator_t), intent(in)    :: tri

    PUSH_SUB(propagator_copy)
    
    call propagator_nullify(tro)
    tro%method = tri%method

    select case(tro%method)
    case(PROP_MAGNUS)
      call loct_pointer_copy(tro%vmagnus, tri%vmagnus)
    case(PROP_CRANK_NICOLSON_SRC_MEM)
      message(1) = 'Internal error at propagator_copy.'
      call messages_fatal(1)
#ifdef HAVE_SPARSKIT
    case(PROP_CRANK_NICOLSON_SPARSKIT)
      SAFE_ALLOCATE(tro%tdsk)
      tro%tdsk_size = tri%tdsk_size
      call sparskit_solver_init(tro%tdsk_size, tro%tdsk, is_complex = .true.)
#endif
#ifdef HAVE_SPARSKIT
    case(PROP_RUNGE_KUTTA4)
      SAFE_ALLOCATE(tro%tdsk)
      tro%tdsk_size = tri%tdsk_size
      call sparskit_solver_init(tro%tdsk_size, tro%tdsk, is_complex = .true.)
#endif
#ifdef HAVE_SPARSKIT
    case(PROP_RUNGE_KUTTA2)
      SAFE_ALLOCATE(tro%tdsk)
      tro%tdsk_size = tri%tdsk_size
      call sparskit_solver_init(tro%tdsk_size, tro%tdsk, is_complex = .true.)
#endif
    end select

    call vksinterp_copy(tro%vksold, tri%vksold)

    call exponential_copy(tro%te, tri%te)
    tro%scf_propagation_steps = tri%scf_propagation_steps

    tro%scf_threshold = tri%scf_threshold
    POP_SUB(propagator_copy)
  end subroutine propagator_copy
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_init(gr, st, hm, tr, dt, max_iter, have_fields)
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(in)    :: st
    type(hamiltonian_t), intent(in)    :: hm
    type(propagator_t),  intent(inout) :: tr
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: max_iter
    !> whether there is an associated "field"
    !! that must be propagated (currently ions
    !! or a gauge field).
    logical,             intent(in)    :: have_fields 

    integer :: default_propagator
    logical :: cmplxscl
    
    PUSH_SUB(propagator_init)
    
    cmplxscl = st%cmplxscl%space
    
    call propagator_nullify(tr)

    !%Variable TDPropagator
    !%Type integer
    !%Default etrs
    !%Section Time-Dependent::Propagation
    !%Description
    !% This variable determines which algorithm will be used to approximate
    !% the evolution operator <math>U(t+\delta t, t)</math>. That is, given
    !% <math>\psi(\tau)</math> and <math>H(\tau)</math> for <math>tau \le t</math>,
    !% calculate <math>t+\delta t</math>. Note that in general the Hamiltonian
    !% is not known at times in the interior of the interval <math>[t,t+\delta t]</math>.
    !% This is due to the self-consistent nature of the time-dependent Kohn-Sham problem:
    !% the Hamiltonian at a given time <math>\tau</math> is built from the
    !% "solution" wavefunctions at that time.
    !%
    !% Some methods, however, do require the knowledge of the Hamiltonian at some
    !% point of the interval <math>[t,t+\delta t]</math>. This problem is solved by making
    !% use of extrapolation: given a number <math>l</math> of time steps previous to time
    !% <math>t</math>, this information is used to build the Hamiltonian at arbitrary times
    !% within <math>[t,t+\delta t]</math>. To be fully precise, one should then proceed
    !% <i>self-consistently</i>: the obtained Hamiltonian at time <math>t+\delta t</math>
    !% may then be used to interpolate the Hamiltonian, and repeat the evolution
    !% algorithm with this new information. Whenever iterating the procedure does
    !% not change the solution wavefunctions, the cycle is stopped. In practice,
    !% in <tt>Octopus</tt> we perform a second-order extrapolation without a
    !% self-consistency check, except for the first two iterations, where obviously
    !% the extrapolation is not reliable.
    !%
    !% The proliferation of methods is certainly excessive. The reason for it is that 
    !% the propagation algorithm is currently a topic of active development. We
    !% hope that in the future the optimal schemes are clearly identified. In the
    !% mean time, if you do not feel like testing, use the default choices and
    !% make sure the time step is small enough.
    !%Option etrs 2
    !% The idea is to make use of time-reversal symmetry from the beginning:
    !%
    !% <MATH>
    !%   \exp \left(-i\delta t/2 H_{n}\right)\psi_n = exp \left(i\delta t/2 H_{n+1}\right)\psi_{n+1},
    !% </MATH>
    !%
    !% and then invert to obtain:
    !%
    !% <MATH>
    !%   \psi_{n+1} = \exp \left(-i\delta t/2 H_{n+1}\right) exp \left(-i\delta t/2 H_{n}\right)\psi_{n}.
    !% </MATH>
    !%
    !% But we need to know <math>H_{n+1}</math>, which can only be known exactly through the solution
    !% <math>\psi_{n+1}</math>. What we do is to estimate it by performing a single exponential:
    !% <math>\psi^{*}_{n+1}=\exp \left( -i\delta t H_{n} \right) \psi_n</math>, and then
    !% <math>H_{n+1} = H[\psi^{*}_{n+1}]</math>. Thus no extrapolation is performed in this case.
    !%Option aetrs 3
    !% Approximated Enforced Time-Reversal Symmetry (AETRS).
    !% A modification of previous method to make it faster.
    !% It is based on extrapolation of the time-dependent potentials. It is faster
    !% by about 40%.
    !%
    !% The only difference is the procedure to estimate @math{H_{n+1}}: in this case
    !% it is extrapolated via a second-order polynomial by making use of the
    !% Hamiltonian at time @math{t-2\delta t}, @math{t-\delta t} and @math{t}.
    !%Option caetrs 12
    !% (experimental) Corrected Approximated Enforced Time-Reversal
    !% Symmetry (AETRS), this is the previous propagator but including
    !% a correction step to the exponential.
    !%Option exp_mid 4
    !% Exponential Midpoint Rule (EM).
    !% This is maybe the simplest method, but it is very well grounded theoretically:
    !% it is unitary (if the exponential is performed correctly) and preserves
    !% time-reversal symmetry (if the self-consistency problem is dealt with correctly).
    !% It is defined as:
    !%
    !% <MATH>
    !%   U_{\rm EM}(t+\delta t, t) = \exp \left( -i\delta t H_{t+\delta t/2}\right)\,.
    !% </MATH>
    !%Option crank_nicholson 5
    !%Option crank_nicolson 5
    !% Classical Crank-Nicolson propagator.
    !%
    !% <MATH>
    !%  (1 + i\delta t/2 H_{n+1/2}) \psi_{n+1} = (1 - i\delta t/2 H_{n+1/2}) \psi_{n}  
    !% </MATH>
    !%Option crank_nicholson_sparskit 6
    !%Option crank_nicolson_sparskit 6
    !% Classical Crank-Nicolson propagator. Requires the SPARSKIT library.
    !%
    !% <MATH>
    !%  (1 + i\delta t/2 H_{n+1/2}) \psi_{n+1} = (1 - i\delta t/2 H_{n+1/2}) \psi_{n}  
    !% </MATH>
    !%Option magnus 7
    !% Magnus Expansion (M4).
    !% This is the most sophisticated approach. It is a fourth-order scheme (a feature
    !% which it shares with the ST scheme; the other schemes are in principle second-order).
    !% It is tailored for making use of very large time steps, or equivalently,
    !% dealing with problem with very high-frequency time-dependence.
    !% It is still in a experimental state; we are not yet sure of when it is
    !% advantageous.
    !%Option crank_nicholson_src_mem 8
    !%Option crank_nicolson_src_mem 8
    !% Crank-Nicolson propagator with source and memory term for transport
    !% calculations.
    !%Option qoct_tddft_propagator 10
    !% WARNING: EXPERIMENTAL
    !%Option runge_kutta4 13
    !% WARNING: EXPERIMENTAL. Implicit Gauss-Legendre 4th order Runge-Kutta.
    !%Option runge_kutta2 14
    !% WARNING: EXPERIMENTAL. Implicit 2nd order Runge-Kutta (trapezoidal rule).
    !% Similar, but not identical, to Crank-Nicolson method.
    !%Option expl_runge_kutta4 15
    !% WARNING: EXPERIMENTAL. Explicit RK4 method.
    !%End
    call messages_obsolete_variable('TDEvolutionMethod', 'TDPropagator')

    default_propagator = PROP_ETRS
    if(gr%ob_grid%open_boundaries) default_propagator = PROP_CRANK_NICOLSON_SRC_MEM

    call parse_integer(datasets_check('TDPropagator'), default_propagator, tr%method)
    if(.not.varinfo_valid_option('TDPropagator', tr%method)) call input_error('TDPropagator')

    if(gr%ob_grid%open_boundaries.and.tr%method /= PROP_CRANK_NICOLSON_SRC_MEM) then
      message(1) = 'The time-evolution method for time-dependent run cannot'
      message(2) = 'be chosen freely. The Crank-Nicolson propagator'
      message(3) = 'with source and memory term has to be used. Either set'
      message(4) = ''
      message(5) = '  TDPropagator = crank_nicolson_src_mem'
      message(6) = ''
      message(7) = 'in your input or remove the TDPropagator variable.'
      call messages_fatal(7)
    end if

    select case(tr%method)
    case(PROP_ETRS)
    case(PROP_AETRS, PROP_CAETRS)
    case(PROP_EXPONENTIAL_MIDPOINT)
    case(PROP_CRANK_NICOLSON)
      ! set up pointer for zmf_dotu_aux, zmf_nrm2_aux
      call mesh_init_mesh_aux(gr%mesh)
    case(PROP_RUNGE_KUTTA4)
#ifdef HAVE_SPARSKIT
      ! set up pointer for zmf_dotu_aux, zmf_nrm2_aux
      call mesh_init_mesh_aux(gr%mesh)
      sp_distdot_mode = 3
      tr%tdsk_size = 2 * st%d%dim * gr%mesh%np * (st%st_end - st%st_start + 1) * (st%d%kpt%end - st%d%kpt%start + 1)
      SAFE_ALLOCATE(tr%tdsk)
      call sparskit_solver_init(tr%tdsk_size, tr%tdsk, is_complex = .true.)
#else
      message(1) = 'Octopus was not compiled with support for the SPARSKIT library. This'
      message(2) = 'library is required if the "runge_kutta4" propagator is selected.'
      message(3) = 'Try using a different propagation scheme or recompile with SPARSKIT support.'
      call messages_fatal(3)
#endif
    case(PROP_RUNGE_KUTTA2)
#ifdef HAVE_SPARSKIT
      ! set up pointer for zmf_dotu_aux, zmf_nrm2_aux
      call mesh_init_mesh_aux(gr%mesh)
      sp_distdot_mode = 2
      tr%tdsk_size = st%d%dim * gr%mesh%np * (st%st_end - st%st_start + 1) * (st%d%kpt%end - st%d%kpt%start + 1)
      SAFE_ALLOCATE(tr%tdsk)
      call sparskit_solver_init(tr%tdsk_size, tr%tdsk, is_complex = .true.)
#else
      message(1) = 'Octopus was not compiled with support for the SPARSKIT library. This'
      message(2) = 'library is required if the "runge_kutta2" propagator is selected.'
      message(3) = 'Try using a different propagation scheme or recompile with SPARSKIT support.'
      call messages_fatal(3)
#endif
    case(PROP_CRANK_NICOLSON_SPARSKIT)
#ifdef HAVE_SPARSKIT
      ! set up pointer for zmf_dotu_aux
      call mesh_init_mesh_aux(gr%mesh)
      sp_distdot_mode = 1
      tr%tdsk_size = st%d%dim*gr%mesh%np
      SAFE_ALLOCATE(tr%tdsk)
      call sparskit_solver_init(st%d%dim*gr%mesh%np, tr%tdsk, is_complex = .true.)
#else
      message(1) = 'Octopus was not compiled with support for the SPARSKIT library. This'
      message(2) = 'library is required if the "crank_nicolson_sparskit" propagator is selected.'
      message(3) = 'Try using a different propagation scheme or recompile with SPARSKIT support.'
      call messages_fatal(3)
#endif
    case(PROP_MAGNUS)
      SAFE_ALLOCATE(tr%vmagnus(1:gr%mesh%np, 1:st%d%nspin, 1:2))
    case(PROP_CRANK_NICOLSON_SRC_MEM)
      call ob_propagator_init(st, gr, hm, tr%ob, dt, max_iter)
    case(PROP_QOCT_TDDFT_PROPAGATOR)
    case(PROP_EXPLICIT_RUNGE_KUTTA4)
    case default
      call input_error('TDPropagator')
    end select
    call messages_print_var_option(stdout, 'TDPropagator', tr%method)

    if(have_fields) then
      if(tr%method /= PROP_ETRS .and.    &
         tr%method /= PROP_AETRS .and. &
         tr%method /= PROP_EXPONENTIAL_MIDPOINT .and. &
         tr%method /= PROP_QOCT_TDDFT_PROPAGATOR .and. &
         tr%method /= PROP_CRANK_NICOLSON .and. &
         tr%method /= PROP_RUNGE_KUTTA4 .and. &
         tr%method /= PROP_EXPLICIT_RUNGE_KUTTA4 .and. &
         tr%method /= PROP_RUNGE_KUTTA2 .and. &
         tr%method /= PROP_CRANK_NICOLSON_SPARSKIT ) then
        message(1) = "To move the ions or put in a gauge field, use the etrs, aetrs or exp_mid propagators." 
        call messages_fatal(1)
      end if
    end if

    call vksinterp_init(tr%vksold, cmplxscl, gr%mesh%np, st%d%nspin)

    call exponential_init(tr%te) ! initialize propagator

    call messages_obsolete_variable('TDSelfConsistentSteps', 'TDStepsWithSelfConsistency')

    !%Variable TDStepsWithSelfConsistency
    !%Type integer
    !%Default 3
    !%Section Time-Dependent::Propagation
    !%Description
    !% Since the KS propagator is non-linear, each propagation step
    !% should be performed self-consistently.  In practice, for most
    !% purposes this is not necessary, except perhaps in the first
    !% iterations. This variable holds the number of propagation steps
    !% for which the propagation is done self-consistently. 
    !%
    !% The special value <tt>all_steps</tt> forces self-consistency to
    !% be imposed on all propagation steps. A value of 0 means that
    !% self-consistency will not be imposed.  The default is 3, which
    !% means that self-consistency is only enforced during the first three
    !% steps.
    !%Option all_steps -1
    !% Self-consistency is imposed for all propagation steps.
    !%End

    call parse_integer(datasets_check('TDStepsWithSelfConsistency'), 3, tr%scf_propagation_steps)
    if(tr%scf_propagation_steps == -1) tr%scf_propagation_steps = HUGE(tr%scf_propagation_steps)
    if(tr%scf_propagation_steps < 0) call input_error('TDStepsWithSelfConsistency')

    !%Variable TDSCFThreshold
    !%Type float
    !%Default 1.0e-3
    !%Section Time-Dependent::Propagation
    !%Description
    !% Since the KS propagator is non-linear, each propagation stepp
    !% should be performed self-consistently.  In practice, for most
    !% purposes this is not necessary, except perhaps in the first
    !% iterations. This variable holds the number of propagation steps
    !% for which the propagation is done self-consistently. 
    !%
    !% The self consistency has to be measured against some accuracy 
    !% threshold. This variable controls the value of that threshold.
    !%End
    call parse_float(datasets_check('TDSCFThreshold'), CNST(1.0e-3), tr%scf_threshold)

    POP_SUB(propagator_init)
  end subroutine propagator_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_set_scf_prop(tr, threshold)
    type(propagator_t), intent(inout) :: tr
    FLOAT, intent(in), optional :: threshold

    PUSH_SUB(propagator_set_scf_prop)

    tr%scf_propagation_steps = huge(1)
    if(present(threshold)) then
      tr%scf_threshold = threshold
    end if

    POP_SUB(propagator_set_scf_prop)
  end subroutine propagator_set_scf_prop
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_remove_scf_prop(tr)
    type(propagator_t), intent(inout) :: tr

    PUSH_SUB(propagator_remove_scf_prop)

    tr%scf_propagation_steps = -1

    POP_SUB(propagator_remove_scf_prop)
  end subroutine propagator_remove_scf_prop
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_end(tr)
    type(propagator_t), intent(inout) :: tr

    PUSH_SUB(propagator_end)

    call vksinterp_end(tr%vksold)

    select case(tr%method)
    case(PROP_MAGNUS)
      ASSERT(associated(tr%vmagnus))
      SAFE_DEALLOCATE_P(tr%vmagnus)
#ifdef HAVE_SPARSKIT
    case(PROP_RUNGE_KUTTA4, PROP_RUNGE_KUTTA2, PROP_CRANK_NICOLSON_SPARSKIT)
      call sparskit_solver_end(tr%tdsk)
      SAFE_DEALLOCATE_P(tr%tdsk)
#endif
    case(PROP_CRANK_NICOLSON_SRC_MEM)
      call ob_propagator_end(tr%ob)
    end select
    
    call exponential_end(tr%te)       ! clean propagator method

    POP_SUB(propagator_end)
  end subroutine propagator_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_run_zero_iter(hm, gr, tr)
    type(hamiltonian_t),  intent(in)    :: hm
    type(grid_t),         intent(in)    :: gr
    type(propagator_t),   intent(inout) :: tr

    integer :: ip, ispin, idim
    PUSH_SUB(propagator_run_zero_iter)

    if(hm%cmplxscl%space) then
      call vksinterp_run_zero_iter(tr%vksold, gr%mesh%np, hm%d%nspin, hm%vhxc, hm%imvhxc)   
    else
      call vksinterp_run_zero_iter(tr%vksold, gr%mesh%np, hm%d%nspin, hm%vhxc)   
    end if

    POP_SUB(propagator_run_zero_iter)
  end subroutine propagator_run_zero_iter


  ! ---------------------------------------------------------
  !> Propagates st from time - dt to t.
  !! If dt<0, it propagates *backwards* from t+|dt| to t
  ! ---------------------------------------------------------
  subroutine propagator_dt(ks, hm, gr, st, tr, time, dt, mu, max_iter, nt, ions, geo, &
    gauge_force, scsteps, update_energy, qcchi)
    type(v_ks_t), target,            intent(inout) :: ks
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    type(propagator_t),  target,     intent(inout) :: tr
    FLOAT,                           intent(in)    :: time
    FLOAT,                           intent(in)    :: dt
    FLOAT,                           intent(in)    :: mu
    integer,                         intent(in)    :: max_iter
    integer,                         intent(in)    :: nt
    type(ion_dynamics_t),            intent(inout) :: ions
    type(geometry_t),                intent(inout) :: geo
    type(gauge_force_t),  optional,  intent(inout) :: gauge_force
    integer,              optional,  intent(out)   :: scsteps
    logical,              optional,  intent(in)    :: update_energy
    type(opt_control_state_t), optional, target, intent(inout) :: qcchi

    integer :: is, iter, ik, ist, ispin
    FLOAT   :: d, d_max
    logical :: self_consistent, cmplxscl, generate, update_energy_
    CMPLX, allocatable :: zpsi1(:, :, :, :)
    FLOAT, allocatable :: dtmp(:), vaux(:, :), imvaux(:, :)
    type(profile_t), save :: prof

    call profiling_in(prof, "TD_PROPAGATOR")
    PUSH_SUB(propagator_dt)

    cmplxscl = hm%cmplxscl%space

    if(gauge_field_is_applied(hm%ep%gfield)) then
      ASSERT(present(gauge_force))
    end if

    self_consistent = self_consistent_step()
    if(self_consistent) then
      SAFE_ALLOCATE(zpsi1(1:gr%mesh%np, 1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
           call states_get_state(st, gr%mesh, ist, ik, zpsi1(:, :, ist, ik))
        end do
      end do
      SAFE_ALLOCATE(vaux(1:gr%mesh%np, 1:st%d%nspin))
      vaux(1:gr%mesh%np, 1:st%d%nspin) = hm%vhxc(1:gr%mesh%np, 1:st%d%nspin)
      if (cmplxscl) then
        SAFE_ALLOCATE(Imvaux(1:gr%mesh%np, 1:st%d%nspin))
        Imvaux(1:gr%mesh%np, 1:st%d%nspin) = hm%Imvhxc(1:gr%mesh%np, 1:st%d%nspin)
      end if 
    end if

    if(cmplxscl) then
      call vksinterp_new(tr%vksold, gr%mesh%np, st%d%nspin, time, dt, hm%vhxc, hm%imvhxc)
    else
      call vksinterp_new(tr%vksold, gr%mesh%np, st%d%nspin, time, dt, hm%vhxc)
    end if

    select case(tr%method)
    case(PROP_ETRS)
      call td_etrs(ks, hm, gr, st, tr, time, dt, mu, ions, geo, gauge_force)
    case(PROP_AETRS, PROP_CAETRS)
      call td_aetrs(ks, hm, gr, st, tr, time, dt, mu, ions, geo, gauge_force)
    case(PROP_EXPONENTIAL_MIDPOINT)
      call exponential_midpoint(ks, hm, gr, st, tr, time, dt, mu, ions, geo, gauge_force)
    case(PROP_CRANK_NICOLSON)
      call td_crank_nicolson(ks, hm, gr, st, tr, time, dt, ions, geo, .false.)
    case(PROP_RUNGE_KUTTA4)
      call td_runge_kutta4(ks, hm, gr, st, tr, time, dt, ions, geo)
    case(PROP_RUNGE_KUTTA2)
      call td_runge_kutta2(ks, hm, gr, st, tr, time, dt, ions, geo)
    case(PROP_CRANK_NICOLSON_SPARSKIT)
      call td_crank_nicolson(ks, hm, gr, st, tr, time, dt, ions, geo, .true.)
    case(PROP_MAGNUS)
      call td_magnus(ks, hm, gr, st, tr, time, dt)
    case(PROP_CRANK_NICOLSON_SRC_MEM)
      call td_crank_nicolson_src_mem(hm, gr, st, tr, max_iter, nt, time, dt)
    case(PROP_QOCT_TDDFT_PROPAGATOR)
      call td_qoct_tddft_propagator(hm, ks%xc, gr, st, tr, time, dt, ions, geo)
    case(PROP_EXPLICIT_RUNGE_KUTTA4)
      if(present(qcchi)) then
        call td_explicit_runge_kutta4(ks, hm, gr, st, tr, time, dt, ions, geo, qcchi)
      else
        call td_explicit_runge_kutta4(ks, hm, gr, st, tr, time, dt, ions, geo)
      end if
    end select

    if(present(scsteps)) scsteps = 1

    if(self_consistent) then

      ! First, compare the new potential to the extrapolated one.
      call v_ks_calc(ks, hm, st, geo, time = time - dt)
      call vksinterp_diff(tr%vksold, gr, st%d%nspin, hm%vhxc, 0, d_max)

      if(d_max > tr%scf_threshold) then

        ! We do a maximum of 10 iterations. If it still not converged, probably the propagation
        ! will not be good anyways.
        do iter = 1, 10
          if(present(scsteps)) INCR(scsteps, 1)

          do ik = st%d%kpt%start, st%d%kpt%end
            do ist = st%st_start, st%st_end
              call states_set_state(st, gr%mesh, ist, ik, zpsi1(:, :, ist, ik))
            end do
          end do
          call vksinterp_set(tr%vksold, gr%mesh%np, st%d%nspin, 0, hm%vhxc)
          vaux(:, :) = hm%vhxc(:, :)

          select case(tr%method)
          case(PROP_ETRS)
            call td_etrs(ks, hm, gr, st, tr, time, dt, mu, ions, geo, gauge_force)
          case(PROP_AETRS, PROP_CAETRS)
            call td_aetrs(ks, hm, gr, st, tr, time, dt, mu, ions, geo, gauge_force)
          case(PROP_EXPONENTIAL_MIDPOINT)
            call exponential_midpoint(ks, hm, gr, st, tr, time, dt, mu, ions, geo, gauge_force)
          case(PROP_CRANK_NICOLSON)
            call td_crank_nicolson(ks, hm, gr, st, tr, time, dt, ions, geo, .false.)
          case(PROP_RUNGE_KUTTA4)
            call td_runge_kutta4(ks, hm, gr, st, tr, time, dt, ions, geo)
          case(PROP_RUNGE_KUTTA2)
            call td_runge_kutta2(ks, hm, gr, st, tr, time, dt, ions, geo)
          case(PROP_CRANK_NICOLSON_SPARSKIT)
            call td_crank_nicolson(ks, hm, gr, st, tr, time, dt, ions, geo, .true.)
          case(PROP_MAGNUS)
            call td_magnus(ks, hm, gr, st, tr, time, dt)
          case(PROP_CRANK_NICOLSON_SRC_MEM)
            call td_crank_nicolson_src_mem(hm, gr, st, tr, max_iter, nt, time, dt)
          case(PROP_QOCT_TDDFT_PROPAGATOR)
            call td_qoct_tddft_propagator(hm, ks%xc, gr, st, tr, time, dt, ions, geo)
          case(PROP_EXPLICIT_RUNGE_KUTTA4)
            call td_explicit_runge_kutta4(ks, hm, gr, st, tr, time, dt, ions, geo)
          end select

          call v_ks_calc(ks, hm, st, geo, time = time - dt)
          SAFE_ALLOCATE(dtmp(1:gr%mesh%np))
          d_max = M_ZERO
          do is = 1, st%d%nspin
            dtmp(1:gr%mesh%np) = hm%vhxc(1:gr%mesh%np, is) - vaux(1:gr%mesh%np, is)
            d = dmf_nrm2(gr%mesh, dtmp)
            if(d > d_max) d_max = d
          end do
          SAFE_DEALLOCATE_A(dtmp)

          if(d_max < tr%scf_threshold) exit
        end do
        
      end if

      SAFE_DEALLOCATE_A(zpsi1)
      SAFE_DEALLOCATE_A(vaux)
      SAFE_DEALLOCATE_A(Imvaux)
    end if
    
    generate = .false.
    if(ion_dynamics_ions_move(ions)) then
      if(.not. propagator_ions_are_propagated(tr)) then
        call ion_dynamics_propagate(ions, gr%sb, geo, abs(nt*dt), dt)
        generate = .true.
      end if
    end if

    if(gauge_field_is_applied(hm%ep%gfield) .and. .not. propagator_ions_are_propagated(tr)) then
      call gauge_field_propagate(hm%ep%gfield, gauge_force, dt)
    end if

    if(generate .or. geometry_species_time_dependent(geo)) then
      call hamiltonian_epot_generate(hm, gr, geo, st, time = abs(nt*dt))
    end if

    update_energy_ = optional_default(update_energy, .false.)
    call v_ks_calc(ks, hm, st, geo, calc_eigenval = update_energy_, time = abs(nt*dt), calc_energy = update_energy_)
    if(update_energy_) call energy_calc_total(hm, gr, st, iunit = -1)

    ! Recalculate forces, update velocities...
    if(ion_dynamics_ions_move(ions) .and. tr%method .ne. PROP_EXPLICIT_RUNGE_KUTTA4) then
      call forces_calculate(gr, geo, hm, st, abs(nt*dt), dt)
      call ion_dynamics_propagate_vel(ions, geo, atoms_moved = generate)
      if(generate) call hamiltonian_epot_generate(hm, gr, geo, st, time = abs(nt*dt))
      geo%kinetic_energy = ion_dynamics_kinetic_energy(geo)
    end if

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_get_force(gr, geo, hm%ep%proj, hm%phase, st, gauge_force)
      call gauge_field_propagate_vel(hm%ep%gfield, gauge_force, dt)
    end if

    POP_SUB(propagator_dt)
    call profiling_out(prof)

  contains

    ! ---------------------------------------------------------
    logical function self_consistent_step() result(scs)
      scs = .false.
      if( hm%theory_level /= INDEPENDENT_PARTICLES .and. &
          tr%method /= PROP_CAETRS .and. &
          tr%method /= PROP_RUNGE_KUTTA4 .and. &
          tr%method /= PROP_EXPLICIT_RUNGE_KUTTA4 .and. &
          tr%method /= PROP_RUNGE_KUTTA2 ) then
        if(time <= tr%scf_propagation_steps*abs(dt) + M_EPSILON) then
          scs = .true.
        end if
      end if
    end function self_consistent_step
    ! ---------------------------------------------------------

  end subroutine propagator_dt
  ! ---------------------------------------------------------




  ! ---------------------------------------------------------
  logical pure function propagator_ions_are_propagated(tr) result(propagated)
    type(propagator_t), intent(in) :: tr

    select case(tr%method)
    case(PROP_ETRS, PROP_AETRS, PROP_CAETRS, PROP_EXPLICIT_RUNGE_KUTTA4)
      propagated = .true.
    case default
      propagated = .false.
    end select

  end function propagator_ions_are_propagated
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------

  subroutine propagator_dt_bo(scf, gr, ks, st, hm, gauge_force, geo, mc, outp, iter, dt, ions, scsteps)
    type(scf_t), intent(inout)         :: scf
    type(grid_t), intent(inout)        :: gr
    type(v_ks_t), intent(inout)        :: ks
    type(states_t), intent(inout)      :: st
    type(hamiltonian_t), intent(inout) :: hm
    type(gauge_force_t), intent(inout) :: gauge_force
    type(geometry_t), intent(inout)    :: geo
    type(multicomm_t), intent(inout)   :: mc    !< index and domain communicators
    type(output_t), intent(inout)      :: outp
    integer, intent(in)                :: iter
    FLOAT, intent(in)                  :: dt
    type(ion_dynamics_t), intent(inout) :: ions
    integer, intent(inout)             :: scsteps

    PUSH_SUB(propagator_dt_bo)

    ! move the hamiltonian to time t
    call ion_dynamics_propagate(ions, gr%sb, geo, iter*dt, dt)
    call hamiltonian_epot_generate(hm, gr, geo, st, time = iter*dt)
    ! now calculate the eigenfunctions
    call scf_run(scf, mc, gr, geo, st, ks, hm, outp, &
      gs_run = .false., verbosity = VERB_COMPACT, iters_done = scsteps)

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_propagate(hm%ep%gfield, gauge_force, dt)
    end if

    call hamiltonian_epot_generate(hm, gr, geo, st, time = iter*dt)

    ! update Hamiltonian and eigenvalues (fermi is *not* called)
    call v_ks_calc(ks, hm, st, geo, calc_eigenval = .true., time = iter*dt, calc_energy = .true.)

    ! Get the energies.
    call energy_calc_total(hm, gr, st, iunit = -1)

    call ion_dynamics_propagate_vel(ions, geo)
    call hamiltonian_epot_generate(hm, gr, geo, st, time = iter*dt)
     geo%kinetic_energy = ion_dynamics_kinetic_energy(geo)

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_propagate_vel(hm%ep%gfield, gauge_force, dt)
    end if

    POP_SUB(propagator_dt_bo)
  end subroutine propagator_dt_bo

#include "propagator_cn_inc.F90"
#include "propagator_rk_inc.F90"
#include "propagator_magnus_inc.F90"
#include "propagator_qoct_inc.F90"
#include "propagator_expmid_inc.F90"
#include "propagator_etrs_inc.F90"

end module propagator_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

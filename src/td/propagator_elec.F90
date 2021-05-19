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

module propagator_elec_oct_m
  use energy_calc_oct_m
  use exponential_oct_m
  use forces_oct_m
  use gauge_field_oct_m
  use grid_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use ion_dynamics_oct_m
  use ions_oct_m
  use lda_u_oct_m
  use parser_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use opt_control_state_oct_m
  use output_oct_m
  use potential_interpolation_oct_m
  use profiling_oct_m
  use propagator_base_oct_m
  use propagator_cn_oct_m
  use propagator_etrs_oct_m
  use propagator_expmid_oct_m
  use propagator_magnus_oct_m
  use propagator_qoct_oct_m
  use propagator_rk_oct_m
  use propagator_verlet_oct_m
  use scf_oct_m
  use sparskit_oct_m
  use space_oct_m
  use states_elec_oct_m
  use v_ks_oct_m
  use varinfo_oct_m
  use xc_oct_m

  implicit none

  private
  public ::                         &
    propagator_elec_init,                &
    propagator_elec_end,                 &
    propagator_elec_copy,                &
    propagator_elec_run_zero_iter,       &
    propagator_elec_dt,                  &
    propagator_elec_set_scf_prop,        &
    propagator_elec_remove_scf_prop,     &
    propagator_elec_ions_are_propagated, &
    propagator_elec_dt_bo

contains

  ! ---------------------------------------------------------
  subroutine propagator_elec_copy(tro, tri)
    type(propagator_base_t), intent(inout) :: tro
    type(propagator_base_t), intent(in)    :: tri

    PUSH_SUB(propagator_elec_copy)
    
    tro%method = tri%method

    select case(tro%method)
    case(PROP_MAGNUS)
      SAFE_ALLOCATE_SOURCE_A(tro%vmagnus, tri%vmagnus)

    case(PROP_CRANK_NICOLSON_SPARSKIT)
      tro%tdsk_size = tri%tdsk_size
      call sparskit_solver_copy(tro%tdsk, tri%tdsk)

    case(PROP_RUNGE_KUTTA4)
      tro%tdsk_size = tri%tdsk_size
      call sparskit_solver_copy(tro%tdsk, tri%tdsk)

    case(PROP_RUNGE_KUTTA2)
      tro%tdsk_size = tri%tdsk_size
      call sparskit_solver_copy(tro%tdsk, tri%tdsk)

    end select

    call potential_interpolation_copy(tro%vksold, tri%vksold)

    call exponential_copy(tro%te, tri%te)
    tro%scf_propagation_steps = tri%scf_propagation_steps

    tro%scf_threshold = tri%scf_threshold
    POP_SUB(propagator_elec_copy)
  end subroutine propagator_elec_copy
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_elec_init(gr, namespace, st, tr, have_fields, family_is_mgga_with_exc)
    type(grid_t),            intent(in)    :: gr
    type(namespace_t),       intent(in)    :: namespace
    type(states_elec_t),     intent(in)    :: st
    type(propagator_base_t), intent(inout) :: tr
    !> whether there is an associated "field"
    !! that must be propagated (currently ions
    !! or a gauge field).
    logical,                 intent(in)    :: have_fields 
    logical,                 intent(in)    :: family_is_mgga_with_exc

    PUSH_SUB(propagator_elec_init)
    
    !%Variable TDPropagator
    !%Type integer
    !%Default etrs
    !%Section Time-Dependent::Propagation
    !%Description
    !% This variable determines which algorithm will be used to approximate
    !% the evolution operator <math>U(t+\delta t, t)</math>. That is, given
    !% <math>\psi(\tau)</math> and <math>H(\tau)</math> for <math>\tau \le t</math>,
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
    !% <math>
    !%   \exp \left(-i\delta t H_{n} / 2 \right)\psi_n = \exp \left(i\delta t H_{n+1} / 2 \right)\psi_{n+1},
    !% </math>
    !%
    !% and then invert to obtain:
    !%
    !% <math>
    !%   \psi_{n+1} = \exp \left(-i\delta t H_{n+1} / 2 \right) \exp \left(-i\delta t H_{n} / 2 \right)\psi_{n}.
    !% </math>
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
    !% The only difference is the procedure to estimate <math>H_{n+1}</math>: in this case
    !% it is extrapolated via a second-order polynomial by making use of the
    !% Hamiltonian at time <math>t-2\delta t</math>, <math>t-\delta t</math> and <math>t</math>.
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
    !% <math>
    !%   U_{\rm EM}(t+\delta t, t) = \exp \left( -i\delta t H_{t+\delta t/2}\right)\,.
    !% </math>
    !%Option crank_nicholson 5
    !%Option crank_nicolson 5
    !% Classical Crank-Nicolson propagator.
    !% <math>
    !%  (1 + i\delta t H_{n+1/2} / 2) \psi_{n+1} = (1 - i\delta t H_{n+1/2} / 2) \psi_{n}  
    !% </math>
    !%Option crank_nicholson_sparskit 6
    !%Option crank_nicolson_sparskit 6
    !% Classical Crank-Nicolson propagator. Requires the SPARSKIT library.
    !% <math>
    !%  (1 + i\delta t H_{n+1/2} / 2) \psi_{n+1} = (1 - i\delta t H_{n+1/2} / 2) \psi_{n}  
    !% </math>
    !%Option magnus 7
    !% Magnus Expansion (M4).
    !% This is the most sophisticated approach. It is a fourth-order scheme (a feature
    !% which it shares with the ST scheme; the other schemes are in principle second-order).
    !% It is tailored for making use of very large time steps, or equivalently,
    !% dealing with problem with very high-frequency time-dependence.
    !% It is still in a experimental state; we are not yet sure of when it is
    !% advantageous.
    !%Option qoct_tddft_propagator 10
    !% WARNING: EXPERIMENTAL
    !%Option runge_kutta4 13
    !% WARNING: EXPERIMENTAL. Implicit Gauss-Legendre 4th order Runge-Kutta.
    !%Option runge_kutta2 14
    !% WARNING: EXPERIMENTAL. Implicit 2nd order Runge-Kutta (trapezoidal rule).
    !% Similar, but not identical, to Crank-Nicolson method.
    !%Option expl_runge_kutta4 15
    !% WARNING: EXPERIMENTAL. Explicit RK4 method.
    !%Option cfmagnus4 16
    !% WARNING EXPERIMENTAL
    !%End
    call messages_obsolete_variable(namespace, 'TDEvolutionMethod', 'TDPropagator')

    call parse_variable(namespace, 'TDPropagator', PROP_ETRS, tr%method)
    if(.not.varinfo_valid_option('TDPropagator', tr%method)) call messages_input_error(namespace, 'TDPropagator')

    select case(tr%method)
    case(PROP_ETRS)
    case(PROP_AETRS)
    case(PROP_CAETRS)
      call messages_experimental("CAETRS propagator")
    case(PROP_EXPONENTIAL_MIDPOINT)
    case(PROP_CRANK_NICOLSON)
      ! set up pointer for zmf_dotu_aux, zmf_nrm2_aux
      call mesh_init_mesh_aux(gr%mesh)
    case(PROP_RUNGE_KUTTA4)
      ! set up pointer for zmf_dotu_aux, zmf_nrm2_aux
      call mesh_init_mesh_aux(gr%mesh)
      sp_distdot_mode = 3
      tr%tdsk_size = 2 * st%d%dim * gr%mesh%np * (st%st_end - st%st_start + 1) * (st%d%kpt%end - st%d%kpt%start + 1)
      call sparskit_solver_init(namespace, tr%tdsk_size, tr%tdsk, .true.)

#ifndef HAVE_SPARSKIT
      message(1) = 'Octopus was not compiled with support for the SPARSKIT library. This'
      message(2) = 'library is required if the "runge_kutta4" propagator is selected.'
      message(3) = 'Try using a different propagation scheme or recompile with SPARSKIT support.'
      call messages_fatal(3, namespace=namespace)
#endif

      call messages_experimental("Runge-Kutta 4 propagator")
    case(PROP_RUNGE_KUTTA2)

      ! set up pointer for zmf_dotu_aux, zmf_nrm2_aux
      call mesh_init_mesh_aux(gr%mesh)
      sp_distdot_mode = 2
      tr%tdsk_size = st%d%dim * gr%mesh%np * (st%st_end - st%st_start + 1) * (st%d%kpt%end - st%d%kpt%start + 1)
      call sparskit_solver_init(namespace, tr%tdsk_size, tr%tdsk, .true.)

#ifndef HAVE_SPARSKIT
      message(1) = 'Octopus was not compiled with support for the SPARSKIT library. This'
      message(2) = 'library is required if the "runge_kutta2" propagator is selected.'
      message(3) = 'Try using a different propagation scheme or recompile with SPARSKIT support.'
      call messages_fatal(3, namespace=namespace)
#endif

      call messages_experimental("Runge-Kutta 2 propagator")
    case(PROP_CRANK_NICOLSON_SPARSKIT)
      ! set up pointer for zmf_dotu_aux
      call mesh_init_mesh_aux(gr%mesh)
      sp_distdot_mode = 1
      tr%tdsk_size = st%d%dim*gr%mesh%np
      call sparskit_solver_init(namespace, st%d%dim*gr%mesh%np, tr%tdsk, .true.)

#ifndef HAVE_SPARSKIT
      message(1) = 'Octopus was not compiled with support for the SPARSKIT library. This'
      message(2) = 'library is required if the "crank_nicolson_sparskit" propagator is selected.'
      message(3) = 'Try using a different propagation scheme or recompile with SPARSKIT support.'
      call messages_fatal(3, namespace=namespace)
#endif
    case(PROP_MAGNUS)
      call messages_experimental("Magnus propagator")
      if (family_is_mgga_with_exc) then
        message(1) = "Magnus propagator with MGGA"
        call messages_fatal(1, namespace=namespace)
      end if
      SAFE_ALLOCATE(tr%vmagnus(1:gr%mesh%np, 1:st%d%nspin, 1:2))
    case(PROP_QOCT_TDDFT_PROPAGATOR)
      call messages_experimental("QOCT+TDDFT propagator")
    case(PROP_EXPLICIT_RUNGE_KUTTA4)
      call messages_experimental("explicit Runge-Kutta 4 propagator")
    case(PROP_CFMAGNUS4)
      call messages_experimental("Commutator-Free Magnus propagator")
    case default
      call messages_input_error(namespace, 'TDPropagator')
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
        call messages_fatal(1, namespace=namespace)
      end if
    end if

    select case(tr%method)
    case(PROP_CFMAGNUS4)
      call potential_interpolation_init(tr%vksold, gr%mesh%np, st%d%nspin, family_is_mgga_with_exc, order = 4)
    case default
      call potential_interpolation_init(tr%vksold, gr%mesh%np, st%d%nspin, family_is_mgga_with_exc)
    end select

    call exponential_init(tr%te, namespace) ! initialize propagator

    call messages_obsolete_variable(namespace, 'TDSelfConsistentSteps', 'TDStepsWithSelfConsistency')

    !%Variable TDStepsWithSelfConsistency
    !%Type integer
    !%Default 0
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
    !% self-consistency will not be imposed.  The default is 0.
    !%Option all_steps -1
    !% Self-consistency is imposed for all propagation steps.
    !%End

    call parse_variable(namespace, 'TDStepsWithSelfConsistency', 0, tr%scf_propagation_steps)

    if(tr%scf_propagation_steps == -1) tr%scf_propagation_steps = HUGE(tr%scf_propagation_steps)
    if(tr%scf_propagation_steps < 0) call messages_input_error(namespace, 'TDStepsWithSelfConsistency', 'Cannot be negative')

    if(tr%scf_propagation_steps /= 0) then
      call messages_experimental('TDStepsWithSelfConsistency')

      if(tr%method /= PROP_ETRS) then
        call messages_write('TDStepsWithSelfConsistency only works with the ETRS propagator')
        call messages_fatal(namespace=namespace)
      end if
    end if

    !%Variable TDSCFThreshold
    !%Type float
    !%Default 1.0e-6
    !%Section Time-Dependent::Propagation
    !%Description
    !% Since the KS propagator is non-linear, each propagation step
    !% should be performed self-consistently.  In practice, for most
    !% purposes this is not necessary, except perhaps in the first
    !% iterations. This variable holds the number of propagation steps
    !% for which the propagation is done self-consistently. 
    !%
    !% The self consistency has to be measured against some accuracy 
    !% threshold. This variable controls the value of that threshold.
    !%End
    call parse_variable(namespace, 'TDSCFThreshold', CNST(1.0e-6), tr%scf_threshold)

    POP_SUB(propagator_elec_init)
  end subroutine propagator_elec_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_elec_set_scf_prop(tr, threshold)
    type(propagator_base_t), intent(inout) :: tr
    FLOAT, optional,         intent(in)    :: threshold

    PUSH_SUB(propagator_elec_set_scf_prop)

    tr%scf_propagation_steps = huge(1)
    if(present(threshold)) then
      tr%scf_threshold = threshold
    end if

    POP_SUB(propagator_elec_set_scf_prop)
  end subroutine propagator_elec_set_scf_prop
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_elec_remove_scf_prop(tr)
    type(propagator_base_t), intent(inout) :: tr

    PUSH_SUB(propagator_elec_remove_scf_prop)

    tr%scf_propagation_steps = -1

    POP_SUB(propagator_elec_remove_scf_prop)
  end subroutine propagator_elec_remove_scf_prop
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_elec_end(tr)
    type(propagator_base_t), intent(inout) :: tr

    PUSH_SUB(propagator_elec_end)

    call potential_interpolation_end(tr%vksold)

    select case(tr%method)
    case(PROP_MAGNUS)
      ASSERT(allocated(tr%vmagnus))
      SAFE_DEALLOCATE_A(tr%vmagnus)

    case(PROP_RUNGE_KUTTA4, PROP_RUNGE_KUTTA2, PROP_CRANK_NICOLSON_SPARSKIT)
      call sparskit_solver_end(tr%tdsk)
      
    end select

    POP_SUB(propagator_elec_end)
  end subroutine propagator_elec_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_elec_run_zero_iter(hm, gr, tr)
    type(hamiltonian_elec_t),  intent(in)    :: hm
    type(grid_t),              intent(in)    :: gr
    type(propagator_base_t),   intent(inout) :: tr

    PUSH_SUB(propagator_elec_run_zero_iter)

    if (family_is_mgga_with_exc(hm%xc)) then
      call potential_interpolation_run_zero_iter(tr%vksold, gr%mesh%np, hm%d%nspin, &
              hm%vhxc, vtau = hm%vtau)   
    else
      call potential_interpolation_run_zero_iter(tr%vksold, gr%mesh%np, hm%d%nspin, &
                hm%vhxc)
    end if

    POP_SUB(propagator_elec_run_zero_iter)
  end subroutine propagator_elec_run_zero_iter


  ! ---------------------------------------------------------
  !> Propagates st from time - dt to t.
  !! If dt<0, it propagates *backwards* from t+|dt| to t
  ! ---------------------------------------------------------
  subroutine propagator_elec_dt(ks, namespace, space, hm, gr, st, tr, time, dt, ionic_scale, nt, ions_dyn, ions, outp, &
    scsteps, update_energy, qcchi)
    type(v_ks_t),                        target, intent(inout) :: ks
    type(namespace_t),                           intent(in)    :: namespace
    type(space_t),                               intent(in)    :: space
    type(hamiltonian_elec_t),            target, intent(inout) :: hm
    type(grid_t),                        target, intent(inout) :: gr
    type(states_elec_t),                 target, intent(inout) :: st
    type(propagator_base_t),             target, intent(inout) :: tr
    FLOAT,                                       intent(in)    :: time
    FLOAT,                                       intent(in)    :: dt
    FLOAT,                                       intent(in)    :: ionic_scale
    integer,                                     intent(in)    :: nt
    type(ion_dynamics_t),                        intent(inout) :: ions_dyn
    type(ions_t),                                intent(inout) :: ions
    type(output_t),                              intent(in)    :: outp
    integer,                   optional,         intent(out)   :: scsteps
    logical,                   optional,         intent(in)    :: update_energy
    type(opt_control_state_t), optional, target, intent(inout) :: qcchi

    logical :: generate, update_energy_, move_ions
    type(profile_t), save :: prof

    call profiling_in(prof, "TD_PROPAGATOR")
    PUSH_SUB(propagator_elec_dt)

    move_ions = ion_dynamics_ions_move(ions_dyn)

    update_energy_ = optional_default(update_energy, .true.)

    if (family_is_mgga_with_exc(hm%xc)) then
      call potential_interpolation_new(tr%vksold, gr%mesh%np, st%d%nspin, time, dt, &
                hm%vhxc, vtau = hm%vtau)
    else
      call potential_interpolation_new(tr%vksold, gr%mesh%np, st%d%nspin, time, dt, &
                hm%vhxc)
    end if

    if(present(scsteps)) scsteps = 1
   
    select case(tr%method)
    case(PROP_ETRS)
      if(self_consistent_step()) then
        call td_etrs_sc(ks, namespace, space, hm, gr, st, tr, time, dt, ionic_scale, ions_dyn, ions, move_ions, &
          tr%scf_threshold, scsteps)
      else
        call td_etrs(ks, namespace, space, hm, gr, st, tr, time, dt, ionic_scale, ions_dyn, ions, move_ions)
      end if
    case(PROP_AETRS)
      call td_aetrs(namespace, space, hm, gr, st, tr, time, dt, ionic_scale, ions_dyn, ions, move_ions)
    case(PROP_CAETRS)
      call td_caetrs(ks, namespace, space, hm, gr, st, tr, time, dt, ionic_scale, ions_dyn, ions, move_ions)
    case(PROP_EXPONENTIAL_MIDPOINT)
      call exponential_midpoint(hm, namespace, space, gr, st, tr, time, dt, ionic_scale, ions_dyn, ions, move_ions)
    case(PROP_CRANK_NICOLSON)
      call td_crank_nicolson(hm, namespace, space, gr, st, tr, time, dt, ions_dyn, ions, .false.)
    case(PROP_RUNGE_KUTTA4)
      call td_runge_kutta4(ks, namespace, space, hm, gr, st, tr, time, dt, ions_dyn, ions)
    case(PROP_RUNGE_KUTTA2)
      call td_runge_kutta2(ks, namespace, space, hm, gr, st, tr, time, dt, ions_dyn, ions)
    case(PROP_CRANK_NICOLSON_SPARSKIT)
      call td_crank_nicolson(hm, namespace, space, gr, st, tr, time, dt, ions_dyn, ions, .true.)
    case(PROP_MAGNUS)
      call td_magnus(hm, gr, st, tr, namespace, time, dt)
    case(PROP_QOCT_TDDFT_PROPAGATOR)
      call td_qoct_tddft_propagator(hm, namespace, space, gr, st, tr, time, dt, ions_dyn, ions)
    case(PROP_EXPLICIT_RUNGE_KUTTA4)
      if(present(qcchi)) then
        call td_explicit_runge_kutta4(ks, namespace, space, hm, gr, st, time, dt, ions_dyn, ions, qcchi)
      else
        call td_explicit_runge_kutta4(ks, namespace, space, hm, gr, st, time, dt, ions_dyn, ions)
      end if
    case(PROP_CFMAGNUS4)
      call td_cfmagnus4(ks, namespace, space, hm, gr, st, tr, time, dt, ions_dyn, ions, nt)
    end select

    generate = .false.
    if(move_ions) then
      if(.not. propagator_elec_ions_are_propagated(tr)) then
        call ion_dynamics_propagate(ions_dyn, ions, abs(nt*dt), ionic_scale*dt, namespace)
        generate = .true.
      end if
    end if

    if(gauge_field_is_applied(hm%ep%gfield) .and. .not. propagator_elec_ions_are_propagated(tr)) then
      call gauge_field_do_td(hm%ep%gfield, OP_VERLET_COMPUTE_ACC, dt, time, namespace)
    end if

    if(generate .or. ions%has_time_dependent_species()) then
      call hamiltonian_elec_epot_generate(hm, namespace,  space, gr, ions, st, time = abs(nt*dt))
    end if

    call v_ks_calc(ks, namespace, space, hm, st, ions, calc_eigenval = update_energy_, time = abs(nt*dt), &
      calc_energy = update_energy_)
    if(update_energy_) call energy_calc_total(namespace, space, hm, gr, st, iunit = -1)

    ! Recalculate forces, update velocities...
    if(move_ions .and. tr%method .ne. PROP_EXPLICIT_RUNGE_KUTTA4) then
      call forces_calculate(gr, namespace, ions, hm, st, ks, t = abs(nt*dt), dt = dt)
      call ion_dynamics_propagate_vel(ions_dyn, ions, atoms_moved = generate)
      if(generate) call hamiltonian_elec_epot_generate(hm, namespace, space, gr, ions, st, time = abs(nt*dt))
      ions%kinetic_energy = ion_dynamics_kinetic_energy(ions)
    else
      if (outp%what(OPTION__OUTPUT__FORCES)) then
        call forces_calculate(gr, namespace, ions, hm, st, ks, t = abs(nt*dt), dt = dt)
      end if
    end if

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_get_force(hm%ep%gfield, gr, st)
      call gauge_field_do_td(hm%ep%gfield, OP_VERLET_COMPUTE_VEL, dt, time, namespace)
    end if

    !We update the occupation matrices
    call lda_u_update_occ_matrices(hm%lda_u, namespace, gr%mesh, st, hm%hm_base, hm%energy)

    POP_SUB(propagator_elec_dt)
    call profiling_out(prof)

  contains

    ! ---------------------------------------------------------
    logical pure function self_consistent_step() result(scs)
      scs = (hm%theory_level /= INDEPENDENT_PARTICLES) .and. (time <= tr%scf_propagation_steps*abs(dt) + M_EPSILON)
    end function self_consistent_step
    ! ---------------------------------------------------------

  end subroutine propagator_elec_dt
  ! ---------------------------------------------------------




  ! ---------------------------------------------------------
  logical pure function propagator_elec_ions_are_propagated(tr) result(propagated)
    type(propagator_base_t), intent(in) :: tr

    select case(tr%method)
    case(PROP_ETRS, PROP_AETRS, PROP_CAETRS, PROP_EXPLICIT_RUNGE_KUTTA4)
      propagated = .true.
    case default
      propagated = .false.
    end select

  end function propagator_elec_ions_are_propagated
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------

  subroutine propagator_elec_dt_bo(scf, namespace, space, gr, ks, st, hm, ions, mc, outp, iter, dt, ions_dyn, scsteps)
    type(scf_t),              intent(inout) :: scf
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(grid_t),             intent(inout) :: gr
    type(v_ks_t),             intent(inout) :: ks
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(ions_t),             intent(inout) :: ions
    type(multicomm_t),        intent(inout) :: mc    !< index and domain communicators
    type(output_t),           intent(inout) :: outp
    integer,                  intent(in)    :: iter
    FLOAT,                    intent(in)    :: dt
    type(ion_dynamics_t),     intent(inout) :: ions_dyn
    integer,                  intent(inout) :: scsteps

    PUSH_SUB(propagator_elec_dt_bo)

    ! move the hamiltonian to time t
    call ion_dynamics_propagate(ions_dyn, ions, iter*dt, dt, namespace)
    call hamiltonian_elec_epot_generate(hm, namespace, space, gr, ions, st, time = iter*dt)
    ! now calculate the eigenfunctions
    call scf_run(scf, namespace, space, mc, gr, ions, st, ks, hm, outp, &
      gs_run = .false., verbosity = VERB_COMPACT, iters_done = scsteps)

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_do_td(hm%ep%gfield, OP_VERLET_COMPUTE_ACC, dt, iter*dt, namespace)
    end if

    !TODO: we should update the occupation matrices here 
    if(hm%lda_u_level /= DFT_U_NONE) then
      call messages_not_implemented("DFT+U with propagator_elec_dt_bo", namespace=namespace)
    end if

    call hamiltonian_elec_epot_generate(hm, namespace, space, gr, ions, st, time = iter*dt)

    ! update Hamiltonian and eigenvalues (fermi is *not* called)
    call v_ks_calc(ks, namespace, space, hm, st, ions, calc_eigenval = .true., time = iter*dt, calc_energy = .true.)

    ! Get the energies.
    call energy_calc_total(namespace, space, hm, gr, st, iunit = -1)

    call ion_dynamics_propagate_vel(ions_dyn, ions)
    call hamiltonian_elec_epot_generate(hm, namespace, space, gr, ions, st, time = iter*dt)
     ions%kinetic_energy = ion_dynamics_kinetic_energy(ions)

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_do_td(hm%ep%gfield, OP_VERLET_COMPUTE_VEL, dt, iter*dt, namespace)
    end if

    POP_SUB(propagator_elec_dt_bo)
  end subroutine propagator_elec_dt_bo

end module propagator_elec_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

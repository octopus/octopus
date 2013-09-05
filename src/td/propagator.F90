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
  use cpmd_m
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
  use output_m
  use profiling_m
  use scf_m
  use states_dim_m
  use solvers_m
  use sparskit_m
  use states_m
  use types_m
  use v_ks_m
  use varinfo_m

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
    propagator_dens_is_propagated,  &
    propagator_requires_vks,        &
    propagator_dt_cpmd,             &
    propagator_dt_bo


  integer, public, parameter ::        &
    PROP_ETRS                    = 2,  &
    PROP_AETRS                   = 3,  &
    PROP_EXPONENTIAL_MIDPOINT    = 4,  &
    PROP_CRANK_NICOLSON          = 5,  &
    PROP_CRANK_NICOLSON_SPARSKIT = 6,  &
    PROP_MAGNUS                  = 7,  &
    PROP_CRANK_NICOLSON_SRC_MEM  = 8,  &
    PROP_QOCT_TDDFT_PROPAGATOR   = 10, &
    PROP_CAETRS                  = 12

  FLOAT :: scf_threshold = CNST(1.0e-3)

  type propagator_t
    integer             :: method           !< Which evolution method to use.
    type(exponential_t) :: te               !< How to apply the propagator \f$ e^{-i H \Delta t} \f$.
    !> Storage of the KS potential of previous iterations.
    FLOAT, pointer      :: v_old(:, :, :) => null()
    FLOAT, pointer      :: Imv_old(:, :, :) => null()
    !> Auxiliary function to store the Magnus potentials.
    FLOAT, pointer      :: vmagnus(:, :, :) => null() 
    type(ob_terms_t)    :: ob               !< For open boundaries: leads, memory
    integer             :: scf_propagation_steps 
    logical             :: first
  end type propagator_t

#ifdef HAVE_SPARSKIT
  type(sparskit_solver_t), pointer, private :: tdsk
#endif
  type(grid_t),            pointer, private :: grid_p
  type(hamiltonian_t),     pointer, private :: hm_p
  type(propagator_t),      pointer, private :: tr_p
  integer,                 private :: ik_op, ist_op, idim_op, dim_op, nst_op
  type(states_t),          private :: st_op
  FLOAT,                   private :: t_op, dt_op

contains

  ! ---------------------------------------------------------
  subroutine propagator_nullify(this)
    type(propagator_t), intent(out) :: this
    !
    PUSH_SUB(propagator_nullify)
    !this%method
    !call exponential_nullify(this%te)
    this%v_old=>null()
    this%Imv_old=>null()
    this%vmagnus=>null() 
    !call ob_terms_nullify(this%ob)
    !this%scf_propagation_steps 
    !this%first
    POP_SUB(propagator_nullify)
    return
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
    end select

    call loct_pointer_copy(tro%v_old, tri%v_old)
    call loct_pointer_copy(tro%Imv_old, tri%Imv_old)
    call exponential_copy(tro%te, tri%te)
    tro%scf_propagation_steps = tri%scf_propagation_steps

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
    case(PROP_CRANK_NICOLSON_SPARSKIT)
#ifdef HAVE_SPARSKIT
      ! set up pointer for zmf_dotu_aux
      call mesh_init_mesh_aux(gr%mesh)
      SAFE_ALLOCATE(tdsk)
      call zsparskit_solver_init(st%d%dim*gr%mesh%np, tdsk)
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
    case default
      call input_error('TDPropagator')
    end select
    call messages_print_var_option(stdout, 'TDPropagator', tr%method)

    if(have_fields) then
      if(tr%method /= PROP_ETRS .and.    &
         tr%method /= PROP_AETRS .and. &
         tr%method /= PROP_EXPONENTIAL_MIDPOINT) then
        message(1) = "To move the ions or put in a gauge field, use the etrs, aetrs or exp_mid propagators." 
        call messages_fatal(1)
      end if
    end if

    ! Allocate memory to store the old KS potentials
    SAFE_ALLOCATE(tr%v_old(1:gr%mesh%np, 1:st%d%nspin, 0:3))
    tr%v_old(:, :, :) = M_ZERO
    if(cmplxscl) then
      SAFE_ALLOCATE(tr%Imv_old(1:gr%mesh%np, 1:st%d%nspin, 0:3))
      tr%Imv_old(:, :, :) = M_ZERO
    end if

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
      scf_threshold = threshold
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

    ! sanity check
    ASSERT(associated(tr%v_old)) 
    SAFE_DEALLOCATE_P(tr%v_old)         ! clean old KS potentials
    SAFE_DEALLOCATE_P(tr%Imv_old) 

    select case(tr%method)
    case(PROP_MAGNUS)
      ASSERT(associated(tr%vmagnus))
      SAFE_DEALLOCATE_P(tr%vmagnus)
    case(PROP_CRANK_NICOLSON_SPARSKIT)
#ifdef HAVE_SPARSKIT
      call zsparskit_solver_end(tdsk)
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

    forall(idim = 1:3, ispin = 1:hm%d%nspin, ip = 1:gr%mesh%np)
      tr%v_old(ip, ispin, idim) = hm%vhxc(ip, ispin)
    end forall
    if(hm%cmplxscl%space) then
      forall(idim = 1:3, ispin = 1:hm%d%nspin, ip = 1:gr%mesh%np)
        tr%Imv_old(ip, ispin, idim) = hm%Imvhxc(ip, ispin)
      end forall
    end if

    POP_SUB(propagator_run_zero_iter)
  end subroutine propagator_run_zero_iter


  ! ---------------------------------------------------------
  !> Propagates st from time - dt to t.
  !! If dt<0, it propagates *backwards* from t+|dt| to t
  ! ---------------------------------------------------------
  subroutine propagator_dt(ks, hm, gr, st, tr, time, dt, mu, max_iter, nt, ions, geo, gauge_force, scsteps, update_energy)
    type(v_ks_t),                    intent(inout) :: ks
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

    integer :: is, iter, ik, ist, ispin
    FLOAT   :: d, d_max
    logical :: self_consistent, cmplxscl, generate, update_energy_
    CMPLX, allocatable :: zpsi1(:, :, :, :)
    FLOAT, allocatable :: dtmp(:), vaux(:, :), vold(:, :)
    FLOAT, allocatable :: Imvaux(:, :), Imvold(:, :)
    type(profile_t), save :: prof

    call profiling_in(prof, "TD_PROPAGATOR")
    PUSH_SUB(propagator_dt)

    cmplxscl = hm%cmplxscl%space

    if(gauge_field_is_applied(hm%ep%gfield)) then
      ASSERT(present(gauge_force))
    end if

    self_consistent = .false.
    if(hm%theory_level /= INDEPENDENT_PARTICLES .and. tr%method /= PROP_CAETRS) then
      if(time <= tr%scf_propagation_steps*abs(dt) + M_EPSILON) then
        self_consistent = .true.
        SAFE_ALLOCATE(zpsi1(1:gr%mesh%np, 1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))

        do ik = st%d%kpt%start, st%d%kpt%end
          do ist = st%st_start, st%st_end
            call states_get_state(st, gr%mesh, ist, ik, zpsi1(:, :, ist, ik))
          end do
        end do

      end if
    end if

    SAFE_ALLOCATE(vaux(1:gr%mesh%np, 1:st%d%nspin))
    vaux(1:gr%mesh%np, 1:st%d%nspin) = hm%vhxc(1:gr%mesh%np, 1:st%d%nspin)
    if (cmplxscl) then
      SAFE_ALLOCATE(Imvaux(1:gr%mesh%np, 1:st%d%nspin))
      Imvaux(1:gr%mesh%np, 1:st%d%nspin) = hm%Imvhxc(1:gr%mesh%np, 1:st%d%nspin)
    end if 

    if(.not. propagator_requires_vks(tr)) then
      SAFE_ALLOCATE(vold(1:gr%mesh%np, 1:st%d%nspin))
      call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 1), vold)
      if(cmplxscl) then
        SAFE_ALLOCATE(Imvold(1:gr%mesh%np, 1:st%d%nspin))
        call lalg_copy(gr%mesh%np, st%d%nspin, tr%Imv_old(:, :, 1), Imvold)
      end if
    else
      call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 2), tr%v_old(:, :, 3))
      call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 1), tr%v_old(:, :, 2))
      call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc(:, :),     tr%v_old(:, :, 1))
      call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), tr%v_old(:, :, 1:3), time, tr%v_old(:, :, 0))
      if(cmplxscl) then
        call lalg_copy(gr%mesh%np, st%d%nspin, tr%Imv_old(:, :, 2), tr%Imv_old(:, :, 3))
        call lalg_copy(gr%mesh%np, st%d%nspin, tr%Imv_old(:, :, 1), tr%Imv_old(:, :, 2))
        call lalg_copy(gr%mesh%np, st%d%nspin, hm%Imvhxc(:, :),     tr%Imv_old(:, :, 1))
        call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), tr%Imv_old(:, :, 1:3), time, tr%Imv_old(:, :, 0))
      end if
    end if

    select case(tr%method)
    case(PROP_ETRS);                     call td_etrs()
    case(PROP_AETRS, PROP_CAETRS);       call td_aetrs()
    case(PROP_EXPONENTIAL_MIDPOINT);     call exponential_midpoint()
    case(PROP_CRANK_NICOLSON);           call td_crank_nicolson(.false.)
    case(PROP_CRANK_NICOLSON_SPARSKIT);  call td_crank_nicolson(.true.)
    case(PROP_MAGNUS);                   call td_magnus()
    case(PROP_CRANK_NICOLSON_SRC_MEM);   call td_crank_nicolson_src_mem()
    case(PROP_QOCT_TDDFT_PROPAGATOR)
      call td_qoct_tddft_propagator(hm, gr, st, tr, time, dt)
    end select

    if(present(scsteps)) scsteps = 1

    if(self_consistent) then

      ! First, compare the new potential to the extrapolated one.
      if(.not. cmplxscl) then
        call density_calc(st, gr, st%rho)
      else
        call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
      end if
      call v_ks_calc(ks, hm, st, geo, time = time - dt)
      SAFE_ALLOCATE(dtmp(1:gr%mesh%np))
      d_max = M_ZERO
      do is = 1, st%d%nspin
        dtmp(1:gr%mesh%np) = hm%vhxc(1:gr%mesh%np, is) - tr%v_old(1:gr%mesh%np, is, 0)
        d = dmf_nrm2(gr%mesh, dtmp)
        if(d > d_max) d_max = d
      end do
      SAFE_DEALLOCATE_A(dtmp)

      if(d_max > scf_threshold) then

        ! We do a maximum of 10 iterations. If it still not converged, probably the propagation
        ! will not be good anyways.
        do iter = 1, 10
          if(present(scsteps)) INCR(scsteps, 1)

          do ik = st%d%kpt%start, st%d%kpt%end
            do ist = st%st_start, st%st_end
              call states_set_state(st, gr%mesh, ist, ik, zpsi1(:, :, ist, ik))
            end do
          end do
          
          tr%v_old(:, :, 0) = hm%vhxc(:, :)
          vaux(:, :) = hm%vhxc(:, :)
          if(cmplxscl) then
            tr%v_old(:, :, 0) = hm%vhxc(:, :)
            vaux(:, :) = hm%vhxc(:, :)
          end if
          select case(tr%method)
          case(PROP_ETRS);                     call td_etrs()
          case(PROP_AETRS, PROP_CAETRS);       call td_aetrs()
          case(PROP_EXPONENTIAL_MIDPOINT);     call exponential_midpoint()
          case(PROP_CRANK_NICOLSON);           call td_crank_nicolson(.false.)
          case(PROP_CRANK_NICOLSON_SPARSKIT);  call td_crank_nicolson(.true.)
          case(PROP_MAGNUS);                   call td_magnus()
          case(PROP_CRANK_NICOLSON_SRC_MEM);   call td_crank_nicolson_src_mem()
          case(PROP_QOCT_TDDFT_PROPAGATOR)
            call td_qoct_tddft_propagator(hm, gr, st, tr, time, dt)
          end select

          if(.not. cmplxscl) then
            call density_calc(st, gr, st%rho)
          else
            call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
          end if
          call v_ks_calc(ks, hm, st, geo, time = time - dt)
          SAFE_ALLOCATE(dtmp(1:gr%mesh%np))
          d_max = M_ZERO
          do is = 1, st%d%nspin
            dtmp(1:gr%mesh%np) = hm%vhxc(1:gr%mesh%np, is) - vaux(1:gr%mesh%np, is)
            d = dmf_nrm2(gr%mesh, dtmp)
            if(d > d_max) d_max = d
          end do
          SAFE_DEALLOCATE_A(dtmp)

          if(d_max < scf_threshold) exit
        end do
        
      end if

      SAFE_DEALLOCATE_A(zpsi1)
    end if
    
    SAFE_DEALLOCATE_A(vold)
    SAFE_DEALLOCATE_A(vaux)
    
    SAFE_DEALLOCATE_A(Imvold)
    SAFE_DEALLOCATE_A(Imvaux)

    ! update density
    if(.not. propagator_dens_is_propagated(tr)) then 
      if(.not. cmplxscl) then
        call density_calc(st, gr, st%rho)
      else
        call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
      end if  
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

    if(update_energy_ .or. propagator_requires_vks(tr)) then

      ! save the vhxc potential for later
      if(.not. propagator_requires_vks(tr)) then
        SAFE_ALLOCATE(vold(1:gr%mesh%np, 1:st%d%nspin))
        do ispin = 1, st%d%nspin
          call lalg_copy(gr%mesh%np, hm%vhxc(:, ispin), vold(:, ispin))
        end do
        if(cmplxscl) then
          SAFE_ALLOCATE(Imvold(1:gr%mesh%np, 1:st%d%nspin))
          do ispin = 1, st%d%nspin
            call lalg_copy(gr%mesh%np, hm%Imvhxc(:, ispin), Imvold(:, ispin))
          end do
        end if
      end if

      ! update Hamiltonian and eigenvalues (fermi is *not* called)
      call v_ks_calc(ks, hm, st, geo, calc_eigenval = update_energy_, time = abs(nt*dt), calc_energy = update_energy_)

      ! Get the energies.
      if(update_energy_) call energy_calc_total(hm, gr, st, iunit = -1)

      ! restore the vhxc
      if(.not. propagator_requires_vks(tr)) then
        do ispin = 1, st%d%nspin
          call lalg_copy(gr%mesh%np, vold(:, ispin), hm%vhxc(:, ispin))
        end do
        SAFE_DEALLOCATE_A(vold)
        if(cmplxscl) then
          do ispin = 1, st%d%nspin
            call lalg_copy(gr%mesh%np, Imvold(:, ispin), hm%Imvhxc(:, ispin))
          end do
          SAFE_DEALLOCATE_A(Imvold)
        end if
      end if
    end if

    ! Recalculate forces, update velocities...
    if(ion_dynamics_ions_move(ions)) then
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
    !> Propagator with enforced time-reversal symmetry
    subroutine td_etrs()

      FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
      integer :: ik, ib
      type(batch_t) :: zpsib_save
      type(density_calc_t) :: dens_calc

      PUSH_SUB(propagator_dt.td_etrs)

      if(hm%theory_level /= INDEPENDENT_PARTICLES) then

        SAFE_ALLOCATE(vhxc_t1(1:gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE(vhxc_t2(1:gr%mesh%np, 1:st%d%nspin))
        call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t1)

        call density_calc_init(dens_calc, st, gr, st%rho)

        do ik = st%d%kpt%start, st%d%kpt%end
          do ib = st%group%block_start, st%group%block_end

            !save the state
            call batch_copy(st%group%psib(ib, ik), zpsib_save, reference = .false.)
            if(batch_is_packed(st%group%psib(ib, ik))) call batch_pack(zpsib_save, copy = .false.)
            call batch_copy_data(gr%der%mesh%np, st%group%psib(ib, ik), zpsib_save)

            !propagate the state dt with H(time - dt)
            call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, dt/mu, time - dt)
            call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik))

            !restore the saved state
            call batch_copy_data(gr%der%mesh%np, zpsib_save, st%group%psib(ib, ik))

            call batch_end(zpsib_save)

          end do
        end do

        call density_calc_end(dens_calc)

        call v_ks_calc(ks, hm, st, geo)

        call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t2)
        call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t1, hm%vhxc)
        call hamiltonian_update(hm, gr%mesh, time = time - dt)
      end if

      ! propagate dt/2 with H(time - dt)
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, dt/(mu*M_TWO), time - dt)
        end do
      end do

      ! propagate dt/2 with H(t)

      ! first move the ions to time t
      if(ion_dynamics_ions_move(ions)) then
        call ion_dynamics_propagate(ions, gr%sb, geo, time, dt)
        call hamiltonian_epot_generate(hm, gr, geo, st, time = time)
      end if

      if(gauge_field_is_applied(hm%ep%gfield)) then
        call gauge_field_propagate(hm%ep%gfield, gauge_force, dt)
      end if

      if(hm%theory_level /= INDEPENDENT_PARTICLES) then
        call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t2, hm%vhxc)
      end if
      call hamiltonian_update(hm, gr%mesh, time = time)

      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, dt/(M_TWO*mu), time)
        end do
      end do

      if(hm%theory_level /= INDEPENDENT_PARTICLES) then
        SAFE_DEALLOCATE_A(vhxc_t1)
        SAFE_DEALLOCATE_A(vhxc_t2)
      end if

      POP_SUB(propagator_dt.td_etrs)
    end subroutine td_etrs


    ! ---------------------------------------------------------
    !> Propagator with approximate enforced time-reversal symmetry
    subroutine td_aetrs()

      integer :: ik, ispin, ip, ist, ib
      FLOAT :: vv
      CMPLX :: phase
      type(density_calc_t)  :: dens_calc
      type(profile_t), save :: phase_prof
      integer               :: pnp, iprange
      type(opencl_mem_t)    :: phase_buff

      PUSH_SUB(propagator_dt.td_aetrs)

      if(tr%method == PROP_CAETRS) then
        call lalg_copy(gr%mesh%np, st%d%nspin, vold, hm%vhxc)
        if(cmplxscl) call lalg_copy(gr%mesh%np, st%d%nspin, Imvold, hm%Imvhxc)
        call hamiltonian_update(hm, gr%mesh, time = time - dt)
        call v_ks_calc_start(ks, hm, st, geo, time = time - dt, calc_energy = .false.)
      end if

      ! propagate half of the time step with H(time - dt)
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, dt/(M_TWO*mu), time - dt)
        end do
      end do

      if(tr%method == PROP_CAETRS) then
        call v_ks_calc_finish(ks, hm)
        call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 2), tr%v_old(:, :, 3))
        call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 1), tr%v_old(:, :, 2))
        call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc(:, :),     tr%v_old(:, :, 1))
        call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), tr%v_old(:, :, 1:3), time, tr%v_old(:, :, 0))
        forall(ispin = 1:st%d%nspin, ip = 1:gr%mesh%np) 
          vold(ip, ispin) =  dt/(M_TWO*mu)*(hm%vhxc(ip, ispin) - vold(ip, ispin))
        end forall

        ! copy vold to a cl buffer
        if(opencl_is_enabled() .and. hamiltonian_apply_packed(hm, gr%mesh)) then
#ifdef HAVE_OPENCL
          pnp = opencl_padded_size(gr%mesh%np)
          call opencl_create_buffer(phase_buff, CL_MEM_READ_ONLY, TYPE_FLOAT, pnp*st%d%nspin)
          ASSERT(ubound(vold, dim = 1) == gr%mesh%np)
          do ispin = 1, st%d%nspin
            call opencl_write_buffer(phase_buff, gr%mesh%np, vold(:, ispin), offset = (ispin - 1)*pnp)
          end do
#endif
        end if

      end if

      ! interpolate the Hamiltonian to time t
      call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 0), hm%vhxc)

      ! move the ions to time t
      if(ion_dynamics_ions_move(ions)) then
        call ion_dynamics_propagate(ions, gr%sb, geo, time, dt)
        call hamiltonian_epot_generate(hm, gr, geo, st, time = time)
      end if

      if(gauge_field_is_applied(hm%ep%gfield)) then
        call gauge_field_propagate(hm%ep%gfield, gauge_force, dt)
      end if

      call hamiltonian_update(hm, gr%mesh, time = time)

      call density_calc_init(dens_calc, st, gr, st%rho)

      ! propagate the other half with H(t)
      do ik = st%d%kpt%start, st%d%kpt%end
        ispin = states_dim_get_spin_index(st%d, ik)

        do ib = st%group%block_start, st%group%block_end
          if(hamiltonian_apply_packed(hm, gr%mesh)) call batch_pack(st%group%psib(ib, ik))
          
          if(tr%method == PROP_CAETRS) then
            call profiling_in(phase_prof, "CAETRS_PHASE")
            select case(batch_status(st%group%psib(ib, ik)))
            case(BATCH_NOT_PACKED)
              do ip = 1, gr%mesh%np
                vv = vold(ip, ispin)
                phase = TOCMPLX(cos(vv), -sin(vv))
                forall(ist = 1:st%group%psib(ib, ik)%nst_linear)
                  st%group%psib(ib, ik)%states_linear(ist)%zpsi(ip) = st%group%psib(ib, ik)%states_linear(ist)%zpsi(ip)*phase
                end forall
              end do
            case(BATCH_PACKED)
              do ip = 1, gr%mesh%np
                vv = vold(ip, ispin)
                phase = TOCMPLX(cos(vv), -sin(vv))
                forall(ist = 1:st%group%psib(ib, ik)%nst_linear)
                  st%group%psib(ib, ik)%pack%zpsi(ist, ip) = st%group%psib(ib, ik)%pack%zpsi(ist, ip)*phase
                end forall
              end do
            case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
              call opencl_set_kernel_arg(kernel_phase, 0, pnp*(ispin - 1))
              call opencl_set_kernel_arg(kernel_phase, 1, phase_buff)
              call opencl_set_kernel_arg(kernel_phase, 2, st%group%psib(ib, ik)%pack%buffer)
              call opencl_set_kernel_arg(kernel_phase, 3, log2(st%group%psib(ib, ik)%pack%size(1)))

              iprange = opencl_max_workgroup_size()/st%group%psib(ib, ik)%pack%size(1)

              call opencl_kernel_run(kernel_phase, (/st%group%psib(ib, ik)%pack%size(1), pnp/), &
                (/st%group%psib(ib, ik)%pack%size(1), iprange/))
#endif
            end select
            call profiling_out(phase_prof)
          end if

          call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, dt/(M_TWO*mu), time)
          call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik))

          if(hamiltonian_apply_packed(hm, gr%mesh)) call batch_unpack(st%group%psib(ib, ik))
        end do
      end do

#ifdef HAVE_OPENCL
      if(tr%method == PROP_CAETRS .and. opencl_is_enabled() .and. hamiltonian_apply_packed(hm, gr%mesh)) then
        call opencl_release_buffer(phase_buff)
      end if
#endif
      
      call density_calc_end(dens_calc)

      POP_SUB(propagator_dt.td_aetrs)
    end subroutine td_aetrs


    ! ---------------------------------------------------------
    !> Exponential midpoint
    subroutine exponential_midpoint()
      integer :: ib, ik
      type(ion_state_t) :: ions_state
      FLOAT :: vecpot(1:MAX_DIM), vecpot_vel(1:MAX_DIM)
      CMPLX :: zt, zdt

      PUSH_SUB(propagator_dt.exponential_midpoint)

      
      if(hm%cmplxscl%time) then
        zt =  TOCMPLX(time, M_ZERO) *exp(M_zI * TOCMPLX(hm%cmplxscl%alphaR, M_ZERO))
        zdt = TOCMPLX(dt,   M_ZERO) *exp(M_zI * TOCMPLX(hm%cmplxscl%alphaR, M_ZERO))
        
        !FIXME: not adapted yet
        if(hm%theory_level /= INDEPENDENT_PARTICLES) then
          call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%v_old(:, :, 0:2), time - dt/M_TWO, hm%vhxc(:, :))
          if(hm%cmplxscl%space) &
            call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%Imv_old(:, :, 0:2), time - dt/M_TWO, hm%Imvhxc(:, :))
        end if

        !FIXME: not implemented yet
        !move the ions to time 'time - dt/2'
        if(ion_dynamics_ions_move(ions)) then
          call ion_dynamics_save_state(ions, geo, ions_state)
          call ion_dynamics_propagate(ions, gr%sb, geo, time - dt/M_TWO, M_HALF*dt)
          call hamiltonian_epot_generate(hm, gr, geo, st, time = time - dt/M_TWO)
        end if
        
        !FIXME: not implemented yet      
        if(gauge_field_is_applied(hm%ep%gfield)) then
          vecpot = gauge_field_get_vec_pot(hm%ep%gfield)
          vecpot_vel = gauge_field_get_vec_pot_vel(hm%ep%gfield)
          call gauge_field_propagate(hm%ep%gfield, gauge_force, M_HALF*dt)
        end if

        call hamiltonian_update(hm, gr%mesh, time = real(zt - zdt/M_z2, REAL_PRECISION), Imtime = aimag(zt - zdt/M_z2  ))

        do ik = st%d%kpt%start, st%d%kpt%end
          do ib = st%group%block_start, st%group%block_end

            call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, &
              real(zdt/mu,REAL_PRECISION), real(zt - zdt/M_z2,REAL_PRECISION), &
              Imdeltat = aimag(zdt/mu), Imtime = aimag(zt -  zdt / M_z2 ) )

          end do
        end do
        
      else

        if(hm%theory_level /= INDEPENDENT_PARTICLES) then
          call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%v_old(:, :, 0:2), time - dt/M_TWO, hm%vhxc(:, :))
          if(hm%cmplxscl%space) &
            call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%Imv_old(:, :, 0:2), time - dt/M_TWO, hm%Imvhxc(:, :))
        end if

        !move the ions to time 'time - dt/2'
        if(ion_dynamics_ions_move(ions)) then
          call ion_dynamics_save_state(ions, geo, ions_state)
          call ion_dynamics_propagate(ions, gr%sb, geo, time - dt/M_TWO, M_HALF*dt)
          call hamiltonian_epot_generate(hm, gr, geo, st, time = time - dt/M_TWO)
        end if
      
        if(gauge_field_is_applied(hm%ep%gfield)) then
          vecpot = gauge_field_get_vec_pot(hm%ep%gfield)
          vecpot_vel = gauge_field_get_vec_pot_vel(hm%ep%gfield)
          call gauge_field_propagate(hm%ep%gfield, gauge_force, M_HALF*dt)
        end if

        call hamiltonian_update(hm, gr%mesh, time = time - M_HALF*dt)

        do ik = st%d%kpt%start, st%d%kpt%end
          do ib = st%group%block_start, st%group%block_end

            call exponential_apply_batch(tr%te, gr%der, hm, st%group%psib(ib, ik), ik, dt/mu, time - dt/M_TWO)

          end do
        end do
      end if

      if(hm%cmplxscl%space) then ! Propagate the left state
        ! (L(t+dt)| = (L|U(t-dt) = (L|e^{i H_\theta(t-dt/2) (dt)}
        
        if(hm%cmplxscl%time) then
          zt =  TOCMPLX(time,M_ZERO) *exp(M_zI * TOCMPLX(hm%cmplxscl%alphaL, M_ZERO))
          zdt = TOCMPLX(dt,  M_ZERO) *exp(M_zI * TOCMPLX(hm%cmplxscl%alphaL, M_ZERO))

          ! FIXME: check this interpolation!! 
          ! probably need some rethinking 
          if(hm%theory_level /= INDEPENDENT_PARTICLES) then
            call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%v_old(:, :, 0:2), time + dt/M_TWO, hm%vhxc(:, :))
            if(hm%cmplxscl%space) &
              call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%Imv_old(:, :, 0:2), time + dt/M_TWO, hm%Imvhxc(:, :))
          end if
        
          call hamiltonian_update(hm, gr%mesh, time = real(zt + zdt/M_z2, REAL_PRECISION), Imtime = aimag(zt + zdt/M_z2  ))
          
          do ik = st%d%kpt%start, st%d%kpt%end
            do ib = st%group%block_start, st%group%block_end
              call exponential_apply_batch(tr%te, gr%der, hm, st%psibL(ib, ik), ik,&
                real(-zdt/mu, REAL_PRECISION), real(zt + zdt/M_z2, REAL_PRECISION), &
                Imdeltat = aimag(-zdt/mu), Imtime = aimag(zt +  zdt / M_z2 ) )

            end do
          end do
        
        else
          
          ! FIXME: check this interpolation!! 
          ! probably need some rethinking 
          if(hm%theory_level /= INDEPENDENT_PARTICLES) then
            call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%v_old(:, :, 0:2), time + dt/M_TWO, hm%vhxc(:, :))
            if(hm%cmplxscl%space) &
              call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%Imv_old(:, :, 0:2), time + dt/M_TWO, hm%Imvhxc(:, :))
          end if

          call hamiltonian_update(hm, gr%mesh, time = time + M_HALF*dt)

          do ik = st%d%kpt%start, st%d%kpt%end
            do ib = st%group%block_start, st%group%block_end
              call exponential_apply_batch(tr%te, gr%der, hm, st%psibL(ib, ik), ik, -dt/mu, time + dt/M_TWO)

            end do
          end do
        
        end if
                
      end if

      !restore to time 'time - dt'
      if(ion_dynamics_ions_move(ions)) call ion_dynamics_restore_state(ions, geo, ions_state)

      if(gauge_field_is_applied(hm%ep%gfield)) then
        call gauge_field_set_vec_pot(hm%ep%gfield, vecpot)
        call gauge_field_set_vec_pot_vel(hm%ep%gfield, vecpot_vel)
        call hamiltonian_update(hm, gr%mesh)
      end if

      POP_SUB(propagator_dt.exponential_midpoint)
    end subroutine exponential_midpoint


    ! ---------------------------------------------------------
    !> Crank-Nicolson propagator
    subroutine td_crank_nicolson(use_sparskit)
      logical, intent(in) :: use_sparskit

      CMPLX, allocatable :: zpsi_rhs(:,:), zpsi(:), rhs(:), inhpsi(:)
      integer :: ik, ist, idim, ip, isize, np_part, np, iter
      FLOAT :: dres
      FLOAT :: cgtol = CNST(1.0e-12)
      logical :: converged

      PUSH_SUB(propagator_dt.td_crank_nicolson)

#ifndef HAVE_SPARSKIT
      if(use_sparskit) then
        message(1) = "Cannot use SPARSKIT in Crank-Nicolson propagator: not compiled with SPARSKIT support."
        call messages_fatal(1)
      endif
#endif

      np_part = gr%mesh%np_part
      np = gr%mesh%np
      isize = np_part*st%lnst*st%d%kpt%nlocal*st%d%dim

      ! define pointer and variables for usage in td_zop, td_zopt routines
      grid_p    => gr
      hm_p      => hm
      tr_p      => tr
      dt_op = dt
      t_op  = time - dt/M_TWO
      dim_op = st%d%dim

      ! we (ab)use exponential_apply to compute (1-i\delta t/2 H_n)\psi^n
      ! exponential order needs to be only 1
      tr%te%exp_method = EXP_TAYLOR
      tr%te%exp_order  = 1

      SAFE_ALLOCATE(zpsi_rhs(1:np_part, 1:st%d%dim))
      SAFE_ALLOCATE(zpsi(1:np*st%d%dim))
      SAFE_ALLOCATE(rhs(1:np*st%d%dim))
        
      call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%v_old(:, :, 0:2), time - dt/M_TWO, hm%vhxc(:, :))
      if(hm%cmplxscl%space) & 
        call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%Imv_old(:, :, 0:2), time - dt/M_TWO, hm%Imvhxc(:, :))
    
      call hamiltonian_update(hm, gr%mesh, time = time - dt/M_TWO)
      
      
      ! solve (1+i\delta t/2 H_n)\psi^{predictor}_{n+1} = (1-i\delta t/2 H_n)\psi^n
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end

          call states_get_state(st, gr%mesh, ist, ik, zpsi_rhs)
          call exponential_apply(tr%te, gr%der, hm, zpsi_rhs, ist, ik, dt/M_TWO, time - dt/M_TWO)

          if(hamiltonian_inh_term(hm)) then
            SAFE_ALLOCATE(inhpsi(1:gr%mesh%np))
            do idim = 1, st%d%dim
              call states_get_state(hm%inh_st, gr%mesh, idim, ist, ik, inhpsi)
              forall(ip = 1:gr%mesh%np) zpsi_rhs(ip, idim) = zpsi_rhs(ip, idim) + dt*inhpsi(ip)
            end do
            SAFE_DEALLOCATE_A(inhpsi)
          end if

          ! put the values in a continuous array
          do idim = 1, st%d%dim
            call states_get_state(st, gr%mesh, idim, ist, ik, zpsi((idim - 1)*np+1:idim*np))
            rhs((idim - 1)*np + 1:idim*np) = zpsi_rhs(1:np, idim)
          end do

          ist_op = ist
          ik_op = ik

          if(use_sparskit) then
#ifdef HAVE_SPARSKIT
            call zsparskit_solver_run(tdsk, td_zop, td_zopt, zpsi, rhs)
#endif
          else
            iter = 2000
            call zqmr_sym(np*st%d%dim, zpsi, rhs, propagator_qmr_op, zmf_dotu_aux, zmf_nrm2_aux, &
            propagator_qmr_prec, iter, dres, cgtol, showprogress = .false., converged = converged)

            if(.not.converged) then
              write(message(1),'(a)')        'The linear solver used for the Crank-Nicolson'
              write(message(2),'(a,es14.4)') 'propagator did not converge: Residual = ', dres
              call messages_warning(2)
            end if

          endif

          do idim = 1, st%d%dim
            call states_set_state(st, gr%mesh, idim, ist, ik, zpsi((idim-1)*np + 1:(idim - 1)*np + np))
          end do          

        end do
      end do

      if(hm%cmplxscl%space) then !Left states
        
        dt_op = - dt !propagate backwards
        t_op  = time + dt/M_TWO
        
        
        call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%v_old(:, :, 0:2), time + dt/M_TWO, hm%vhxc(:, :))
        if(hm%cmplxscl%space) & 
          call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%Imv_old(:, :, 0:2), time + dt/M_TWO, hm%Imvhxc(:, :))
    
        call hamiltonian_update(hm, gr%mesh, time = time + dt/M_TWO)
      
      
        ! solve (1+i\delta t/2 H_n)\psi^{predictor}_{n+1} = (1-i\delta t/2 H_n)\psi^n
        do ik = st%d%kpt%start, st%d%kpt%end
          do ist = st%st_start, st%st_end

            call states_get_state(st, gr%mesh, ist, ik, zpsi_rhs,left = .true. )
            call exponential_apply(tr%te, gr%der, hm, zpsi_rhs, ist, ik, -dt/M_TWO, time + dt/M_TWO)

            if(hamiltonian_inh_term(hm)) then
              SAFE_ALLOCATE(inhpsi(1:gr%mesh%np))
              do idim = 1, st%d%dim
                call states_get_state(hm%inh_st, gr%mesh, idim, ist, ik, inhpsi, left = .true.)
                forall(ip = 1:gr%mesh%np) zpsi_rhs(ip, idim) = zpsi_rhs(ip, idim) - dt*inhpsi(ip)
              end do
              SAFE_DEALLOCATE_A(inhpsi)
            end if

            ! put the values in a continuous array
            do idim = 1, st%d%dim
              call states_get_state(st, gr%mesh, idim, ist, ik, zpsi((idim - 1)*np+1:idim*np), left = .true.)
              rhs((idim - 1)*np + 1:idim*np) = zpsi_rhs(1:np, idim)
            end do

            ist_op = ist
            ik_op = ik
            if(use_sparskit) then
#ifdef HAVE_SPARSKIT
              call zsparskit_solver_run(tdsk, td_zop, td_zopt, zpsi, rhs)
#endif
            else
              iter = 2000
              call zqmr_sym(np*st%d%dim, zpsi, rhs, propagator_qmr_op, zmf_dotu_aux, zmf_nrm2_aux, &
              propagator_qmr_prec, iter, dres, cgtol, showprogress = .false., converged = converged)

              if(.not.converged) then
                write(message(1),'(a)')        'The linear solver used for the Crank-Nicolson'
                write(message(2),'(a,es14.4)') 'propagator did not converge for left states: Residual = ', dres
                call messages_warning(2)
              end if

            endif

            do idim = 1, st%d%dim
              call states_set_state(st, gr%mesh, idim, ist, ik, zpsi((idim-1)*np + 1:(idim - 1)*np + np), left = .true.)
            end do

          end do
        end do
        
        
      end if


      SAFE_DEALLOCATE_A(zpsi_rhs)
      SAFE_DEALLOCATE_A(zpsi)
      SAFE_DEALLOCATE_A(rhs)
      POP_SUB(propagator_dt.td_crank_nicolson)

    end subroutine td_crank_nicolson
    ! ---------------------------------------------------------



    ! ---------------------------------------------------------
    !> Magnus propagator
    subroutine td_magnus()

      integer :: j, is, ist, ik, i
      FLOAT :: atime(2)
      FLOAT, allocatable :: vaux(:, :, :), pot(:)

      PUSH_SUB(propagator_dt.td_magnus)

      SAFE_ALLOCATE(vaux(1:gr%mesh%np, 1:st%d%nspin, 1:2))

      atime(1) = (M_HALF-sqrt(M_THREE)/M_SIX)*dt
      atime(2) = (M_HALF+sqrt(M_THREE)/M_SIX)*dt

      if(hm%theory_level /= INDEPENDENT_PARTICLES) then
        do j = 1, 2
          call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%v_old(:, :, 0:2), atime(j) - dt, hm%vhxc(:, :))
          call hamiltonian_update(hm, gr%mesh)
        end do
      else
        vaux = M_ZERO
      end if

      do j = 1, 2
        ! WARNING: This should be carefully tested, and extended to allow for velocity-gauge laser fields.
        do i = 1, hm%ep%no_lasers
          select case(laser_kind(hm%ep%lasers(i)))
          case(E_FIELD_ELECTRIC)
            SAFE_ALLOCATE(pot(1:gr%mesh%np))
            pot = M_ZERO
            call laser_potential(hm%ep%lasers(i), gr%mesh, pot, time - dt + atime(j))
            do is = 1, st%d%nspin
              vaux(:, is, j) = vaux(:, is, j) + pot(:)
            end do
            SAFE_DEALLOCATE_A(pot)
          case(E_FIELD_MAGNETIC, E_FIELD_VECTOR_POTENTIAL)
            write(message(1),'(a)') 'The Magnus propagator cannot be used with magnetic fields, or'
            write(message(2),'(a)') 'with an electric field described in the velocity gauge.'
            call messages_fatal(2)
          end select
        end do
      end do

      tr%vmagnus(:, :, 2)  = M_HALF*(vaux(:, :, 1) + vaux(:, :, 2))
      tr%vmagnus(:, :, 1) = (sqrt(M_THREE)/CNST(12.0))*dt*(vaux(:, :, 2) - vaux(:, :, 1))

     do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          call exponential_apply(tr%te, gr%der, hm, st%zpsi(:,:, ist, ik), ist, ik, dt, M_ZERO, &
            vmagnus = tr%vmagnus)
        end do
      end do

      SAFE_DEALLOCATE_A(vaux)
      POP_SUB(propagator_dt.td_magnus)
    end subroutine td_magnus


    ! ---------------------------------------------------------
    !> Crank-Nicolson scheme with source and memory term.
    subroutine td_crank_nicolson_src_mem()
      PUSH_SUB(propagator_dt.td_crank_nicolson_src_mem)
      
      select case(tr%ob%mem_type)
      case(SAVE_CPU_TIME)
        call cn_src_mem_dt(tr%ob, st, hm, gr, dt, time, nt)
      case(SAVE_RAM_USAGE)
        call cn_src_mem_sp_dt(tr%ob, st, hm, gr, max_iter, dt, time, nt)
      end select

      POP_SUB(propagator_dt.td_crank_nicolson_src_mem)
    end subroutine td_crank_nicolson_src_mem

  end subroutine propagator_dt
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> operators for Crank-Nicolson scheme
  subroutine propagator_qmr_op(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)
    integer :: idim
    CMPLX, allocatable :: zpsi(:, :)

    PUSH_SUB(propagator_qmr_op)

    SAFE_ALLOCATE(zpsi(1:grid_p%mesh%np_part, 1:dim_op))
    zpsi = M_z0
    forall(idim = 1:dim_op)
      zpsi(1:grid_p%mesh%np, idim) = x((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np)
    end forall

    call exponential_apply(tr_p%te, grid_p%der, hm_p, zpsi, ist_op, ik_op, -dt_op/M_TWO, t_op)

    forall(idim = 1:dim_op)
      y((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) = zpsi(1:grid_p%mesh%np, idim)
    end forall

    SAFE_DEALLOCATE_A(zpsi)
    POP_SUB(propagator_qmr_op)
  end subroutine propagator_qmr_op
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_qmr_prec(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    PUSH_SUB(propagator_qmr_prec)
    y = x

    POP_SUB(propagator_qmr_prec)
  end subroutine propagator_qmr_prec
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> operators for Crank-Nicolson scheme
  subroutine td_zop(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:)
    FLOAT, intent(in)  :: xim(:)
    FLOAT, intent(out) :: yre(:)
    FLOAT, intent(out) :: yim(:)
    integer :: idim
    CMPLX, allocatable :: zpsi(:, :)

    PUSH_SUB(td_zop)

    SAFE_ALLOCATE(zpsi(1:grid_p%mesh%np_part, 1:dim_op))
    zpsi = M_z0
    forall(idim = 1:dim_op)
      zpsi(1:grid_p%mesh%np, idim) = xre((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) + &
                                     M_zI * xim((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np)
    end forall

    call exponential_apply(tr_p%te, grid_p%der, hm_p, zpsi, ist_op, ik_op, -dt_op/M_TWO, t_op)

    forall(idim = 1:dim_op)
      yre((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) = real(zpsi(1:grid_p%mesh%np, idim))
      yim((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) = aimag(zpsi(1:grid_p%mesh%np, idim))
    end forall

    SAFE_DEALLOCATE_A(zpsi)

    POP_SUB(td_zop)
  end subroutine td_zop
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Transpose of H (called e.g. by bi-conjugate gradient solver)
  subroutine td_zopt(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:)
    FLOAT, intent(in)  :: xim(:)
    FLOAT, intent(out) :: yre(:)
    FLOAT, intent(out) :: yim(:)
    integer :: idim
    CMPLX, allocatable :: zpsi(:, :)

    PUSH_SUB(td_zopt)

    ! To act with the transpose of H on the wfn we apply H to the conjugate of psi
    ! and conjugate the resulting hpsi (note that H is not a purely real operator
    ! for scattering wavefunctions anymore).
    SAFE_ALLOCATE(zpsi(1:grid_p%mesh%np_part, 1:dim_op))
    zpsi = M_z0
    forall(idim = 1:dim_op)
      zpsi(1:grid_p%mesh%np, idim) = xre((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) - &
                                     M_zI * xim((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np)
    end forall

    call exponential_apply(tr_p%te, grid_p%der, hm_p, zpsi, ist_op, ik_op, -dt_op/M_TWO, t_op)

    forall(idim = 1:dim_op)
      yre((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) =    real(zpsi(1:grid_p%mesh%np, idim))
      yim((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) = - aimag(zpsi(1:grid_p%mesh%np, idim))
    end forall

    SAFE_DEALLOCATE_A(zpsi)

    POP_SUB(td_zopt)
  end subroutine td_zopt
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  logical pure function propagator_ions_are_propagated(tr) result(propagated)
    type(propagator_t), intent(in) :: tr

    select case(tr%method)
    case(PROP_ETRS, PROP_AETRS, PROP_CAETRS)
      propagated = .true.
    case default
      propagated = .false.
    end select

  end function propagator_ions_are_propagated

  ! ---------------------------------------------------------

  logical pure function propagator_dens_is_propagated(tr) result(propagated)
    type(propagator_t), intent(in) :: tr

    select case(tr%method)
    case(PROP_AETRS, PROP_CAETRS)
      propagated = .true.
    case default
      propagated = .false.
    end select

  end function propagator_dens_is_propagated

  ! ---------------------------------------------------------

  logical pure function propagator_requires_vks(tr) result(requires)
    type(propagator_t), intent(in) :: tr
    
    select case(tr%method)
    case(PROP_CAETRS)
      requires = .false.
    case default
      requires = .true.
    end select
    
  end function propagator_requires_vks

  ! ---------------------------------------------------------


  subroutine propagator_dt_cpmd(cp_propagator, gr, ks, st, hm, gauge_force, geo, iter, dt, ions, scsteps, update_energy)
    type(cpmd_t), intent(inout)        :: cp_propagator
    type(grid_t), intent(inout)        :: gr
    type(v_ks_t), intent(inout)        :: ks
    type(states_t), intent(inout)      :: st
    type(hamiltonian_t), intent(inout) :: hm
    type(gauge_force_t), intent(inout) :: gauge_force
    type(geometry_t), intent(inout)    :: geo
    integer, intent(in)                :: iter
    FLOAT, intent(in)                  :: dt
    type(ion_dynamics_t), intent(inout) :: ions
    integer, intent(inout)             :: scsteps
    logical, intent(in)                :: update_energy

    logical :: cmplxscl, generate
    PUSH_SUB(propagator_dt_cpmd)

    cmplxscl = hm%cmplxscl%space

    if(states_are_real(st)) then
      call dcpmd_propagate(cp_propagator, gr, hm, st, iter, dt)
    else
      call zcpmd_propagate(cp_propagator, gr, hm, st, iter, dt)
    end if
    scsteps = 1

    if(.not. cmplxscl) then
        call density_calc(st, gr, st%rho)
    else
        call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
    end if  

    call ion_dynamics_propagate(ions, gr%sb, geo, iter*dt, dt)
    generate = .true.

    if(gauge_field_is_applied(hm%ep%gfield)) then
       call gauge_field_propagate(hm%ep%gfield, gauge_force, dt)
    end if

    call hamiltonian_epot_generate(hm, gr, geo, st, time = iter*dt)

    ! update Hamiltonian and eigenvalues (fermi is *not* called)
    call v_ks_calc(ks, hm, st, geo, calc_eigenval = update_energy, time = iter*dt, calc_energy = update_energy)
    ! Get the energies.
    if(update_energy) call energy_calc_total(hm, gr, st, iunit = -1)

    if(states_are_real(st)) then
      call dcpmd_propagate_vel(cp_propagator, gr, hm, st, dt)
    else
      call zcpmd_propagate_vel(cp_propagator, gr, hm, st, dt)
    end if

    ! Recalculate forces, update velocities...
    call forces_calculate(gr, geo, hm, st, iter*dt, dt)
    call ion_dynamics_propagate_vel(ions, geo, atoms_moved = generate)
    call hamiltonian_epot_generate(hm, gr, geo, st, time = iter*dt)
    geo%kinetic_energy = ion_dynamics_kinetic_energy(geo)

    if(gauge_field_is_applied(hm%ep%gfield)) then
      call gauge_field_get_force(gr, geo, hm%ep%proj, hm%phase, st, gauge_force)
      call gauge_field_propagate_vel(hm%ep%gfield, gauge_force, dt)
    end if

    POP_SUB(propagator_dt_cpmd)
  end subroutine propagator_dt_cpmd


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


#include "propagator_qoct_inc.F90"

end module propagator_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

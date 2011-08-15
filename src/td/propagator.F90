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

module propagator_m
  use batch_m
  use blas_m
  use cube_function_m
  use datasets_m
  use density_m
  use exponential_m
  use gauge_field_m
  use grid_m
  use geometry_m
  use global_m
  use hamiltonian_m
  use ion_dynamics_m
  use lalg_basic_m
  use lasers_m
  use loct_m
  use parser_m
  use math_m
  use mesh_function_m
  use messages_m
  use ob_mem_m
  use ob_propagator_m
  use ob_terms_m
  use opencl_m
  use profiling_m
  use states_dim_m
  use solvers_m
  use sparskit_m
  use states_m
  use types_m
  use v_ks_m
  use varinfo_m

  implicit none

  private
  public ::                   &
    propagator_t,                 &
    propagator_init,              &
    propagator_end,               &
    propagator_copy,              &
    propagator_run_zero_iter,     &
    propagator_dt,                &
    td_zop,                   &
    td_zopt,                  &
    propagator_qmr_op,            &
    propagator_qmr2_op,           &
    propagator_qmr_prec,          &
    propagator_set_scf_prop,      &
    propagator_remove_scf_prop,   &
    propagator_ions_are_propagated, &
    propagator_dens_is_propagated,  &
    propagator_requires_vks

  integer, public, parameter ::        &
    PROP_ETRS                    = 2,  &
    PROP_AETRS                   = 3,  &
    PROP_EXPONENTIAL_MIDPOINT    = 4,  &
    PROP_CRANK_NICHOLSON         = 5,  &
    PROP_CRANK_NICHOLSON_SPARSKIT= 6,  &
    PROP_MAGNUS                  = 7,  &
    PROP_CRANK_NICHOLSON_SRC_MEM = 8,  &
    PROP_QOCT_TDDFT_PROPAGATOR   = 10, &
    PROP_QOCT_TDDFT_PROPAGATOR_2 = 11, &
    PROP_CAETRS                  = 12

  FLOAT, parameter :: scf_threshold = CNST(1.0e-3)

  type propagator_t
    integer             :: method           ! Which evolution method to use.
    type(exponential_t) :: te               ! How to apply the propagator (e^{-i H \Delta t}).
    FLOAT, pointer      :: v_old(:, :, :) => null()
                                            ! Storage of the KS potential of previous iterations.
    FLOAT, pointer      :: vmagnus(:, :, :) => null() 
                                            ! Auxiliary function to store the Magnus potentials.
    type(ob_terms_t)    :: ob               ! For open boundaries: leads, memory
    integer             :: scf_propagation_steps 
    logical             :: first
  end type propagator_t

#ifdef HAVE_SPARSKIT
  type(sparskit_solver_t), pointer, private :: tdsk
#endif
  type(grid_t),            pointer, private :: grid_p
  type(hamiltonian_t),     pointer, private :: hm_p
  type(propagator_t),      pointer, private :: tr_p
  CMPLX, allocatable,      private :: zpsi_tmp(:,:,:,:)
  integer,                 private :: ik_op, ist_op, idim_op, dim_op, nst_op
  type(states_t),          private :: st_op
  FLOAT,                   private :: t_op, dt_op

contains


  ! ---------------------------------------------------------
  subroutine propagator_copy(tro, tri)
    type(propagator_t), intent(inout) :: tro
    type(propagator_t), intent(in)    :: tri
    PUSH_SUB(tr_rti_copy)

    tro%method = tri%method

    select case(tro%method)
    case(PROP_MAGNUS)
      call loct_pointer_copy(tro%vmagnus, tri%vmagnus)
    case(PROP_CRANK_NICHOLSON_SRC_MEM)
      message(1) = 'Internal error at propagator_copy.'
      call messages_fatal(1)
    end select

    call loct_pointer_copy(tro%v_old, tri%v_old)
    call exponential_copy(tro%te, tri%te)
    tro%scf_propagation_steps = tri%scf_propagation_steps

    POP_SUB(tr_rti_copy)
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
    logical,             intent(in)    :: have_fields ! whether there is an associated "field"
                                                      ! that must be propagated (currently ions
                                                      ! or a gauge field).

    integer :: default_propagator

    PUSH_SUB(propagator_init)

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
    !% Classical Crank-Nicholson propagator.
    !%
    !% <MATH>
    !%  (1 + i\delta t/2 H_{n+1/2}) \psi_{n+1} = (1 - i\delta t/2 H_{n+1/2}) \psi_{n}  
    !% </MATH>
    !%Option crank_nicholson_sparskit 6
    !% Classical Crank-Nicholson propagator. Requires the SPARSKIT library.
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
    !% Crank-Nicholson propagator with source and memory term for transport
    !% calculations.
    !%Option qoct_tddft_propagator 10
    !% WARNING: EXPERIMENTAL
    !%Option qoct_tddft_propagator_2 11
    !% WARNING: EXPERIMENTAL
    !%End
    call messages_obsolete_variable('TDEvolutionMethod', 'TDPropagator')

    default_propagator = PROP_ETRS
    if(gr%ob_grid%open_boundaries) default_propagator = PROP_CRANK_NICHOLSON_SRC_MEM

    call parse_integer(datasets_check('TDPropagator'), default_propagator, tr%method)
    if(.not.varinfo_valid_option('TDPropagator', tr%method)) call input_error('TDPropagator')

    if(gr%ob_grid%open_boundaries.and.tr%method.ne.PROP_CRANK_NICHOLSON_SRC_MEM) then
      message(1) = 'The time-evolution method for time-dependent run cannot'
      message(2) = 'be chosen freely. The Crank-Nicholson propagator'
      message(3) = 'with source and memory term has to be used. Either set'
      message(4) = ''
      message(5) = '  TDPropagator = crank_nicholson_src_mem'
      message(6) = ''
      message(7) = 'in your input or remove the TDPropagator variable.'
      call messages_fatal(7)
    end if

    select case(tr%method)
    case(PROP_ETRS)
    case(PROP_AETRS, PROP_CAETRS)
    case(PROP_EXPONENTIAL_MIDPOINT)
    case(PROP_CRANK_NICHOLSON)
    case(PROP_CRANK_NICHOLSON_SPARSKIT)
#ifdef HAVE_SPARSKIT
      SAFE_ALLOCATE(tdsk)
      call zsparskit_solver_init(gr%mesh%np, tdsk)
      SAFE_ALLOCATE(zpsi_tmp(1:gr%mesh%np_part, 1:st%d%dim, 1:st%nst, st%d%kpt%start:st%d%kpt%end))
#else
      message(1) = 'Octopus was not compiled with support for the SPARSKIT library. This'
      message(2) = 'library is required if the "crank_nicholson_sparskit" propagator is selected.'
      message(3) = 'Try using a different propagation scheme or recompile with SPARSKIT support.'
      call messages_fatal(3)
#endif
    case(PROP_MAGNUS)
      SAFE_ALLOCATE(tr%vmagnus(1:gr%mesh%np, 1:st%d%nspin, 1:2))
    case(PROP_CRANK_NICHOLSON_SRC_MEM)
      call ob_propagator_init(st, gr, hm, tr%ob, dt, max_iter)
    case(PROP_QOCT_TDDFT_PROPAGATOR)
    case(PROP_QOCT_TDDFT_PROPAGATOR_2)
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
    call exponential_init(tr%te, gr%der) ! initialize propagator

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
  subroutine propagator_set_scf_prop(tr)
    type(propagator_t), intent(inout) :: tr
    tr%scf_propagation_steps = huge(1)
  end subroutine propagator_set_scf_prop
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_remove_scf_prop(tr)
    type(propagator_t), intent(inout) :: tr
    tr%scf_propagation_steps = -1
  end subroutine propagator_remove_scf_prop
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine propagator_end(tr)
    type(propagator_t), intent(inout) :: tr

    PUSH_SUB(propagator_end)

    ! sanity check
    ASSERT(associated(tr%v_old)) 
    SAFE_DEALLOCATE_P(tr%v_old)         ! clean old KS potentials

    select case(tr%method)
    case(PROP_MAGNUS)
      ASSERT(associated(tr%vmagnus))
      SAFE_DEALLOCATE_P(tr%vmagnus)
    case(PROP_CRANK_NICHOLSON_SPARSKIT)
#ifdef HAVE_SPARSKIT
      call zsparskit_solver_end()
      SAFE_DEALLOCATE_A(zpsi_tmp)
#endif
    case(PROP_CRANK_NICHOLSON_SRC_MEM)
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

    POP_SUB(propagator_run_zero_iter)
  end subroutine propagator_run_zero_iter


  ! ---------------------------------------------------------
  ! Propagates st from time - dt to t.
  ! If dt<0, it propagates *backwards* from t+|dt| to t
  ! ---------------------------------------------------------
  subroutine propagator_dt(ks, hm, gr, st, tr, time, dt, mu, max_iter, nt, gauge_force, ions, geo, scsteps)
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
    type(gauge_force_t),  optional,  intent(inout) :: gauge_force
    type(ion_dynamics_t), optional,  intent(inout) :: ions
    type(geometry_t),     optional,  intent(inout) :: geo
    integer,              optional,  intent(out)   :: scsteps

    integer :: is, iter, ik, ist, idim
    FLOAT   :: d, d_max
    logical :: self_consistent
    CMPLX, allocatable :: zpsi1(:, :, :, :)
    FLOAT, allocatable :: dtmp(:), vaux(:, :), vold(:, :)
    type(profile_t), save :: prof

    call profiling_in(prof, "TD_PROPAGATOR")
    PUSH_SUB(propagator_dt)

    if(present(ions)) then
      ASSERT(present(geo))
    end if

    if(gauge_field_is_applied(hm%ep%gfield)) then
      ASSERT(present(gauge_force))
    end if

    self_consistent = .false.
    if(hm%theory_level .ne. INDEPENDENT_PARTICLES .and. tr%method /= PROP_CAETRS) then
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

    if(.not. propagator_requires_vks(tr)) then
      SAFE_ALLOCATE(vold(1:gr%mesh%np, 1:st%d%nspin))
      call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 1), vold)
    else
      call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 2), tr%v_old(:, :, 3))
      call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 1), tr%v_old(:, :, 2))
      call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc(:, :),     tr%v_old(:, :, 1))
      call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), tr%v_old(:, :, 1:3), time, tr%v_old(:, :, 0))
    end if

    select case(tr%method)
    case(PROP_ETRS);                     call td_etrs
    case(PROP_AETRS, PROP_CAETRS);       call td_aetrs
    case(PROP_EXPONENTIAL_MIDPOINT);     call exponential_midpoint
    case(PROP_CRANK_NICHOLSON);          call td_crank_nicholson
    case(PROP_CRANK_NICHOLSON_SPARSKIT); call td_crank_nicholson_sparskit
    case(PROP_MAGNUS);                   call td_magnus
    case(PROP_CRANK_NICHOLSON_SRC_MEM);  call td_crank_nicholson_src_mem
    case(PROP_QOCT_TDDFT_PROPAGATOR)
      call td_qoct_tddft_propagator(hm, gr, st, tr, time, dt)
    case(PROP_QOCT_TDDFT_PROPAGATOR_2)
      call td_qoct_tddft_propagator_2(hm, gr, st, tr, time, dt)
    end select

    if(present(scsteps)) scsteps = 1

    if(self_consistent) then
      SAFE_ALLOCATE(vaux(1:gr%mesh%np, 1:st%d%nspin))

      ! First, compare the new potential to the extrapolated one.
      call density_calc(st, gr, st%rho)
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
          select case(tr%method)
          case(PROP_ETRS);                     call td_etrs
          case(PROP_AETRS, PROP_CAETRS);       call td_aetrs
          case(PROP_EXPONENTIAL_MIDPOINT);     call exponential_midpoint
          case(PROP_CRANK_NICHOLSON);          call td_crank_nicholson
          case(PROP_CRANK_NICHOLSON_SPARSKIT); call td_crank_nicholson_sparskit
          case(PROP_MAGNUS);                   call td_magnus
          case(PROP_CRANK_NICHOLSON_SRC_MEM);  call td_crank_nicholson_src_mem
          case(PROP_QOCT_TDDFT_PROPAGATOR)
            call td_qoct_tddft_propagator(hm, gr, st, tr, time, dt)
          case(PROP_QOCT_TDDFT_PROPAGATOR_2)
            call td_qoct_tddft_propagator_2(hm, gr, st, tr, time, dt)
          end select

          call density_calc(st, gr, st%rho)
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
      SAFE_DEALLOCATE_A(vaux)
    end if
    
    SAFE_DEALLOCATE_A(vold)
    
    POP_SUB(propagator_dt)
    call profiling_out(prof)

  contains

    ! ---------------------------------------------------------
    ! Propagator with enforced time-reversal symmetry
    subroutine td_etrs
      FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
      integer :: ik, ist, ib
      type(batch_t) :: zpsib_save
      type(density_calc_t) :: dens_calc

      PUSH_SUB(propagator_dt.td_etrs)

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then

        SAFE_ALLOCATE(vhxc_t1(1:gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE(vhxc_t2(1:gr%mesh%np, 1:st%d%nspin))
        call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t1)

        call density_calc_init(dens_calc, st, gr, st%rho)

        do ik = st%d%kpt%start, st%d%kpt%end
          do ib = st%block_start, st%block_end

            !save the state
            call batch_copy(st%psib(ib, ik), zpsib_save, reference = .false.)
            if(batch_is_packed(st%psib(ib, ik))) call batch_pack(zpsib_save, copy = .false.)
            call batch_copy_data(gr%der%mesh%np, st%psib(ib, ik), zpsib_save)

            !propagate the state dt with H(time - dt)
            call exponential_apply_batch(tr%te, gr%der, hm, st%psib(ib, ik), ik, dt/mu, time - dt)
            call density_calc_accumulate(dens_calc, ik, st%psib(ib, ik))

            !restore the saved state
            call batch_copy_data(gr%der%mesh%np, zpsib_save, st%psib(ib, ik))

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
        do ib = st%block_start, st%block_end
          call exponential_apply_batch(tr%te, gr%der, hm, st%psib(ib, ik), ik, dt/(mu*M_TWO), time - dt)
        end do
      end do

      ! propagate dt/2 with H(t)

      ! first move the ions to time t
      if(present(ions)) then
        call ion_dynamics_propagate(ions, gr%sb, geo, time, dt)
        call hamiltonian_epot_generate(hm, gr, geo, st, time = time)
      end if

      if(gauge_field_is_applied(hm%ep%gfield)) then
        call gauge_field_propagate(hm%ep%gfield, gauge_force, dt)
      end if

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t2, hm%vhxc)
      end if
      call hamiltonian_update(hm, gr%mesh, time = time)

      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%block_start, st%block_end
          call exponential_apply_batch(tr%te, gr%der, hm, st%psib(ib, ik), ik, dt/(M_TWO*mu), time)
        end do
      end do

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        SAFE_DEALLOCATE_A(vhxc_t1)
        SAFE_DEALLOCATE_A(vhxc_t2)
      end if

      POP_SUB(propagator_dt.td_etrs)
    end subroutine td_etrs


    ! ---------------------------------------------------------
    ! Propagator with approximate enforced time-reversal symmetry
    subroutine td_aetrs
      integer :: ik, ispin, ip, ist, ib
      type(batch_t) :: zpsib
      FLOAT :: vv
      CMPLX :: phase
      type(density_calc_t)  :: dens_calc
      type(profile_t), save :: phase_prof
      integer               :: pnp, iprange
      type(opencl_mem_t)    :: phase_buff

      PUSH_SUB(propagator_dt.td_aetrs)

      if(tr%method == PROP_CAETRS) then
        call lalg_copy(gr%mesh%np, st%d%nspin, vold, hm%vhxc)
        call hamiltonian_update(hm, gr%mesh, time = time - dt)
        call v_ks_calc_start(ks, hm, st, geo, time = time - dt, calc_energy = .false.)
      end if

      ! propagate half of the time step with H(time - dt)
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%block_start, st%block_end
          call exponential_apply_batch(tr%te, gr%der, hm, st%psib(ib, ik), ik, dt/(M_TWO*mu), time - dt)
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
          pnp = opencl_padded_size(gr%mesh%np)
          call opencl_create_buffer(phase_buff, CL_MEM_READ_ONLY, TYPE_FLOAT, pnp*st%d%nspin)
          ASSERT(ubound(vold, dim = 1) == gr%mesh%np)
          do ispin = 1, st%d%nspin
            call opencl_write_buffer(phase_buff, gr%mesh%np, vold(:, ispin), offset = (ispin - 1)*pnp)
          end do
        end if

      end if

      ! interpolate the Hamiltonian to time t
      call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 0), hm%vhxc)

      ! move the ions to time t
      if(present(ions)) then      
        call ion_dynamics_propagate(ions, gr%sb, geo, time, dt)
        call hamiltonian_epot_generate(hm, gr, geo, st, time = time)
      end if

      if(gauge_field_is_applied(hm%ep%gfield)) then
        call gauge_field_propagate(hm%ep%gfield, gauge_force, dt)
      end if

      call hamiltonian_update(hm, gr%mesh, time = time)

      call density_calc_init(dens_calc, st, gr, st%rho, packed = hamiltonian_apply_packed(hm, gr%mesh))

      ! propagate the other half with H(t)
      do ik = st%d%kpt%start, st%d%kpt%end
        ispin = states_dim_get_spin_index(st%d, ik)

        do ib = st%block_start, st%block_end
          if(hamiltonian_apply_packed(hm, gr%mesh)) call batch_pack(st%psib(ib, ik))
          
          if(tr%method == PROP_CAETRS) then
            call profiling_in(phase_prof, "CAETRS_PHASE")
            select case(batch_status(st%psib(ib, ik)))
            case(BATCH_NOT_PACKED)
              do ip = 1, gr%mesh%np
                vv = vold(ip, ispin)
                phase = cmplx(cos(vv), -sin(vv), kind = REAL_PRECISION)
                forall(ist = 1:st%psib(ib, ik)%nst_linear)
                  st%psib(ib, ik)%states_linear(ist)%zpsi(ip) = st%psib(ib, ik)%states_linear(ist)%zpsi(ip)*phase
                end forall
              end do
            case(BATCH_PACKED)
              do ip = 1, gr%mesh%np
                vv = vold(ip, ispin)
                phase = cmplx(cos(vv), -sin(vv), kind = REAL_PRECISION)
                forall(ist = 1:st%psib(ib, ik)%nst_linear)
                  st%psib(ib, ik)%pack%zpsi(ist, ip) = st%psib(ib, ik)%pack%zpsi(ist, ip)*phase
                end forall
              end do
            case(BATCH_CL_PACKED)
              call opencl_set_kernel_arg(kernel_phase, 0, pnp*(ispin - 1))
              call opencl_set_kernel_arg(kernel_phase, 1, phase_buff)
              call opencl_set_kernel_arg(kernel_phase, 2, st%psib(ib, ik)%pack%buffer)
              call opencl_set_kernel_arg(kernel_phase, 3, log2(st%psib(ib, ik)%pack%size(1)))

              iprange = opencl_max_workgroup_size()/st%psib(ib, ik)%pack%size(1)

              call opencl_kernel_run(kernel_phase, (/st%psib(ib, ik)%pack%size(1), pnp/), &
                (/st%psib(ib, ik)%pack%size(1), iprange/))

            end select
            call profiling_out(phase_prof)
          end if

          call exponential_apply_batch(tr%te, gr%der, hm, st%psib(ib, ik), ik, dt/(M_TWO*mu), time)
          call density_calc_accumulate(dens_calc, ik, st%psib(ib, ik))

          if(hamiltonian_apply_packed(hm, gr%mesh)) call batch_unpack(st%psib(ib, ik))
        end do
      end do

      if(tr%method == PROP_CAETRS .and. opencl_is_enabled() .and. hamiltonian_apply_packed(hm, gr%mesh)) then
        call opencl_release_buffer(phase_buff)
      end if
      
      call density_calc_end(dens_calc)

      POP_SUB(propagator_dt.td_aetrs)
    end subroutine td_aetrs


    ! ---------------------------------------------------------
    ! Exponential midpoint
    subroutine exponential_midpoint
      integer :: ib, ik
      type(ion_state_t) :: ions_state
      FLOAT :: vecpot(1:MAX_DIM), vecpot_vel(1:MAX_DIM)

      PUSH_SUB(propagator_dt.exponential_midpoint)

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        call interpolate( (/time, time - dt, time - M_TWO*dt/), tr%v_old(:, :, 0:2), time - dt/M_TWO, hm%vhxc(:, :))
      end if

      !move the ions to time 'time - dt/2'
      if(present(ions)) then
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
        do ib = st%block_start, st%block_end
          call exponential_apply_batch(tr%te, gr%der, hm, st%psib(ib, ik), ik, dt/mu, time - dt/M_TWO)
        end do
      end do

      !restore to time 'time - dt'
      if(present(ions)) call ion_dynamics_restore_state(ions, geo, ions_state)

      if(gauge_field_is_applied(hm%ep%gfield)) then
        call gauge_field_set_vec_pot(hm%ep%gfield, vecpot)
        call gauge_field_set_vec_pot_vel(hm%ep%gfield, vecpot_vel)
        call hamiltonian_update(hm, gr%mesh)
      end if

      POP_SUB(propagator_dt.exponential_midpoint)
    end subroutine exponential_midpoint


    ! ---------------------------------------------------------
    ! Crank-Nicholson propagator, linear solver from SPARSKIT.
    subroutine td_crank_nicholson_sparskit
#ifdef HAVE_SPARSKIT
      FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
      CMPLX, allocatable :: zpsi_rhs_pred(:,:,:,:), zpsi_rhs_corr(:,:,:,:)
      integer :: ik, ist, idim, np_part

      PUSH_SUB(propagator_dt.td_crank_nicholson_sparskit)

      np_part = gr%mesh%np_part
      SAFE_ALLOCATE(zpsi_rhs_corr(1:np_part, 1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
      zpsi_rhs_corr = st%zpsi ! store zpsi for corrector step

      ! define pointer and variables for usage in td_zop, td_zopt routines
      grid_p    => gr
      hm_p      => hm
      tr_p      => tr
      dt_op = dt
      t_op  = time

      ! we (ab)use exponential_apply to compute (1-i\delta t/2 H_n)\psi^n
      ! exponential order needs to be only 1
      tr%te%exp_method = EXP_TAYLOR
      tr%te%exp_order  = 1

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        np_part = gr%mesh%np_part
        SAFE_ALLOCATE(zpsi_rhs_pred(1:np_part, 1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
        zpsi_rhs_pred = st%zpsi ! store zpsi for predictor step
        
        SAFE_ALLOCATE(vhxc_t1(1:gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE(vhxc_t2(1:gr%mesh%np, 1:st%d%nspin))
        vhxc_t1 = hm%vhxc

        ! get rhs of CN linear system (rhs1 = (1-i\delta t/2 H_n)\psi^n)
        do ik = st%d%kpt%start, st%d%kpt%end
          do ist = st%st_start, st%st_end
            call exponential_apply(tr%te, gr%der, hm, zpsi_rhs_pred(:, :, ist, ik), ist, ik, dt/M_TWO, time - dt)
            if(hamiltonian_inh_term(hm)) then
              zpsi_rhs_pred(:, :, ist, ik) = zpsi_rhs_pred(:, :, ist, ik) + dt * hm%inh_st%zpsi(:, :, ist, ik)
            end if
          end do
        end do

        ! predictor step: 
        ! solve (1+i\delta t/2 H_n)\psi^{predictor}_{n+1} = (1-i\delta t/2 H_n)\psi^n
        do idim = 1, st%d%dim
          do ik = st%d%kpt%start, st%d%kpt%end
            do ist = st%st_start, st%st_end
              idim_op = idim
              ist_op = ist
              ik_op = ik
              call zsparskit_solver_run(tdsk, td_zop, td_zopt, &
                st%zpsi(1:gr%mesh%np, idim, ist, ik), zpsi_rhs_pred(1:gr%mesh%np, idim, ist, ik))
            end do
          end do
        end do

        call density_calc(st, gr, st%rho)
        call v_ks_calc(ks, hm, st, geo)

        vhxc_t2 = hm%vhxc
        ! compute potential at n+1/2 as average
        hm%vhxc = (vhxc_t1 + vhxc_t2)/M_TWO
        call hamiltonian_update(hm, gr%mesh)
      end if

      ! get rhs of CN linear system (rhs2 = (1-i\delta t H_{n+1/2})\psi^n)
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          call exponential_apply(tr%te, gr%der, hm, zpsi_rhs_corr(:, :, ist, ik), ist, ik, dt/M_TWO, time - dt)
          if(hamiltonian_inh_term(hm)) then
            zpsi_rhs_corr(:, :, ist, ik) = zpsi_rhs_corr(:, :, ist, ik) + dt * hm%inh_st%zpsi(:, :, ist, ik)
          end if
        end do
      end do

      ! corrector step: 
      ! solve (1+i\delta t/2 H_{n+1/2})\psi_{n+1} = (1-i\delta t/2 H_{n+1/2})\psi^n
      do idim = 1, st%d%dim
        do ik = st%d%kpt%start, st%d%kpt%end
          do ist = st%st_start, st%st_end
            idim_op = idim
            ist_op = ist
            ik_op = ik
            call zsparskit_solver_run(tdsk, td_zop, td_zopt, &
              st%zpsi(1:gr%mesh%np, idim, ist, ik), zpsi_rhs_corr(1:gr%mesh%np, idim, ist, ik))
          end do
        end do
      end do
      
      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        SAFE_DEALLOCATE_A(vhxc_t1)
        SAFE_DEALLOCATE_A(vhxc_t2)
        SAFE_DEALLOCATE_A(zpsi_rhs_pred)
      end if
      SAFE_DEALLOCATE_A(zpsi_rhs_corr)

      POP_SUB(propagator_dt.td_crank_nicholson_sparskit)
#endif
    end subroutine td_crank_nicholson_sparskit
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    ! Crank-Nicholson propagator, QMR linear solver.
    subroutine td_crank_nicholson
      CMPLX, allocatable :: zpsi_rhs(:,:), zpsi(:), rhs(:), inhpsi(:)
      integer :: ik, ist, idim, ip, isize, np_part, np, iter
      FLOAT :: dres
      FLOAT :: cgtol = CNST(1.0e-12)
      logical :: converged

      PUSH_SUB(propagator_dt.td_crank_nicholson)

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

          ! put the values is a continuous array
          do idim = 1, st%d%dim
            call states_get_state(st, gr%mesh, idim, ist, ik, zpsi((idim - 1)*np+1:idim*np))
            rhs((idim - 1)*np + 1:idim*np) = zpsi_rhs(1:np, idim)
          end do

          ist_op = ist
          ik_op = ik
          iter = 2000
          call zqmr_sym(np*st%d%dim, zpsi, rhs, propagator_qmr_op, propagator_qmr_prec, iter, dres, cgtol, &
            showprogress = .false., converged = converged)

          do idim = 1, st%d%dim
            call states_set_state(st, gr%mesh, idim, ist, ik, zpsi((idim-1)*np + 1:(idim - 1)*np + np))
          end do

          if(.not.converged) then
            write(message(1),'(a)')        'The linear solver used for the Crank-Nicholson'
            write(message(2),'(a,es14.4)') 'propagator did not converge: Residual = ', dres
            call messages_warning(2)
          end if

        end do
      end do

      SAFE_DEALLOCATE_A(zpsi_rhs)
      SAFE_DEALLOCATE_A(zpsi)
      SAFE_DEALLOCATE_A(rhs)
      POP_SUB(propagator_dt.td_crank_nicholson)
    end subroutine td_crank_nicholson
    ! ---------------------------------------------------------



    ! ---------------------------------------------------------
    ! Magnus propagator
    subroutine td_magnus
      integer :: j, is, ist, ik, i
      FLOAT :: atime(2)
      FLOAT, allocatable :: vaux(:, :, :), pot(:)

      PUSH_SUB(propagator_dt.td_magnus)

      SAFE_ALLOCATE(vaux(1:gr%mesh%np, 1:st%d%nspin, 1:2))

      atime(1) = (M_HALF-sqrt(M_THREE)/M_SIX)*dt
      atime(2) = (M_HALF+sqrt(M_THREE)/M_SIX)*dt

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
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
    ! Crank-Nicholson scheme with source and memory term.
    subroutine td_crank_nicholson_src_mem()
      PUSH_SUB(propagator_dt.td_crank_nicholson_src_mem)
      
      select case(tr%ob%mem_type)
      case(SAVE_CPU_TIME)
        call cn_src_mem_dt(tr%ob, st, ks, hm, gr, max_iter, dt, time, nt)
      case(SAVE_RAM_USAGE)
        call cn_src_mem_sp_dt(tr%ob, st, ks, hm, gr, max_iter, dt, time, nt)
      end select

      POP_SUB(propagator_dt.td_crank_nicholson_src_mem)
    end subroutine td_crank_nicholson_src_mem

  end subroutine propagator_dt
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! operators for Crank-Nicholson scheme
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
  ! operators for Crank-Nicholson scheme
  subroutine td_zop(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:)
    FLOAT, intent(in)  :: xim(:)
    FLOAT, intent(out) :: yre(:)
    FLOAT, intent(out) :: yim(:)

    PUSH_SUB(td_zop)
#ifdef HAVE_SPARSKIT    
    zpsi_tmp(1:grid_p%mesh%np, idim_op, ist_op, ik_op) = &
      xre(1:grid_p%mesh%np) + M_zI * xim(1:grid_p%mesh%np)

    ! propagate backwards
    call exponential_apply(tr_p%te, grid_p%der, hm_p, zpsi_tmp(:, :, ist_op, ik_op), ist_op, ik_op, -dt_op / M_TWO, t_op)

    yre(1:grid_p%mesh%np) =  real(zpsi_tmp(1:grid_p%mesh%np, idim_op, ist_op, ik_op))
    yim(1:grid_p%mesh%np) = aimag(zpsi_tmp(1:grid_p%mesh%np, idim_op, ist_op, ik_op))
#endif    
    POP_SUB(td_zop)
  end subroutine td_zop
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Transpose of H (called e.g. by bi-conjugate gradient solver)
  subroutine td_zopt(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:)
    FLOAT, intent(in)  :: xim(:)
    FLOAT, intent(out) :: yre(:)
    FLOAT, intent(out) :: yim(:)
    
    PUSH_SUB(td_zopt)
#ifdef HAVE_SPARSKIT        
    ! To act with the transpose of H on the wfn we apply H to the conjugate of psi
    ! and conjugate the resulting hpsi (note that H is not a purely real operator
    ! for scattering wavefunctions anymore).
    zpsi_tmp(1:grid_p%mesh%np, idim_op, ist_op, ik_op) = &
      xre(1:grid_p%mesh%np) - M_zI * xim(1:grid_p%mesh%np)
    
    ! propagate backwards
    call exponential_apply(tr_p%te, grid_p%der, hm_p, zpsi_tmp(:, :, ist_op, ik_op), ist_op, ik_op, -dt_op/M_TWO, t_op)

    yre(1:grid_p%mesh%np) =    real(zpsi_tmp(1:grid_p%mesh%np, idim_op, ist_op, ik_op))
    yim(1:grid_p%mesh%np) = - aimag(zpsi_tmp(1:grid_p%mesh%np, idim_op, ist_op, ik_op))
#endif        
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

#include "propagator_qoct_inc.F90"

end module propagator_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

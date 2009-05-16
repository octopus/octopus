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

module td_rti_m
  use loct_m
  use batch_m
  use cube_function_m
  use datasets_m
  use exponential_m
  use exponential_split_m
  use fourier_space_m
  use gauge_field_m
  use geometry_m
  use hamiltonian_m
  use ion_dynamics_m
  use lalg_basic_m
  use lasers_m
  use loct_parser_m
  use math_m
  use mesh_function_m
  use messages_m
  use ob_mem_m
  use ob_rti_m
  use ob_terms_m
  use profiling_m
  use sparskit_m
  use states_m
  use varinfo_m
  use v_ks_m
  use solvers_m

  implicit none

  private
  public ::                   &
    td_rti_t,                 &
    td_rti_init,              &
    td_rti_end,               &
    td_rti_copy,              &
    td_rti_run_zero_iter,     &
    td_rti_dt,                &
    td_zop,                   &
    td_zopt,                  &
    td_rti_qmr_op,            &
    td_rti_qmr_prec,          &
    td_rti_set_scf_prop,      &
    td_rti_remove_scf_prop,   &
    td_rti_ions_are_propagated

  integer, public, parameter ::       &
    PROP_SPLIT_OPERATOR          = 0, &
    PROP_SUZUKI_TROTTER          = 1, &
    PROP_REVERSAL                = 2, &
    PROP_APP_REVERSAL            = 3, &
    PROP_EXPONENTIAL_MIDPOINT    = 4, &
    PROP_CRANK_NICHOLSON         = 5, &
    PROP_CRANK_NICHOLSON_SPARSKIT= 6, &
    PROP_MAGNUS                  = 7, &
    PROP_CRANK_NICHOLSON_SRC_MEM = 8, &
    PROP_VISSCHER                = 9, &
    PROP_QOCT_TDDFT_PROPAGATOR   = 10

  FLOAT, parameter :: scf_threshold = CNST(1.0e-3)

  type td_rti_t
    integer             :: method           ! Which evolution method to use.
    type(exponential_t) :: te               ! How to apply the propagator (e^{-i H \Delta t}).
    FLOAT, pointer      :: v_old(:, :, :) => null()
                                            ! Storage of the KS potential of previous iterations.
    FLOAT, pointer      :: vmagnus(:, :, :) => null() 
                                            ! Auxiliary function to store the Magnus potentials.
    type(zcf_t)         :: cf               ! Auxiliary cube for split operator methods.
    type(ob_terms_t)    :: ob               ! For open boundaries: leads, memory
    logical             :: scf_propagation
    FLOAT, pointer      :: prev_psi(:, :, :, :) => null()
    logical             :: first
  end type td_rti_t

#ifdef HAVE_SPARSKIT
  type(sparskit_solver_t), pointer, private :: tdsk
#endif
  type(grid_t),            pointer, private :: grid_p
  type(hamiltonian_t),     pointer, private :: hm_p
  type(td_rti_t),          pointer, private :: tr_p
  CMPLX, allocatable,      private :: zpsi_tmp(:,:,:,:)
  integer,                 private :: ik_op, ist_op, idim_op, dim_op
  FLOAT,                   private :: t_op, dt_op

contains


  ! ---------------------------------------------------------
  subroutine td_rti_copy(tro, tri)
    type(td_rti_t), intent(inout) :: tro
    type(td_rti_t), intent(in)    :: tri
    call push_sub('tr_rti.tr_rti_copy')

    tro%method = tri%method

    select case(tro%method)
    case(PROP_SPLIT_OPERATOR)
      call zcf_new_from(tro%cf, tri%cf)
    case(PROP_SUZUKI_TROTTER)
      call zcf_new_from(tro%cf, tri%cf)
    case(PROP_MAGNUS)
      call loct_pointer_copy(tro%vmagnus, tri%vmagnus)
    case(PROP_CRANK_NICHOLSON_SRC_MEM)
      message(1) = 'Internal error at td_rti_copy'
      call write_fatal(1)
    end select

    call loct_pointer_copy(tro%v_old, tri%v_old)
    call exponential_copy(tro%te, tri%te)
    tro%scf_propagation = tri%scf_propagation
    call loct_pointer_copy(tro%prev_psi, tri%prev_psi)

    call pop_sub()
  end subroutine td_rti_copy
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine td_rti_init(gr, st, hm, tr, dt, max_iter, have_fields)
    type(grid_t),   intent(in)      :: gr
    type(states_t), intent(in)      :: st
    type(hamiltonian_t), intent(in) :: hm
    type(td_rti_t), intent(inout)   :: tr
    FLOAT,          intent(in)      :: dt
    integer,        intent(in)      :: max_iter
    logical,        intent(in)      :: have_fields ! whether there is an associated "field"
                                                   ! that must be propagated (currently ions
                                                   ! or a gauge field).

    integer :: default_propagator

    call push_sub('td_rti.td_rti_init')

    !%Variable TDEvolutionMethod
    !%Type integer
    !%Default etrs
    !%Section Time Dependent::Propagation
    !%Description
    !% This variable determines which algorithm will be used to approximate
    !% the evolution operator <math>U(t+\delta t, t)</math>. That is, known
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
    !% not change the solution wave-functions, the cycle is stopped. In practice,
    !% in <tt>octopus</tt> we perform a second-order extrapolation without
    !% self-consistente check, except for the first two iterations, where obviously
    !% the extrapolation is not reliable.
    !%
    !% The proliferation of methods is certainly excessive; The reason for it is that 
    !% the propagation algorithm is currently a topic of active development. We
    !% hope that in the future the optimal schemes are clearly identified. In the
    !% mean time, if you do not feel like testing, use the default choices and
    !% make sure the time step is small enough.
    !%Option split 0
    !% Split Operator (SO).
    !% This is one of the most traditional methods. It splits the full Hamiltonian
    !% into a kinetic and a potential part, performing the first in Fourier-space,
    !% and the latter in real space. The necessary transformations are performed
    !% with the FFT algorithm.
    !%
    !% <MATH>
    !%    U_{\rm SO}(t+\delta t, t) = \exp \lbrace - {i \over 2}\delta t T \rbrace
    !%       \exp \lbrace -i\delta t V^* \rbrace
    !%       \exp \lbrace - {i \over 2}\delta t T \rbrace
    !% </MATH>
    !%
    !% Since those three exponentials may be calculated exactly, one does not
    !% need to use any of the methods specified by variable <tt>TDExponentialMethod</tt>
    !% to perform them. 
    !%
    !% The key point is the specification of <math>V^*</math>. Let <math>V(t)</math> be divided into
    !% <math>V_{\rm int}(t)</math>, the "internal" potential which depends self-consistently
    !% on the density, and <math>V_{\rm ext}(t)</math>, the external potential that we know
    !% at all times since it is imposed to the system by us (e.g. a laser field):
    !% <math>V(t)=V_{\rm int}(t)+V_{\rm ext}(t)</math>. Then we define to be <math>V^*</math> to
    !% be the sum of <math>V_{\rm ext}(t+\delta t/2)</math> and the internal potential built
    !% from the wavefunctions <i>after</i> applying the right-most kinetic term
    !% of the equation, <math>\exp \lbrace -i\delta t/2 T \rbrace</math>.
    !%
    !% It may the be demonstrated that the order of the error of the algorithm is the
    !% same that the one that we would have by making use of the Exponential Midpoint Rule
    !% (EM, described below), the SO algorithm to calculate the action of the 
    !% exponential of the Hamiltonian, and a full self-consistent procedure.
    !%Option suzuki_trotter 1
    !% This is the generalization of the Suzuki-Trotter algorithm, described
    !% as one of the choices of the <tt>TDExponentialMethod</tt>,
    !% to time-dependent problem. Consult the paper by O. Sugino and M. Miyamoto,
    !% Phys. Rev. B <b>59</b>, 2579 (1999), for details.
    !%
    !% It requires of Hamiltonian extrapolations.
    !%Option etrs 2
    !% The idea is to make use of the time-reversal symmetry from the beginning:
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
    !% it is extrapolated trough a second-order polynomial by making use of the
    !% Hamiltonian at time @math{t-2\delta t}, @math{t-\delta t} and @math{t}.
    !%Option exp_mid 4
    !% Exponential Midpoint Rule (EM).
    !% This is maybe the simplest method, but is is very well grounded theretically:
    !% it is unitary (if the exponential is performed correctly) and preserves
    !% time symmetry (if the self-consistency problem is dealt with correctly).
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
    !% Classical Crank-Nicholson propagator. Requires the sparskit library.
    !%
    !% <MATH>
    !%  (1 + i\delta t/2 H_{n+1/2}) \psi_{n+1} = (1 - i\delta t/2 H_{n+1/2}) \psi_{n}  
    !% </MATH>
    !%Option magnus 7
    !% Magnus Expansion (M4).
    !% This is the most sophisticated approach. It is a fourth order scheme (feature
    !% that shares with the ST scheme; the other schemes are in principle second order).
    !% It is tailored for making use of very large time steps, or equivalently,
    !% dealing with problem with very high-frequency time dependence.
    !% It is still in a experimental state; we are not yet sure of when it is
    !% advantageous.
    !%Option crank_nicholson_src_mem 8
    !% Crank-Nicholson propagator with source and memory term for transport
    !% calculations.
    !%Option visscher 9
    !% (experimental) Visscher integration scheme. Computational Physics 5 596 (1991).
    !%Option qoct_tddft_propagator 10
    !% WARNING: EXPERIMENTAL
    !%End

    if(gr%sb%open_boundaries) then
      default_propagator = PROP_CRANK_NICHOLSON_SRC_MEM
    else
      default_propagator = PROP_REVERSAL
    end if
    call loct_parse_int(datasets_check('TDEvolutionMethod'), default_propagator, tr%method)
    if(.not.varinfo_valid_option('TDEVolutionMethod', tr%method)) call input_error('TDEvolutionMethod')

    if(gr%sb%open_boundaries.and.tr%method.ne.PROP_CRANK_NICHOLSON_SRC_MEM) then
      message(1) = 'The time evolution method for time dependent cannot'
      message(2) = 'be chosen freely. The Crank-Nicholson propagator'
      message(3) = 'with source and memory term has to be used. Either set'
      message(4) = ''
      message(5) = '  TDEvolutionmethod = crank_nicholson_src_mem'
      message(6) = ''
      message(7) = 'in your input or remove the TDEvolutionMethod variable.'
      call write_fatal(7)
    end if

    if(tr%method .eq. SPLIT_OPERATOR .or. tr%method .eq. SUZUKI_TROTTER) then
      message(1) = "You cannnot use the split operator evolution method, or the"
      message(2) = "Suzuki-Trotter, if the code was compiled without FFTW support."
      call write_fatal(2)
    end if

    select case(tr%method)
    case(PROP_SPLIT_OPERATOR)
      call zcf_new(gr%mesh%idx%ll, tr%cf)
      call zcf_fft_init(tr%cf, gr%sb)
    case(PROP_SUZUKI_TROTTER)
      call zcf_new(gr%mesh%idx%ll, tr%cf)
      call zcf_fft_init(tr%cf, gr%sb)
    case(PROP_REVERSAL)
    case(PROP_APP_REVERSAL)
    case(PROP_EXPONENTIAL_MIDPOINT)
    case(PROP_CRANK_NICHOLSON)
    case(PROP_CRANK_NICHOLSON_SPARSKIT)
#ifdef HAVE_SPARSKIT
      SAFE_ALLOCATE(tdsk)
      call zsparskit_solver_init(gr%mesh%np, tdsk)
      SAFE_ALLOCATE(zpsi_tmp(1:gr%mesh%np_part, 1:st%d%dim, 1:st%nst, st%d%kpt%start:st%d%kpt%end))
#else
      message(1) = 'Octopus was not compiled with support for the sparskit library. This'
      message(2) = 'library is required if the "crank_nicholson_sparskit" propagator is selected.'
      message(3) = 'Try to use a different propagation scheme or recompile with sparskit support.'
      call write_fatal(3)
#endif
    case(PROP_MAGNUS)
      SAFE_ALLOCATE(tr%vmagnus(1:gr%mesh%np, 1:st%d%nspin, 1:2))
    case(PROP_CRANK_NICHOLSON_SRC_MEM)
      call ob_rti_init(st, gr, hm, tr%ob, dt, max_iter)
    case(PROP_VISSCHER)
    case(PROP_QOCT_TDDFT_PROPAGATOR)
    case default
      call input_error('TDEvolutionMethod')
    end select
    call messages_print_var_option(stdout, 'TDEvolutionMethod', tr%method)

    if(have_fields) then
      if(tr%method /= PROP_REVERSAL .and.    &
         tr%method /= PROP_APP_REVERSAL .and. &
         tr%method /= PROP_VISSCHER .and. &
         tr%method /= PROP_EXPONENTIAL_MIDPOINT) then
        message(1) = "To move the ions or put a gauge field use the etrs, aetrs, visscher or exp_mid propagators." 
        call write_fatal(1)
      end if
    end if

    ! Allocate memory to store the old KS potentials
    SAFE_ALLOCATE(tr%v_old(1:gr%mesh%np, 1:st%d%nspin, 0:3))
    tr%v_old(:, :, :) = M_ZERO
    call exponential_init(tr%te, gr) ! initialize propagator

    ! By default, the propagation is only self-consistent in the first iterations
    ! (unless we are doing a QOCT run)
    tr%scf_propagation = .false.

    if(tr%method == PROP_VISSCHER) then
      SAFE_ALLOCATE(tr%prev_psi(1:gr%mesh%np, 1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
      tr%first = .true.
    else
      nullify(tr%prev_psi)
    end if

    call pop_sub()
  end subroutine td_rti_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine td_rti_set_scf_prop(tr)
    type(td_rti_t), intent(inout) :: tr
    tr%scf_propagation = .true.
  end subroutine td_rti_set_scf_prop
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine td_rti_remove_scf_prop(tr)
    type(td_rti_t), intent(inout) :: tr
    tr%scf_propagation = .false.
  end subroutine td_rti_remove_scf_prop
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine td_rti_end(tr)
    type(td_rti_t), intent(inout) :: tr

    SAFE_DEALLOCATE_P(tr%prev_psi)

    ! sanity check
    ASSERT(associated(tr%v_old)) 
    SAFE_DEALLOCATE_P(tr%v_old)         ! clean ols KS potentials
    nullify(tr%v_old)

    select case(tr%method)
    case(PROP_MAGNUS)
      ASSERT(associated(tr%vmagnus))
      SAFE_DEALLOCATE_P(tr%vmagnus); nullify(tr%vmagnus)
    case(PROP_CRANK_NICHOLSON_SPARSKIT)
#ifdef HAVE_SPARSKIT
      call zsparskit_solver_end()
      SAFE_DEALLOCATE_A(zpsi_tmp)
#endif
    case(PROP_SUZUKI_TROTTER, PROP_SPLIT_OPERATOR)
      call zcf_free(tr%cf)
    case(PROP_CRANK_NICHOLSON_SRC_MEM)
      call ob_rti_end(tr%ob)
    end select
    
    call exponential_end(tr%te)       ! clean propagator method
  end subroutine td_rti_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine td_rti_run_zero_iter(hm, tr)
    type(hamiltonian_t), intent(in)    :: hm
    type(td_rti_t),      intent(inout) :: tr

    tr%v_old(:, :, 2) = hm%vhxc(:, :)
    tr%v_old(:, :, 3) = hm%vhxc(:, :)
    tr%v_old(:, :, 1) = hm%vhxc(:, :)
  end subroutine td_rti_run_zero_iter


  ! ---------------------------------------------------------
  ! Propagates st from t-dt to t.
  ! If dt<0, it propagates *backwards* from t+|dt| to t
  ! ---------------------------------------------------------
  subroutine td_rti_dt(ks, hm, gr, st, tr, t, dt, max_iter, nt, gauge_force, ions, geo, ionic_dt)
    type(v_ks_t),                    intent(inout) :: ks
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    type(td_rti_t),      target,     intent(inout) :: tr
    FLOAT,                           intent(in)    :: t, dt
    integer,                         intent(in)    :: max_iter
    integer,                         intent(in)    :: nt
    type(gauge_force_t),  optional,  intent(inout) :: gauge_force
    type(ion_dynamics_t), optional,  intent(inout) :: ions
    type(geometry_t),     optional,  intent(inout) :: geo
    FLOAT,                optional,  intent(in)    :: ionic_dt

    integer :: is, iter
    FLOAT   :: d, d_max
    logical :: self_consistent
    CMPLX, allocatable :: zpsi1(:, :, :, :)
    FLOAT, allocatable :: dtmp(:), vaux(:, :)
    type(profile_t), save :: prof

    call profiling_in(prof, "TD_PROPAGATOR")
    call push_sub('td_rti.td_rti_dt')

    if(present(ions)) then
      ASSERT(present(geo))
      ASSERT(present(ionic_dt))
    end if

    if(gauge_field_is_applied(hm%ep%gfield)) then
      ASSERT(present(gauge_force))
    end if

    self_consistent = .false.
    if(hm%theory_level .ne. INDEPENDENT_PARTICLES) then
      if( (t < 3*dt)  .or.  (tr%scf_propagation) ) then
        self_consistent = .true.
        SAFE_ALLOCATE(zpsi1(1:gr%mesh%np_part, 1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
        zpsi1 = st%zpsi
      end if
    end if

    call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 2), tr%v_old(:, :, 3))
    call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 1), tr%v_old(:, :, 2))
    call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc(:, :),      tr%v_old(:, :, 1))
    call interpolate( (/t-dt, t-2*dt, t-3*dt/), tr%v_old(:, :, 1:3), t, tr%v_old(:, :, 0))

    select case(tr%method)
    case(PROP_SPLIT_OPERATOR);          call td_split_operator
    case(PROP_SUZUKI_TROTTER);          call td_suzuki_trotter
    case(PROP_REVERSAL);                call td_reversal
    case(PROP_APP_REVERSAL);            call td_app_reversal
    case(PROP_EXPONENTIAL_MIDPOINT);    call exponential_midpoint
    case(PROP_CRANK_NICHOLSON);         call td_crank_nicholson
    case(PROP_CRANK_NICHOLSON_SPARSKIT);call td_crank_nicholson_sparskit
    case(PROP_MAGNUS);                  call td_magnus
    case(PROP_CRANK_NICHOLSON_SRC_MEM); call td_crank_nicholson_src_mem
    case(PROP_VISSCHER);                call td_visscher
    case(PROP_QOCT_TDDFT_PROPAGATOR);   call td_qoct_tddft_propagator
    end select

    if(self_consistent) then
      SAFE_ALLOCATE(vaux(1:gr%mesh%np, 1:st%d%nspin))

      ! First, compare the new potential to the extrapolated one.
      call states_calc_dens(st, gr%mesh%np)
      call v_ks_calc(gr, ks, hm, st)
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

          st%zpsi = zpsi1
          tr%v_old(:, :, 0) = hm%vhxc(:, :)
          vaux(:, :) = hm%vhxc(:, :)
          select case(tr%method)
          case(PROP_SPLIT_OPERATOR);          call td_split_operator
          case(PROP_SUZUKI_TROTTER);          call td_suzuki_trotter
          case(PROP_REVERSAL);                call td_reversal
          case(PROP_APP_REVERSAL);            call td_app_reversal
          case(PROP_EXPONENTIAL_MIDPOINT);    call exponential_midpoint
          case(PROP_CRANK_NICHOLSON);         call td_crank_nicholson
          case(PROP_CRANK_NICHOLSON_SPARSKIT);call td_crank_nicholson_sparskit
          case(PROP_MAGNUS);                  call td_magnus
          case(PROP_CRANK_NICHOLSON_SRC_MEM); call td_crank_nicholson_src_mem
          case(PROP_VISSCHER);                call td_visscher
          case(PROP_QOCT_TDDFT_PROPAGATOR);   call td_qoct_tddft_propagator
          end select

          call states_calc_dens(st, gr%mesh%np)
          call v_ks_calc(gr, ks, hm, st)
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

    call pop_sub()
    call profiling_out(prof)

  contains

    ! ---------------------------------------------------------
    ! Split operator.
    subroutine td_split_operator
      integer :: ik, ist
      call push_sub('td_rti.td_split_operator')

      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = 1, st%nst
          call zexp_kinetic(gr, hm, st%zpsi(:, :, ist, ik), tr%cf, -M_HALF*M_zI*dt)
        end do
      end do
      call states_calc_dens(st, gr%mesh%np)
      call v_ks_calc(gr, ks, hm, st)
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = 1, st%nst
          if (hm%ep%non_local) call zexp_vnlpsi (gr%mesh, hm, st%zpsi(:, :, ist, ik), -M_zI*dt, .true.)
          call zexp_vlpsi (gr, hm, st%zpsi(:, :, ist, ik), ik, t-dt*M_HALF, -M_zI*dt)
          if (hm%ep%non_local) call zexp_vnlpsi (gr%mesh, hm, st%zpsi(:, :, ist, ik), -M_zI*dt, .false.)
        end do
      end do
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = 1, st%nst
          call zexp_kinetic(gr, hm, st%zpsi(:, :, ist, ik), tr%cf, -M_HALF*M_zI*dt)
        end do
      end do

      call pop_sub()
    end subroutine td_split_operator


    ! ---------------------------------------------------------
    ! Suzuki-Trotter.
    subroutine td_suzuki_trotter
      FLOAT :: p, pp(5), time(5), dtime(5)
      integer :: ik, ist, k

      call push_sub('td_rti.td_suzuki_trotter')

      p = M_ONE/(M_FOUR - M_FOUR**(M_THIRD))
      pp = (/ p, p, M_ONE-M_FOUR*p, p, p /)
      dtime = pp*dt
      time(1) = t-dt+pp(1)/M_TWO*dt
      time(2) = t-dt+(pp(1)+pp(2)/M_TWO)*dt
      time(3) = t-dt+(pp(1)+pp(2)+pp(3)/M_TWO)*dt
      time(4) = t-dt+(pp(1)+pp(2)+pp(3)+pp(4)/M_TWO)*dt
      time(5) = t-dt+(pp(1)+pp(2)+pp(3)+pp(4)+pp(5)/M_TWO)*dt

      do k = 1, 5
        call interpolate( (/t, t-dt, t-2*dt/), tr%v_old(:, :, 0:2), time(k), hm%vhxc(:, :))
        call hamiltonian_update_potential(hm, gr%mesh)
        do ik = st%d%kpt%start, st%d%kpt%end
          do ist = 1, st%nst
            call zexp_vlpsi (gr, hm, st%zpsi(:, :, ist, ik), ik, time(k), -M_zI*dtime(k)/M_TWO)
            if (hm%ep%non_local) call zexp_vnlpsi (gr%mesh, hm, &
              st%zpsi(:, :, ist, ik), -M_zI*dtime(k)/M_TWO, .true.)

            call zexp_kinetic(gr, hm, st%zpsi(:, :, ist, ik), tr%cf, -M_zI*dtime(k))
            if (hm%ep%non_local) call zexp_vnlpsi (gr%mesh, hm, &
              st%zpsi(:, :, ist, ik), -M_zI*dtime(k)/M_TWO, .false.)
            call zexp_vlpsi (gr, hm, st%zpsi(:, :, ist, ik), ik, time(k), -M_zI*dtime(k)/M_TWO)
          end do
        end do
      end do

      call pop_sub()
    end subroutine td_suzuki_trotter

    ! ---------------------------------------------------------
    ! Propagator with enforced time-reversal symmetry
    subroutine td_reversal
      FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
      CMPLX, allocatable :: zpsi1(:,:)
      integer :: ik, ist, idim

      call push_sub('td_rti.td_reversal')

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then

        SAFE_ALLOCATE(vhxc_t1(1:gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE(vhxc_t2(1:gr%mesh%np, 1:st%d%nspin))
        call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t1)

        SAFE_ALLOCATE(zpsi1(1:gr%mesh%np_part, 1:st%d%dim))

        st%rho(1:gr%mesh%np, 1:st%d%nspin) = M_ZERO

        do ik = st%d%kpt%start, st%d%kpt%end
          do ist = st%st_start, st%st_end
            
            !save the state
            do idim = 1, st%d%dim
              call lalg_copy(gr%mesh%np, st%zpsi(:, idim, ist, ik), zpsi1(:, idim))
            end do
            
            !propagate the state dt with H(t-dt)
            call exponential_apply(tr%te, gr, hm, st%zpsi(:,:, ist, ik), ist, ik, dt, t-dt)
            
            !calculate the contribution to the density
            call states_dens_accumulate(st, gr%mesh%np, st%rho, ist, ik)
            
            !restore the saved state
            do idim = 1, st%d%dim
              call lalg_copy(gr%mesh%np, zpsi1(:, idim), st%zpsi(:, idim, ist, ik))
            end do
          end do
          
        end do

        SAFE_DEALLOCATE_A(zpsi1)
        
        ! finish the calculation of the density
        call states_dens_reduce(st, gr%mesh%np, st%rho)

        call v_ks_calc(gr, ks, hm, st)

        call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t2)
        call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t1, hm%vhxc)
        call hamiltonian_update_potential(hm, gr%mesh)
      end if

      ! propagate dt/2 with H(t-dt)
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          call exponential_apply(tr%te, gr, hm, st%zpsi(:,:, ist, ik), ist, ik, dt/M_TWO, t-dt)
        end do
      end do

      ! propagate dt/2 with H(t)

      ! first move the ions to time t
      if(present(ions)) then
        call ion_dynamics_propagate(ions, gr%sb, geo, t, ionic_dt)
        call hamiltonian_epot_generate(hm, gr, geo, st, time = t)
      end if

      if(gauge_field_is_applied(hm%ep%gfield)) call gauge_field_propagate(hm%ep%gfield, gauge_force, dt)

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t2, hm%vhxc)
        call hamiltonian_update_potential(hm, gr%mesh)
      end if

      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          call exponential_apply(tr%te, gr, hm, st%zpsi(:,:, ist, ik), ist, ik, dt/M_TWO, t)
        end do
      end do

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        SAFE_DEALLOCATE_A(vhxc_t1)
        SAFE_DEALLOCATE_A(vhxc_t2)
      end if

      call pop_sub()
    end subroutine td_reversal


    ! ---------------------------------------------------------
    ! Propagator with approximate enforced time-reversal symmetry
    subroutine td_app_reversal
      integer :: ik, sts, ste
      type(batch_t) :: zpsib

      call push_sub('td_rti.td_app_reversal')

      ! propagate half of the time step with H(t-dt)
      do ik = st%d%kpt%start, st%d%kpt%end
        do sts = st%st_start, st%st_end, st%d%block_size
          ste = min(st%st_end, sts + st%d%block_size - 1)
          call batch_init(zpsib, st%d%dim, sts, ste, st%zpsi(:, :, sts:, ik))
          call exponential_apply_batch(tr%te, gr, hm, zpsib, ik, dt/M_TWO, t-dt)
          call batch_end(zpsib)
        end do
      end do

      ! interpolate the hamiltonian to time t
      call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 0), hm%vhxc)
      call hamiltonian_update_potential(hm, gr%mesh)

      ! move the ions to time t
      if(present(ions)) then      
        call ion_dynamics_propagate(ions, gr%sb, geo, t, ionic_dt)
        call hamiltonian_epot_generate(hm, gr, geo, st, time = t)
      end if

      if(gauge_field_is_applied(hm%ep%gfield)) call gauge_field_propagate(hm%ep%gfield, gauge_force, dt)

      ! propagate the other half with H(t)
      do ik = st%d%kpt%start, st%d%kpt%end
        do sts = st%st_start, st%st_end, st%d%block_size
          ste = min(st%st_end, sts + st%d%block_size - 1)
          call batch_init(zpsib, hm%d%dim, sts, ste, st%zpsi(:, :, sts:, ik))
          call exponential_apply_batch(tr%te, gr, hm, zpsib, ik, dt/M_TWO, t)
          call batch_end(zpsib)
        end do
      end do
      
      call pop_sub()
    end subroutine td_app_reversal


    ! ---------------------------------------------------------
    ! Exponential midpoint
    subroutine exponential_midpoint
      integer :: ist, ik
      type(ion_state_t) :: ions_state
      FLOAT :: vecpot(1:MAX_DIM), vecpot_vel(1:MAX_DIM)

      call push_sub('td_rti.exponential_midpoint')

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        call interpolate( (/t, t-dt, t-2*dt/), tr%v_old(:, :, 0:2), t-dt/M_TWO, hm%vhxc(:, :))
        call hamiltonian_update_potential(hm, gr%mesh)
      end if

      !move the ions to time t - dt/2
      if(present(ions)) then
        call ion_dynamics_save_state(ions, geo, ions_state)
        call ion_dynamics_propagate(ions, gr%sb, geo, t - ionic_dt/M_TWO, M_HALF*ionic_dt)
        call hamiltonian_epot_generate(hm, gr, geo, st, time = t - ionic_dt/M_TWO)
      end if
      
      if(gauge_field_is_applied(hm%ep%gfield)) then
        vecpot = gauge_field_get_vec_pot(hm%ep%gfield)
        vecpot_vel = gauge_field_get_vec_pot_vel(hm%ep%gfield)
        call gauge_field_propagate(hm%ep%gfield, gauge_force, M_HALF*dt)
      end if

      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          call exponential_apply(tr%te, gr, hm, st%zpsi(:,:, ist, ik), ist, ik, dt, t - dt/M_TWO)
        end do
      end do

      if(present(ions)) call ion_dynamics_restore_state(ions, geo, ions_state)

      if(gauge_field_is_applied(hm%ep%gfield)) then
        call gauge_field_set_vec_pot(hm%ep%gfield, vecpot)
        call gauge_field_set_vec_pot_vel(hm%ep%gfield, vecpot_vel)
      end if

      call pop_sub()
    end subroutine exponential_midpoint


    ! ---------------------------------------------------------
    ! Crank-Nicholson propagator, linear solver from sparskit.
    subroutine td_crank_nicholson_sparskit
#ifdef HAVE_SPARSKIT
      FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
      CMPLX, allocatable :: zpsi_rhs_pred(:,:,:,:), zpsi_rhs_corr(:,:,:,:)
      integer :: ik, ist, idim, np_part

      call push_sub('td_rti.td_crank_nicholson_sparskit')

      np_part = gr%mesh%np_part
      SAFE_ALLOCATE(zpsi_rhs_corr(1:np_part, 1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
      zpsi_rhs_corr = st%zpsi ! store zpsi for corrector step

      ! define pointer and variables for usage in td_zop, td_zopt routines
      grid_p    => gr
      hm_p      => hm
      tr_p      => tr
      dt_op = dt
      t_op  = t

      ! we (ab)use exponential_apply to compute (1-i\delta t/2 H_n)\psi^n
      ! exponential order needs to be only 1
      tr%te%exp_method = 3 ! == taylor expansion
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
            call exponential_apply(tr%te, gr, hm, zpsi_rhs_pred(:, :, ist, ik), ist, ik, dt/M_TWO, t-dt)
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
              idim_op = idim; ist_op = ist; ik_op = ik
              call zsparskit_solver_run(tdsk, td_zop, td_zopt, &
                st%zpsi(1:gr%mesh%np, idim, ist, ik), zpsi_rhs_pred(1:gr%mesh%np, idim, ist, ik))
            end do
          end do
        end do

        call states_calc_dens(st, gr%mesh%np)
        call v_ks_calc(gr, ks, hm, st)

        vhxc_t2 = hm%vhxc
        ! compute potential at n+1/2 as average
        hm%vhxc = (vhxc_t1 + vhxc_t2)/M_TWO
        call hamiltonian_update_potential(hm, gr%mesh)
      end if

      ! get rhs of CN linear system (rhs2 = (1-i\delta t H_{n+1/2})\psi^n)
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          call exponential_apply(tr%te, gr, hm, zpsi_rhs_corr(:, :, ist, ik), ist, ik, dt/M_TWO, t-dt)
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
            idim_op = idim; ist_op = ist; ik_op = ik
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

      call pop_sub()
#endif
    end subroutine td_crank_nicholson_sparskit
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    ! Crank-Nicholson propagator, QMR linear solver.
    subroutine td_crank_nicholson
      CMPLX, allocatable :: zpsi_rhs(:,:), zpsi(:), rhs(:)
      integer :: ik, ist, idim, isize, np_part, np, iter
      FLOAT :: dres
      FLOAT :: cgtol = CNST(1.0e-8)
      logical :: converged

      call push_sub('td_rti.td_crank_nicholson')

      np_part = gr%mesh%np_part
      np = gr%mesh%np
      isize = np_part*st%lnst*st%d%kpt%nlocal*st%d%dim

      ! define pointer and variables for usage in td_zop, td_zopt routines
      grid_p    => gr
      hm_p      => hm
      tr_p      => tr
      dt_op = dt
      t_op  = t
      dim_op = st%d%dim

      ! we (ab)use exponential_apply to compute (1-i\delta t/2 H_n)\psi^n
      ! exponential order needs to be only 1
      tr%te%exp_method = TAYLOR
      tr%te%exp_order  = 1

      SAFE_ALLOCATE(zpsi_rhs(1:np_part, 1:st%d%dim))
      SAFE_ALLOCATE(zpsi(1:np*st%d%dim))
      SAFE_ALLOCATE(rhs(1:np*st%d%dim))
        
      call interpolate( (/t, t-dt, t-2*dt/), tr%v_old(:, :, 0:2), t-dt/M_TWO, hm%vhxc(:, :))
      call hamiltonian_update_potential(hm, gr%mesh)

      ! solve (1+i\delta t/2 H_n)\psi^{predictor}_{n+1} = (1-i\delta t/2 H_n)\psi^n
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end

          zpsi_rhs(:, :) = st%zpsi(:, :, ist, ik)
          call exponential_apply(tr%te, gr, hm, zpsi_rhs, ist, ik, dt/M_TWO, t-dt)
          if(hamiltonian_inh_term(hm)) &
            zpsi_rhs(:, :) = zpsi_rhs(:, :) + dt * hm%inh_st%zpsi(:, :, ist, ik)

          forall(idim = 1:st%d%dim)
            zpsi((idim-1)*np+1:idim*np) = st%zpsi(1:np, idim, ist, ik)
            rhs((idim-1)*np+1:idim*np) = zpsi_rhs(1:np, idim)
          end forall

          ist_op = ist; ik_op = ik; iter = 2000
          call zqmr_sym(np*st%d%dim, zpsi, rhs, &
            td_rti_qmr_op, td_rti_qmr_prec, iter, dres, cgtol, &
            showprogress = .false., converged = converged)

          forall(idim = 1:st%d%dim)
            st%zpsi(1:np, idim, ist, ik) = &
              zpsi((idim-1)*np_part+1:(idim-1)*np_part+np)
          end forall

          if(.not.converged) then
            write(message(1),'(a)')        'The linear solver used for the Crank-Nicholson'
            write(message(2),'(a,es14.4)') 'propagator did not converge: Residual = ', dres
            call write_warning(2)
          end if

        end do
      end do

      SAFE_DEALLOCATE_A(zpsi_rhs)
      SAFE_DEALLOCATE_A(zpsi)
      SAFE_DEALLOCATE_A(rhs)
      call pop_sub()
    end subroutine td_crank_nicholson
    ! ---------------------------------------------------------



    ! ---------------------------------------------------------
    ! Magnus propagator
    subroutine td_magnus
      integer :: j, is, ist, ik, i
      FLOAT :: time(2)
      FLOAT, allocatable :: vaux(:, :, :), pot(:)

      call push_sub('td_rti.td_magnus')

      SAFE_ALLOCATE(vaux(1:gr%mesh%np, 1:st%d%nspin, 1:2))

      time(1) = (M_HALF-sqrt(M_THREE)/M_SIX)*dt
      time(2) = (M_HALF+sqrt(M_THREE)/M_SIX)*dt

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        do j = 1, 2
          call interpolate( (/t, t-dt, t-2*dt/), tr%v_old(:, :, 0:2), time(j)-dt, hm%vhxc(:, :))
          call hamiltonian_update_potential(hm, gr%mesh)
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
            call laser_potential(gr%sb, hm%ep%lasers(i), gr%mesh, pot, t-dt+time(j))
            do is = 1, st%d%nspin
              vaux(:, is, j) = vaux(:, is, j) + pot(:)
            end do
            SAFE_DEALLOCATE_A(pot)
          case(E_FIELD_MAGNETIC, E_FIELD_VECTOR_POTENTIAL)
            write(message(1),'(a)') 'The Magnus propagator cannot be used with magnetic fields, or'
            write(message(2),'(a)') 'with an electric field described in the velocity gauge.'
            call write_fatal(2)
          end select
        end do
      end do

      tr%vmagnus(:, :, 2)  = M_HALF*(vaux(:, :, 1) + vaux(:, :, 2))
      tr%vmagnus(:, :, 1) = (sqrt(M_THREE)/CNST(12.0))*dt*(vaux(:, :, 2) - vaux(:, :, 1))

     do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          call exponential_apply(tr%te, gr, hm, st%zpsi(:,:, ist, ik), ist, ik, dt, M_ZERO, &
            vmagnus = tr%vmagnus)
        end do
      end do

      SAFE_DEALLOCATE_A(vaux)
      call pop_sub()
    end subroutine td_magnus


    ! ---------------------------------------------------------
    ! Crank-Nicholson scheme with source and memory term.
    subroutine td_crank_nicholson_src_mem()
      call push_sub('td_rti.td_crank_nicholson_src_mem')
      
      select case(tr%ob%mem_type)
      case(SAVE_CPU_TIME)
        call cn_src_mem_dt(tr%ob, st, ks, hm, gr, max_iter, dt, t, nt)
      case(SAVE_RAM_USAGE)
        call cn_src_mem_sp_dt(tr%ob, st, ks, hm, gr, max_iter, dt, t, nt)
      end select

      call pop_sub()
    end subroutine td_crank_nicholson_src_mem

    ! ---------------------------------------------------------
    ! Propagator with approximate enforced time-reversal symmetry
    subroutine td_visscher
      integer :: ik, ist, idim, ip
      FLOAT, allocatable :: dpsi(:, :), hpsi(:, :)
      CMPLX, allocatable :: zpsi(:,:)

      call push_sub('td_rti.td_app_reversal')
      
      SAFE_ALLOCATE(dpsi(1:gr%mesh%np_part, 1:st%d%dim))
      SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim))

      ! we have to initialize the imaginary part of \psi at dt/2, so we
      ! use the exponential midpoint rule.
      if(tr%first) then

        !propagate the hamiltonian to t - dt/4
        if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
          call interpolate( (/t, t-dt, t-2*dt/), tr%v_old(:, :, 0:2), t - CNST(0.75)*dt, hm%vhxc(:, :))
        end if
        
        if(present(ions)) then
          call ion_dynamics_propagate(ions, gr%sb, geo, t - CNST(0.75)*ionic_dt, CNST(0.25)*ionic_dt)
          call hamiltonian_epot_generate(hm, gr, geo, st, time = t - CNST(0.75)/M_FOUR)
        else
          call hamiltonian_update_potential(hm, gr%mesh)
        end if
        
        if(gauge_field_is_applied(hm%ep%gfield)) call gauge_field_propagate(hm%ep%gfield, gauge_force, CNST(0.25)*dt)
        
        SAFE_ALLOCATE(zpsi(1:gr%mesh%np_part, 1:st%d%dim))
        do ik = st%d%kpt%start, st%d%kpt%end
          do ist = st%st_start, st%st_end
            do idim = 1, st%d%dim
              call lalg_copy(gr%mesh%np, st%zpsi(:, idim, ist, ik), zpsi(:, idim))
            end do
            call exponential_apply(tr%te, gr, hm, zpsi, ist, ik, dt/M_TWO, t - CNST(0.75)*dt)
            forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np) tr%prev_psi(ip, idim, ist, ik) = aimag(zpsi(ip, idim))
          end do
        end do
        SAFE_DEALLOCATE_A(zpsi)
        tr%first = .false.
        
        !finish to propagate the hamiltonian to t - dt/2
        
        if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
          call interpolate( (/t, t-dt, t-2*dt/), tr%v_old(:, :, 0:2), t - CNST(0.5)*dt, hm%vhxc(:, :))
        end if
        
        if(present(ions)) then
          call ion_dynamics_propagate(ions, gr%sb, geo, t - CNST(0.5)*ionic_dt, CNST(0.25)*ionic_dt)
          call hamiltonian_epot_generate(hm, gr, geo, st, time = t - CNST(0.5)/M_FOUR)
        else
          call hamiltonian_update_potential(hm, gr%mesh)
        end if
        
        if(gauge_field_is_applied(hm%ep%gfield)) call gauge_field_propagate(hm%ep%gfield, gauge_force, CNST(0.25)*dt)
        
      else
        
        !directly propagate the hamiltonian to t - dt/2
        
        if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
          call interpolate( (/t, t-dt, t-2*dt/), tr%v_old(:, :, 0:2), t-dt/M_TWO, hm%vhxc(:, :))
          call hamiltonian_update_potential(hm, gr%mesh)
        end if
        
        if(present(ions)) then
          call ion_dynamics_propagate(ions, gr%sb, geo, t - ionic_dt/M_TWO, M_HALF*ionic_dt)
          call hamiltonian_epot_generate(hm, gr, geo, st, time = t - ionic_dt/M_TWO)
        end if
        
        if(gauge_field_is_applied(hm%ep%gfield)) call gauge_field_propagate(hm%ep%gfield, gauge_force, M_HALF*dt)
        
      end if

      ! propagate the real part
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          
          do idim = 1, st%d%dim
            call lalg_copy(gr%mesh%np, tr%prev_psi(:, idim, ist, ik), dpsi(:, idim))
          end do

          call dhamiltonian_apply(hm, gr, dpsi, hpsi, ist, ik, t - dt/2)

          forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np) 
            st%zpsi(ip, idim, ist, ik) = real(st%zpsi(ip, idim, ist, ik)) + cmplx(dt*hpsi(ip, idim), &
                 ! the imaginary part is calculated as the average
                 ! this is the first half
                 CNST(0.5)*tr%prev_psi(ip, idim, ist, ik), REAL_PRECISION)
          end forall

        end do
      end do

      ! propagate the hamiltonian to time t
      call lalg_copy(gr%mesh%np, st%d%nspin, tr%v_old(:, :, 0), hm%vhxc)
      call hamiltonian_update_potential(hm, gr%mesh)

      if(present(ions)) then
        call ion_dynamics_propagate(ions, gr%sb, geo, t, M_HALF*ionic_dt)
        call hamiltonian_epot_generate(hm, gr, geo, st, time = t)
      end if

      if(gauge_field_is_applied(hm%ep%gfield)) call gauge_field_propagate(hm%ep%gfield, gauge_force, M_HALF*dt)

      ! propagate the imaginary part
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end

          forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np) dpsi(ip, idim) = real(st%zpsi(ip, idim, ist, ik))

          call dhamiltonian_apply(hm, gr, dpsi, hpsi, ist, ik, t)

          forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np)
            tr%prev_psi(ip, idim, ist, ik) = tr%prev_psi(ip, idim, ist, ik) - dt*hpsi(ip, idim)
            ! this is the second half of the average
            st%zpsi(ip, idim, ist, ik) = st%zpsi(ip, idim, ist, ik) + &
                 cmplx(M_ZERO, CNST(0.5)*tr%prev_psi(ip, idim, ist, ik), REAL_PRECISION)
          end forall

        end do
      end do

      call pop_sub()
    end subroutine td_visscher


    ! ---------------------------------------------------------
    ! Propagator specifically designed for the QOCT+TDDFT problem
    subroutine td_qoct_tddft_propagator
      type(ion_state_t) :: ions_state
      FLOAT :: vecpot(1:MAX_DIM), vecpot_vel(1:MAX_DIM)
      call push_sub('td_rti.td_qoct_tddft_propagator')

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        call interpolate( (/t, t-dt, t-2*dt/), tr%v_old(:, :, 0:2), t-dt/M_TWO, hm%vhxc(:, :))
        call hamiltonian_update_potential(hm, gr%mesh)
      end if

      !move the ions to time t - dt/2
      if(present(ions)) then
        call ion_dynamics_save_state(ions, geo, ions_state)
        call ion_dynamics_propagate(ions, gr%sb, geo, t - ionic_dt/M_TWO, M_HALF*ionic_dt)
        call hamiltonian_epot_generate(hm, gr, geo, st, time = t - ionic_dt/M_TWO)
      end if
      
      if(gauge_field_is_applied(hm%ep%gfield)) then
        vecpot = gauge_field_get_vec_pot(hm%ep%gfield)
        vecpot_vel = gauge_field_get_vec_pot_vel(hm%ep%gfield)
        call gauge_field_propagate(hm%ep%gfield, gauge_force, M_HALF*dt)
      end if

      call exponential_apply_all(tr%te, gr, hm, st, dt, t - dt/M_TWO)

      if(present(ions)) call ion_dynamics_restore_state(ions, geo, ions_state)

      if(gauge_field_is_applied(hm%ep%gfield)) then
        call gauge_field_set_vec_pot(hm%ep%gfield, vecpot)
        call gauge_field_set_vec_pot_vel(hm%ep%gfield, vecpot_vel)
      end if

      call pop_sub()
    end subroutine td_qoct_tddft_propagator

  end subroutine td_rti_dt
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! operators for Crank-Nicholson scheme
  subroutine td_rti_qmr_op(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)
    integer :: idim
    CMPLX, allocatable :: zpsi(:, :)

    SAFE_ALLOCATE(zpsi(1:grid_p%mesh%np_part, 1:dim_op))
    zpsi = M_z0
    forall(idim = 1:dim_op) &
      zpsi(1:grid_p%mesh%np, idim) = x((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np)
    call exponential_apply(tr_p%te, grid_p, hm_p, zpsi, ist_op, ik_op, -dt_op/M_TWO, t_op)
    forall(idim = 1:dim_op) &
      y((idim-1)*grid_p%mesh%np+1:idim*grid_p%mesh%np) = zpsi(1:grid_p%mesh%np, idim)
    SAFE_DEALLOCATE_A(zpsi)

  end subroutine td_rti_qmr_op
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine td_rti_qmr_prec(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)
    y = x
  end subroutine td_rti_qmr_prec
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! operators for Crank-Nicholson scheme
  subroutine td_zop(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:), xim(:)
    FLOAT, intent(out) :: yre(:), yim(:)

    call push_sub('td_rti.td_zop')
#ifdef HAVE_SPARSKIT    
    zpsi_tmp(1:grid_p%mesh%np, idim_op, ist_op, ik_op) = &
      xre(1:grid_p%mesh%np) + M_zI*xim(1:grid_p%mesh%np)

    ! propagate backwards
    call exponential_apply(tr_p%te, grid_p, hm_p, zpsi_tmp(:, :, ist_op, ik_op), ist_op, ik_op, -dt_op/M_TWO, t_op)

    yre(1:grid_p%mesh%np) =  real(zpsi_tmp(1:grid_p%mesh%np, idim_op, ist_op, ik_op))
    yim(1:grid_p%mesh%np) = aimag(zpsi_tmp(1:grid_p%mesh%np, idim_op, ist_op, ik_op))
#endif    
    call pop_sub()
  end subroutine td_zop
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Transpose of H (called e.g. by bi-conjugate gradient solver)
  subroutine td_zopt(xre, xim, yre, yim)
    FLOAT, intent(in)  :: xre(:), xim(:)
    FLOAT, intent(out) :: yre(:), yim(:)
    
    call push_sub('td_rti.td_zopt')
#ifdef HAVE_SPARSKIT        
    ! To act with the transpose of H on the wfn we apply H to the conjugate of psi
    ! and conjugate the resulting hpsi (note that H is not a purely real operator
    ! for scattering wavefunctions anymore).
    zpsi_tmp(1:grid_p%mesh%np, idim_op, ist_op, ik_op) = &
      xre(1:grid_p%mesh%np) - M_zI*xim(1:grid_p%mesh%np)
    
    ! propagate backwards
    call exponential_apply(tr_p%te, grid_p, hm_p, zpsi_tmp(:, :, ist_op, ik_op), ist_op, ik_op, -dt_op/M_TWO, t_op)

    yre(1:grid_p%mesh%np) =    real(zpsi_tmp(1:grid_p%mesh%np, idim_op, ist_op, ik_op))
    yim(1:grid_p%mesh%np) = - aimag(zpsi_tmp(1:grid_p%mesh%np, idim_op, ist_op, ik_op))
#endif        
    call pop_sub()
  end subroutine td_zopt
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  logical pure function td_rti_ions_are_propagated(tr) result(propagated)
    type(td_rti_t), intent(in) :: tr

    select case(tr%method)
    case(PROP_REVERSAL, PROP_APP_REVERSAL, PROP_VISSCHER)
      propagated = .true.
    case default
      propagated = .false.
    end select

  end function td_rti_ions_are_propagated
  ! ---------------------------------------------------------

end module td_rti_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

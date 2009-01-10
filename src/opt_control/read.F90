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
!! $Id: opt_control.F90 2862 2007-04-25 23:34:20Z acastro $


  ! ---------------------------------------------------------
  ! read the parameters for the optimal control run     
  subroutine oct_read_inp
    call push_sub('opt_control_read.oct_read_inp')  

    !%Variable OCTScheme
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default oct_algorithm_zbr98
    !%Description
    !% In order to find the optimal laser field for a given task, e.g., the excitation from an
    !% initial state to a predefined final state at the final time, optimal control theory can 
    !% be applied to quantum mechanics. The mathematical derivation leads a set of equations 
    !% which require the propagation of the wavefunction and a lagrange multiplier (sometimes
    !% comparable to a wavefunction). Several schemes have been sought to solve these control 
    !% equations which boils down to forward and backward propagations. However, the order in
    !% which these equations are solved makes a huge difference. Some schemes can be proven 
    !% to increase the value of the target functional (merit function) in each step. (In 
    !% practice this can be violated if the accuracy of the numerical time propagation is
    !% small. Most likely in 3D.)
    !%Option oct_algorithm_zbr98 1 
    !% Backward-Forward-Backward scheme described in JCP 108, 1953 (1998).
    !% Only possible if targetoperator is a projection operator
    !% Provides the fastest and most stable convergence.
    !% Monotonic convergence.
    !%Option oct_algorithm_zr98  2
    !% Forward-Backward-Forward scheme described in JCP 109,385 (1998).
    !% Works for projection and local target operators
    !% Convergence is stable but slower than ZBR98. 
    !% Note that local operators show an extremely slow convergence.
    !% Monotonic convergence.
    !%Option oct_algorithm_wg05  3
    !% Forward-Backward scheme described in J. Opt. B. 7 300 (2005).
    !% Works for all kind target operators and 
    !% can be used with all kind of filters and allows a fixed fluence.
    !% The price is a rather instable convergence. 
    !% If the restrictions set by the filter and fluence are reasonable, a good overlap can be 
    !% expected with 20 iterations.
    !% No monotonic convergence.
    !%Option oct_algorithm_mt03 4
    !% Basically an improved and generalized scheme. 
    !% Comparable to ZBR98/ZR98. See [Y. Maday and G. Turinici, J. Chem. Phys. 118, 
    !% 8191 (2003)].
    !%Option oct_algorithm_krotov 5
    !% The procedure reported in [D. Tannor, V. Kazakov and V.
    !% Orlov, in "Time Dependent Quantum Molecular Dynamics", edited by J. Broeckhove
    !% and L. Lathouweres (Plenum, New York, 1992), pp. 347-360].
    !%Option oct_algorithm_straight_iteration 6
    !% Straight iteration: one forward and one backward propagation is performed at each
    !% iteration, both with the same control field. An output field is calculated with the
    !% resulting wave functions. Note that this scheme typically does not converge, unless
    !% some mixing ("OCTMixing = yes") is used.
    !%Option oct_algorithm_direct 7
    !% Direct optimization (experimental)
    !%Option oct_algorithm_newuoa 8
    !% Direct optimization with the newuoa algorithm (experimental)
    !%End
    call loct_parse_int(datasets_check('OCTScheme'), oct_algorithm_zr98, oct%algorithm)
    if(.not.varinfo_valid_option('OCTScheme', oct%algorithm)) call input_error('OCTScheme')
    select case(oct%algorithm)
    case(oct_algorithm_mt03)
      oct%delta = M_TWO; oct%eta = M_ZERO
    case(oct_algorithm_zr98)
      oct%delta = M_ONE; oct%eta = M_ONE
    case(oct_algorithm_krotov)
      oct%delta = M_ONE; oct%eta = M_ZERO
    case(oct_algorithm_str_iter)
      oct%delta = M_ZERO; oct%eta = M_ONE
    case(oct_algorithm_direct)
      call loct_parse_float(datasets_check('OCTEta'), M_ONE, oct%eta)
      call loct_parse_float(datasets_check('OCTDelta'), M_ZERO, oct%delta)
    case(oct_algorithm_newuoa)
#if defined(HAVE_NEWUOA)
      call loct_parse_float(datasets_check('OCTEta'), M_ONE, oct%eta)
      call loct_parse_float(datasets_check('OCTDelta'), M_ZERO, oct%delta)
#else
      write(message(1), '(a)') '"OCTScheme = oct_algorithm_newuoa" is only possible if the newuoa'
      write(message(2), '(a)') 'code has been compiled. You must configure octopus passing the'
      write(message(3), '(a)') 'the "--enable-newuoa" switch.'
      call write_fatal(3)
#endif
    case default
      oct%delta = M_ONE; oct%eta = M_ONE
    end select

    !%Variable OCTDoubleCheck
    !%Type logical
    !%Section Calculation Modes::Optimal Control
    !%Default true
    !%Description 
    !% Run a normal propagation after the optimization using the optimized field.
    !%End
    call loct_parse_logical(datasets_check('OCTDoubleCheck'), .true., oct%oct_double_check)

    !%Variable OCTMixing
    !%Type logical
    !%Section Calculation Modes::Optimal Control
    !%Default false
    !%Description 
    !% Use mixing algorithms to create the input fields in the iterative OCT schemes.
    !% Note that this idea is still a little bit experimental, and depending on the
    !% kind of mixing that you use, and the parameters that you set, it may or may
    !% not accelerate the convergence, or even spoil the convergence.
    !%
    !% Using "TypeOfMixing = broyden", "Mixing = 0.1" and "MixNumberSteps = 3" seems
    !% to work in many cases, but your mileage may vary.
    !%End
    call loct_parse_logical(datasets_check('OCTMixing'), .false., oct%use_mixing)

    !%Variable OCTDirectStep
    !%Type float
    !%Section Calculation Modes::Optimal Control
    !%Default 0.25
    !%Description 
    !%
    !%End
    call loct_parse_float(datasets_check('OCTDirectStep'), CNST(0.25), oct%direct_step)

    !%Variable OCTDumpIntermediate
    !%Type logical
    !%Section Calculation Modes::Optimal Control
    !%Default true
    !%Description 
    !% Writes to disk some data during the OCT algorithm at intermediate steps.
    !% This is rather technical and it should be considered only for debugging
    !% purposes. Nevertheless, since the whole OCT infrastructure is at a very
    !% preliminary developing stage, it is set to true by default.
    !%End
    call loct_parse_logical(datasets_check('OCTDumpIntermediate'), .true., oct%dump_intermediate)

    !%Variable OCTNumberCheckPoints
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default 0
    !%Description 
    !% During an OCT propagation, the code may write down at some time steps (the
    !% "check points", the wavefunctions. When the inverse backward or forward propagation
    !% is performed in a following step, the wave function should reverse its path
    !% (almost) exactly. This can be checked to make sure that it is the case -- otherwise
    !% one should try reducing the time-step, or altering in some other way the
    !% variables that control the propagation.
    !%
    !% If the backward (or forward) propagation is not retracing the steps of the previous
    !% forward (or backward) propation, the code will emit a warning.
    !%End
    call loct_parse_int(datasets_check('OCTNumberCheckPoints'), 0, oct%number_checkpoints)

    !%Variable OCTRandomInitialGuess
    !%Type logical
    !%Section Calculation Modes::Optimal Control
    !%Default false
    !%Description 
    !%
    !%End
    call loct_parse_logical(datasets_check('OCTRandomInitialGuess'), .false., oct%random_initial_guess)

    call pop_sub()
  end subroutine oct_read_inp


  ! ---------------------------------------------------------
  ! Tries to avoid ill defined combinations of run modes.
  ! ---------------------------------------------------------
  subroutine check_faulty_runmodes(sys, h, tr)
    type(system_t),                 intent(in)    :: sys
    type(hamiltonian_t),            intent(in)    :: h
    type(td_rti_t),                 intent(in)    :: tr

    integer :: no_electrons, n_filled, n_partially_filled, n_half_filled
    call push_sub('read.check_faulty_runmodes')

    ! Only dipole approximation in length gauge.
    if(h%gauge.ne.LENGTH) then
      write(message(1),'(a)') "So far only length gauge is supported in optimal control runs."
      call write_fatal(1)
    end if

    ! No QOCT runs with periodic boundary conditions.
    if(simul_box_is_periodic(sys%gr%sb)) then
      write(message(1), '(a)') 'No QOCT runs with periodic boundary conditions. '
      call write_fatal(1)
    end if

    ! This should check that we really have occupation one for
    ! one of the spin-orbitals, and occupation zero for all the others.
    ! Otherwise the algorithms are bound to fail.
    select case(sys%st%d%ispin)
    case(UNPOLARIZED)
      call occupied_states(sys%st, 1, n_filled, n_partially_filled, n_half_filled)
      no_electrons = 2*n_filled + n_half_filled
      if(n_partially_filled > 0 ) then
        write(message(1),'(a)') 'No partially filled orbitals are allowd in OCT calculations'
        call write_fatal(1)
      end if
    case(SPIN_POLARIZED)
      call occupied_states(sys%st, 1, n_filled, n_partially_filled, n_half_filled)
      if(n_partially_filled > 0 .or. n_half_filled > 0) then
        write(message(1),'(a)') 'No partially filled orbitals are allowd in OCT calculations'
        call write_fatal(1)
      end if
      no_electrons = n_filled
      call occupied_states(sys%st, 2, n_filled, n_partially_filled, n_half_filled)
      no_electrons = n_filled + no_electrons
      if(n_partially_filled > 0 .or. n_half_filled > 0) then
        write(message(1),'(a)') 'No partially filled orbitals are allowd in OCT calculations'
        call write_fatal(1)
      end if
    case(SPINORS)
      call occupied_states(sys%st, 1, n_filled, n_partially_filled, n_half_filled)
      no_electrons = n_filled
      if(n_partially_filled > 0 .or. n_half_filled > 0) then
        write(message(1),'(a)') 'No partially filled orbitals are allowd in OCT calculations'
        call write_fatal(1)
      end if
    end select

    if(abs(sys%st%qtot - real(no_electrons, REAL_PRECISION) ) > CNST(1.0e-8)) then
      write(message(1), '(a)') 'Error inv check_faulty_runmodes'
      call write_fatal(1)
    end if

    if((oct%mode_fixed_fluence)) then
      if( (oct%algorithm.ne.oct_algorithm_wg05)     .and. &
          (oct%algorithm.ne.oct_algorithm_str_iter) .and. &
          (oct%algorithm.ne.oct_algorithm_direct)   .and. &
          (oct%algorithm.ne.oct_algorithm_newuoa) ) then
        write(message(1),'(a)') "Cannot optimize to a given fluence with the chosen algorithm."
        write(message(2),'(a)') "Switching to scheme WG05."         
        call write_info(2)
        oct%algorithm = oct_algorithm_wg05
      end if
    else
      if( oct%algorithm.eq.oct_algorithm_direct) then
        write(message(1),'(a)') 'The direct QOCT optimization can be only done in fixed-fluence '
        write(message(2),'(a)') 'mode (i.e. use "OCTFixFluenceTo" input variable).'
        call write_fatal(2)
      end if
    end if

    if(oct%algorithm .eq. oct_algorithm_zbr98) then
      if( (target_type(target) .ne. oct_tg_groundstate) .and. &
          (target_type(target) .ne. oct_tg_gstransformation) ) then
        write(message(1), '(a)') 'The scheme "OCTScheme = oct_algorithm_zbr98 can only be used if'
        write(message(2), '(a)') 'the target state is "OCTTargetOperator = oct_tg_gstransformation"'
        write(message(3), '(a)') 'or "OCTTargetOperator = oct_tg_groundstate".'
        call write_fatal(3)
      end if
    end if
      
    ! Filters only with the WG05 scheme.
    if(filter_number(filter).ne.0) then
      if(oct%algorithm .ne. oct_algorithm_wg05) then
        write(message(1), '(a)') 'Filters can only be used with the WG05 QOCT algorithm.'
        call write_fatal(1)
      end if
    end if
      
    ! local targets only in ZR98 and WG05
    if(target_type(target) .eq. oct_tg_local .or. &
       target_type(target) .eq. oct_tg_density .or. &
       target_type(target) .eq. oct_tg_td_local) then
      if(oct%algorithm .eq. oct_algorithm_zbr98) then
        write(message(1), '(a)') 'Cannot use ZBR98 OCT scheme if the target is oct_tg_density,'
        write(message(2), '(a)') 'oct_tg_local or oct_tg_td_local.'
        call write_fatal(2)
      end if
    end if

    if(target_type(target) .eq. oct_tg_td_local) then
      if(tr%method .ne. PROP_CRANK_NICHOLSON) then
        write(message(1), '(a)') 'If OCTTargetMode = oct_tg_td_local, then you must set'
        write(message(2), '(a)') 'TDEvolutionMethod = crank_nicholson'
        call write_fatal(2)
      end if
    end if

    if(target_type(target) .eq. oct_tg_excited) then
      if(sys%st%d%ispin .eq. UNPOLARIZED) then
        write(message(1), '(a)') 'If OCTTargetMode = oct_tg_excited, then you must run either with'
        write(message(1), '(a)') 'SpinComponents = spin_polarized or SpinComponents = spinors.'
        call write_fatal(2)
      end if
    end if

    if( h%theory_level.ne.INDEPENDENT_PARTICLES ) then
      if(h%theory_level.ne.KOHN_SHAM_DFT) then
        write(message(1), '(a)') 'In optimal control theory mode, you can only use either independent'
        write(message(2), '(a)') 'particles "TheoryLevel = independent_particles", or Kohn-Sham DFT'
        write(message(3), '(a)') '"TheoryLevel = dft".'
        call write_fatal(3)
      end if
      if( tr%method .ne. PROP_EXPONENTIAL_MIDPOINT ) then
        write(message(1), '(a)') 'When doing QOCT with interacting electronsl, then you must set'
        write(message(2), '(a)') 'TDEvolutionMethod = exp_mid'
        call write_fatal(2)
      end if
    end if

    if( h%ab .eq. MASK_ABSORBING) then
      write(message(1), '(a)') 'Cannot do QOCT with mask absorbing boundaries. Use either'
      write(message(2), '(a)') '"AbsorbingBoudaries = sin2" or "AbsorbingBoundaries = no".'
      call write_fatal(2)
    end if

    if(target_type(target) .eq. oct_tg_exclude_state ) then
      if(no_electrons > 1) then
        write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_exclude_state", you can only do'
        write(message(2), '(a)') 'one-electron runs.'
        call write_fatal(2)
      end if
      if(sys%st%d%ispin .eq. SPIN_POLARIZED) then
        write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_exclude_state", you can only do'
        write(message(2), '(a)') 'runs in spin restricted, or in spinors mode (spin-polarized is'
        write(message(3), '(a)') 'is not allowed.'
        call write_fatal(3)
      end if
    end if

    if(oct%mode_basis_set) then
      if( (oct%algorithm .ne. oct_algorithm_str_iter)  .and.  &
          (oct%algorithm .ne. oct_algorithm_direct)    .and.  &
          (oct%algorithm .ne. oct_algorithm_newuoa) ) then
        write(message(1), '(a)') 'If the control parameters are to be represented with a basis set'
        write(message(2), '(a)') '(i.e. "OCTParameterRepresentation = control_parameters_fourier_space"),'
        write(message(3), '(a)') 'the only acceptable algorithms are:'
        write(message(4), '(a)') ' (1) "OCTScheme = oct_algorithm_straight_iteration".'
        write(message(5), '(a)') ' (2) "OCTScheme = oct_algorithm_direct".'
        write(message(6), '(a)') ' (3) "OCTScheme = oct_algorithm_newuoa".'
        call write_fatal(6)
      endif
    end if

    if(oct%algorithm .eq. oct_algorithm_direct) then
      if(.not.oct%mode_basis_set) then
        write(message(1), '(a)') 'If you want to use "OCTScheme = oct_algorithm_direct", then you'
        write(message(2), '(a)') 'must represent the control parameters with a basis set (i.e.'
        write(message(3), '(a)') '"OCTParameterRepresentation = control_parameters_fourier_space"'
        call write_fatal(3)
      end if
    end if
      
    call pop_sub()      
  end subroutine check_faulty_runmodes


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

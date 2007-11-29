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
  subroutine oct_read_inp(oct)
    type(oct_t), intent(inout) :: oct


    call push_sub('opt_control_read.oct_read_inp')  

    !%Variable OCTFixFluenceTo
    !%Type float
    !%Section Optimal Control
    !%Default -1.0
    !%Description
    !% The algorithm tries to obtain the specified fluence for the laser field. 
    !% This works only in conjunction with the WG05 scheme.
    !%End
    oct%mode_fixed_fluence = .false.
    call loct_parse_float(check_inp('OCTFixFluenceTo'), -M_ONE, oct%targetfluence)
    if (oct%targetfluence.ne.-M_ONE) oct%mode_fixed_fluence = .true.


    if(oct%mode_fixed_fluence) then
      !%Variable OCTFixInitialFluence
      !%Type logical
      !%Section Optimal Control
      !%Default yes
      !%Description
      !% By default, when asking for a fixed-fluence optimization ("OCTFixFluenceTo = whatever"), 
      !% the initial laser guess provided in the input file is scaled to match this
      !% fluence. However, you can force the program to use that initial laser as the initial
      !% guess, no matter the fluence, by setting "OCTFixInitialFluence = no".
      !%End
      call loct_parse_logical(check_inp('OCTFixInitialFluence'), .true., oct%fix_initial_fluence)
    else
      oct%fix_initial_fluence = .false.
    end if

      
    !%Variable OCTScheme
    !%Type integer
    !%Section Optimal Control
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
    !% Yet to be implemented and tested. Basically an improved and generalized scheme. 
    !% Comparable to ZBR98/ZR98. See [Y. Maday and G. Turinici, J. Chem. Phys. 118, 
    !% 8191 (2003)].
    !%Option oct_algorithm_krotov 5
    !% Yet to be implemented and tested. Reported in [D. Tannor, V. Kazakov and V.
    !% Orlov, in "Time Dependent Quantum Molecular Dynamics", edited by J. Broeckhove
    !% and L. Lathouweres (Plenum, New York, 1992), pp. 347-360].
    !%End
    call loct_parse_int(check_inp('OCTScheme'), oct_algorithm_zr98, oct%algorithm_type)
    if(.not.varinfo_valid_option('OCTScheme', oct%algorithm_type)) call input_error('OCTScheme')
    select case(oct%algorithm_type)
    case(oct_algorithm_mt03)
      oct%delta = M_TWO; oct%eta = M_ZERO
    case(oct_algorithm_zr98)
      oct%delta = M_ONE; oct%eta = M_ONE
    case(oct_algorithm_krotov)
      oct%delta = M_ONE; oct%eta = M_ZERO
    case default
      oct%delta = M_ONE; oct%eta = M_ONE
    end select

    !%Variable OCTZBR98zero_iteration
    !%Type logical
    !%Section Optimal Control
    !%Default true
    !%Description 
    !% No documentation is available.
    !%End
    call loct_parse_logical(check_inp('OCTZBR98zero_iteration'),.true., oct%zbr98_zero_iteration)

    !%Variable OCTDoubleCheck
    !%Type logical
    !%Section Optimal Control
    !%Default true
    !%Description 
    !% Run a normal propagation after the optimization using the optimized field.
    !%End
    call loct_parse_logical(check_inp('OCTDoubleCheck'), .true., oct%oct_double_check)

    !%Variable OCTMixing
    !%Type logical
    !%Section Optimal Control
    !%Default false
    !%Description 
    !% Use mixing algorithms to create the input fields in the iterative OCT schemes.
    !%
    !% WARNING: Very experimental.
    !%
    !% WARNING: The oct_algorithm_wg05 scheme is not affected by this.
    !%End
    call loct_parse_logical(check_inp('OCTMixing'), .false., oct%use_mixing)

    !%Variable OCTDumpIntermediate
    !%Type logical
    !%Section Optimal Control
    !%Default true
    !%Description 
    !% Writes to disk some data during the OCT algorithm at intermediate steps.
    !% This is rather technical and it should be considered only for debugging
    !% purposes. Nevertheless, since the whole OCT infrastructure is at a very
    !% preliminary developing stage, it is set to true by default, and in fact
    !% all the intermediate information is printed always.
    !%End
    call loct_parse_logical(check_inp('OCTDumpIntermediate'), .true., oct%dump_intermediate)

    !%Variable OCTNumberCheckPoints
    !%Type integer
    !%Section Optimal Control
    !%Default 0
    !%Description 
    !% WARNING: EXPERIMENTAL, description missing.
    !%End
    call loct_parse_int(check_inp('OCTNumberCheckPoints'), 0, oct%number_checkpoints)

    call pop_sub()
  end subroutine oct_read_inp


  ! ---------------------------------------------------------
  ! Tries to avoid ill defined combinations of run modes.
  ! ---------------------------------------------------------
  subroutine check_faulty_runmodes(oct, sys, h, target, tr)
    type(oct_t), intent(inout) :: oct
    type(system_t), target, intent(in) :: sys
    type(hamiltonian_t),    intent(in) :: h
    type(target_t),         intent(in) :: target
    type(td_rti_t),         intent(in) :: tr

    integer :: no_electrons, n_filled, n_partially_filled, n_half_filled
    call push_sub('read.check_faulty_runmodes')

    ! Only dipole approximation in length gauge.
    if(h%gauge.ne.LENGTH) then
      write(message(1),'(a)') "So far only length gauge is supported in optimal control runs."
      call write_fatal(1)
    end if

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

    ! FixedFluence and Filter only work with WG05
    if((oct%mode_fixed_fluence).and.(oct%algorithm_type.ne.oct_algorithm_wg05)) then
      write(message(1),'(a)') "Cannot optimize to a given fluence with the chosen algorithm."
      write(message(2),'(a)') "Switching to scheme WG05."         
      call write_info(2)
      oct%algorithm_type = oct_algorithm_wg05
    end if

    if(oct%algorithm_type .eq. oct_algorithm_zbr98) then
      if( (target%type .ne. oct_tg_groundstate) .and. &
          (target%type .ne. oct_tg_gstransformation) ) then
        write(message(1), '(a)') 'The scheme "OCTScheme = oct_algorithm_zbr98 can only be used if'
        write(message(2), '(a)') 'the target state is "OCTTargetOperator = oct_tg_gstransformation"'
        write(message(3), '(a)') 'or "OCTTargetOperator = oct_tg_groundstate".'
        call write_fatal(3)
      end if
    end if
      
    ! WARNING filters can only be used with WG05, and this is not checked any more.
      
    ! local targets only in ZR98 and WG05
    if(target%type .eq. oct_tg_local .or. &
       target%type .eq. oct_tg_density .or. &
       target%type .eq. oct_tg_td_local) then
      if(oct%algorithm_type .eq. oct_algorithm_zbr98) then
        write(message(1), '(a)') 'Cannot use ZBR98 OCT scheme if the target is oct_tg_density,'
        write(message(2), '(a)') 'oct_tg_local or oct_tg_td_local.'
        call write_fatal(2)
      end if
    end if

    if(target%type .eq. oct_tg_td_local) then
      if(tr%method .ne. PROP_CRANK_NICHOLSON) then
        write(message(1), '(a)') 'If OCTTargetMode = oct_tg_td_local, then you must set'
        write(message(2), '(a)') 'TDEvolutionMethod = crank_nicholson'
        call write_fatal(2)
      end if
    end if

    if(target%type .eq. oct_tg_excited) then
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

    if(target%type .eq. oct_tg_exclude_state ) then
      if(no_electrons > 1) then
        write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_exclude_state", you can only do'
        write(message(2), '(a)') 'one-electron runs.'
        call write_fatal(2)
      end if
      if(sys%st%d%ispin .eq. SPIN_POLARIZED) then
        write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_exclude_state", you can only do'
        write(message(2), '(a)') 'runs in spin restricted, or in spinors mode (spin-polarized is'
        writE(message(3), '(a)') 'is not allowed.'
        call write_fatal(3)
      end if
    end if
      
    call pop_sub()      
  end subroutine check_faulty_runmodes


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

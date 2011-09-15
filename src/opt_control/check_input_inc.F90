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
  ! Tries to avoid ill-defined combinations of run modes.
  ! ---------------------------------------------------------
  subroutine check_faulty_runmodes(sys, hm, tr)
    type(system_t),      intent(in) :: sys
    type(hamiltonian_t), intent(in) :: hm
    type(propagator_t),  intent(in) :: tr

    integer :: no_electrons, n_filled, n_partially_filled, n_half_filled

    PUSH_SUB(check_faulty_runmodes)

    ! No QOCT runs with periodic boundary conditions.
    if(simul_box_is_periodic(sys%gr%sb)) then
      write(message(1), '(a)') 'No QOCT runs with periodic boundary conditions. '
      call messages_fatal(1)
    end if

    ! This should check that we really have occupation one for
    ! one of the spin-orbitals, and occupation zero for all the others.
    ! Otherwise the algorithms are bound to fail.
    select case(sys%st%d%ispin)
    case(UNPOLARIZED)
      call occupied_states(sys%st, 1, n_filled, n_partially_filled, n_half_filled)
      no_electrons = 2*n_filled + n_half_filled
      if(n_partially_filled > 0 ) then
        write(message(1),'(a)') 'No partially filled orbitals are allowed in OCT calculations.'
        call messages_fatal(1)
      end if
    case(SPIN_POLARIZED)
      call occupied_states(sys%st, 1, n_filled, n_partially_filled, n_half_filled)
      if(n_partially_filled > 0 .or. n_half_filled > 0) then
        write(message(1),'(a)') 'No partially filled orbitals are allowed in OCT calculations.'
        call messages_fatal(1)
      end if
      no_electrons = n_filled
      call occupied_states(sys%st, 2, n_filled, n_partially_filled, n_half_filled)
      no_electrons = n_filled + no_electrons
      if(n_partially_filled > 0 .or. n_half_filled > 0) then
        write(message(1),'(a)') 'No partially filled orbitals are allowed in OCT calculations.'
        call messages_fatal(1)
      end if
    case(SPINORS)
      call occupied_states(sys%st, 1, n_filled, n_partially_filled, n_half_filled)
      no_electrons = n_filled
      if(n_partially_filled > 0 .or. n_half_filled > 0) then
        write(message(1),'(a)') 'No partially filled orbitals are allowed in OCT calculations.'
        call messages_fatal(1)
      end if
    end select

    if(abs(sys%st%qtot - real(no_electrons, REAL_PRECISION) ) > CNST(1.0e-8)) then
      write(message(1), '(a)') 'Error inv check_faulty_runmodes'
      call messages_fatal(1)
    end if

    if(oct%algorithm .eq. oct_algorithm_zbr98) then
      select case(target_type(target))
      case(oct_tg_groundstate, oct_tg_gstransformation, &
           oct_tg_userdefined)
      case default
        write(message(1), '(a)') 'The scheme "OCTScheme = oct_algorithm_zbr98 can only be used if'
        write(message(2), '(a)') 'the target state is "OCTTargetOperator = oct_tg_gstransformation"'
        write(message(3), '(a)') 'or "OCTTargetOperator = oct_tg_groundstate"'
        write(message(4), '(a)') 'or "OCTTargetOperator = oct_tg_userdefined".'
        call messages_fatal(4)
      end select 
    end if
      
    ! Filters only with the WG05 scheme.
    if(filter_number(filter).ne.0) then
      if(oct%algorithm .ne. oct_algorithm_wg05) then
        write(message(1), '(a)') 'Filters can only be used with the WG05 QOCT algorithm.'
        call messages_fatal(1)
      end if
    end if
      
    ! local targets only in ZR98 and WG05
    if(target_type(target) .eq. oct_tg_local .or. &
       target_type(target) .eq. oct_tg_density .or. &
       target_type(target) .eq. oct_tg_td_local) then
      if(oct%algorithm .eq. oct_algorithm_zbr98) then
        write(message(1), '(a)') 'Cannot use ZBR98 OCT scheme if the target is oct_tg_density,'
        write(message(2), '(a)') 'oct_tg_local or oct_tg_td_local.'
        call messages_fatal(2)
      end if
    end if
    
    ! the inh term in the bwd evolution of chi is taken into
    ! consideration only for certain propagators
    if(.not.oct_algorithm_is_direct(oct)) then
      if(target_mode(target) .eq. oct_targetmode_td) then
        select case(tr%method)
        case(PROP_CRANK_NICHOLSON)
       ! for the moment exp mid point with lanczos is broken. needs to be fixed.  
       ! case(PROP_EXPONENTIAL_MIDPOINT)
       !   if(tr%te%exp_method .ne. EXP_LANCZOS) then
       !     WRITE(message(1), '(a)') 'If you use time-dependent target, and you set'
       !     write(message(2), '(a)') '"TDPropagator = exp_mid", then you must set'
       !     write(message(3), '(a)') '"TDExponentialMethod = lanczos".'
       !     call messages_fatal(3)
       !   end if
        case default
          write(message(1), '(a)') 'If you use time-dependent target, then you must set'
          write(message(2), '(a)') '"TDPropagator = crank_nicholson"'
       !   write(message(3), '(a)') '"TDPropagator = exp_mid".'
          call messages_fatal(3)
        end select
      end if
    end if

    if(target_type(target) .eq. oct_tg_excited) then
      if(sys%st%d%ispin .eq. UNPOLARIZED) then
        write(message(1), '(a)') 'If OCTTargetMode = oct_tg_excited, then you must run either with'
        write(message(1), '(a)') 'SpinComponents = spin_polarized or SpinComponents = spinors.'
        call messages_fatal(2)
      end if
    end if

    if( hm%theory_level.ne.INDEPENDENT_PARTICLES ) then
      if(hm%theory_level.ne.KOHN_SHAM_DFT) then
        write(message(1), '(a)') 'In optimal control theory mode, you can only use either independent'
        write(message(2), '(a)') 'particles "TheoryLevel = independent_particles", or Kohn-Sham DFT'
        write(message(3), '(a)') '"TheoryLevel = dft".'
        call messages_fatal(3)
      end if
      if( (tr%method .ne. PROP_QOCT_TDDFT_PROPAGATOR) .and. &
          (tr%method .ne. PROP_QOCT_TDDFT_PROPAGATOR_2) )  then
        if( .not. oct_algorithm_is_direct(oct) ) then
          write(message(1), '(a)') 'When doing QOCT with interacting electrons, then you must set'
          write(message(2), '(a)') 'TDPropagator = qoct_tddft_propagator'
          call messages_fatal(2)
        end if
      end if
    end if

    if( hm%ab .eq. MASK_ABSORBING) then
      if( (oct%algorithm .ne. oct_algorithm_direct) .and. &
          (oct%algorithm .ne. oct_algorithm_newuoa) ) then
        write(message(1), '(a)') 'Cannot do QOCT with mask absorbing boundaries. Use either'
        write(message(2), '(a)') '"AbsorbingBoudaries = sin2" or "AbsorbingBoundaries = no".'
        call messages_fatal(2)
      end if
    end if

    if(target_type(target) .eq. oct_tg_exclude_state ) then
      if(no_electrons > 1) then
        write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_exclude_state", you can only do'
        write(message(2), '(a)') 'one-electron runs.'
        call messages_fatal(2)
      end if
      if(sys%st%d%ispin .eq. SPIN_POLARIZED) then
        write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_exclude_state", you can only do'
        write(message(2), '(a)') 'runs in spin restricted, or in spinors mode (spin-polarized is'
        write(message(3), '(a)') 'is not allowed.'
        call messages_fatal(3)
      end if
    end if
    
    if(target_type(target) .eq. oct_tg_velocity) then
       if( (oct%algorithm .ne. oct_algorithm_direct) .and. &
            (oct%algorithm .ne. oct_algorithm_newuoa) .and. & 
            (oct%algorithm .ne. oct_algorithm_cg) ) then
          write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_velocity", you can only use'
          write(message(2), '(a)') '"OCTScheme = oct_algorithm_direct" or'
          write(message(3), '(a)') '"OCTScheme = oct_algorithm_newuoa" or'
          write(message(4), '(a)') '"OCTScheme = oct_algorithm_cg" for the optimization.'
          call messages_fatal(4)
       end if
       if((oct%algorithm .eq. oct_algorithm_cg) .and. target_move_ions(target)) then
          write(message(1), '(a)') 'If "OCTTargetOperator = oct_tg_velocity", and'
          write(message(2), '(a)') '"OCTScheme = oct_algorithm_cg", then you have to'
          write(message(3), '(a)') 'set "OCTMoveIons = false"'
          call messages_fatal(3)
       end if
    end if
    
    if(target_curr_functional(target) .ne. oct_no_curr) then
      select case(sys%st%d%ispin)
      case(UNPOLARIZED)
      case(SPIN_POLARIZED)
        message(1) = 'Spin_polarized! Do not use OCT current functionals.'
        call messages_fatal(1)
      case(SPINORS)
        message(1) = 'Spinors! Do not use OCT current functionals.'
        call messages_fatal(1)
      end select
    end if

    select case(controlfunction_mode())
    case(controlfunction_mode_f, controlfunction_mode_phi)
      if(.not. oct_algorithm_is_direct(oct)) then
        message(1) = 'If you attempt an envelope-only or phase-only optimization, then'
        message(2) = 'you must use a gradient-free algorithm.'
        call messages_fatal(2)
      end if
    end select
      
    POP_SUB(check_faulty_runmodes)
  end subroutine check_faulty_runmodes


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

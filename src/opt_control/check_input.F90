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
  ! Tries to avoid ill defined combinations of run modes.
  ! ---------------------------------------------------------
  subroutine check_faulty_runmodes(sys, hm, tr)
    type(system_t),                 intent(in)    :: sys
    type(hamiltonian_t),            intent(in)    :: hm
    type(td_rti_t),                 intent(in)    :: tr

    integer :: no_electrons, n_filled, n_partially_filled, n_half_filled
    call push_sub('read.check_faulty_runmodes')

    ! Only dipole approximation in length gauge.
    if(hm%gauge.ne.LENGTH) then
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

    if( hm%theory_level.ne.INDEPENDENT_PARTICLES ) then
      if(hm%theory_level.ne.KOHN_SHAM_DFT) then
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

    if( hm%ab .eq. MASK_ABSORBING) then
      if( (oct%algorithm .ne. oct_algorithm_direct) .and. &
          (oct%algorithm .ne. oct_algorithm_newuoa) ) then
        write(message(1), '(a)') 'Cannot do QOCT with mask absorbing boundaries. Use either'
        write(message(2), '(a)') '"AbsorbingBoudaries = sin2" or "AbsorbingBoundaries = no".'
        call write_fatal(2)
      end if
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

    call pop_sub()      
  end subroutine check_faulty_runmodes


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

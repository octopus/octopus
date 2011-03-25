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

! ---------------------------------------------------------
subroutine td_init(td, sys, hm)
  type(td_t),            intent(inout) :: td
  type(system_t),        intent(inout) :: sys
  type(hamiltonian_t),   intent(inout) :: hm

  integer :: dummy

  PUSH_SUB(td_init)

  call ion_dynamics_init(td%ions, sys%geo)

  td%iter = 0

  !%Variable TDTimeStep
  !%Type float
  !%Default 0.07 a.u.
  !%Section Time-Dependent::Propagation
  !%Description
  !% Time-step for the time propagation. The default is 0.07
  !% [hbar/Hartree].
  !%End
  call parse_float(datasets_check('TDTimeStep'), units_from_atomic(units_inp%time, CNST(0.07)), td%dt)
  td%dt = units_to_atomic(units_inp%time, td%dt)

  if (td%dt <= M_ZERO) then
    write(message(1),'(a,f14.6,a)') "Input: '", td%dt, "' is not a valid TDTimeStep"
    message(2) = '(0 < TDTimeStep)'
    call messages_fatal(2)
  end if

  !%Variable TDIonicTimeScale
  !%Type float
  !%Default 1.0
  !%Section Time-Dependent::Propagation
  !%Description
  !% This variable defines the factor between the timescale of ionic
  !% and electronic movement. It allows reasonably fast
  !% Born-Oppenheimer molecular-dynamics simulations based on
  !% Ehrenfest dynamics. The value of this variable is equivalent to
  !% the role of <math>\mu</math> in Car-Parrinello. Increasing it
  !% linearly accelerates the time step of the ion
  !% dynamics, but also increases the deviation of the system from the
  !% Born-Oppenheimer surface. The default is 1, which means that both
  !% timescales are the same. Note that a value different than 1
  !% implies that the electrons will not follow physical behaviour.
  !%
  !% According to our tests, values around 10 are reasonable, but it
  !% will depend on your system, mainly on the width of the gap.
  !%
  !% Important: The electronic time step will be the value of
  !% <tt>TDTimeStep</tt> divided by this variable, so if you have determined an
  !% optimal electronic time step (that we can call <i>dte</i>), it is
  !% recommended that you define your time step as:
  !%
  !% <tt>TDTimeStep</tt> = <i>dte</i> * <tt>TDIonicTimeScale</tt>
  !%
  !% so you will always use the optimal electronic time step.
  !%
  !% For more details see: <tt>http://arxiv.org/abs/0710.3321</tt>
  !%
  !%End

  call parse_float(datasets_check('TDIonicTimeScale'), CNST(1.0), td%mu)

  !%Variable TDMaximumIter
  !%Type integer
  !%Default 1500
  !%Section Time-Dependent::Propagation
  !%Description
  !% Number of time-propagation steps that will be performed. By default 1500.
  !%
  !% Tip: If you would like to specify the real time of the
  !% propagation, rather than the number of steps, just use something
  !% like:
  !%
  !% <tt>TDMaximumIter</tt> = 1000.0 / <tt>TDTimeStep</tt>
  !%
  !%End 

  call parse_integer(datasets_check('TDMaximumIter'), 1500, td%max_iter)

  if(td%max_iter < 1) then
    write(message(1), '(a,i6,a)') "Input: '", td%max_iter, "' is not a valid TDMaximumIter"
    message(2) = '(TDMaximumIter <= 1)'
    call messages_fatal(2)
  end if
  
  ! Initialize the kick (if optical-spectrum calculations are to be performed)
  call kick_init(td%kick, sys%st%d%nspin, sys%gr%mesh%sb%dim)

  ! now the photoelectron stuff
  call parse_integer(datasets_check('AbsorbingBoundaries'), 0, dummy)
  call PES_init(td%PESv, sys%gr%mesh, sys%gr%sb, sys%st, dummy, sys%outp%iter,hm,td%max_iter,td%dt)

  !%Variable TDDynamics
  !%Type integer
  !%Default ehrenfest
  !%Section Time-Dependent::Propagation
  !%Description
  !% Type of dynamics to follow during a time propagation. By default
  !% it is Ehrenfest TDDFT.
  !%Option ehrenfest 1
  !% Ehrenfest dynamics.
  !%Option bo 2
  !% Born-Oppenheimer (Experimental).
  !%Option cp 3
  !% Car-Parrinello molecular dynamics.
  !%End

  call parse_integer(datasets_check('TDDynamics'), EHRENFEST, td%dynamics)
  if(.not.varinfo_valid_option('TDDynamics', td%dynamics)) call input_error('TDDynamics')
  call messages_print_var_option(stdout, 'TDDynamics', td%dynamics)

  !%Variable RecalculateGSDuringEvolution
  !%Type logical
  !%Default no
  !%Section Time-Dependent::Propagation
  !%Description
  !% In order to calculate some information about the system along the
  !% evolution (e.g. projection onto the ground-state KS determinant,
  !% projection of the TDKS spin-orbitals onto the ground-state KS
  !% spin-orbitals), the ground-state KS orbitals are needed. If the
  !% ionic potential changes -- that is, the ions move -- one may want
  !% to recalculate the ground state. You may do this by setting this
  !% variable.
  !%
  !% The recalculation is not done every time step, but only every
  !% OutputEvery time steps.
  !%End
  call parse_logical(datasets_check("RecalculateGSDuringEvolution"), .false., td%recalculate_gs)

  call messages_obsolete_variable(datasets_check('TDScissor'))

  call propagator_init(sys%gr, sys%st, hm, td%tr, td%dt, td%max_iter, &
       ion_dynamics_ions_move(td%ions) .or. gauge_field_is_applied(hm%ep%gfield))
  if(td%dynamics == BO)  call scf_init(td%scf, sys%gr, sys%geo, sys%st, hm)
  if(hm%ep%no_lasers>0.and.mpi_grp_is_root(mpi_world)) then
    call messages_print_stress(stdout, "Time-dependent external fields")
    call laser_messages_info(hm%ep%lasers, stdout, td%dt, td%max_iter)
    call messages_print_stress(stdout)
  end if

  !%Variable TDEnergyUpdateIter
  !%Type integer
  !%Default 10
  !%Section Time-Dependent::Propagation
  !%Description
  !% This variable controls how often Octopus updates the total energy
  !% during a time propagation run. The default is every 10
  !% iterations. For iterations where the energy is not updated, the
  !% last calculated value is reported. If you set this variable to 1,
  !% the energy will be calculated in each step.
  !%End 

  call parse_integer(datasets_check('TDEnergyUpdateIter'), 10, td%energy_update_iter)

  POP_SUB(td_init)
end subroutine td_init


! ---------------------------------------------------------
subroutine td_end(td)
  type(td_t), intent(inout) :: td

  PUSH_SUB(td_end)

  call PES_end(td%PESv)
  call propagator_end(td%tr)  ! clean the evolution method

  if(td%dynamics == BO) call scf_end(td%scf)
  
  POP_SUB(td_end)
end subroutine td_end

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

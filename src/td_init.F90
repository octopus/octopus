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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

! ---------------------------------------------------------
subroutine td_init(gr, td, st, outp)
  type(td_t),     intent(out)   :: td
  type(grid_t),   intent(inout) :: gr
  type(states_t), intent(inout) :: st
  type(output_t), intent(in)    :: outp

  integer :: dummy

  call push_sub('td_init.td_init')

  td%iter = 0

  !%Variable TDTimeStep
  !%Type float
  !%Default 0.07 a.u.
  !%Section Time Dependent::Propagation
  !%Description
  !% Time-step for the propagation;
  !% in previous notation, <math>\delta t</math>.
  !%End
  call loct_parse_float(check_inp('TDTimeStep'), CNST(0.07)/units_inp%time%factor, td%dt)
  td%dt = td%dt * units_inp%time%factor
  if (td%dt <= M_ZERO) then
    write(message(1),'(a,f14.6,a)') "Input: '", td%dt, "' is not a valid TDTimeStep"
    message(2) = '(0 < TDTimeStep)'
    call write_fatal(2)
  end if

  !%Variable TDMaximumIter
  !%Type integer
  !%Default 1500
  !%Section Time Dependent::Propagation
  !%Description
  !% Number of time steps in which the total integration time is divided;
  !% in previous notation, <i>N</i>.
  !%End
  call loct_parse_int(check_inp('TDMaximumIter'), 1500, td%max_iter)
  if(td%max_iter < 1) then
    write(message(1), '(a,i6,a)') "Input: '", td%max_iter, "' is not a valid TDMaximumIter"
    message(2) = '(1 <= TDMaximumIter)'
    call write_fatal(2)
  end if

  ! Initialize the kick (if optical spectrum calculations are to be performed)
  call kick_init(td%kick, st%d%nspin)

  ! now the photoelectron stuff
#if !defined(DISABLE_PES) && defined(HAVE_FFT)
  call loct_parse_int(check_inp('AbsorbingBoundaries'), 0, dummy)
  call PES_init(td%PESv, gr%m, gr%sb, st, dummy, outp%iter)
#endif

  !%Variable MoveIons
  !%Type integer
  !%Default static_ions
  !%Section Time Dependent::Propagation
  !%Description
  !% What kind of simulation to perform.
  !%Option static_ions 0
  !% Do not move the ions.
  !%Option verlet 3
  !% Newtonian dynamics using Verlet.
  !%Option vel_verlet 4
  !% Newtonian dynamics using velocity Verlet.
  !%End
  call loct_parse_int(check_inp('MoveIons'), STATIC_IONS, td%move_ions)
  if(.not.varinfo_valid_option('MoveIons', td%move_ions)) call input_error('MoveIons')
  call messages_print_var_option(stdout, 'MoveIons', td%move_ions)

  if( td%move_ions .eq. NORMAL_VERLET) then
    write(message(1),'(a)') "Normal Verlet algorithm temporarily disabled."
    write(message(2),'(a)') "Using Velocity Verlet (MoveIons = 4) instead."
    call write_warning(2)
    td%move_ions = VELOCITY_VERLET
  end if

  !%Variable RecalculateGSDuringEvolution
  !%Type logical
  !%Default no
  !%Section Time Dependent::Propagation
  !%Description
  !% In order to calculate some information about the system along the
  !% evolution (e.g. projection onto the ground-state KS determinant,
  !% projection of the TDKS spin-orbitals onto the grond-state KS
  !% spin-orbitals), the ground-state KS orbitals are needed. If the
  !% ionic potential changes -- that is, the ions move --, one may want
  !% to recalculate the ground state. You may do this by setting this
  !% variable.
  !%
  !% The recalculation is not done every time step, but only every
  !% OutputEvery time steps.
  !%End
  call loct_parse_logical(check_inp("RecalculateGSDuringEvolution"), .false., td%recalculate_gs)

  call loct_parse_int(check_inp('TDFastEpotGeneration'), 1, td%epot_regenerate)
  if(td%epot_regenerate < 1) td%epot_regenerate = 1

  call td_rti_init(gr, st, td%tr)

  call pop_sub()
end subroutine td_init


! ---------------------------------------------------------
subroutine td_end(td)
  type(td_t), intent(inout) :: td

  call push_sub('td_init.td_end')

#if !defined(DISABLE_PES) && defined(HAVE_FFT)
  call PES_end(td%PESv)
#endif
  call td_rti_end(td%tr)  ! clean the evolution method

  call pop_sub()
end subroutine td_end

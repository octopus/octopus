!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
subroutine td_init(gr, td, st, outp)
  type(td_type),          intent(out)   :: td
  type(grid_type),        intent(inout) :: gr
  type(states_type),      intent(inout) :: st
  type(output_type),      intent(in)    :: outp

  integer :: dummy

  call push_sub('td_init.td_init')

  td%iter = 0

  !%Variable TDTimeStep
  !%Type float
  !%Default 0.07 a.u.
  !%Section Time Dependent
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
  !%Section Time Dependent
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
  !%Type logical
  !%Default static_ions
  !%Section Time Dependent
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

  if( td%move_ions .eq. NORMAL_VERLET) then
    write(message(1),'(a)') "Normal Verlet algorithm temporarily disabled."
    write(message(2),'(a)') "Using Velocity Verlet (MoveIons = 4) instead."
    call write_warning(2)
    td%move_ions = VELOCITY_VERLET
  end if

  call loct_parse_int(check_inp('TDFastEpotGeneration'), 1, td%epot_regenerate)
  if(td%epot_regenerate < 1) td%epot_regenerate = 1

  call td_rti_init(gr, st, td%tr)

  call pop_sub()
end subroutine td_init


! ---------------------------------------------------------
subroutine td_end(td)
  type(td_type), intent(inout) :: td

  call push_sub('td_init.td_end')

#if !defined(DISABLE_PES) && defined(HAVE_FFT)
  call PES_end(td%PESv)
#endif
  call td_rti_end(td%tr)  ! clean the evolution method

  call pop_sub()
end subroutine td_end

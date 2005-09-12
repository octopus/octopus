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
subroutine td_init(gr, td, st, h, outp)
  type(td_type),          intent(out)   :: td
  type(grid_type),        intent(inout) :: gr
  type(states_type),      intent(inout) :: st
  type(hamiltonian_type), intent(in)    :: h
  type(output_type),      intent(in)    :: outp

  integer :: i, j, n, dummy
  integer(POINTER_SIZE) :: blk

  call push_sub('td_init.td_init')

  td%iter = 0

  call loct_parse_float(check_inp('TDTimeStep'), CNST(0.07)/units_inp%time%factor, td%dt)
  td%dt = td%dt * units_inp%time%factor
  if (td%dt <= M_ZERO) then
    write(message(1),'(a,f14.6,a)') "Input: '", td%dt, "' is not a valid TDTimeStep"
    message(2) = '(0 < TDTimeStep)'
    call write_fatal(2)
  end if

  call loct_parse_int(check_inp('TDMaximumIter'), 1500, td%max_iter)
  if(td%max_iter < 1) then
    write(message(1), '(a,i6,a)') "Input: '", td%max_iter, "' is not a valid TDMaximumIter"
    message(2) = '(1 <= TDMaximumIter)'
    call write_fatal(2)
  end if

  call loct_parse_float(check_inp('TDDeltaStrength'), M_ZERO, td%delta_strength)
  ! units are 1/length
  td%delta_strength = td%delta_strength / units_inp%length%factor

  !%Variable TDDeltaStrengthMode
  !%Type integer
  !%Section 10 Time Dependent
  !%Description
  !% When calculating the linear response of the density via the propagation
  !% in real time, one needs to perfrom an initical kick on the KS system, at 
  !% time zero. Depending on what kind response property one wants to obtain,
  !% this kick may be done in several modes.
  !%Option kick_density 0
  !% The total density of the system is perturbed.
  !%Option kick_spin 1
  !% The individual spin densities are perturbed differently. Note that this mode
  !% is only possible if the run is done in spin polarized mode, or with spinors.
  !%Option kick_spin_and_density 2
  !% A combination of the two above. Note that this mode
  !% is only possible if the run is done in spin polarized mode, or with spinors.
  !%End
  call loct_parse_int(check_inp('TDDeltaStrengthMode'), KICK_DENSITY_MODE, td%delta_strength_mode)
  select case (td%delta_strength_mode)
    case (KICK_DENSITY_MODE)
    case (KICK_SPIN_MODE, KICK_SPIN_DENSITY_MODE)
      if (st%d%ispin == UNPOLARIZED) call input_error('TDDeltaStrengthMode')
    case default
      call input_error('TDDeltaStrengthMode')
  end select

  !!! read in the the polarization directions -- or eventually learn them from symmetry.
  if(td%delta_strength > M_ZERO) then


     ! Find out how many equivalent axis we have...
     ! WARNING: TODO: document this variable.
     call loct_parse_int(check_inp('TDPolarizationEquivAxis'), 0, td%pol_equiv_axis)

     call loct_parse_int(check_inp('TDPolarizationDirection'), 1, td%pol_dir)

     td%pol(:, :) = M_ZERO
     if(loct_parse_block(check_inp('TDPolarization'), blk)==0) then
        n = loct_parse_block_n(blk)
        do j = 1, n
           do i = 1, NDIM
              call loct_parse_block_float(blk, j-1, i-1, td%pol(i, j))
           end do
        enddo
        if(n==2) td%pol(1:3, 3) = (/ M_ZERO, M_ZERO, M_ONE /)
        if(n==1) td%pol(1:3, 2) = (/ M_ZERO, M_ONE, M_ZERO /)
        call loct_parse_block_end(blk)
     else
        ! Here the symmetry of the system should be analized, and the polarization
        ! basis, built accordingly.
        td%pol_dir = 1
        td%pol(1:3, 1) = (/ M_ONE, M_ZERO, M_ZERO /)
        td%pol(1:3, 2) = (/ M_ZERO, M_ONE, M_ZERO /)
        td%pol(1:3, 3) = (/ M_ZERO, M_ZERO, M_ONE /)
     endif
     ! Normalize:
     do i = 1, 3
        td%pol(1:3, i) = td%pol(1:3, i)/sqrt(sum(td%pol(1:3, i)**2))
     enddo

  else

     td%pol(1:3, 1) = (/ M_ONE, M_ZERO, M_ZERO /)
     td%pol(1:3, 2) = (/ M_ZERO, M_ONE, M_ZERO /)
     td%pol(1:3, 3) = (/ M_ZERO, M_ZERO, M_ONE /)

  endif

  ! now the photoelectron stuff
#if !defined(DISABLE_PES) && defined(HAVE_FFT)
  call loct_parse_int(check_inp('AbsorbingBoundaries'), 0, dummy)
  call PES_init(td%PESv, gr%m, gr%sb, st, dummy, outp%iter)
#endif

  ! should we move the ions during the simulation?
  call loct_parse_int(check_inp('MoveIons'), 0, td%move_ions)
  if( (td%move_ions .ne. STATIC_IONS) .and.   &
      (td%move_ions .ne. NORMAL_VERLET) .and. &
      (td%move_ions .ne. VELOCITY_VERLET) ) then
    write(message(1),'(a,i4,a)') "Input: '", td%move_ions, "' is not a valid MoveIons"
    message(2) = '  MoveIons = 0 <= do not move'
    message(3) = '  MoveIons = 3 <= verlet'
    message(4) = '  MoveIons = 4 <= velocity verlet'
    call write_fatal(4)
  endif
  if( td%move_ions .eq. NORMAL_VERLET) then
    write(message(1),'(a)') "Normal Verlet algorithm temporarily disabled."
    write(message(2),'(a)') "Using Velocity Verlet (MoveIons = 4) instead."
    call write_warning(2)
    td%move_ions = VELOCITY_VERLET
  endif

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

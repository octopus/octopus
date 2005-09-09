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

  integer :: i, dummy
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

  !!! read in the default direction for the polarization
  td%pol(:) = M_ZERO
  if(loct_parse_block(check_inp('TDPolarization'), blk)==0) then
    do i = 1, NDIM
      call loct_parse_block_float(blk, 0, i-1, td%pol(i))
    end do
    call loct_parse_block_end(blk)
  else  !default along the x-direction
    td%pol(1) = M_ONE
  endif
  ! Normalize:
  td%pol(1:3) = td%pol(1:3)/sqrt(sum(td%pol(1:3)**2))

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

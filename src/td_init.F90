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

subroutine td_init(td, sys, m, st, h)
  type(td_type), intent(out) :: td
  type(system_type), intent(IN) :: sys
  type(mesh_type), intent(inout) :: m
  type(states_type), intent(inout) :: st
  type(hamiltonian_type), intent(IN) :: h

  integer :: i, iunit, dummy

  call push_sub('td_init')

  td%iter = 0

  call oct_parse_double("TDTimeStep", 0.07_r8/units_inp%time%factor, td%dt)
  td%dt = td%dt * units_inp%time%factor
  if (td%dt <= 0._r8) then
    write(message(1),'(a,f14.6,a)') "Input: '", td%dt, "' is not a valid TDTimeStep"
    message(2) = '(0 < TDTimeStep)'
    call write_fatal(2)
  end if

  call oct_parse_int("TDMaximumIter", 1500, td%max_iter)
  if(td%max_iter < 1) then
    write(message(1), '(a,i6,a)') "Input: '", td%max_iter, "' is not a valid TDMaximumIter"
    message(2) = '(1 <= TDMaximumIter)'
    call write_fatal(2)
  end if
    
  call oct_parse_int("TDDipoleLmax", 1, td%lmax)
  if (td%lmax < 0 .or. td%lmax > 4) then
    write(message(1), '(a,i6,a)') "Input: '", td%lmax, "' is not a valid TDDipoleLmax"
    message(2) = '(0 <= TDDipoleLmax <= 4 )'
    call write_fatal(2)
  end if

  !!! read in the default direction for the polarization
  td%pol(:) = M_ZERO
  if(oct_parse_isdef('TDPolarization') .ne. 0) then
    do i = 1, conf%dim
      call oct_parse_block_double('TDPolarization', 0, i-1, td%pol(i))
    end do
  else  !default along the x-direction
    td%pol(1) = M_ONE
  endif

  ! now the photoelectron stuff
#if !defined(DISABLE_PES) && defined(HAVE_FFT)
  call oct_parse_int("AbsorbingBoundaries", 0, dummy)
  call PES_init(td%PESv, m, sys%st, dummy, sys%outp%iter)
#endif

  ! should we move the ions during the simulation?
  call oct_parse_int("MoveIons", 0, td%move_ions)
  if(td%move_ions.ne.0 .and. td%move_ions<3 .and. td%move_ions>4) then
    write(message(1),'(a,i4)') "Input: '", td%move_ions, &
         "' is not a valid MoveIons"
    message(2) = '  MoveIons = 0 <= do not move'
    message(3) = '  MoveIons = 3 <= verlet'
    message(4) = '  MoveIons = 4 <= velocity verlet'
    call write_fatal(4)
  endif
  
  ! Check what should be output
  call oct_parse_logical("TDOutputMultipoles", .true., td%out_multip)
  if(td%move_ions>0) then
    call oct_parse_logical("TDOutputCoordinates", .true., td%out_coords)
  else
    td%out_coords = .false.
  end if
  call oct_parse_logical("TDOutputAngularMomentum", .false., td%out_angular)
  call oct_parse_logical("TDOutputGSProjection", .false., td%out_gsp)
  call oct_parse_logical("TDOutputAcceleration", .false., td%out_acc)
  if(td%out_acc.and.td%move_ions>0) then
    message(1) = 'Error. If harmonic spectrum is to be calculated'
    message(2) = 'Atoms should not be allowed to move'
    call write_fatal(2)
  endif
  call oct_parse_logical("TDOutputLaser", h%ep%no_lasers>0, td%out_laser)
  call oct_parse_logical("TDOutputElEnergy", .false., td%out_energy)
  call oct_parse_logical("TDOutputOccAnalysis", .false., td%out_proj)

  call td_rti_init(m, st, td%tr)
  call td_init_states()

  call pop_sub()
contains
  
  subroutine td_init_states()
    integer :: i, ix

    ! MPI stuff
#if defined(HAVE_MPI) && defined(MPI_TD)
    if(st%nst < mpiv%numprocs) then
      message(1) = "Have more processors than necessary"
      write(message(2),'(i4,a,i4,a)') mpiv%numprocs, " processors and ", st%nst, " states."
      call write_fatal(2)
    end if

    i = st%nst / mpiv%numprocs
    ix = st%nst - i*mpiv%numprocs
    if(ix > 0 .and. mpiv%node < ix) then
      i = i + 1
      st%st_start = mpiv%node*i + 1
      st%st_end = st%st_start + i - 1
    else
      st%st_end = st%nst - (mpiv%numprocs - mpiv%node - 1)*i
      st%st_start = st%st_end - i + 1
    end if
    write(stdout, '(a,i4,a,i4,a,i4)') "Info: Node ", mpiv%node, " will propagate state ", &
         st%st_start, " - ", st%st_end

    ! syncronize to get the output properly
    call MPI_Barrier(MPI_COMM_WORLD, i)
  
#else
    st%st_start = 1
    st%st_end   = st%nst
#endif

    ! allocate memory
    allocate(st%zpsi(m%np, st%dim, st%st_start:st%st_end, st%nik))
  end subroutine td_init_states

end subroutine td_init

subroutine td_end(td)
  type(td_type), intent(inout) :: td

  call push_sub('td_end')

#if !defined(DISABLE_PES) && defined(HAVE_FFT)
  call PES_end(td%PESv)
#endif

  call td_rti_end(td%tr)  ! clean the evolution method

  call pop_sub()
end subroutine td_end

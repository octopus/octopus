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

subroutine td_init(td, sys, m, st)
  type(td_type), intent(out) :: td
  type(system_type), intent(IN) :: sys
  type(mesh_type), intent(inout) :: m
  type(states_type), intent(inout) :: st

  integer :: i, iunit, dummy

  sub_name = 'td_init'; call push_sub()

  td%iter = 0
  allocate(td%pol(conf%dim))

  call oct_parse_int("TDMaximumIter", 1500, td%max_iter)
  if(td%max_iter < 1) then
    write(message(1), '(a,i6,a)') "Input: '", td%max_iter, "' is not a valid TDMaximumIter"
    message(2) = '(1 <= TDMaximumIter)'
    call write_fatal(2)
  end if
  
  call oct_parse_double("TDTimeStep", 0.07_r8/units_inp%time%factor, td%dt)
  td%dt = td%dt * units_inp%time%factor
  if (td%dt <= 0._r8) then
    write(message(1),'(a,f14.6,a)') "Input: '", td%dt, "' is not a valid TDTimeStep"
    message(2) = '(0 < TDTimeStep)'
    call write_fatal(2)
  end if
  
  call oct_parse_int("TDEvolutionMethod", REVERSAL, td%evolution_method)
  select case(td%evolution_method)
    case(SIMPLE_EXP);           message(1) = 'Info: Evolution method:  Simple Exponential Method.'
    case(OLD_REVERSAL);         message(1) = 'Info: Evolution method:  Old-Style.'
    case(REVERSAL);             message(1) = 'Info: Evolution method:  Enforced Time-Reversal Symmetry'
    case(APP_REVERSAL);         message(1) = 'Info: Evolution method:  Approx.Enforced Time-Reversal Symmetry' 
    case(EXPONENTIAL_MIDPOINT); message(1) = 'Info: Evolution method:  Exponential Midpoint Rule.'
    case default
     write(message(1), '(a,i6,a)') "Input: '", td%evolution_method, "' is not a valid TDEvolutionMethod"
     message(2) = '(1 <= TDEvolutionMethod <= 4)'
     call write_fatal(2)
  end select
  call write_info(1)

  call oct_parse_int("TDExponentialMethod", FOURTH_ORDER, td%exp_method)
  select case(td%exp_method)
    case(FOURTH_ORDER);         message(1) = 'Info: Exponential method: 4th order expansion.'
    case(LANCZOS_EXPANSION);    message(1) = 'Info: Exponential method: Lanczos subspace approximation.'
    case(SPLIT_OPERATOR);       message(1) = 'Info: Exponential method: Split-Operator.'
      call mesh_alloc_ffts(m, 1)
    case(SUZUKI_TROTTER);       message(1) = 'Info: Exponential method: Suzuki-Trotter.'
      call mesh_alloc_ffts(m, 1)
    case(CHEBYSHEV);            message(1) = 'Info: Exponential method: Chebyshev.'
    case default
     write(message(1), '(a,i6,a)') "Input: '", td%exp_method, "' is not a valid TDEvolutionMethod"
     message(2) = '(1 <= TDExponentialMethod <= 4)'
     call write_fatal(2)
  end select
  call write_info(1)

  call oct_parse_double("TDLanczosTol", 5e-4_r8, td%lanczos_tol)
  if (td%lanczos_tol <= 0._r8) then
    write(message(1),'(a,f14.6,a)') "Input: '", td%lanczos_tol, "' is not a valid TDLanczosTol"
    message(2) = '(0 < TDLanczosTol)'
    call write_fatal(2)
  end if

  call oct_parse_int("TDExpOrder", 4, td%exp_order)
  if (td%exp_order < 2) then
    write(message(1), '(a,i6,a)') "Input: '", td%exp_order, "' is not a valid TDExpOrder"
    message(2) = '(2 <= TDExpOrder)'
    call write_fatal(2)
  end if

  call oct_parse_int("TDDipoleLmax", 1, td%lmax)
  if (td%lmax < 0 .or. td%lmax > 4) then
    write(message(1), '(a,i6,a)') "Input: '", td%lmax, "' is not a valid TDDipoleLmax"
    message(2) = '(0 <= TDDipoleLmax <= 4 )'
    call write_fatal(2)
  end if

  ! delta impulse used to calculate optical spectrum
  ! units are 1/length
  call oct_parse_double("TDDeltaStrength", 0._r8, td%delta_strength)
  td%delta_strength = td%delta_strength / units_inp%length%factor
  if(oct_parse_isdef('TDPolarization') .ne. 0) then
    do i = 1, conf%dim
      call oct_parse_block_double('TDPolarization', 0, i-1, td%pol(i))
    end do
  else  !default along the x-direction
    td%pol(:) = M_ZERO
    td%pol(1) = M_ONE
  endif

  ! now the photoelectron stuff
  call oct_parse_int("AbsorbingBoundaries", 0, dummy)
  call PES_init(td%PESv, m, sys%st, dummy, sys%outp%iter)

  ! occupational analysis stuff
  call oct_parse_logical("TDOccupationalAnalysis", .false., td%occ_analysis)

  ! harmonic spectrum or not
  call oct_parse_logical("TDWriteHarmonicSpectrum", .false., td%harmonic_spectrum)
  if(td%harmonic_spectrum) then
    message(1) = 'Warning: The harmonic spectrum, calculated from Ehrenfest theorem, '
    message(2) = '  is not yet well calculated if the ions move... Sorry!'
    call write_warning(2)
  endif

  ! Ground state component..
  call oct_parse_logical("TDWriteGScomponent", .false., td%gs_projection)

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
  
  call td_init_states()

  call pop_sub(); return
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
    allocate(td%v_old(m%np, st%nspin, 3))
    allocate(st%zpsi(0:m%np, st%dim, st%st_start:st%st_end, st%nik))
  end subroutine td_init_states

end subroutine td_init

subroutine td_end(td)
  type(td_type), intent(inout) :: td

  sub_name = 'td_end'; call push_sub()

  deallocate(td%pol)

#ifndef NO_PES
  call PES_end(td%PESv)
#endif

  if(associated(td%v_old)) then
    deallocate(td%v_old);  nullify(td%v_old)
  end if

  call pop_sub(); return
end subroutine td_end

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

#include "config_F90.h"

program excitations
  use global
  use liboct
  use states
  use system
  use mix
  use linear

  implicit none

  type(system_type) :: sys
  type(states_type), target :: st
  type(hartree_type) :: hart
  integer :: ierr, n_occ, n_unocc, flags(32)
  character(len=100) :: ch
  logical :: l

  mpiv%node = 0
  mpiv%numprocs = 1

  ! init some of the stuff
  ierr = oct_parse_init(C_string('inp'), C_string('out.oct'))
  if(ierr .ne. 0) then
    ierr = oct_parse_init(C_string("-"), C_string('out.oct'))
    if(ierr .ne. 0) then
      message(1) = "Error initializing liboct"
      call write_fatal(1)
    end if
  end if

  call oct_parse_int(C_string("verbose"), 30, conf%verbose)

  ! This only works for three dimensions
  conf%dim     = 3

  ! Initialize stuff
  call units_init()
  call system_init(sys)

  ! how many states do we have?
  call oct_parse_int(C_string("UnoccNumberStates"), 0, n_unocc)
  call states_dup(sys%st, st)
  n_occ = sys%st%nst
  st%nst = n_occ + n_unocc
  st%st_end   = st%nst
  allocate(st%R_FUNC(psi) (0:sys%m%np, st%dim, st%nst, st%nik))
  allocate(st%eigenval(st%nst, st%nik), st%occ(st%nst, st%nik), st%rho(sys%m%np, st%nspin))
  st%occ = 0._r8
  st%occ(1:n_occ,:) = sys%st%occ(1:n_occ,:)
  call states_choose_kpoints(st, sys%m)
  deallocate(sys%st)
  sys%st => st

  allocate(sys%st%R_FUNC(psi) (0:sys%m%np, sys%st%dim, sys%st%nst, sys%st%nik))
  if(R_FUNC(states_load_restart) ("restart.occ", sys%m, sys%st)) then
    call R_FUNC(calcdens)(sys%st, sys%m%np, st%rho)
  else
    message(1) = "Error opening 'restart.occ' file"
    call write_fatal(1)
  endif

  ! which states to take into account
  call oct_parse_str("ExciteStates", "1-1024", ch)
  call oct_wfs_list(C_string(ch), flags)

  ! initialize Hartree potential
  call hartree_init(hart, sys%m)

  ! calculate resonances
  message(1) = "Info: Eigenvalue differences"
  call write_info(1)
  call oct_parse_logical("LinEigenvalues", .true., l)
  if(l) call calc_petersilka(0, sys%st, sys%m, hart, n_occ, n_unocc, flags, 'linear', 'eps-diff')
  

  message(1) = "Info: Calculating resonance energies a la Petersilka"
  call write_info(1)
  call oct_parse_logical("LinPetersilka", .true., l)
  if(l) call calc_petersilka(1, sys%st, sys%m, hart, n_occ, n_unocc, flags, 'linear', 'petersilka')

  message(1) = "Info: Calculating resonance energies a la Casida"
  call write_info(1)
  call oct_parse_logical("LinCasida", .true., l)
  if(l)call calc_petersilka(2, sys%st, sys%m, hart, n_occ, n_unocc, flags, 'linear', 'full')

  call hartree_end(hart)

end program excitations

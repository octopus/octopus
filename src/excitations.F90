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
  call oct_parse_int(C_string('Dimensions'), 3, conf%dim)
  if(conf%dim<1 .or. conf%dim>3) then
    message(1) = 'Dimensions must be either 1, 2, or 3'
    call write_fatal(1)
  end if
  write(message(1), '(a,i1,a)') 'Octupus will run in ', conf%dim, ' dimension(s)'

  ! Initialize stuff
  call units_init()
  call system_init(sys)

  ! how many states do we have?
  n_occ = sys%st%nst
  call oct_parse_int(C_string("UnoccNumberStates"), 5, n_unocc)
  call states_init(st, sys%m, sys%val_charge)
  st%nst = n_unocc + n_occ
  st%st_end = st%nst
  allocate(st%R_FUNC(psi) (0:sys%m%np, st%dim, st%nst, st%nik))
  allocate(st%eigenval(st%nst, st%nik), st%occ(st%nst, st%nik))
  st%eigenval = M_ZERO
  st%occ      = M_ZERO
  st%occ(1:sys%st%nst,:) = sys%st%occ(1:sys%st%nst,:)
  call states_end(sys%st)

  if(R_FUNC(states_load_restart) ("tmp/restart.occ", sys%m, st)) then
    call R_FUNC(calcdens)(st, sys%m%np, st%rho)
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
  if(l) call calc_petersilka(0, st, sys%m, hart, n_occ, n_unocc, flags, 'linear', 'eps-diff')
  

  message(1) = "Info: Calculating resonance energies a la Petersilka"
  call write_info(1)
  call oct_parse_logical("LinPetersilka", .true., l)
  if(l) call calc_petersilka(1, st, sys%m, hart, n_occ, n_unocc, flags, 'linear', 'petersilka')

  message(1) = "Info: Calculating resonance energies a la Casida"
  call write_info(1)
  call oct_parse_logical("LinCasida", .true., l)
  if(l)call calc_petersilka(2, st, sys%m, hart, n_occ, n_unocc, flags, 'linear', 'casida')

  call hartree_end(hart)

end program excitations

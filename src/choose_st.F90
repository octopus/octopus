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

#include "global.h"

program choose_st
  use global
  use lib_oct_parser
  use lib_oct
  use units
  use states
  use system

  implicit none

  type(system_type) :: sys
  type(states_type) :: st1
  integer :: ierr, i, n_st, n_unocc, flags(32)
  character(len=100) :: ch

  ! Initialize stuff
  call global_init()
  call units_init()
  call system_init(sys)

  deallocate(sys%st%rho, sys%st%occ, sys%st%eigenval)
  call states_init(st1, sys%m, sys%val_charge)
  deallocate(st1%rho, st1%occ, st1%eigenval)

  ! how many states do we have?
  call loct_parse_int("UnoccNumberStates", 5, n_unocc)

  ! setup variables
  sys%st%nst = sys%st%nst + n_unocc
  sys%st%st_end = sys%st%nst
  allocate(sys%st%X(psi) (sys%m%np, sys%st%dim, sys%st%nst, sys%st%nik), &
       sys%st%eigenval(sys%st%nst, sys%st%nik))
  if(.not.X(states_load_restart) ("tmp/restart.occ", sys%m, sys%st)) then
    message(1) = "Error opening 'restart.occ' file"
    call write_fatal(1)
  endif

  ! which states to take into account
  call loct_parse_string("ChooseStates", "1-1024", ch)
  call loct_wfs_list(ch, flags)

  ! create directory if it does not exist
  call loct_mkdir('opt-control')

  n_st = 0
  do i = 1, sys%st%nst
    if(iand(flags((i-1)/32 + 1), 2**(modulo(i-1, 32))).ne.0) n_st = n_st + 1
  end do

  st1%nst = n_st
  st1%st_end = n_st
  allocate(st1%X(psi) (sys%m%np, st1%dim, st1%nst, st1%nik), &
       st1%eigenval(st1%nst, st1%nik))
  n_st = 1
  do i = 1, sys%st%nst
    if(iand(flags((i-1)/32 + 1), 2**(modulo(i-1, 32))).ne.0) then
      write(stdout, '(a,i4,a,i4)') "Including state ", i, " => ", n_st
      st1%X(psi) (:, 1:st1%dim, n_st, 1:st1%nik) = &
           sys%st%X(psi) (:, 1:st1%dim, i, 1:st1%nik)
      st1%eigenval(n_st, 1:st1%nik) = sys%st%eigenval(i, 1:st1%nik)      
      n_st = n_st + 1
    end if
  end do

  call loct_parse_string("ChooseStatesFilename", "wf.initial", ch)
  call X(states_write_restart) ("opt-control/"+trim(ch), sys%m, st1)

  stop
end program choose_st

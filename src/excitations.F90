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

program excitations
  use global
  use lib_oct_parser
  use lib_oct
#ifdef HAVE_FFT
  use fft
#endif
  use states
  use restart
  use system
  use linear

  implicit none

  type(system_type) :: sys
  type(states_type), target :: st
  integer :: n_occ, n_unocc, flags(32), ierr
  character(len=100) :: ch
  logical :: from_scratch, l

  ! Initialize stuff
  call global_init()

  ! initialize ffts
#ifdef HAVE_FFT
  call fft_all_init()
#endif

  call units_init()
  call system_init(sys)

  ! how many states do we have?
  n_occ = sys%st%nst
  call loct_parse_int("UnoccNumberStates", 5, n_unocc)
  call states_init(st, sys%m, sys%geo, sys%val_charge)
  st%nst = n_unocc + n_occ
  st%st_end = st%nst
  allocate(st%X(psi) (sys%m%np, st%dim, st%nst, st%nik))
  allocate(st%eigenval(st%nst, st%nik), st%occ(st%nst, st%nik))
  st%eigenval = M_ZERO
  st%occ      = M_ZERO
  st%occ(1:sys%st%nst,:) = sys%st%occ(1:sys%st%nst,:)
  call states_end(sys%st)

  call restart_load("tmp/restart_occ", st, sys%m, ierr)
  if(ierr.ne.0) then
    message(1) = "Error opening 'restart.occ' file"
    call write_fatal(1)
  endif

  ! Get the density...
  call X(calcdens)(st, sys%m%np, st%rho)

  ! which states to take into account
  call loct_parse_string("ExciteStates", "1-1024", ch)
  call loct_wfs_list(ch, flags)

  ! initialize Poisson solver
  call poisson_init(sys%m)

  ! should we start from scratch, or use restart file
  call loct_parse_logical("FromScratch", .false., from_scratch)

  ! calculate resonances
  message(1) = "Info: Eigenvalue differences"
  call write_info(1)
  call loct_parse_logical("LinEigenvalues", .true., l)
  if(l) call LR_el_excitations(0, st, sys%m, n_occ, n_unocc, flags, 'linear', 'eps-diff', from_scratch)
  

  message(1) = "Info: Calculating resonance energies a la Petersilka"
  call write_info(1)
  call loct_parse_logical("LinPetersilka", .true., l)
  if(l) call LR_el_excitations(1, st, sys%m, n_occ, n_unocc, flags, 'linear', 'petersilka', from_scratch)

  message(1) = "Info: Calculating resonance energies a la Casida"
  call write_info(1)
  call loct_parse_logical("LinCasida", .true., l)
  if(l)call LR_el_excitations(2, st, sys%m, n_occ, n_unocc, flags, 'linear', 'casida', from_scratch)

  call poisson_end()
  call global_end()
end program excitations

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
  use system
  use linear

  implicit none

  type(system_type) :: sys
  type(states_type), target :: st
  integer :: n_occ, n_unocc, flags(32)
  character(len=100) :: ch
  logical :: l

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
  call states_init(st, sys%m, sys%val_charge)
  st%nst = n_unocc + n_occ
  st%st_end = st%nst
  allocate(st%X(psi) (sys%m%np, st%dim, st%nst, st%nik))
  allocate(st%eigenval(st%nst, st%nik), st%occ(st%nst, st%nik))
  st%eigenval = M_ZERO
  st%occ      = M_ZERO
  st%occ(1:sys%st%nst,:) = sys%st%occ(1:sys%st%nst,:)
  call states_end(sys%st)

  if(X(states_load_restart) ("tmp/restart.occ", sys%m, st)) then
    call X(calcdens)(st, sys%m%np, st%rho)
  else
    message(1) = "Error opening 'restart.occ' file"
    call write_fatal(1)
  endif

  ! which states to take into account
  call loct_parse_string("ExciteStates", "1-1024", ch)
  call loct_wfs_list(ch, flags)

  ! initialize Poisson solver
  call poisson_init(sys%m, sys%geo)

  ! calculate resonances
  message(1) = "Info: Eigenvalue differences"
  call write_info(1)
  call loct_parse_logical("LinEigenvalues", .true., l)
  if(l) call calc_petersilka(0, st, sys%m, n_occ, n_unocc, flags, 'linear', 'eps-diff')
  

  message(1) = "Info: Calculating resonance energies a la Petersilka"
  call write_info(1)
  call loct_parse_logical("LinPetersilka", .true., l)
  if(l) call calc_petersilka(1, st, sys%m, n_occ, n_unocc, flags, 'linear', 'petersilka')

  message(1) = "Info: Calculating resonance energies a la Casida"
  call write_info(1)
  call loct_parse_logical("LinCasida", .true., l)
  if(l)call calc_petersilka(2, st, sys%m, n_occ, n_unocc, flags, 'linear', 'casida')

  call poisson_end()

end program excitations

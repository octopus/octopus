!!! Copyright (C) 2008-2009 David Strubbe
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
!! $Id: kdotp_calc.F90 2548 2006-11-06 21:42:27Z xavier $

#include "global.h"

module kdotp_calc_m
  use global_m
  use hamiltonian_m
  use linear_response_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use pert_m
  use profiling_m
  use states_m
  use states_calc_m
  use sternheimer_m
  use system_m

  implicit none

  private
  public ::                        &
    dcalc_eff_mass_inv,            &  
    zcalc_eff_mass_inv,            &  
    zcalc_dipole_periodic,         &
    kdotp_wfs_tag,                 &
    kdotp_rho_tag

contains

! ---------------------------------------------------------
  character(len=100) function kdotp_rho_tag(dir) result(str)
    integer, intent(in) :: dir

    call push_sub('kdotp_calc.kdotp_rho_tag')

    write(str, '(a,i1)') 'rho_', dir

    call pop_sub()

  end function kdotp_rho_tag
  

! ---------------------------------------------------------
  character(len=100) function kdotp_wfs_tag(dir) result(str)
    integer, intent(in) :: dir 

    call push_sub('kdotp_calc.kdotp_wfs_tag')

    write(str, '(a,i1)') "wfs_", dir

    call pop_sub()

  end function kdotp_wfs_tag

! ---------------------------------------------------------
! This routine cannot be used with d/dk wavefunctions calculated
! by perturbation theory.
subroutine zcalc_dipole_periodic(sys, lr, dipole)
  type(system_t),         intent(inout) :: sys
  type(lr_t),             intent(in)    :: lr(:,:)
  FLOAT,                  intent(out)   :: dipole(:)

  integer idir, ist, ik, idim
  type(mesh_t), pointer :: m
  CMPLX :: term, moment
  m => sys%gr%mesh

  ! mu_i = sum(m occ, k) <u_mk(0)|(-id/dk_i|u_mk(0)>)
  !      = Im sum(m occ, k) <u_mk(0)|(d/dk_i|u_mk(0)>)

  call push_sub('kdotp_calc.zcalc_dipole_periodic')

  do idir = 1, sys%gr%sb%dim
    moment = M_ZERO

    do ik = 1, sys%st%d%nik
      term = M_ZERO

      do ist = 1, sys%st%nst
        do idim = 1, sys%st%d%dim
          term = term + zmf_dotp(m, sys%st%zpsi(1:m%np, idim, ist, ik), lr(idir, 1)%zdl_psi(1:m%np, idim, ist, ik))
        enddo
      enddo

      moment = moment + term * sys%st%d%kweights(ik) * sys%st%smear%el_per_state
    enddo

    dipole(idir) = -moment
  enddo

  call pop_sub()

end subroutine zcalc_dipole_periodic

#include "undef.F90"
#include "real.F90"
#include "kdotp_calc_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "kdotp_calc_inc.F90"

end module kdotp_calc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

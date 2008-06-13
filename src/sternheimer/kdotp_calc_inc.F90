!! Copyright (C) 2004 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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
!! $Id: kdotp_calc_inc.F90 2548 2006-11-06 21:42:27Z xavier $

subroutine X(lr_calc_eff_mass_inv)(sys, h, lr, perturbation, eff_mass_inv)
  type(system_t),         intent(inout) :: sys
  type(hamiltonian_t),    intent(in)    :: h
!  type(kdotp_t),          intent(inout) :: kdotp_vars
  type(lr_t),             intent(in)    :: lr(:,:)
  type(pert_t),           intent(inout) :: perturbation
  FLOAT,                  intent(out)   :: eff_mass_inv(:, :, :, :)
!  integer, optional,      intent(in)    :: ndir
!
!  integer :: dir1, dir2, ndir_
!

! m^-1[ij] = delta[ij] + 2*Re<psi0|H'i|psi'j>
! for each state, spin, and k-point
! This routine is not set up for spinors.

  integer ik, ist, idir1, idir2
  R_TYPE, allocatable :: pertpsi(:)

  ALLOCATE(pertpsi(1:sys%gr%m%np), sys%gr%m%np)

  eff_mass_inv(:, :, :, :) = 0

  do ik = 1, sys%st%d%nik
    do ist = 1, sys%st%nst
      do idir1 = 1, sys%gr%sb%dim

        eff_mass_inv(ik, ist, idir1, idir1) = 1
        call pert_setup_dir(perturbation, idir1)
        call X(pert_apply) (perturbation, sys%gr, sys%geo, h, ik, &
          sys%st%X(psi)(:, 1, ist, ik), pertpsi)

        do idir2 = 1, sys%gr%sb%dim
          eff_mass_inv(ik, ist, idir1, idir2) = eff_mass_inv(ik, ist, idir1, idir2) + &
             M_TWO*R_REAL(X(mf_integrate)(sys%gr%m, lr(idir2, 1)%X(dl_psi)(1:sys%gr%m%np, 1, ist, ik)*pertpsi(1:sys%gr%m%np)))
        enddo
      enddo
    enddo
  enddo

end subroutine X(lr_calc_eff_mass_inv)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

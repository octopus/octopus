!!! Copyright (C) 2004 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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
!  use kdotp_lr_m
  use linear_response_m
  use mesh_function_m
  use messages_m
  use pert_m
  use profiling_m
  use states_m
  use sternheimer_m
  use system_m
!  use units_m

  implicit none

  private
  public ::                       &
    zlr_calc_eff_mass_inv,        &  
    kdotp_wfs_tag,                &
    kdotp_rho_tag

contains

  character(len=100) function kdotp_rho_tag(dir) result(str)
    integer, intent(in) :: dir

    call push_sub('kdotp_calc.kdotp_rho_tag')

    write(str, '(a,i1)') 'rho_', dir

    call pop_sub()

  end function kdotp_rho_tag
  
  character(len=100) function kdotp_wfs_tag(dir) result(str)
    integer, intent(in) :: dir 

    call push_sub('kdotp_calc.kdotp_wfs_tag')

    write(str, '(a,i1)') "wfs_", dir

    call pop_sub()

  end function kdotp_wfs_tag

#include "complex.F90"

subroutine zlr_calc_eff_mass_inv(sys, h, lr, perturbation, eff_mass_inv, occ_solution_method, degen_thres)
  type(system_t),         intent(inout) :: sys
  type(hamiltonian_t),    intent(in)    :: h
  type(lr_t),             intent(in)    :: lr(:,:)
  type(pert_t),           intent(inout) :: perturbation
  FLOAT,                  intent(out)   :: eff_mass_inv(:, :, :, :)
  integer,                intent(in)    :: occ_solution_method
  FLOAT,                  intent(in)    :: degen_thres

! m^-1[ij] = delta[ij] + 2*Re<psi0|H'i|psi'j>
! for each state, spin, and k-point
! This routine is not set up for spinors.
! The off-diagonal elements are not correct in a degenerate subspace

  integer ik, ist, ist2, idir1, idir2
  R_TYPE term
  R_TYPE, allocatable :: pertpsi(:,:)       ! H'i|psi0>
  R_TYPE, allocatable :: proj_dl_psi(:)   ! (1-Pn')|psi'j>

  ALLOCATE(pertpsi(sys%gr%sb%dim, sys%gr%m%np), sys%gr%sb%dim * sys%gr%m%np)
  ALLOCATE(proj_dl_psi(sys%gr%m%np), sys%gr%m%np)

  eff_mass_inv(:, :, :, :) = 0

  do ik = 1, sys%st%d%nik

    do ist = 1, sys%st%nst

      do idir1 = 1, sys%gr%sb%dim
        call pert_setup_dir(perturbation, idir1)
        call zpert_apply (perturbation, sys%gr, sys%geo, h, ik, &
          sys%st%zpsi(:, 1, ist, ik), pertpsi(idir1,:))
      enddo

      do idir1 = 1, sys%gr%sb%dim
        eff_mass_inv(ik, ist, idir1, idir1) = 1

        do idir2 = 1, sys%gr%sb%dim

          if (idir2 < idir1) then
            eff_mass_inv(ik, ist, idir1, idir2) = eff_mass_inv(ik, ist, idir2, idir1)
            ! utilizing symmetry of inverse effective mass tensor
            cycle
          end if

          proj_dl_psi(1:sys%gr%m%np) = lr(idir2, 1)%zdl_psi(1:sys%gr%m%np, 1, ist, ik)
          
          if (occ_solution_method == 0) then
          ! project out components of other states in degenerate subspace
             do ist2 = 1, sys%st%nst
                if (ist2 == ist) cycle ! projection on unperturbed wfn already removed in Sternheimer eqn
                if (abs(sys%st%eigenval(ist2, ik) - sys%st%eigenval(ist, ik)) < degen_thres) then
                   proj_dl_psi(1:sys%gr%m%np) = proj_dl_psi(1:sys%gr%m%np) &
                     - sys%st%zpsi(1:sys%gr%m%np, 1, ist2, ik) * &
                     zmf_integrate(sys%gr%m, sys%st%zpsi(1:sys%gr%m%np, 1, ist2, ik) * proj_dl_psi(1:sys%gr%m%np))
                endif
             enddo
          endif

          ! contribution from Sternheimer equation
          term = zmf_integrate(sys%gr%m, proj_dl_psi(1:sys%gr%m%np) * conjg(pertpsi(idir1, 1:sys%gr%m%np)))
          eff_mass_inv(ik, ist, idir1, idir2) = eff_mass_inv(ik, ist, idir1, idir2) + M_TWO * term
!          write(*,*) "unocc: ", eff_mass_inv(ik, ist, idir1, idir2)

          if (occ_solution_method == 1) then
          ! contribution from linear-response projection onto occupied states, by sum over states and perturbation theory
          ! this could be sped up yet more by storing each (ist,ist2) term, which will be used again as the complex
          !   conjugate as the (ist2,ist) term.  same will apply to hyperpolarizability
             do ist2 = 1, sys%st%nst
                if (ist2 == ist .or. abs(sys%st%eigenval(ist2, ik) - sys%st%eigenval(ist, ik)) < degen_thres) cycle
                term = zmf_integrate(sys%gr%m, conjg(pertpsi(idir1, 1:sys%gr%m%np)) * sys%st%zpsi(1:sys%gr%m%np, 1, ist2, ik)) * &
                     zmf_integrate(sys%gr%m, conjg(sys%st%zpsi(1:sys%gr%m%np, 1, ist2, ik)) * pertpsi(idir2, 1:sys%gr%m%np)) / &
                     (sys%st%eigenval(ist, ik) - sys%st%eigenval(ist2, ik))
                eff_mass_inv(ik, ist, idir1, idir2) = eff_mass_inv(ik, ist, idir1, idir2) + M_TWO * term
!                write(*,'(a,i2,a,f20.6)') "occ(", ist2, "): ", eff_mass_inv(ik, ist, idir1, idir2)
             enddo
          endif
!          write(*,'(a,i1,a,i1,a,f20.6)') "eff_mass_inv(", idir1, ", ", idir2, ") = ", eff_mass_inv(ik, ist, idir1, idir2)

        enddo !idir2
      enddo !idir1
    enddo !ist
  enddo !ik

end subroutine zlr_calc_eff_mass_inv
  
end module kdotp_calc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

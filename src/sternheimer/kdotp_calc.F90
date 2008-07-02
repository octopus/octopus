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
  use units_m

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

! m^-1[ij] = delta[ij] + 4*Re<psi0|H'i|psi'j>
! for each state, spin, and k-point
! This routine is not set up for spinors.

  integer ik, ist, ist2, idir1, idir2
  R_TYPE term
  R_TYPE, allocatable :: pertpsi(:,:)

  ALLOCATE(pertpsi(sys%gr%sb%dim, sys%gr%m%np), sys%gr%sb%dim * sys%gr%m%np)

  eff_mass_inv(:, :, :, :) = 0

  do ik = 1, sys%st%d%nik

    ist = 1
    do while (ist < sys%st%nst)

      ! test for degeneracies
      call messages_print_stress(stdout, 'Degenerate subspace')
      write(*,'(a, i3, a, f12.8, a, a)') 'State #', ist, ', Energy = ', &
           sys%st%eigenval(ist, ik)/units_out%energy%factor, ' ', units_out%energy%abbrev              
      ist2 = ist + 1
      do while (ist2 <= sys%st%nst .and. sys%st%eigenval(ist2, ik) - sys%st%eigenval(ist, ik) < degen_thres)
         ! eigenvalues are supposed to be in ascending order; if they are not, it is a sign
         ! of being in a degenerate subspace; hence no absolute value here
!         write(*,*) ist2, ist, sys%st%eigenval(ist2, ik) - sys%st%eigenval(ist, ik)
         write(*,'(a, i3, a, f12.8, a, a)') 'State #', ist2, ', Energy = ', &
            sys%st%eigenval(ist2, ik)/units_out%energy%factor, ' ', units_out%energy%abbrev              
         ist2 = ist2 + 1
      enddo

      ist = ist2
    enddo

    do ist = 1, sys%st%nst

      do idir1 = 1, sys%gr%sb%dim
        call pert_setup_dir(perturbation, idir1)
        call zpert_apply (perturbation, sys%gr, sys%geo, h, ik, &
          sys%st%zpsi(:, 1, ist, ik), pertpsi(idir1,:))
      enddo

      do idir1 = 1, sys%gr%sb%dim
        eff_mass_inv(ik, ist, idir1, idir1) = 1

        do idir2 = 1, sys%gr%sb%dim
          ! contribution from Sternheimer equation
          if (idir2 < idir1) then
            eff_mass_inv(ik, ist, idir1, idir2) = eff_mass_inv(ik, ist, idir2, idir1)
            ! utilizing symmetry of inverse effective mass tensor
            cycle
          end if

          term = zmf_integrate(sys%gr%m, lr(idir2, 1)%zdl_psi(1:sys%gr%m%np, 1, ist, ik) &
            * conjg(pertpsi(idir1, 1:sys%gr%m%np)))
          eff_mass_inv(ik, ist, idir1, idir2) = eff_mass_inv(ik, ist, idir1, idir2) + M_FOUR * real(term)
!          write(*,*) "unocc: ", eff_mass_inv(ik, ist, idir1, idir2)

          if (occ_solution_method == 1) then
          ! contribution from linear-response projection onto occupied states, by sum over states and perturbation theory
          ! this could be sped up yet more by storing each (ist,ist2) term, which will be used again as the complex
          !   conjugate as the (ist2,ist) term.  same will apply to hyperpolarizability
             do ist2 = 1, sys%st%nst
                if (ist2 == ist) cycle
                term = zmf_integrate(sys%gr%m, conjg(pertpsi(idir1, 1:sys%gr%m%np)) * sys%st%zpsi(1:sys%gr%m%np, 1, ist2, ik)) * &
                     zmf_integrate(sys%gr%m, conjg(sys%st%zpsi(1:sys%gr%m%np, 1, ist2, ik)) * pertpsi(idir2, 1:sys%gr%m%np)) / &
                     (sys%st%eigenval(ist, ik) - sys%st%eigenval(ist2, ik))
                eff_mass_inv(ik, ist, idir1, idir2) = eff_mass_inv(ik, ist, idir1, idir2) + M_FOUR * real(term)
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

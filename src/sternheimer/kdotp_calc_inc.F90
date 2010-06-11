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

! ---------------------------------------------------------
subroutine X(calc_eff_mass_inv)(sys, hm, lr, perturbation, eff_mass_inv, &
  occ_solution_method, degen_thres)
  type(system_t),         intent(inout) :: sys
  type(hamiltonian_t),    intent(in)    :: hm
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
  R_TYPE, allocatable   :: pertpsi(:, :, :)       ! H`i|psi0>
  R_TYPE, allocatable   :: proj_dl_psi(:, :)   ! (1-Pn`)|psi`j>
  type(mesh_t), pointer :: mesh
  logical, allocatable  :: orth_mask(:)

  call push_sub('kdotp_calc_inc.Xcalc_eff_mass_inv')

  mesh => sys%gr%mesh

  SAFE_ALLOCATE(pertpsi(1:mesh%np, 1:hm%d%dim, 1:sys%gr%sb%dim))
  SAFE_ALLOCATE(proj_dl_psi(1:mesh%np, 1)) ! second index should be sys%st%d%dim, i.e. spinors
  SAFE_ALLOCATE(orth_mask(1:sys%st%nst))

  eff_mass_inv(:, :, :, :) = 0

  do ik = 1, sys%st%d%nik

    do ist = 1, sys%st%nst

      ! start by computing all the wavefunctions acted on by perturbation
      do idir1 = 1, sys%gr%sb%dim
        call pert_setup_dir(perturbation, idir1)
        call X(pert_apply)(perturbation, sys%gr, sys%geo, hm, ik, &
          sys%st%X(psi)(:, :, ist, ik), pertpsi(:, :, idir1))
      enddo

      do idir1 = 1, sys%gr%sb%dim
        eff_mass_inv(ik, ist, idir1, idir1) = 1

        do idir2 = 1, sys%gr%sb%dim

          if (idir2 < idir1) then
            eff_mass_inv(ik, ist, idir1, idir2) = eff_mass_inv(ik, ist, idir2, idir1)
            ! utilizing symmetry of inverse effective mass tensor
            cycle
          end if

          proj_dl_psi(1:mesh%np, 1) = lr(idir2, 1)%X(dl_psi)(1:mesh%np, 1, ist, ik)
          
          if (occ_solution_method == 0) then
          ! project out components of other states in degenerate subspace
            do ist2 = 1, sys%st%nst
!              alternate direct method
!              if (abs(sys%st%eigenval(ist2, ik) - sys%st%eigenval(ist, ik)) < degen_thres) then
!                   proj_dl_psi(1:mesh%np) = proj_dl_psi(1:mesh%np) - sys%st%X(psi)(1:mesh%np, 1, ist2, ik) * &
!                     X(mf_dotp)(m, sys%st%X(psi)(1:mesh%np, 1, ist2, ik), proj_dl_psi(1:mesh%np))
              orth_mask(ist2) = .not. (abs(sys%st%eigenval(ist2, ik) - sys%st%eigenval(ist, ik)) < degen_thres)
              ! mask == .false. means do projection; .true. means do not
            enddo

!            orth_mask(ist) = .true. ! projection on unperturbed wfn already removed in Sternheimer eqn

            call X(states_orthogonalization)(mesh, sys%st%nst, sys%st%d%dim, sys%st%X(psi)(1:mesh%np, 1:1, 1:sys%st%nst, ik), &
              proj_dl_psi(1:mesh%np, 1:1), mask = orth_mask(1:sys%st%nst))
          endif

          ! contribution from Sternheimer equation
          term = X(mf_dotp)(mesh, sys%st%d%dim, proj_dl_psi, pertpsi(:, :, idir1))
          eff_mass_inv(ik, ist, idir1, idir2) = eff_mass_inv(ik, ist, idir1, idir2) + M_TWO * term

          if (occ_solution_method == 1) then
          ! contribution from linear-response projection onto occupied states, by sum over states and perturbation theory
          ! this could be sped up yet more by storing each (ist,ist2) term, which will be used again as the complex
          !   conjugate as the (ist2,ist) term.  same will apply to hyperpolarizability
             do ist2 = 1, sys%st%nst
                if (ist2 == ist .or. abs(sys%st%eigenval(ist2, ik) - sys%st%eigenval(ist, ik)) < degen_thres) cycle

                term = X(mf_dotp)(mesh, sys%st%d%dim, pertpsi(:, :, idir1), sys%st%X(psi)(:, :, ist2, ik)) * &
                     X(mf_dotp)(mesh, sys%st%d%dim, sys%st%X(psi)(:, :, ist2, ik), pertpsi(:, :, idir2)) / &
                     (sys%st%eigenval(ist, ik) - sys%st%eigenval(ist2, ik))
                eff_mass_inv(ik, ist, idir1, idir2) = eff_mass_inv(ik, ist, idir1, idir2) + M_TWO * term
             enddo
          endif

        enddo !idir2
      enddo !idir1
    enddo !ist
  enddo !ik

  call pop_sub('kdotp_calc_inc.Xcalc_eff_mass_inv')

end subroutine X(calc_eff_mass_inv)
  
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

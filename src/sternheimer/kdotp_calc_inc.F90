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
  type(hamiltonian_t),    intent(inout) :: hm
  type(lr_t),             intent(in)    :: lr(:,:)
  type(pert_t),           intent(inout) :: perturbation
  FLOAT,                  intent(out)   :: eff_mass_inv(:,:,:,:)
  integer,                intent(in)    :: occ_solution_method
  FLOAT,                  intent(in)    :: degen_thres

! m^-1[ij] = <psi0|H2ij|psi0> + 2*Re<psi0|H'i|psi'j>
! for each state, spin, and k-point
! The off-diagonal elements are not correct in a degenerate subspace

  integer :: ik, ist, ist2, idir1, idir2, pdim
  R_TYPE :: term
  R_TYPE, allocatable   :: pertpsi(:,:,:)     ! H`i|psi0>
  R_TYPE, allocatable   :: pertpsi2(:,:)      ! H2i|psi0>
  R_TYPE, allocatable   :: proj_dl_psi(:,:)   ! (1-Pn`)|psi`j>
  type(mesh_t), pointer :: mesh
  logical, allocatable  :: orth_mask(:)
#ifdef HAVE_MPI
  FLOAT, allocatable :: eff_mass_inv_temp(:,:,:,:)
#endif

  PUSH_SUB(X(calc_eff_mass_inv))

  mesh => sys%gr%mesh
  pdim = sys%gr%sb%periodic_dim

  SAFE_ALLOCATE(pertpsi(1:mesh%np, 1:hm%d%dim, 1:pdim))
  SAFE_ALLOCATE(pertpsi2(1:mesh%np, 1:hm%d%dim))
  SAFE_ALLOCATE(proj_dl_psi(1:mesh%np, 1:hm%d%dim))
  SAFE_ALLOCATE(orth_mask(1:sys%st%nst))
#ifdef HAVE_MPI
  SAFE_ALLOCATE(eff_mass_inv_temp(1:sys%st%d%nik, 1:sys%st%nst, 1:pdim, 1:pdim))
#endif

  eff_mass_inv(:,:,:,:) = M_ZERO

  do ik = sys%st%d%kpt%start, sys%st%d%kpt%end

    do ist = sys%st%st_start, sys%st%st_end

      ! start by computing all the wavefunctions acted on by perturbation
      do idir1 = 1, pdim
        call pert_setup_dir(perturbation, idir1)
        call X(pert_apply)(perturbation, sys%gr, sys%geo, hm, ik, &
          sys%st%X(psi)(:, :, ist, ik), pertpsi(:, :, idir1))
      enddo

      do idir2 = 1, pdim
        do idir1 = 1, pdim

          if (idir2 < idir1) then
            eff_mass_inv(idir1, idir2, ist, ik) = eff_mass_inv(idir1, idir2, ist, ik)
            ! utilizing symmetry of inverse effective mass tensor
            cycle
          end if

          proj_dl_psi(1:mesh%np, 1:hm%d%dim) = lr(idir2, 1)%X(dl_psi)(1:mesh%np, 1:hm%d%dim, ist, ik)
          
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

            call X(states_orthogonalize_single)(sys%st, mesh, sys%st%nst, ik, proj_dl_psi, mask = orth_mask)
          endif

          ! contribution from Sternheimer equation
          term = X(mf_dotp)(mesh, sys%st%d%dim, proj_dl_psi, pertpsi(:, :, idir1))
          eff_mass_inv(idir1, idir2, ist, ik) = eff_mass_inv(idir1, idir2, ist, ik) + M_TWO * term

          if (occ_solution_method == 1) then
          ! contribution from linear-response projection onto occupied states, by sum over states and perturbation theory
             do ist2 = 1, sys%st%nst
                if (ist2 == ist .or. abs(sys%st%eigenval(ist2, ik) - sys%st%eigenval(ist, ik)) < degen_thres) cycle

                term = X(mf_dotp)(mesh, sys%st%d%dim, pertpsi(:, :, idir1), sys%st%X(psi)(:, :, ist2, ik)) * &
                     X(mf_dotp)(mesh, sys%st%d%dim, sys%st%X(psi)(:, :, ist2, ik), pertpsi(:, :, idir2)) / &
                     (sys%st%eigenval(ist, ik) - sys%st%eigenval(ist2, ik))
                eff_mass_inv(idir1, idir2, ist, ik) = eff_mass_inv(idir1, idir2, ist, ik) + M_TWO * term
             enddo
          endif

          call pert_setup_dir(perturbation, idir1, idir2)
          call X(pert_apply_order_2)(perturbation, sys%gr, sys%geo, hm, ik, &
            sys%st%X(psi)(1:mesh%np, 1:hm%d%dim, ist, ik), pertpsi2(1:mesh%np, 1:hm%d%dim))
          eff_mass_inv(idir1, idir2, ist, ik) = eff_mass_inv(idir1, idir2, ist, ik) + &
            X(mf_dotp)(mesh, hm%d%dim, sys%st%X(psi)(1:mesh%np, 1:hm%d%dim, ist, ik), &
            pertpsi2(1:mesh%np, 1:hm%d%dim))

        enddo !idir2
      enddo !idir1
    enddo !ist
  enddo !ik

#ifdef HAVE_MPI
  if(sys%st%parallel_in_states) then
    call MPI_Allreduce(eff_mass_inv, eff_mass_inv_temp, pdim**2 * sys%st%nst * sys%st%d%nik, &
      MPI_FLOAT, MPI_SUM, sys%st%mpi_grp%comm, mpi_err)
    eff_mass_inv(:,:,:,:) = eff_mass_inv_temp(:,:,:,:)
  endif
  if(sys%st%d%kpt%parallel) then
    call MPI_Allreduce(eff_mass_inv, eff_mass_inv_temp, pdim**2 * sys%st%nst * sys%st%d%nik, &
      MPI_FLOAT, MPI_SUM, sys%st%d%kpt%mpi_grp%comm, mpi_err)
    eff_mass_inv(:,:,:,:) = eff_mass_inv_temp(:,:,:,:)
  endif
  SAFE_DEALLOCATE_A(eff_mass_inv_temp)
#endif

  SAFE_DEALLOCATE_A(pertpsi)
  SAFE_DEALLOCATE_A(pertpsi2)
  SAFE_DEALLOCATE_A(proj_dl_psi)
  SAFE_DEALLOCATE_A(orth_mask)

  POP_SUB(X(calc_eff_mass_inv))

end subroutine X(calc_eff_mass_inv)
  

! ---------------------------------------------------------
subroutine X(kdotp_add_diagonal)(sys, hm, em_pert, kdotp_lr)
  type(system_t),      intent(inout) :: sys
  type(hamiltonian_t), intent(inout) :: hm
  type(pert_t),        intent(inout) :: em_pert
  type(lr_t),          intent(inout) :: kdotp_lr(:)

  integer :: ik, ist, idir
  R_TYPE, allocatable :: ppsi(:,:)
  FLOAT :: expectation

  PUSH_SUB(X(kdotp_add_diagonal))

  SAFE_ALLOCATE(ppsi(sys%gr%mesh%np, sys%st%d%dim))

  do idir = 1, sys%gr%sb%periodic_dim
    call pert_setup_dir(em_pert, idir)
    do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
      do ist = 1, sys%st%nst
        call X(pert_apply)(em_pert, sys%gr, sys%geo, hm, ik, sys%st%X(psi)(:, :, ist, ik), ppsi)
        expectation = X(mf_dotp)(sys%gr%mesh, sys%st%d%dim, sys%st%X(psi)(:, :, ist, ik), ppsi)
        kdotp_lr(idir)%X(dl_psi)(1:sys%gr%mesh%np, 1:sys%st%d%dim, ist, ik) = &
          kdotp_lr(idir)%X(dl_psi)(1:sys%gr%mesh%np, 1:sys%st%d%dim, ist, ik) + &
          expectation * sys%st%X(psi)(1:sys%gr%mesh%np, 1:sys%st%d%dim, ist, ik)
      enddo
    end do
  enddo

  SAFE_DEALLOCATE_A(ppsi) 
  POP_SUB(X(kdotp_add_diagonal))
end subroutine X(kdotp_add_diagonal)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

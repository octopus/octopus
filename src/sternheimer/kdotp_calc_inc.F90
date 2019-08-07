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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

! ---------------------------------------------------------
!> m^-1[ij] = <psi0|H2ij|psi0> + 2*Re<psi0|H'i|psi'j>
!! for each state, spin, and k-point
!! The off-diagonal elements are not correct in a degenerate subspace
subroutine X(calc_eff_mass_inv)(sys, lr, perturbation, eff_mass_inv, degen_thres)
  type(system_t), target, intent(inout) :: sys
  type(lr_t),             intent(in)    :: lr(:,:) !< (1, pdim)
  type(pert_t),           intent(inout) :: perturbation
  FLOAT,                  intent(out)   :: eff_mass_inv(:,:,:,:) !< (pdim, pdim, nik, nst)
  FLOAT,                  intent(in)    :: degen_thres

  integer :: ik, ist, ist2, idir1, idir2, pdim
  R_TYPE :: term
  R_TYPE, allocatable   :: psi(:, :)
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

  SAFE_ALLOCATE(psi(1:mesh%np_part, 1:sys%hm%d%dim))
  SAFE_ALLOCATE(pertpsi(1:mesh%np, 1:sys%hm%d%dim, 1:pdim))
  SAFE_ALLOCATE(pertpsi2(1:mesh%np, 1:sys%hm%d%dim))
  SAFE_ALLOCATE(proj_dl_psi(1:mesh%np, 1:sys%hm%d%dim))
  SAFE_ALLOCATE(orth_mask(1:sys%st%nst))
#ifdef HAVE_MPI
  SAFE_ALLOCATE(eff_mass_inv_temp(1:pdim, 1:pdim, 1:sys%st%nst, 1:sys%st%d%nik))
#endif

  eff_mass_inv(:,:,:,:) = M_ZERO

  do ik = sys%st%d%kpt%start, sys%st%d%kpt%end

    do ist = sys%st%st_start, sys%st%st_end

      call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi)
      
      ! start by computing all the wavefunctions acted on by perturbation
      do idir1 = 1, pdim
        call pert_setup_dir(perturbation, idir1)
        call X(pert_apply)(perturbation, sys%namespace, sys%gr, sys%geo, sys%hm, sys%psolver, ik, psi, pertpsi(:, :, idir1))
      end do

      do idir2 = 1, pdim
        do idir1 = 1, pdim

          if (idir1 < idir2) then
            eff_mass_inv(idir1, idir2, ist, ik) = eff_mass_inv(idir2, idir1, ist, ik)
            ! utilizing symmetry of inverse effective mass tensor
            cycle
          end if

          proj_dl_psi(1:mesh%np, 1:sys%hm%d%dim) = lr(1, idir2)%X(dl_psi)(1:mesh%np, 1:sys%hm%d%dim, ist, ik)
          
          ! project out components of other states in degenerate subspace
          do ist2 = 1, sys%st%nst
!              alternate direct method
!              if (abs(sys%st%eigenval(ist2, ik) - sys%st%eigenval(ist, ik)) < degen_thres) then
!                   proj_dl_psi(1:mesh%np) = proj_dl_psi(1:mesh%np) - sys%st%X(psi)(1:mesh%np, 1, ist2, ik) * &
!                     X(mf_dotp)(m, sys%st%X(psi)(1:mesh%np, 1, ist2, ik), proj_dl_psi(1:mesh%np))
            orth_mask(ist2) = .not. (abs(sys%st%eigenval(ist2, ik) - sys%st%eigenval(ist, ik)) < degen_thres)
            ! mask == .false. means do projection; .true. means do not
          end do

!            orth_mask(ist) = .true. ! projection on unperturbed wfn already removed in Sternheimer eqn

          call X(states_elec_orthogonalize_single)(sys%st, mesh, sys%st%nst, ik, proj_dl_psi, mask = orth_mask)

          ! contribution from Sternheimer equation
          term = X(mf_dotp)(mesh, sys%st%d%dim, proj_dl_psi, pertpsi(:, :, idir1))
          eff_mass_inv(idir1, idir2, ist, ik) = M_TWO * R_REAL(term)

          call pert_setup_dir(perturbation, idir1, idir2)
          call X(pert_apply_order_2)(perturbation, sys%namespace, sys%gr, sys%geo, sys%hm, sys%psolver, ik, psi, &
            pertpsi2(1:mesh%np, 1:sys%hm%d%dim))
          eff_mass_inv(idir1, idir2, ist, ik) = &
            eff_mass_inv(idir1, idir2, ist, ik) - R_REAL(X(mf_dotp)(mesh, sys%hm%d%dim, psi, pertpsi2))

        end do !idir2
      end do !idir1
    end do !ist
  end do !ik

#ifdef HAVE_MPI
  if(sys%st%parallel_in_states) then
    call MPI_Allreduce(eff_mass_inv, eff_mass_inv_temp, pdim**2 * sys%st%nst * sys%st%d%nik, &
      MPI_FLOAT, MPI_SUM, sys%st%mpi_grp%comm, mpi_err)
    eff_mass_inv(:,:,:,:) = eff_mass_inv_temp(:,:,:,:)
  end if
  if(sys%st%d%kpt%parallel) then
    call MPI_Allreduce(eff_mass_inv, eff_mass_inv_temp, pdim**2 * sys%st%nst * sys%st%d%nik, &
      MPI_FLOAT, MPI_SUM, sys%st%d%kpt%mpi_grp%comm, mpi_err)
    eff_mass_inv(:,:,:,:) = eff_mass_inv_temp(:,:,:,:)
  end if
  SAFE_DEALLOCATE_A(eff_mass_inv_temp)
#endif

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(pertpsi)
  SAFE_DEALLOCATE_A(pertpsi2)
  SAFE_DEALLOCATE_A(proj_dl_psi)
  SAFE_DEALLOCATE_A(orth_mask)

  POP_SUB(X(calc_eff_mass_inv))

end subroutine X(calc_eff_mass_inv)
  

! ---------------------------------------------------------
!> add projection onto occupied states, by sum over states
subroutine X(kdotp_add_occ)(sys, pert, kdotp_lr, degen_thres)
  type(system_t),      intent(inout) :: sys
  type(pert_t),        intent(in)    :: pert
  type(lr_t),          intent(inout) :: kdotp_lr
  FLOAT,               intent(in)    :: degen_thres

  integer :: ip, idim, ist, ist2, ik
  R_TYPE :: mtxel
  R_TYPE, allocatable :: pertpsi(:, :), psi1(:, :), psi2(:, :)

  PUSH_SUB(X(kdotp_add_occ))

  if(sys%st%parallel_in_states) then
    call messages_not_implemented("kdotp_add_occ parallel in states")
  end if

  SAFE_ALLOCATE(psi1(1:sys%gr%mesh%np_part, 1:sys%st%d%dim))
  SAFE_ALLOCATE(psi2(1:sys%gr%mesh%np_part, 1:sys%st%d%dim))
  SAFE_ALLOCATE(pertpsi(1:sys%gr%mesh%np, 1:sys%st%d%dim))

  do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
    do ist = 1, sys%st%nst

      call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi1)
      
      call X(pert_apply)(pert, sys%namespace, sys%gr, sys%geo, sys%hm, sys%psolver, ik, psi1, pertpsi)
      
      do ist2 = ist + 1, sys%st%nst

        call states_elec_get_state(sys%st, sys%gr%mesh, ist2, ik, psi2)
        
        ! avoid dividing by zero below; these contributions are arbitrary anyway
        if (abs(sys%st%eigenval(ist2, ik) - sys%st%eigenval(ist, ik)) < degen_thres) cycle

        ! the unoccupied subspace was handled by the Sternheimer equation
        if(sys%st%occ(ist2, ik) < M_HALF) cycle
        
        mtxel = X(mf_dotp)(sys%gr%mesh, sys%st%d%dim, psi2, pertpsi(:, :))

        do idim = 1, sys%st%d%dim
          do ip = 1, sys%gr%mesh%np
            kdotp_lr%X(dl_psi)(ip, idim, ist, ik) = kdotp_lr%X(dl_psi)(ip, idim, ist, ik) + &
              psi2(ip, idim)*mtxel/(sys%st%eigenval(ist, ik) - sys%st%eigenval(ist2, ik))
            
            ! note: there is a minus sign here, because the perturbation is an anti-Hermitian operator
            kdotp_lr%X(dl_psi)(ip, idim, ist2, ik) = kdotp_lr%X(dl_psi)(ip, idim, ist2, ik) + &
              psi1(ip, idim)*R_CONJ(-mtxel)/(sys%st%eigenval(ist2, ik) - sys%st%eigenval(ist, ik))
          end do
        end do
        
      end do
    end do
  end do

  SAFE_DEALLOCATE_A(psi1)
  SAFE_DEALLOCATE_A(psi2)
  SAFE_DEALLOCATE_A(pertpsi)

  POP_SUB(X(kdotp_add_occ))
end subroutine X(kdotp_add_occ)


! ---------------------------------------------------------
subroutine X(kdotp_add_diagonal)(sys, em_pert, kdotp_lr)
  type(system_t),      intent(inout) :: sys
  type(pert_t),        intent(inout) :: em_pert
  type(lr_t),          intent(inout) :: kdotp_lr(:)

  integer :: ik, ist, idir
  R_TYPE, allocatable :: psi(:, :), ppsi(:,:)
  R_TYPE :: expectation

  PUSH_SUB(X(kdotp_add_diagonal))

  SAFE_ALLOCATE(ppsi(1:sys%gr%mesh%np, 1:sys%st%d%dim))
  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np_part, 1:sys%st%d%dim))

  do idir = 1, sys%gr%sb%periodic_dim
    call pert_setup_dir(em_pert, idir)
    do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
      do ist = 1, sys%st%nst

        call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi)
        
        call X(pert_apply)(em_pert, sys%namespace, sys%gr, sys%geo, sys%hm, sys%psolver, ik, psi, ppsi)
        
        expectation = X(mf_dotp)(sys%gr%mesh, sys%st%d%dim, psi, ppsi)
        
        kdotp_lr(idir)%X(dl_psi)(1:sys%gr%mesh%np, 1:sys%st%d%dim, ist, ik) = &
          kdotp_lr(idir)%X(dl_psi)(1:sys%gr%mesh%np, 1:sys%st%d%dim, ist, ik) + expectation*psi(1:sys%gr%mesh%np, 1:sys%st%d%dim)
        
      end do
    end do
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(ppsi)
  
  POP_SUB(X(kdotp_add_diagonal))
end subroutine X(kdotp_add_diagonal)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

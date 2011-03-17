!! Copyright (C) 2004 E.S. Kadantsev, M. Marques
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
!! $Id$


! ---------------------------------------------------------
! Orthogonalizes vec against all the occupied states.
! For details on the metallic part, take a look at
! de Gironcoli, PRB 51, 6773 (1995).
! ---------------------------------------------------------
subroutine X(lr_orth_vector) (mesh, st, vec, ist, ik, omega)
  type(mesh_t),        intent(in)    :: mesh
  type(states_t),      intent(in)    :: st
  R_TYPE,              intent(inout) :: vec(:,:)
  integer,             intent(in)    :: ist, ik
  R_TYPE,              intent(in)    :: omega

  integer :: jst
  FLOAT :: xx, theta, theta_ij, theta_ji, alpha_j, dsmear
  R_TYPE :: delta_e
  FLOAT, allocatable :: theta_Fi(:)
  R_TYPE, allocatable :: beta_ij(:)

  PUSH_SUB(X(lr_orth_vector))

  dsmear = max(CNST(1e-14), st%smear%dsmear)

  SAFE_ALLOCATE(beta_ij(1:st%nst))

  if(smear_is_semiconducting(st%smear) .or. st%smear%integral_occs) then
    theta = st%occ(ist, ik) / st%smear%el_per_state
    do jst = 1, st%nst
      if(abs(st%occ(ist, ik) * st%occ(jst, ik)) .gt. M_EPSILON) then
        beta_ij(jst) = st%occ(jst, ik) / st%smear%el_per_state
      else
        beta_ij(jst) = M_ZERO
      end if
    end do
  else
    SAFE_ALLOCATE(theta_Fi(1:st%nst))
    theta_Fi(1:st%nst) = st%occ(1:st%nst, ik) / st%smear%el_per_state

    do jst = 1, st%nst
      if(st%smear%method .eq. SMEAR_FIXED_OCC) then
        theta_ij = theta_Fi(ist) / (theta_Fi(ist) + theta_Fi(jst))
        theta_ji = theta_Fi(jst) / (theta_Fi(ist) + theta_Fi(jst))
      else
        xx = (st%eigenval(ist, ik) - st%eigenval(jst, ik))/dsmear
        theta_ij = smear_step_function(st%smear,  xx)
        theta_ji = smear_step_function(st%smear, -xx)
        ! In principle, theta_ji = 1 - theta_ji as approximations to a step function,
        ! but in practice, only fermi_dirac, semiconducting, and fixed_occ formulations
        ! satisfy this.
      endif
      
      if(abs(st%occ(ist, ik) * st%occ(jst, ik)) .lt. M_EPSILON) then
        ! Supposedly beta_ij = 0 if i or j is unoccupied. In practice, this may not be
        ! true with the general formula, and must be explicitly enforced.
        beta_ij(jst) = M_ZERO
      else  
        beta_ij(jst) = theta_Fi(ist)*Theta_ij + Theta_Fi(jst)*Theta_ji
          
        alpha_j = lr_alpha_j(st, jst, ik)
        delta_e = st%eigenval(ist, ik) - st%eigenval(jst, ik) - omega
        
        if(abs(delta_e) >= CNST(1e-5)) then
          beta_ij(jst) = beta_ij(jst) + alpha_j*Theta_ji*(Theta_Fi(ist) - Theta_Fi(jst))/delta_e
        else
          if(st%smear%method .ne. SMEAR_FIXED_OCC) then 
            xx = (st%smear%e_fermi - st%eigenval(ist, ik) + CNST(1e-14))/dsmear
            beta_ij(jst) = beta_ij(jst) + alpha_j*Theta_ji*(smear_delta_function(st%smear, xx)/dsmear)
          endif
        end if
      endif

    end do

    theta = theta_Fi(ist)
    SAFE_DEALLOCATE_A(theta_Fi)
  end if

  call X(states_orthogonalization)(mesh, st%nst, st%d%dim, st%X(psi)(:, :, :, ik), vec(:, :), &
    Theta_Fi=theta, beta_ij=beta_ij)

  SAFE_DEALLOCATE_A(beta_ij)

  POP_SUB(X(lr_orth_vector))

end subroutine X(lr_orth_vector)


! --------------------------------------------------------------------
subroutine X(lr_build_dl_rho) (mesh, st, lr, nsigma)
  type(mesh_t),   intent(in)    :: mesh
  type(states_t), intent(in)    :: st
  type(lr_t),     intent(inout) :: lr(:)
  integer,        intent(in)    :: nsigma

  integer :: ip, ist, ik, ispin, isigma
  FLOAT   :: weight
  CMPLX   :: cc
  R_TYPE  :: dd

  PUSH_SUB(X(lr_build_dl_rho))

  if(st%d%ispin == SPINORS) then
    message(1) = "Not yet implemented - please fix me"
    call write_fatal(1)
  end if
  
  ! initialize density
  do isigma = 1, nsigma
    lr(isigma)%X(dl_rho)(:, :) = M_ZERO
  end do

  ! calculate density
  do ik = st%d%kpt%start, st%d%kpt%end
    ispin = states_dim_get_spin_index(st%d, ik)
    do ist  = st%st_start, st%st_end
      weight = st%d%kweights(ik)*st%smear%el_per_state
      
      if(nsigma == 1) then  ! either omega purely real or purely imaginary
        do ip = 1, mesh%np
          dd = weight*st%X(psi)(ip, 1, ist, ik)*R_CONJ(lr(1)%X(dl_psi)(ip, 1, ist, ik))
          lr(1)%X(dl_rho)(ip, ispin) = lr(1)%X(dl_rho)(ip, ispin) + dd + R_CONJ(dd)
        end do
      else
        do ip = 1, mesh%np
          cc = weight*(R_CONJ(st%X(psi)(ip, 1, ist, ik))*lr(1)%X(dl_psi)(ip, 1, ist, ik) + &
            st%X(psi)(ip, 1, ist, ik)*R_CONJ(lr(2)%X(dl_psi)(ip, 1, ist, ik)))
          lr(1)%X(dl_rho)(ip, ispin) = lr(1)%X(dl_rho)(ip, ispin) + cc
          lr(2)%X(dl_rho)(ip, ispin) = lr(2)%X(dl_rho)(ip, ispin) + R_CONJ(cc)
        end do
      end if

    end do
  end do

  ! reduce
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    do isigma = 1, nsigma
      call comm_allreduce(st%st_kpt_mpi_grp%comm, lr(isigma)%X(dl_rho), dim = (/mesh%np, st%d%nspin/))
    end do
  end if
      
  POP_SUB(X(lr_build_dl_rho))
end subroutine X(lr_build_dl_rho)


! ---------------------------------------------------------
! Orthogonalizes response of \alpha KS orbital to all occupied
! \alpha KS orbitals.
! ---------------------------------------------------------
subroutine X(lr_orth_response)(mesh, st, lr, omega)
  type(mesh_t),   intent(in)    :: mesh
  type(states_t), intent(in)    :: st
  type(lr_t),     intent(inout) :: lr
  R_TYPE,         intent(in)    :: omega
  
  integer :: ist, ik
  PUSH_SUB(X(lr_orth_response))
  
  do ik = 1, st%d%nik
    do ist = 1, st%nst
      call X(lr_orth_vector) (mesh, st, lr%X(dl_psi)(:,:, ist, ik), ist, ik, omega)
    end do
  end do
  
  POP_SUB(X(lr_orth_response))
end subroutine X(lr_orth_response)


! ---------------------------------------------------------
subroutine X(lr_swap_sigma)(st, mesh, plus, minus)
  type(states_t), intent(in)    :: st
  type(mesh_t),   intent(in)    :: mesh
  type(lr_t),     intent(inout) :: plus
  type(lr_t),     intent(inout) :: minus

  integer :: ik, idim, ist
  R_TYPE, allocatable :: tmp(:)

  PUSH_SUB(X(lr_swap_sigma))

  SAFE_ALLOCATE(tmp(1:mesh%np))

  do ik = 1, st%d%nspin
    call lalg_copy(mesh%np, plus%X(dl_rho)(:, ik), tmp(:))
    call lalg_copy(mesh%np, minus%X(dl_rho)(:, ik), plus%X(dl_rho)(:, ik))
    call lalg_copy(mesh%np, tmp(:), minus%X(dl_rho)(:, ik))
  enddo

  do ik = 1, st%d%nik
    do ist = 1, st%nst
      do idim = 1, st%d%dim
        call lalg_copy(mesh%np_part, plus%X(dl_psi)(:, idim, ist, ik), tmp(:))
        call lalg_copy(mesh%np_part, minus%X(dl_psi)(:, idim, ist, ik), plus%X(dl_psi)(:, idim, ist, ik))
        call lalg_copy(mesh%np_part, tmp(:), minus%X(dl_psi)(:, idim, ist, ik))
      end do
    end do
  end do

  SAFE_DEALLOCATE_A(tmp)
  POP_SUB(X(lr_swap_sigma))

end subroutine X(lr_swap_sigma)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

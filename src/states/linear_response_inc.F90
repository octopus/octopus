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
subroutine X(lr_orth_vector) (mesh, st, vec, ist, ik)
  type(mesh_t),        intent(in)    :: mesh
  type(states_t),      intent(in)    :: st
  R_TYPE,              intent(inout) :: vec(:,:)
  integer,             intent(in)    :: ist, ik

  integer :: jst
  FLOAT :: xx, theta_ij, theta_ji, alpha_j, delta_e, dsmear
  FLOAT, allocatable :: theta_Fi(:), beta_ij(:)

  call push_sub('linear_response_inc.Xlr_orth_vector')

  dsmear = max(CNST(1e-14), st%smear%dsmear)

  SAFE_ALLOCATE(theta_Fi(1:st%nst))
  do jst = 1, st%nst
    ! epsilon has to be added or we have problem with semiconducting smearing
    xx = (st%smear%e_fermi - st%eigenval(jst, ik) + CNST(1e-14))/dsmear
    theta_Fi(jst) = smear_step_function(st%smear, xx)
  end do

  SAFE_ALLOCATE(beta_ij(1:st%nst))

  if(st%smear%fixed_occ .or. st%smear%method == SMEAR_SEMICONDUCTOR) then
    do jst = 1, st%nst
      if(Theta_Fi(ist).ne.M_ZERO .and. Theta_Fi(jst).ne.M_ZERO) then
        beta_ij(jst) = Theta_Fi(jst)
      else
        beta_ij(jst) = M_ZERO
      end if
    end do
  else
    beta_ij = M_ZERO
    do jst = 1, st%nst
      xx = (st%eigenval(ist, ik) - st%eigenval(jst, ik))/dsmear
      theta_ij = smear_step_function(st%smear,  xx)
      theta_ji = smear_step_function(st%smear, -xx)
      
      beta_ij(jst) = theta_Fi(ist)*Theta_ij + Theta_Fi(jst)*Theta_ji
        
      alpha_j = max(st%smear%e_fermi + M_THREE*dsmear - st%eigenval(jst, ik), M_ZERO)
      delta_e = st%eigenval(ist, ik) - st%eigenval(jst, ik)
      
      if(abs(delta_e) >= CNST(1e-5)) then
        beta_ij(jst) = beta_ij(jst) + alpha_j*Theta_ji*(Theta_Fi(ist) - Theta_Fi(jst))/delta_e
      else
        xx = (st%smear%e_fermi - st%eigenval(ist, ik) + CNST(1e-14))/dsmear
        beta_ij(jst) = beta_ij(jst) + alpha_j*Theta_ji*  &
          (-smear_delta_function(st%smear,  xx)/dsmear)
      end if
    end do
  end if

  call X(states_gram_schmidt)(mesh, st%nst, st%d%dim, st%X(psi)(:, :, :, ik), vec(:, :), &
    Theta_Fi=Theta_Fi(ist), beta_ij=beta_ij)

  SAFE_DEALLOCATE_A(beta_ij)
  SAFE_DEALLOCATE_A(Theta_Fi)

  call pop_sub('linear_response_inc.Xlr_orth_vector')

end subroutine X(lr_orth_vector)


! --------------------------------------------------------------------
subroutine X(lr_build_dl_rho) (mesh, st, lr, nsigma)
  type(mesh_t),   intent(in)    :: mesh
  type(states_t), intent(in)    :: st
  type(lr_t),     intent(inout) :: lr(:)
  integer,        intent(in)    :: nsigma

  integer :: ip, ist, ik, ik2, sp, isigma
  CMPLX   :: cc
  R_TYPE  :: dd

  call push_sub('linear_response_inc.Xlr_build_dl_rho')

  if(st%d%ispin == SPINORS) then
    message(1) = "Not yet implemented - please fix me"
    call write_fatal(1)
  end if
  
  ! initialize density
  do isigma = 1, nsigma
    lr(isigma)%X(dl_rho)(:, :) = M_ZERO
  end do

  sp = 1
  if(st%d%ispin == SPIN_POLARIZED) sp = 2

  ! calculate density
  do ik = 1, st%d%nik, sp
    do ist  = st%st_start, st%st_end
      do ip = 1, mesh%np

        do ik2 = ik, ik+sp-1 ! this loop takes care of the SPIN_POLARIZED case
          dd = st%d%kweights(ik2) * st%smear%el_per_state

          if(nsigma == 1) then  ! either omega purely real or purely imaginary
            dd = dd * st%X(psi)(ip, 1, ist, ik2)*R_CONJ(lr(1)%X(dl_psi)(ip, 1, ist, ik2))
            lr(1)%X(dl_rho)(ip, 1) = lr(1)%X(dl_rho)(ip, 1) + dd + R_CONJ(dd)
          else
            cc = dd * (                                                             &
              R_CONJ(st%X(psi)(ip, 1, ist, ik2))*lr(1)%X(dl_psi)(ip, 1, ist, ik2) + &
              st%X(psi)(ip, 1, ist, ik2)*R_CONJ(lr(2)%X(dl_psi)(ip, 1, ist, ik2)))
            lr(1)%X(dl_rho)(ip, 1) = lr(1)%X(dl_rho)(ip, 1) + cc
            lr(2)%X(dl_rho)(ip, 1) = lr(2)%X(dl_rho)(ip, 1) + R_CONJ(cc)
          end if
        end do

      end do
    end do
  end do

  call pop_sub('linear_response_inc.Xlr_build_dl_rho')
end subroutine X(lr_build_dl_rho)


! ---------------------------------------------------------
! Orthogonalizes response of \alpha KS orbital to all occupied
! \alpha KS orbitals.
! ---------------------------------------------------------
subroutine X(lr_orth_response)(mesh, st, lr)
  type(mesh_t),   intent(in)    :: mesh
  type(states_t), intent(in)    :: st
  type(lr_t),     intent(inout) :: lr
  
  integer :: ist, ik
  call push_sub('linear_response_inc.Xlr_orth_response')
  
  do ik = 1, st%d%nik
    do ist = 1, st%nst
      call X(lr_orth_vector) (mesh, st, lr%X(dl_psi)(:,:, ist, ik), ist, ik)
    end do
  end do
  
  call pop_sub('linear_response_inc.Xlr_orth_response')
end subroutine X(lr_orth_response)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

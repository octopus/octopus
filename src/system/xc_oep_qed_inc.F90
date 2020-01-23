!! Copyright (C) 2017 Johannes Flick
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
!> This file handles the evaluation of the photon-OEP potential,
!! as described in J. Flick et al. ACS Photonics 2018, 5, 3, 992-1005
! ---------------------------------------------------------
subroutine X(xc_oep_pt_phi) (namespace, gr, hm, st, is, oep, phi1)
  type(namespace_t),        intent(in)    :: namespace
  type(grid_t),             intent(in)    :: gr
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t),      intent(in)    :: st
  integer,                  intent(in)    :: is
  type(xc_oep_t),           intent(inout) :: oep
  R_TYPE,                   intent(inout) :: phi1(:,:,:)

  integer :: ist, kst, iter_used
  FLOAT :: rhs_kkbar, residue, kkopii
  R_TYPE, allocatable :: rhs(:,:), psiii(:, :), psikk(:, :)
  FLOAT, allocatable :: rhs2(:)

  PUSH_SUB(X(xc_oep_pt_phi))

  SAFE_ALLOCATE(rhs(1:gr%mesh%np, 1:1))
  SAFE_ALLOCATE(psiii(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psikk(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(rhs2(1:gr%mesh%np))
  
  if (is == 1) then !only initialize for first spin channel
    oep%pt%pt_number = M_ZERO
    oep%pt%correlator = M_ZERO
    oep%pt%ex = M_ZERO
  end if

  do ist = st%st_start, oep%noccst
    call states_elec_get_state(st, gr%mesh, ist, is, psiii)
    rhs(:,1) = -oep%pt%lambda_array(1)*sqrt(M_HALF*oep%pt%omega_array(1))*oep%pt%pol_dipole_array(:,1)*psiii(:,1)

    do kst = st%st_start, oep%noccst
      call states_elec_get_state(st, gr%mesh, kst, is, psikk)
      rhs_kkbar = X(mf_dotp)(gr%mesh, R_CONJ(psikk(:,1)), rhs(:,1))
      call lalg_axpy(gr%mesh%np, -rhs_kkbar, psikk(:,1), rhs(:,1))
    end do

    call X(states_elec_orthogonalize_single)(st, gr%mesh, st%nst, is, rhs, normalize = .false.)

    call X(linear_solver_solve_HXeY)(oep%solver, namespace, hm, gr, st, ist, is, oep%pt%lr%X(dl_psi)(:,:, ist, is), rhs, &
           R_TOTYPE(-st%eigenval(ist, is)) + R_REAL(oep%pt%omega_array(1)), oep%scftol%final_tol, residue, iter_used)

    call X(states_elec_orthogonalize_single)(st, gr%mesh, st%nst, is, &
      oep%pt%lr%X(dl_psi)(:,:, ist, is), normalize = .false.)

    phi1(:, 1:st%d%dim, ist) = oep%pt%lr%X(dl_psi)(:, 1:st%d%dim, ist, is)

    oep%pt%pt_number = oep%pt%pt_number + st%occ(ist, is)*X(mf_dotp)(gr%mesh, R_CONJ(phi1(:, 1, ist)), &
                phi1(:, 1, ist))

    oep%pt%ex = oep%pt%ex + st%occ(ist, is)*oep%pt%lambda_array(1)*sqrt(M_HALF*oep%pt%omega_array(1))* &
              X(mf_dotp)(gr%mesh, R_CONJ(phi1(:, 1, ist)), &
                oep%pt%pol_dipole_array(:,1)*psiii(:,1))

    oep%pt%ex = oep%pt%ex + st%occ(ist, is)*M_HALF*oep%pt%lambda_array(1)**2* &
              X(mf_dotp)(gr%mesh, &
                R_CONJ(oep%pt%pol_dipole_array(:,1)*psiii(:,1)), &
                oep%pt%pol_dipole_array(:,1)*psiii(:,1))

    do kst = st%st_start, oep%noccst
      call states_elec_get_state(st, gr%mesh, kst, is, psikk)
      kkopii = oep%pt%lambda_array(1)*X(mf_dotp)(gr%mesh, R_CONJ(psiii(:,1)), &
        oep%pt%pol_dipole_array(:,1)*psikk(:,1))
      oep%pt%ex = oep%pt%ex - st%occ(kst, is)*M_HALF*(kkopii)**2
    end do

    ! calculate correlator function
    rhs2(:) = psiii(:,1)*phi1(:, 1, ist)
    call lalg_axpy(gr%mesh%np, st%occ(ist, is), rhs2(:), oep%pt%correlator(:,1))

  end do

  SAFE_DEALLOCATE_A(rhs)
  SAFE_DEALLOCATE_A(rhs2)
  SAFE_DEALLOCATE_A(psiii)
  SAFE_DEALLOCATE_A(psikk)
  POP_SUB(X(xc_oep_pt_phi))
end subroutine X(xc_oep_pt_phi)

! ---------------------------------------------------------
subroutine X(xc_oep_pt_rhs) (gr, st, is, oep, phi1, ist, rhs)
  type(grid_t),        intent(in)    :: gr
  type(states_elec_t),      intent(in)    :: st
  integer,             intent(in)    :: is
  type(xc_oep_t),      intent(inout) :: oep
  R_TYPE,              intent(in)    :: phi1(:,:,:)
  integer,             intent(in)    :: ist
  R_TYPE,              intent(inout) :: rhs(:,:)

  integer :: kst
  FLOAT :: abar, kkopii
  R_TYPE, allocatable :: aa(:,:), psiii(:, :), psikk(:, :)

  PUSH_SUB(X(xc_oep_pt_rhs))

  SAFE_ALLOCATE(aa(1:gr%mesh%np, 1:1))
  SAFE_ALLOCATE(psiii(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psikk(1:gr%mesh%np, 1:st%d%dim))

  call states_elec_get_state(st, gr%mesh, ist, is, psiii)

  aa(:,1) = M_HALF*oep%pt%lambda_array(1)**2*oep%pt%pol_dipole_array(:,1)**2 &
            *R_CONJ(psiii(:,1))
  call lalg_axpy(gr%mesh%np, sqrt(M_HALF*oep%pt%omega_array(1))*oep%pt%lambda_array(1), &
    oep%pt%pol_dipole_array(:,1)*R_CONJ(phi1(:, 1, ist)), aa(:,1))

  do kst = st%st_start, oep%noccst
    call states_elec_get_state(st, gr%mesh, kst, is, psikk)
    kkopii = oep%pt%lambda_array(1)*X(mf_dotp)(gr%mesh, R_CONJ(psiii(:,1)), &
      oep%pt%pol_dipole_array(:,1)*psikk(:,1))
    call lalg_axpy(gr%mesh%np, - oep%pt%lambda_array(1)*kkopii, oep%pt%pol_dipole_array(:,1)*R_CONJ(psikk(:,1)), aa(:,1))
    call lalg_axpy(gr%mesh%np, - sqrt(M_HALF*oep%pt%omega_array(1))*kkopii, R_CONJ(phi1(:, 1, kst)), aa(:,1))
  end do

  if (ist/=(oep%eigen_n + 1) .or. (oep%level == XC_OEP_FULL)) then
    abar = X(mf_dotp)(gr%mesh,  aa(:,1), psiii(:,1))
    call lalg_axpy(gr%mesh%np, -abar, R_CONJ(psiii(:,1)), aa(:,1))
  end if

  call lalg_axpy(gr%mesh%np, M_ONE, aa(:,1), rhs(:,1))

  SAFE_DEALLOCATE_A(aa)
  SAFE_DEALLOCATE_A(psiii)
  SAFE_DEALLOCATE_A(psikk)
  POP_SUB(X(xc_oep_pt_rhs))
end subroutine X(xc_oep_pt_rhs)


! ---------------------------------------------------------
subroutine X(xc_oep_pt_inhomog) (gr, st, is, oep, phi1, ist, ss)
  type(grid_t),        intent(in)    :: gr
  type(states_elec_t),      intent(in)    :: st
  integer,             intent(in)    :: is
  type(xc_oep_t),      intent(in)    :: oep
  R_TYPE,              intent(in)    :: phi1(:,:,:)
  integer,             intent(in)    :: ist
  FLOAT,               intent(inout) :: ss(:)

  FLOAT :: phi_bar
  R_TYPE, allocatable :: psiii(:, :)
  FLOAT, allocatable  :: rhs(:)

  PUSH_SUB(X(xc_oep_pt_inhomog))

  SAFE_ALLOCATE(psiii(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(rhs(1:gr%mesh%np))

  call states_elec_get_state(st, gr%mesh, ist, is, psiii)

  phi_bar = X(mf_dotp)(gr%mesh, R_CONJ(phi1(:, 1, ist)), &
                phi1(:, 1, ist))
  rhs = phi1(:, 1, ist)*phi1(:, 1, ist) - phi_bar*(R_CONJ(psiii(:,1))*psiii(:,1))
  call lalg_axpy(gr%mesh%np, -M_ONE, rhs, ss)
  
  SAFE_DEALLOCATE_A(psiii)
  SAFE_DEALLOCATE_A(rhs)
  POP_SUB(X(xc_oep_pt_inhomog))
end subroutine X(xc_oep_pt_inhomog)

! ---------------------------------------------------------
subroutine X(xc_oep_pt_uxcbar) (gr, st, is, oep, phi1, ist, vxbar)
  type(grid_t),        intent(in)    :: gr
  type(states_elec_t),      intent(in)    :: st
  integer,             intent(in)    :: is
  type(xc_oep_t),      intent(in)    :: oep
  R_TYPE,              intent(in)    :: phi1(:,:,:)
  integer,             intent(in)    :: ist
  FLOAT,               intent(inout) :: vxbar

  integer :: kst
  FLOAT :: kkopii, result1
  R_TYPE, allocatable :: aa(:,:), psiii(:,:), psikk(:,:)

  PUSH_SUB(X(xc_oep_pt_uxcbar))

  SAFE_ALLOCATE(aa(1:gr%mesh%np, 1:1))
  SAFE_ALLOCATE(psiii(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psikk(1:gr%mesh%np, 1:st%d%dim))

  call states_elec_get_state(st, gr%mesh, ist, is, psiii)

  aa(:,1) = oep%pt%lambda_array(1)*oep%pt%pol_dipole_array(:,1)*psiii(:,1)
  result1 = sqrt(M_HALF*oep%pt%omega_array(1))*X(mf_dotp)(gr%mesh, R_CONJ(phi1(:,1,ist)), aa(:,1))
  result1 = result1 + M_HALF*dmf_dotp(gr%mesh, R_REAL(aa(:,1)), R_REAL(aa(:,1)))

   do kst = st%st_start, oep%noccst
     call states_elec_get_state(st, gr%mesh, kst, is, psikk)
     kkopii = oep%pt%lambda_array(1)*X(mf_dotp)(gr%mesh, R_CONJ(psiii(:,1)), &
     oep%pt%pol_dipole_array(:,1)*psikk(:,1))
     result1 = result1 - M_HALF*kkopii**2
   end do

  vxbar = vxbar - result1

  SAFE_DEALLOCATE_A(aa)
  SAFE_DEALLOCATE_A(psiii)
  SAFE_DEALLOCATE_A(psikk)
  POP_SUB(X(xc_oep_pt_uxcbar))
end subroutine X(xc_oep_pt_uxcbar)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

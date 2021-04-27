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
subroutine xc_oep_pt_phi(namespace, mesh, hm, st, is, oep, phi1)
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t),      intent(in)    :: st
  integer,                  intent(in)    :: is
  type(xc_oep_t),           intent(inout) :: oep
  FLOAT,                    intent(inout) :: phi1(:,:,:)

  integer :: ist, kst, iter_used
  FLOAT :: rhs_kkbar, residue, kkopii
  FLOAT, allocatable :: rhs(:,:), rhs2(:), psiii(:, :), psikk(:, :), pol_dip_psiii(:, :)

  PUSH_SUB(xc_oep_pt_phi)

  SAFE_ALLOCATE(rhs(1:mesh%np, 1:1))
  SAFE_ALLOCATE(psiii(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(pol_dip_psiii(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psikk(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(rhs2(1:mesh%np))
  
  if (is == 1) then !only initialize for first spin channel
    oep%pt%number = M_ZERO
    oep%pt%correlator = M_ZERO
    oep%pt%ex = M_ZERO
  end if

  do ist = st%st_start, oep%noccst
    call states_elec_get_state(st, mesh, ist, is, psiii)
    pol_dip_psiii(:, 1) = oep%pt%pol_dipole(:, 1)*psiii(:, 1)
    rhs(:,1) = -oep%pt%lambda(1)*sqrt(M_HALF*oep%pt%omega(1))*pol_dip_psiii(:, 1)

    do kst = st%st_start, oep%noccst
      call states_elec_get_state(st, mesh, kst, is, psikk)
      rhs_kkbar = dmf_dotp(mesh, psikk(:, 1), rhs(:, 1))
      call lalg_axpy(mesh%np, -rhs_kkbar, psikk(:, 1), rhs(:, 1))
    end do

    call dstates_elec_orthogonalize_single(st, mesh, st%nst, is, rhs, normalize = .false.)

    call dlinear_solver_solve_HXeY(oep%solver, namespace, hm, mesh, st, ist, is, oep%photon_lr%ddl_psi(:, :, ist, is), rhs, &
           -st%eigenval(ist, is) + real(oep%pt%omega(1)), oep%scftol%final_tol, residue, iter_used)

    call dstates_elec_orthogonalize_single(st, mesh, st%nst, is, oep%photon_lr%ddl_psi(:, :, ist, is), normalize = .false.)

    phi1(1:mesh%np, 1:st%d%dim, ist) = oep%photon_lr%ddl_psi(1:mesh%np, 1:st%d%dim, ist, is)

    oep%pt%number(1) = oep%pt%number(1) + st%occ(ist, is)*dmf_dotp(mesh, phi1(:, 1, ist), phi1(:, 1, ist))

    oep%pt%ex = oep%pt%ex + st%occ(ist, is)*oep%pt%lambda(1)*sqrt(M_HALF*oep%pt%omega(1))* &
      dmf_dotp(mesh, phi1(:, 1, ist), pol_dip_psiii(:, 1))

    oep%pt%ex = oep%pt%ex + st%occ(ist, is)*M_HALF*oep%pt%lambda(1)**2*dmf_dotp(mesh, pol_dip_psiii(:, 1), pol_dip_psiii(:, 1))

    do kst = st%st_start, oep%noccst
      call states_elec_get_state(st, mesh, kst, is, psikk)
      kkopii = oep%pt%lambda(1)*dmf_dotp(mesh, psiii(:, 1), oep%pt%pol_dipole(:, 1)*psikk(:, 1))
      oep%pt%ex = oep%pt%ex - st%occ(kst, is)*M_HALF*kkopii**2
    end do

    ! calculate correlator function
    rhs2(:) = psiii(:,1)*phi1(:, 1, ist)
    call lalg_axpy(mesh%np, st%occ(ist, is), rhs2(:), oep%pt%correlator(:, 1))

  end do

  SAFE_DEALLOCATE_A(rhs)
  SAFE_DEALLOCATE_A(rhs2)
  SAFE_DEALLOCATE_A(psiii)
  SAFE_DEALLOCATE_A(psikk)

  POP_SUB(xc_oep_pt_phi)
end subroutine xc_oep_pt_phi

! ---------------------------------------------------------
subroutine xc_oep_pt_rhs(mesh, st, is, oep, phi1, ist, rhs)
  type(mesh_t),        intent(in)    :: mesh
  type(states_elec_t), intent(in)    :: st
  integer,             intent(in)    :: is
  type(xc_oep_t),      intent(inout) :: oep
  FLOAT,               intent(in)    :: phi1(:,:,:)
  integer,             intent(in)    :: ist
  FLOAT,               intent(inout) :: rhs(:,:)

  integer :: kst
  FLOAT :: abar, kkopii
  FLOAT, allocatable :: aa(:,:), psiii(:, :), psikk(:, :)

  PUSH_SUB(xc_oep_pt_rhs)

  SAFE_ALLOCATE(aa(1:mesh%np, 1:1))
  SAFE_ALLOCATE(psiii(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psikk(1:mesh%np, 1:st%d%dim))

  call states_elec_get_state(st, mesh, ist, is, psiii)

  aa(:,1) = M_HALF*oep%pt%lambda(1)**2*oep%pt%pol_dipole(:, 1)**2 * psiii(:, 1)
  call lalg_axpy(mesh%np, sqrt(M_HALF*oep%pt%omega(1))*oep%pt%lambda(1), oep%pt%pol_dipole(:, 1)*phi1(:, 1, ist), aa(:, 1))

  do kst = st%st_start, oep%noccst
    call states_elec_get_state(st, mesh, kst, is, psikk)
    kkopii = oep%pt%lambda(1)*dmf_dotp(mesh, psiii(:, 1), oep%pt%pol_dipole(:, 1)*psikk(:, 1))
    call lalg_axpy(mesh%np, - oep%pt%lambda(1)*kkopii, oep%pt%pol_dipole(:, 1)*psikk(:, 1), aa(:, 1))
    call lalg_axpy(mesh%np, - sqrt(M_HALF*oep%pt%omega(1))*kkopii, phi1(:, 1, kst), aa(:,1))
  end do

  if (ist /= oep%eigen_n + 1 .or. oep%level == XC_OEP_FULL) then
    abar = dmf_dotp(mesh,  aa(:, 1), psiii(:, 1))
    call lalg_axpy(mesh%np, -abar, psiii(:, 1), aa(:, 1))
  end if

  call lalg_axpy(mesh%np, M_ONE, aa(:, 1), rhs(:, 1))

  SAFE_DEALLOCATE_A(aa)
  SAFE_DEALLOCATE_A(psiii)
  SAFE_DEALLOCATE_A(psikk)

  POP_SUB(xc_oep_pt_rhs)
end subroutine xc_oep_pt_rhs


! ---------------------------------------------------------
subroutine xc_oep_pt_inhomog(mesh, st, is, phi1, ist, ss)
  type(mesh_t),        intent(in)    :: mesh
  type(states_elec_t), intent(in)    :: st
  integer,             intent(in)    :: is
  FLOAT,               intent(in)    :: phi1(:,:,:)
  integer,             intent(in)    :: ist
  FLOAT,               intent(inout) :: ss(:)

  FLOAT :: phi_bar
  FLOAT, allocatable  :: psiii(:, :), rhs(:)

  PUSH_SUB(xc_oep_pt_inhomog)

  SAFE_ALLOCATE(psiii(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(rhs(1:mesh%np))

  call states_elec_get_state(st, mesh, ist, is, psiii)

  phi_bar = dmf_dotp(mesh, phi1(:, 1, ist), phi1(:, 1, ist))
  rhs = phi1(:, 1, ist)*phi1(:, 1, ist) - phi_bar*psiii(:, 1)**2
  call lalg_axpy(mesh%np, -M_ONE, rhs, ss)
  
  SAFE_DEALLOCATE_A(psiii)
  SAFE_DEALLOCATE_A(rhs)

  POP_SUB(xc_oep_pt_inhomog)
end subroutine xc_oep_pt_inhomog

! ---------------------------------------------------------
subroutine xc_oep_pt_uxcbar(mesh, st, is, oep, phi1, ist, vxbar)
  type(mesh_t),        intent(in)    :: mesh
  type(states_elec_t), intent(in)    :: st
  integer,             intent(in)    :: is
  type(xc_oep_t),      intent(in)    :: oep
  FLOAT,               intent(in)    :: phi1(:,:,:)
  integer,             intent(in)    :: ist
  FLOAT,               intent(inout) :: vxbar

  integer :: kst
  FLOAT :: kkopii, result1
  FLOAT, allocatable :: aa(:,:), psiii(:,:), psikk(:,:)

  PUSH_SUB(xc_oep_pt_uxcbar)

  SAFE_ALLOCATE(aa(1:mesh%np, 1:1))
  SAFE_ALLOCATE(psiii(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psikk(1:mesh%np, 1:st%d%dim))

  call states_elec_get_state(st, mesh, ist, is, psiii)

  aa(:, 1) = oep%pt%lambda(1)*oep%pt%pol_dipole(:, 1)*psiii(:, 1)
  result1 = sqrt(M_HALF*oep%pt%omega(1))*dmf_dotp(mesh, phi1(:, 1, ist), aa(:, 1))
  result1 = result1 + M_HALF*dmf_dotp(mesh, aa(:, 1), aa(:, 1))

  do kst = st%st_start, oep%noccst
    call states_elec_get_state(st, mesh, kst, is, psikk)
    kkopii = oep%pt%lambda(1)*dmf_dotp(mesh, psiii(:, 1), oep%pt%pol_dipole(:, 1)*psikk(:, 1))
    result1 = result1 - M_HALF*kkopii**2
  end do

  vxbar = vxbar - result1

  SAFE_DEALLOCATE_A(aa)
  SAFE_DEALLOCATE_A(psiii)
  SAFE_DEALLOCATE_A(psikk)

  POP_SUB(xc_oep_pt_uxcbar)
end subroutine xc_oep_pt_uxcbar

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

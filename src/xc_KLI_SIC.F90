!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

!!! This routine calculates the SIC exchange functional. note that the LDA
!!! part of the functional is already included by the LDA routines
subroutine R_FUNC(kli_x_sic) (xcs, m, st, is, socc, sfact, lx_, ex)
  type(xc_type),     intent(in)  :: xcs
  type(mesh_type),   intent(in)  :: m
  type(states_type), intent(inout)  :: st
  integer,           intent(in)  :: is
  real(r8),          intent(in)  :: socc, sfact
  R_TYPE,            intent(out) :: lx_(m%np, st%nst)
  real(r8),          intent(out) :: ex

  integer  :: i, k, i1, i2
  real(r8) :: ex2, edummy, r
  real(r8), allocatable :: vx2(:, :)
  real(r8), pointer :: rho(:,:), rho_save(:,:)
  type(xc_type) :: xcs2

  allocate(rho(m%np, 2), Vx2(m%np, 2))
  i1 = st%spin_channels; st%spin_channels = 2
  i2 = st%nspin;         st%nspin         = 2
  rho_save => st%rho;    st%rho => rho

! calculate the u_sij using poissons equation
  rho(:, 2) = M_ZERO
  
  xcs2 = xcs
  xcs2%family = XC_FAMILY_LDA
  xcs2%functl = X_FUNC_LDA_NREL
  do i = 1, st%nst
    if(st%occ(i, is) .gt. small) then ! we only need the occupied states
      vx2 = M_ZERO
      ex2 = M_ZERO

      st%rho(:, 1) = socc*st%occ(i, is)*R_ABS(st%X(psi)(:, 1, i, is))**2

      call xc_lda (xcs2, m, st, vx2, ex2, edummy)

      ex = ex - sfact*ex2
      
      vx2(:, 2) = M_ZERO
      call poisson_solve(m, vx2(:, 2), rho(:, 1))
      lx_(:, i) = lx_(:, i) - (vx2(:, 1) + vx2(:,2))*R_CONJ(st%X(psi) (:, 1, i, is))

      ex = ex - M_HALF*sfact*socc*st%occ(i, is)* &
           sum(vx2(:, 2)* R_ABS(st%X(psi)(:, 1, i, is))**2)*m%vol_pp
    end if
  end do

  st%spin_channels = i1
  st%nspin         = i2    
  st%rho => rho_save
  deallocate(rho, Vx2)

end subroutine X(kli_x_sic)

subroutine X(kli_c_sic) (xcs, m, st, is, socc, sfact, lc_, ec)
  type(xc_type),     intent(in)  :: xcs
  type(mesh_type),   intent(in)  :: m
  type(states_type), intent(inout)  :: st
  integer,           intent(in)  :: is
  real(r8),          intent(in)  :: socc, sfact
  R_TYPE,            intent(out) :: lc_(m%np, st%nst)
  real(r8),          intent(out) :: ec

  integer  :: i, k, i1, i2
  real(r8) :: ec2, edummy, r
  real(r8), allocatable :: vc2(:, :)
  real(r8), pointer :: rho(:,:), rho_save(:,:)
  type(xc_type) :: xcs2

  allocate(rho(m%np, 2), Vc2(m%np, 2))
  i1 = st%spin_channels; st%spin_channels = 2
  i2 = st%nspin;         st%nspin         = 2
  rho_save => st%rho;    st%rho => rho

! calculate the u_sij using poissons equation
  rho(:, 2) = M_ZERO
  
  xcs2 = xcs
  xcs2%family = XC_FAMILY_LDA
  xcs2%functl = C_FUNC_LDA_PZ
  do i = 1, st%nst
    if(st%occ(i, is) .gt. small) then ! we only need the occupied states
      vc2 = M_ZERO
      ec2 = M_ZERO

      st%rho(:, 1) = socc*st%occ(i, is)*R_ABS(st%X(psi)(:, 1, i, is))**2

      call xc_lda (xcs2, m, st, vc2, edummy, ec2)

      ec = ec - sfact*ec2
      lc_(:, i) = lc_(:, i) - vc2(:, 1)*R_CONJ(st%X(psi) (:, 1, i, is))
    end if
  end do

  st%spin_channels = i1
  st%nspin         = i2    
  st%rho => rho_save
  deallocate(rho, Vc2)

end subroutine X(kli_c_sic)

!!$subroutine R_FUNC(kli_c_sic) (nlcc, m, st, psi, Vc, ec)
!!$  logical, intent(in) :: nlcc
!!$  type(mesh_type), intent(IN) :: m
!!$  type(states_type), intent(inout) :: st
!!$  R_TYPE, intent(IN) :: psi(m%np, st%nst, st%nspin)
!!$  real(r8), intent(out) :: Vc(m%np, st%nspin), ec
!!$
!!$  integer :: is, i, k, i1, i2
!!$  real(r8) :: socc, sfact, ec2
!!$
!!$  real(r8), allocatable :: rho(:,:), Vc2(:, :)
!!$  real(r8), allocatable :: u_xc(:,:), u_bar_xc(:)
!!$
!!$  allocate(u_xc(m%np, st%nst), u_bar_xc(st%nst))
!!$  allocate(rho(m%np, 2), Vc2(m%np, 2))
!!$
!!$  call getSpinFactor(st%nspin, socc, sfact)
!!$
!!$! first the LDA part
!!$  do is = 1, st%nspin
!!$     do k = 1, m%np
!!$        ! rho_sigma now contains the "total" spin density
!!$        rho(k, is) = sum(st%occ(:, is)*R_ABS(psi(k, :, is))**2)
!!$     end do
!!$  end do
!!$
!!$  call R_FUNC(xc_lda) (C_FUNC_LDA_PZ, nlcc, m, st, Vc, ec)
!!$
!!$! calculate the u_sij using poissons equation
!!$  rho(1:m%np, 2) = 0.0_r8
!!$  do is = 1, st%nspin
!!$     u_xc     = 0.0_r8
!!$     u_bar_xc = 0.0_r8
!!$
!!$     do i = 1, st%nst
!!$        if(st%occ(i, is) .gt. small) then ! we only need the occupied states
!!$           Vc2 = 0.0_r8; Ec2 = 0.0_r8;
!!$           rho(1:m%np, 1) = st%occ(i,is)*socc*R_ABS(psi(1:m%np, i, is))**2
!!$
!!$           i1 = st%ispin; st%ispin = 2
!!$           i2 = st%nspin; st%nspin = 2
!!$           call R_FUNC(xc_lda) (C_FUNC_LDA_PZ, nlcc, m, st, Vc2, ec2)
!!$           st%ispin = i1; st%nspin = i2;
!!$
!!$           u_xc(1:m%np, i) = - Vc2(1:m%np, 1)
!!$           u_bar_xc(i) = sum(u_xc(1:m%np, i) * R_ABS(psi(1:m%np, i, is))**2)*m%vol_pp
!!$           Ec = Ec - sfact*ec2
!!$        end if
!!$     end do
!!$
!!$     Vc2 = 0.0_r8
!!$     call R_FUNC(solve_KLI) (m, st%nspin, is, st%nst, socc, st%occ, st%eigenval, psi, &
!!$         u_xc, u_bar_xc, Vc2(:, 1:st%nspin))
!!$
!!$     Vc(1:m%np, is) = Vc(1:m%np, is) + Vc2(1:m%np, is)
!!$
!!$  end do !spin cycle
!!$
!!$  ! deallocate the rest of the vars..
!!$  deallocate(u_xc, u_bar_xc)
!!$  deallocate(rho, Vc2)
!!$
!!$  return
!!$end subroutine R_FUNC(kli_c_sic)

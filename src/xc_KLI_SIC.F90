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
subroutine X(kli_x_sic) (xcs, m, st, is, socc, sfact, lx_, ex)
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

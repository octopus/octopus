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

subroutine R_FUNC(kli_x_sic) (nlcc, m, st, psi, Vx, ex)
  logical, intent(in) :: nlcc
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(inout) :: st
  R_TYPE, intent(IN) :: psi(m%np, st%nst, st%nspin)
  real(r8), intent(out) :: Vx(m%np, st%nspin), ex

  integer :: is, i, k, i1, i2
  real(r8) :: socc, sfact, ex2

  real(r8), allocatable :: rho(:,:), Vx2(:, :)
  real(r8), allocatable :: u_xc(:,:), u_bar_xc(:)

  allocate(u_xc(m%np, st%nst), u_bar_xc(st%nst))
  allocate(rho(m%np, 2), Vx2(m%np, 2))

  call getSpinFactor(st%nspin, socc, sfact)

! first the LDA part
  do is = 1, st%nspin
     do k=1, m%np
        ! rho_sigma now contains the "total" spin density
        rho(k, is) = sum(st%occ(:, is)*R_ABS(psi(k, :, is))**2)
     end do
  end do

  call R_FUNC(xc_lda) (X_FUNC_LDA_NREL, nlcc, m, st, Vx, ex)

! calculate the u_sij using poissons equation
  rho(1:m%np, 2) = M_ZERO
  do is = 1, st%nspin
     u_xc     = M_ZERO
     u_bar_xc = M_ZERO

     do i = 1, st%nst
        if(st%occ(i, is) .gt. small) then ! we only need the occupied states
           Vx2 = 0.0_r8; Ex2 = 0.0_r8;
           rho(1:m%np, 1) = st%occ(i,is)*socc*R_ABS(psi(1:m%np, i, is))**2

           i1 = st%ispin; st%ispin = 2
           i2 = st%nspin; st%nspin = 2
           call R_FUNC(xc_lda) (X_FUNC_LDA_NREL, nlcc, m, st, Vx2, ex2)
           st%ispin = i1; st%nspin = i2;

           u_xc(1:m%np, i) = - Vx2(1:m%np, 1)
           Ex = Ex - sfact*ex2

           call poisson_solve(m, Vx2(:, 1), rho(:, 1))
           u_xc(1:m%np, i) = u_xc(1:m%np, i) - Vx2(1:m%np, 1)

           u_bar_xc(i) = sum(u_xc(1:m%np, i) * R_ABS(psi(1:m%np, i, is))**2)*m%vol_pp
           Ex = Ex - 0.5_r8*sfact*st%occ(i, is)*socc* &
               sum(Vx2(1:m%np, 1)* R_ABS(psi(1:m%np, i, is))**2)*m%vol_pp
        end if
     end do

     Vx2 = 0.0_r8
     call R_FUNC(solve_KLI) (m, st%nspin, is, st%nst, socc, st%occ, st%eigenval, psi, &
         u_xc, u_bar_xc, Vx2(:, 1:st%nspin))

     Vx(1:m%np, is) = Vx(1:m%np, is) + Vx2(1:m%np, is)

  end do !spin cycle

  ! deallocate the rest of the vars..
  deallocate(u_xc, u_bar_xc)
  deallocate(rho, Vx2)

  return
end subroutine R_FUNC(kli_x_sic)

subroutine R_FUNC(kli_c_sic) (nlcc, m, st, psi, Vc, ec)
  logical, intent(in) :: nlcc
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(inout) :: st
  R_TYPE, intent(IN) :: psi(m%np, st%nst, st%nspin)
  real(r8), intent(out) :: Vc(m%np, st%nspin), ec

  integer :: is, i, k, i1, i2
  real(r8) :: socc, sfact, ec2

  real(r8), allocatable :: rho(:,:), Vc2(:, :)
  real(r8), allocatable :: u_xc(:,:), u_bar_xc(:)

  allocate(u_xc(m%np, st%nst), u_bar_xc(st%nst))
  allocate(rho(m%np, 2), Vc2(m%np, 2))

  call getSpinFactor(st%nspin, socc, sfact)

! first the LDA part
  do is = 1, st%nspin
     do k = 1, m%np
        ! rho_sigma now contains the "total" spin density
        rho(k, is) = sum(st%occ(:, is)*R_ABS(psi(k, :, is))**2)
     end do
  end do

  call R_FUNC(xc_lda) (C_FUNC_LDA_PZ, nlcc, m, st, Vc, ec)

! calculate the u_sij using poissons equation
  rho(1:m%np, 2) = 0.0_r8
  do is = 1, st%nspin
     u_xc     = 0.0_r8
     u_bar_xc = 0.0_r8

     do i = 1, st%nst
        if(st%occ(i, is) .gt. small) then ! we only need the occupied states
           Vc2 = 0.0_r8; Ec2 = 0.0_r8;
           rho(1:m%np, 1) = st%occ(i,is)*socc*R_ABS(psi(1:m%np, i, is))**2

           i1 = st%ispin; st%ispin = 2
           i2 = st%nspin; st%nspin = 2
           call R_FUNC(xc_lda) (C_FUNC_LDA_PZ, nlcc, m, st, Vc2, ec2)
           st%ispin = i1; st%nspin = i2;

           u_xc(1:m%np, i) = - Vc2(1:m%np, 1)
           u_bar_xc(i) = sum(u_xc(1:m%np, i) * R_ABS(psi(1:m%np, i, is))**2)*m%vol_pp
           Ec = Ec - sfact*ec2
        end if
     end do

     Vc2 = 0.0_r8
     call R_FUNC(solve_KLI) (m, st%nspin, is, st%nst, socc, st%occ, st%eigenval, psi, &
         u_xc, u_bar_xc, Vc2(:, 1:st%nspin))

     Vc(1:m%np, is) = Vc(1:m%np, is) + Vc2(1:m%np, is)

  end do !spin cycle

  ! deallocate the rest of the vars..
  deallocate(u_xc, u_bar_xc)
  deallocate(rho, Vc2)

  return
end subroutine R_FUNC(kli_c_sic)

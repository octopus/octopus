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

#ifdef HAVE_LAPACK
subroutine R_FUNC(xc_kli)(func, nlcc, m, st, pot, energy)
  integer, intent(in) :: func
  logical, intent(in) :: nlcc
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(inout) :: st
  real(r8), intent(out) :: pot(m%np, st%nspin), energy
  
  call push_sub('xc_kli')

  ! this routine is only prepared for finite systems, and ispin = 1, 2
  if(st%ispin > 2 .or. st%nik>st%ispin) then
    message(1) = "KLI only works for finite systems and collinear spin!"
    call write_fatal(1)
  end if

  select case(func)
    case(X_FUNC_KLI_X)
      call R_FUNC(kli_x) (m, st%nspin, st%nst, st%occ, st%eigenval, &
           st%R_FUNC(psi) (:,1,:,:), pot, energy, st%rho)

    case(X_FUNC_KLI_SIC)
      call R_FUNC(kli_x_sic) (nlcc, m, st, st%R_FUNC(psi) (:,1,:,:), pot, energy)
    case(C_FUNC_KLI_SIC)
      call R_FUNC(kli_c_sic) (nlcc, m, st, st%R_FUNC(psi) (:,1,:,:), pot, energy)

    case(X_FUNC_KLI_HJU)
      call dkli_hju(m, st, nlcc, X_FUNC_LDA_NREL, pot, energy)
    case(C_FUNC_KLI_HJU)
      call dkli_hju(m, st, nlcc, C_FUNC_LDA_PZ,   pot, energy)
  end select

  call pop_sub()
  return
end subroutine R_FUNC(xc_kli)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates the potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine R_FUNC(solve_KLI) (m, nspin, is, nst, socc, occ, eigenval, &
    psi, u_xc, u_bar_xc, v_xc)
  
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: nspin, is, nst
  real(r8), intent(in) :: socc
  real(r8), intent(IN) :: occ(nst, nspin), eigenval(nst, nspin), &
      u_xc(m%np, nst), u_bar_xc(nst)
  R_TYPE, intent(IN) :: psi(m%np, nst, nspin)
  real(r8), intent(out) :: v_xc(m%np, nspin)

  integer i, j, k, eigen_n
  real(r8), allocatable :: Ma(:,:), x(:), y(:)
  integer,  allocatable :: eigen_type(:), eigen_index(:)
  real(r8) :: max_eigen, rho_sigma(m%np), v_bar_S(nst)

  ! variables needed by LAPACK
  character(len=1) :: la_EQUED
  real(r8), allocatable :: la_AF(:,:), la_R(:), la_C(:), la_work(:)
  real(r8) :: la_R_cond, la_Ferr(1), la_Berr(1)
  integer, allocatable :: la_IPIV(:), la_iwork(:)
  integer :: info

  ! some intermediate quantities
  ! v_xc contains the Slater part!
  do k=1, m%np
    rho_sigma(k) = sum(occ(:, is)*R_ABS(psi(k, :, is))**2)*socc
    v_xc(k, is)  = sum(u_xc(k, :) * occ(:, is)*R_ABS(psi(k, :, is))**2)*socc &
        /(rho_sigma(k) + denom_eps)
  end do

  do i=1, NST
    if(occ(i, is) .gt. small) then
      v_bar_S(i) = sum(R_ABS(psi(1:m%np, i, is))**2 * v_xc(1:m%np, is))*m%vol_pp
    end if
  end do

  ! find out the top occupied state, to correct for the assymptotics
  ! of the potential
  allocate(eigen_type(NST), eigen_index(NST))
  max_eigen = -1e30_r8
  do i=1, nst
    if((occ(i, is) .gt. small).and.(eigenval(i,is).gt.max_eigen)) then
      max_eigen = eigenval(i, is)
    end if
  end do
  
  eigen_n = 1
  do i=1, nst
    if(occ(i, is) .gt. small) then
      ! criterium for degeneracy
      if(abs(eigenval(i,is)-max_eigen).le.1e-3_r8) then
        eigen_type(i) = 2
      else
        eigen_type(i) = 1
        eigen_index(eigen_n) = i
        eigen_n = eigen_n +1
      end if
    else
      eigen_type(i) = 0
    end if
  end do
  eigen_n = eigen_n - 1

  if(eigen_n > 0) then ! there is more than one state, so solve linear equation
    allocate(x(eigen_n))
    x = 0.0_r8
    allocate(Ma(eigen_n, eigen_n), y(eigen_n))
    do i=1,eigen_n
      do j=i,eigen_n
        Ma(i,j) = -sum(                                 &
            R_ABS(psi(1:m%np, eigen_index(i), is))**2 * &
            R_ABS(psi(1:m%np, eigen_index(j), is))**2 / &
            (rho_sigma(1:m%np) + denom_eps))*m%vol_pp
        Ma(j,i) = Ma(i,j)
      end do
      Ma(i,i) = 1 + Ma(i,i)
      y(i) = v_bar_S(eigen_index(i)) - u_bar_xc(eigen_index(i))
    end do

    ! setup lapack arrays
    allocate(la_AF(eigen_n,eigen_n), la_IPIV(eigen_n), la_R(eigen_n), &
        la_C(eigen_n), la_work(4*eigen_n), la_iwork(eigen_n))
    
    call DGESVX('N', 'N', eigen_n, 1, Ma, eigen_n, la_AF, eigen_n,     &
        la_IPIV, la_EQUED, la_R, la_C, y, eigen_n, x, eigen_n,       &
        la_R_cond, la_Ferr, la_Berr, la_work, la_iwork, info)
    if(info.ne.0) then
      write(6,'(a,I5,a)') 'KLI:: error in lapack (info=',info,')'
    end if

    deallocate(la_AF, la_IPIV, la_R, la_C, la_work, la_iwork)
    deallocate(Ma, y)

    ! add contribution of low lying states
    do i=1, eigen_n
      v_xc(1:m%np, is) = v_xc(1:m%np, is) + &
          occ(eigen_index(i),is)*socc* x(i) *  &
          (R_ABS(psi(1:m%np, eigen_index(i), is))**2 / (rho_sigma(1:m%np) + denom_eps))
    end do
    deallocate(x)
  end if
  deallocate(eigen_type, eigen_index)

  return
end subroutine R_FUNC(solve_KLI)

#endif

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
subroutine X(xc_kli)(xcs, m, st, vxc, ex, ec)
  type(xc_type),     intent(in)    :: xcs
  type(mesh_type),   intent(in)    :: m
  type(states_type), intent(inout) :: st
  real(r8),          intent(out)   :: vxc(m%np, st%nspin), ex, ec
  
  R_TYPE,   allocatable :: lxc(:,:)
  real(r8) :: socc, sfact, e
  integer :: is, ixc, ifunc

  call push_sub('xc_kli')

  ! this routine is only prepared for finite systems, and ispin = 1, 2
  if(st%ispin > 2 .or. st%nik>st%ispin) then
    message(1) = "KLI only works for finite systems and collinear spin!"
    call write_fatal(1)
  end if

  call getSpinFactor(st%nspin, socc, sfact)
  allocate(lxc(m%np, st%nst))

  spin: do is = 1, min(st%nspin, 2)
    lxc     = M_ZERO
    e       = M_ZERO

    ! get lxc
    functl_loop: do ixc = 0, N_X_FUNCTL+N_C_FUNCTL-1
      if(.not.btest(xcs%functl, ixc)) cycle
      ifunc = ibset(0, ixc)
      if(.not.( &
           ifunc == X_FUNC_KLI_X    .or. &
           ifunc == X_FUNC_KLI_SIC  .or. &
           ifunc == C_FUNC_KLI_SIC)) cycle
      
      select case(ifunc)
      case(X_FUNC_KLI_X)
        call X(kli_x) (m, st, is, socc, sfact, lxc, e)
        
      case(X_FUNC_KLI_SIC)
        call X(kli_x_sic) (xcs, m, st, is, socc, sfact, lxc, e)
      case(C_FUNC_KLI_SIC)
        call X(kli_c_sic) (xcs, m, st, is, socc, sfact, lxc, e)
      end select
      
      if(ixc < N_X_FUNCTL) then
        ex = ex + e
      else
        ec = ec + e
      end if
    end do functl_loop
  
    ! solve the KLI equation
    call X(solve_KLI) (m, st, is, socc, lxc, vxc)

  end do spin

  deallocate(lxc)

  call pop_sub()
end subroutine X(xc_kli)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates the potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine X(solve_KLI) (m, st, is, socc, lxc, vxc_)
  type(mesh_type),   intent(in)  :: m
  type(states_type), intent(in)  :: st
  integer,           intent(in)  :: is
  real(r8),          intent(in)  :: socc
  R_TYPE,            intent(in)  :: lxc(m%np, st%nst)
  real(r8),          intent(out) :: vxc_(m%np, st%nspin)

  integer i, j, k, eigen_n
  real(r8), allocatable :: rho_sigma(:), vxc(:), v_bar_S(:)
  real(r8), allocatable :: Ma(:,:), x(:), y(:)
  integer,  allocatable :: eigen_type(:), eigen_index(:)
  real(r8) :: max_eigen, uxc_bar

  ! variables needed by LAPACK
  character(len=1) :: la_EQUED
  real(r8), allocatable :: la_AF(:,:), la_R(:), la_C(:), la_work(:)
  real(r8) :: la_R_cond, la_Ferr(1), la_Berr(1)
  integer, allocatable :: la_IPIV(:), la_iwork(:)
  integer :: info

  ! some intermediate quantities
  ! vxc contains the Slater part!
  allocate(vxc(m%np), rho_sigma(m%np), v_bar_S(st%nst))
  do k = 1, m%np
    rho_sigma(k) = max(sum(socc*st%occ(:, is)*R_ABS(st%X(psi)(k, 1, :, is))**2), 1e-20_r8)
    vxc(k)   = socc*sum(st%occ(:, is)*R_REAL(lxc(k, :)*st%X(psi)(k, 1, :, is))) &
         /rho_sigma(k)
  end do

  do i = 1, st%nst
    if(st%occ(i, is) .gt. small) then
      v_bar_S(i) = sum(R_ABS(st%X(psi)(:, 1, i, is))**2 * vxc(:))*m%vol_pp
    end if
  end do

  ! find out the top occupied state, to correct for the assymptotics
  ! of the potential
  allocate(eigen_type(st%nst), eigen_index(st%nst))
  max_eigen = -1e30_r8
  do i = 1, st%nst
    if((st%occ(i, is) .gt. small).and.(st%eigenval(i, is).gt.max_eigen)) then
      max_eigen = st%eigenval(i, is)
    end if
  end do
  
  eigen_n = 1
  do i = 1, st%nst
    if(st%occ(i, is) .gt. small) then
      ! criterium for degeneracy
      if(abs(st%eigenval(i,is)-max_eigen).le.1e-3_r8) then
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
            R_ABS(st%X(psi)(:, 1, eigen_index(i), is))**2 * &
            R_ABS(st%X(psi)(:, 1, eigen_index(j), is))**2 / &
            rho_sigma(:))*m%vol_pp
        Ma(j,i) = Ma(i,j)
      end do
      Ma(i,i) = 1 + Ma(i,i)

      ! y(i) = v_bar_S - uxc_bar
      uxc_bar = sum(R_REAL(st%X(psi)(:, 1, i, is) * lxc(:, i)))*m%vol_pp
      y(i) = v_bar_S(eigen_index(i)) - uxc_bar
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
    do i = 1, eigen_n
      vxc(:) = vxc(:) + &
           socc*st%occ(eigen_index(i),is)* x(i) *  &
           R_ABS(st%X(psi)(:, 1, eigen_index(i), is))**2 / rho_sigma(:)
    end do
    deallocate(x)
  end if

  vxc_(:,is) = vxc_(:,is) + vxc(:)

  deallocate(eigen_type, eigen_index)
  deallocate(vxc, rho_sigma, v_bar_S)

end subroutine X(solve_KLI)

#endif

subroutine R_FUNC(kli_x_sic) (m, nspin, nst, occ, eigenval, psi, rho_core, hartr, Vx, ex)

  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: nspin, nst
  real(r8), intent(IN) :: occ(nst, nspin), eigenval(nst, nspin), &
      rho_core(m%np)
  R_TYPE, intent(IN) :: psi(0:m%np, nst, nspin)
  type(hartree_type), intent(inout) :: hartr
  real(r8), intent(out) :: Vx(m%np, nspin), ex

  integer :: is, i, k
  real(r8) :: socc, sfact, ex2

  real(r8), allocatable :: rho(:,:), Vx2(:, :)
  real(r8), allocatable :: u_xc(:,:), u_bar_xc(:)

  allocate(u_xc(m%np, nst), u_bar_xc(nst))
  allocate(rho(m%np, 2), Vx2(m%np, 2))

  call getSpinFactor(nspin, socc, sfact)

! first the LDA part
  do is = 1, min(nspin, 2)
     do k=1, m%np
        ! rho_sigma now contains the "total" spin density
        rho(k, is) = sum(occ(:, is)*R_ABS(psi(k, :, is))**2)
     end do
  end do

  call xc_lda(X_FUNC_LDA_NREL, m, nspin, nspin, rho(:, 1:nspin), &
      rho_core, Vx, ex)

! calculate the u_sij using poissons equation
  rho(1:m%np, 2) = 0.0_r8
  do is = 1, min(nspin, 2)
     u_xc     = 0.0_r8
     u_bar_xc = 0.0_r8

     do i = 1, NST
        if(occ(i, is) .gt. small) then ! we only need the occupied states
           Vx2 = 0.0_r8; Ex2 = 0.0_r8;
           rho(1:m%np, 1) = occ(i,is)*socc*R_ABS(psi(1:m%np, i, is))**2

           call xc_lda(X_FUNC_LDA_NREL, m, 2, 2, rho, rho_core, Vx2, ex2)
           u_xc(1:m%np, i) = - Vx2(1:m%np, 1)
           Ex = Ex - sfact*ex2

           call hartree_solve(hartr, m, Vx2(:, 1), rho(:, 1:1))
           u_xc(1:m%np, i) = u_xc(1:m%np, i) - Vx2(1:m%np, 1)

           u_bar_xc(i) = sum(u_xc(1:m%np, i) * R_ABS(psi(1:m%np, i, is))**2)*m%vol_pp
           Ex = Ex - 0.5_r8*sfact*occ(i, is)*socc* &
               sum(Vx2(1:m%np, 1)* R_ABS(psi(1:m%np, i, is))**2)*m%vol_pp
        end if
     end do

     Vx2 = 0.0_r8
     call R_FUNC(solve_KLI) (m, nspin, is, nst, socc, occ, eigenval, psi, &
         u_xc, u_bar_xc, Vx2(:, 1:nspin))

     Vx(1:m%np, is) = Vx(1:m%np, is) + Vx2(1:m%np, is)

  end do !spin cycle

  ! deallocate the rest of the vars..
  deallocate(u_xc, u_bar_xc)
  deallocate(rho, Vx2)

  return
end subroutine R_FUNC(kli_x_sic)

subroutine R_FUNC(kli_c_sic) (m, nspin, nst, occ, eigenval, psi, rho_core, hartr, Vc, ec)

  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: nspin, nst
  real(r8), intent(IN) :: occ(nst, nspin), eigenval(nst, nspin), &
      rho_core(m%np)
  R_TYPE, intent(IN) :: psi(0:m%np, nst, nspin)
  type(hartree_type), intent(IN) :: hartr
  real(r8), intent(out) :: Vc(m%np, nspin), ec

  integer :: is, i, k
  real(r8) :: socc, sfact, ec2

  real(r8), allocatable :: rho(:,:), Vc2(:, :)
  real(r8), allocatable :: u_xc(:,:), u_bar_xc(:)

  allocate(u_xc(m%np, nst), u_bar_xc(nst))
  allocate(rho(m%np, 2), Vc2(m%np, 2))

  call getSpinFactor(nspin, socc, sfact)

! first the LDA part
  do is = 1, min(nspin, 2)
     do k=1, m%np
        ! rho_sigma now contains the "total" spin density
        rho(k, is) = sum(occ(:, is)*R_ABS(psi(k, :, is))**2)
     end do
  end do

  call xc_lda(C_FUNC_LDA_PZ, m, nspin, nspin, rho(:, 1:nspin), &
      rho_core, Vc, ec)

! calculate the u_sij using poissons equation
  rho(1:m%np, 2) = 0.0_r8
  do is = 1, min(nspin, 2)
     u_xc     = 0.0_r8
     u_bar_xc = 0.0_r8

     do i = 1, NST
        if(occ(i, is) .gt. small) then ! we only need the occupied states
           Vc2 = 0.0_r8; Ec2 = 0.0_r8;
           rho(1:m%np, 1) = occ(i,is)*socc*R_ABS(psi(1:m%np, i, is))**2

           call xc_lda(C_FUNC_LDA_PZ, m, 2, 2, rho(:, 1:2), rho_core, Vc2, ec2)
           u_xc(1:m%np, i) = - Vc2(1:m%np, 1)
           u_bar_xc(i) = sum(u_xc(1:m%np, i) * R_ABS(psi(1:m%np, i, is))**2)*m%vol_pp
           Ec = Ec - sfact*ec2
        end if
     end do

     Vc2 = 0.0_r8
     call R_FUNC(solve_KLI) (m, nspin, is, nst, socc, occ, eigenval, psi, &
         u_xc, u_bar_xc, Vc2(:, 1:nspin))

     Vc(1:m%np, is) = Vc(1:m%np, is) + Vc2(1:m%np, is)

  end do !spin cycle

  ! deallocate the rest of the vars..
  deallocate(u_xc, u_bar_xc)
  deallocate(rho, Vc2)

  return
end subroutine R_FUNC(kli_c_sic)

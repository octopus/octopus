subroutine R_FUNC(kli_x) (m, nspin, nst, occ, eigenval, psi, hartr, Vx, ex, rho)

  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: nspin, nst
  real(r8), intent(IN) :: occ(nst, nspin), eigenval(nst, nspin)
  R_TYPE, intent(IN) :: psi(0:m%np, nst, nspin)
  type(hartree_type), intent(inout) :: hartr
  real(r8), intent(out) :: Vx(m%np, nspin), ex
  real(r8), intent(in)  :: rho(m%np, nspin)

  integer :: is, i, j, k
  real(r8) :: charge, socc, sfact
  real(r8), allocatable :: rho_ij(:,:)
  real(r8), allocatable :: pot(:)             ! For real
  real(r8), allocatable :: pot_r(:), pot_i(:) ! For complex
  real(r8), allocatable :: u_xc(:,:), u_bar_xc(:)
  
  allocate(u_xc(m%np, nst), u_bar_xc(nst))
  call getSpinFactor(nspin, socc, sfact)

  allocate(pot(m%np))
  if(nst==1 .and. nspin==1) then
    call hartree_solve(hartr, m, pot, rho)
    if(occ(1, 1) == 1.0_r8) then
      vx(:, 1) = -pot(:); return
    elseif(occ(1, 1) == 2.0_r8) then
      vx(:, 1) = -0.5_r8*pot(:); Ex = 0.5_r8*dmesh_dotp(m, vx(:, 1), rho(:, 1)); return
    else
      write(*,*) 'Warning in kli_x', occ(1, 1)
    endif
  endif

  ! calculate the u_sij using poissons equation
  Ex = 0.0_r8
  do is = 1, min(nspin, 2)
    u_xc = 0.0_r8
    
    do i=1, nst
      if(occ(i, is) .gt. small) then ! we only need the occupied states
        do j=1, nst
          if(occ(j, is) .gt. small) then
            allocate(rho_ij(1:m%np,1:1))

#ifdef R_TREAL
            !allocate(pot(m%np))
            pot = 0._r8

            rho_ij(1:m%np, 1) = psi(1:m%np, i, is)*psi(1:m%np, j, is)
            call hartree_solve(hartr, m, pot, rho_ij)
            deallocate(rho_ij)
            do k=1, m%np
              u_xc(k, i) = u_xc(k, i) - occ(j,is)*socc*pot(k)*    &
                  (psi(k,j,is)/(psi(k,i,is) + my_sign(psi(k,i,is))*denom_eps))

            end do
            deallocate(pot)
#else
            allocate(pot_r(m%np), pot_i(m%np))
            pot_r = 0._r8
            pot_i = 0._r8

            rho_ij(1:m%np, 1) = real(psi(1:m%np, i, is))*real(psi(1:m%np, j, is)) + &
                aimag(psi(1:m%np, i, is))*aimag(psi(1:m%np, j, is))
            call hartree_solve(hartr, m, pot_r, rho_ij(1:m%np, 1:1))
            ! and now the imaginary part
            rho_ij(1:m%np, 1) = real(psi(1:m%np, i, is))*aimag(psi(1:m%np, j, is)) - &
                aimag(psi(1:m%np, i, is))*real(psi(1:m%np, j, is))
            call hartree_solve(hartr, m, pot_i, rho_ij(1:m%np, 1:1))
            deallocate(rho_ij)
            do k=1, m%np
              u_xc(k, i) = u_xc(k, i) - occ(j,is)*socc*real((pot_r(k) + M_zI*pot_i(k))*  &
                  (conjg(psi(k,j,is))/ &
                  (conjg(psi(k,i,is)) + my_sign(real(psi(k,i,is)))* denom_eps)))
            end do
            deallocate(pot_r, pot_i)
#endif
            
          end if
        end do
        u_bar_xc(i) = sum(R_ABS(psi(1:m%np, i, is))**2 * u_xc(1:m%np, i))*m%vol_pp
        Ex = Ex + 0.5_r8 * occ(i, is)*socc*u_bar_xc(i)        
      end if
    end do
    
    call R_FUNC(solve_KLI) (m, nspin, is, nst, socc, occ, eigenval, psi, &
        u_xc, u_bar_xc, vx)
    !vx = u_xc

  end do !spin cycle
  
  ! adjust energy
  Ex = Ex * sfact
  
  ! deallocate the rest of the vars..
  deallocate(u_xc, u_bar_xc)

  return
end subroutine R_FUNC(kli_x)

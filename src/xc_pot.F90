subroutine R_FUNC(xc_pot) (xcs, m, st, hart, rho_core, vx, vc, ex, ec)
  type(xc_type), intent(inout) :: xcs
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  type(hartree_type), intent(inout) :: hart
  real(r8), intent(IN) :: rho_core(m%np)
  real(r8), intent(out) :: vx(m%np, st%ispin), vc(m%np, st%ispin), ex, ec

  ! for fxc != vxc...
  ! fxc is always LDA!!!!
!!$  R_TYPE, allocatable, save :: save_vxc(:,:)
!!$  logical, save :: first_time = .true.

  sub_name = 'xc_pot'; call push_sub()

  vx = 0.0_r8; vc = 0.0_r8
  ex = 0.0_r8; ec = 0.0_r8

  select case(xcs%x_family)
  case(XC_FAMILY_ZER)
  case(XC_FAMILY_LDA)
    call xc_lda(xcs%x_func, m, st%ispin, st%rho, rho_core, vx, ex)
  case(XC_FAMILY_GGA)
!    call xc_gga(xcs%x_func, xcs, m, st%ispin, st%rho, rho_core, vx, ex)
#ifdef HAVE_LAPACK
!  case(XC_FAMILY_MGGA)
!    call R_FUNC(xc_mgga) (xcs%x_func, xcs, m, nst, st%ispin, psi, occ, eigenval, &
!        rho, rho_core, vx, ex)
!  case(XC_FAMILY_KLI)
!    call R_FUNC(xc_kli) (xcs%x_func, m, st%ispin, nst, occ, eigenval, &
!        psi, rho_core, hartr, vx, ex)
#endif
  end select

  select case(xcs%c_family)
  case(XC_FAMILY_ZER)
  case(XC_FAMILY_LDA)
    call xc_lda(xcs%c_func, m, st%ispin, st%rho, rho_core, vc, ec)
  case(XC_FAMILY_GGA)
!    call xc_gga(xcs%c_func, xcs, m, st%ispin, st%rho, rho_core, vc, ec)
#ifdef HAVE_LAPACK
!  case(XC_FAMILY_MGGA)
!    call R_FUNC(xc_mgga) (xcs%c_func, xcs, m, nst, st%ispin, psi, occ, eigenval, &
!        rho, rho_core, vc, ec)
!  case(XC_FAMILY_KLI)
!    call R_FUNC(xc_kli) (xcs%c_func, m, st%ispin, nst, occ, eigenval, &
!        psi, rho_core, hartr, vc, ec)
#endif
  end select

  ! Warning: For vxc != vxc
!!$  if(first_time) then
!!$    first_time = .false.
!!$    allocate(save_vxc(m%np, st%ispin))
!!$    save_vxc = vx + vc
!!$
!!$    ! now, we get the LDA xc
!!$    xcs%x_family = XC_FAMILY_LDA
!!$    xcs%x_func   = X_FUNC_LDA_NREL
!!$    call xc_lda(xcs%x_func, m, st%ispin, rho, rho_core, vx, ex)
!!$
!!$    xcs%c_family = XC_FAMILY_LDA
!!$    xcs%c_func   = C_FUNC_LDA_PZ
!!$    call xc_lda(xcs%c_func, m, st%ispin, rho, rho_core, vc, ec)
!!$
!!$    save_vxc = save_vxc - vx - vc
!!$  end if
!!$  vx = vx + save_vxc

  call pop_sub()
  return
end subroutine R_FUNC(xc_pot)

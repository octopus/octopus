!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!!
!! $Id$

subroutine dxc_get_vxc(der, xcs, st, rho, ispin, ioniz_pot, qtot, ex, ec, deltaxc, vxc, vtau)
  type(derivatives_t),  intent(in)    :: der             !< Discretization and the derivative operators and details
  type(xc_t), target,   intent(in)    :: xcs             !< Details about the xc functional used
  type(states_t),       intent(in)    :: st              !< State of the system (wavefunction,eigenvalues...)
  FLOAT,                intent(in)    :: rho(:, :)       !< Electronic density 
  integer,              intent(in)    :: ispin           !< Number of spin channels 
  FLOAT,                intent(in)    :: ioniz_pot
  FLOAT,                intent(in)    :: qtot 
  FLOAT, optional,      intent(inout) :: ex              !< Exchange energy.
  FLOAT, optional,      intent(inout) :: ec              !< Correlation energy.
  FLOAT, optional,      intent(inout) :: deltaxc         !< The XC derivative descontinuity
  FLOAT, optional,      intent(inout) :: vxc(:,:)        !< XC potential
  FLOAT, optional,      intent(inout) :: vtau(:,:)       !< Derivative wrt (two times kinetic energy density)

  integer :: n_block

  FLOAT, allocatable :: l_zk(:)        ! Local block of the energy functional (with the correct memory order for libxc)
  FLOAT, allocatable :: l_dens(:,:)    ! Local block for the density 
  FLOAT, allocatable :: l_dedd(:,:)    ! Local block of the exchange or correl. potential(with the correct memory order for libxc)
  FLOAT, allocatable :: l_sigma(:,:)   
  FLOAT, allocatable :: l_vsigma(:,:)  
  FLOAT, allocatable :: l_tau(:,:)
  FLOAT, allocatable :: l_ldens(:,:)
  FLOAT, allocatable :: l_dedtau(:,:)
  FLOAT, allocatable :: l_dedldens(:,:)

  FLOAT, allocatable :: dens(:,:)      ! Density
  FLOAT, allocatable :: dedd(:,:)      ! (Functional) Derivative of the exchange or correlation energy with
  ! respect to the density (vector used to store the exchange or the correlation potential)
  FLOAT, allocatable :: ex_per_vol(:)  ! Exchange energy per unit volume 
  FLOAT, allocatable :: ec_per_vol(:)  ! Correlation energy per unit volume 
  FLOAT, allocatable :: gdens(:,:,:)   ! Gradient of the density
  FLOAT, allocatable :: dedgd(:,:,:)   ! (Functional) Derivative of the exchange or correlation energy with
  !respect to the gradient of the density.
  FLOAT, allocatable :: current(:,:,:) ! Paramagnetic or total current
  FLOAT, allocatable :: ldens(:,:)     ! Laplacian of the density
  FLOAT, allocatable :: tau(:,:)       ! Kinetic energy density
  FLOAT, allocatable :: dedldens(:,:)  ! (Functional) Derivative of the exchange or correlation energy with
  !respect to the laplacian of the density.
  FLOAT, allocatable :: symmtmp(:, :)  ! Temporary vector for the symmetrizer
  FLOAT, allocatable :: vx(:)
  FLOAT, allocatable :: unp_dens(:), unp_dedd(:)

  integer :: ib, ib2, ip, isp, families, ixc, spin_channels, is
  integer, save :: xc_get_vxc_counter = 0
  FLOAT   :: rr,ipot_to_pass
  logical :: gga, mgga
  type(profile_t), save :: prof
  logical :: calc_energy
  type(xc_functl_t), pointer :: functl(:)
  type(symmetrizer_t) :: symmetrizer

  PUSH_SUB(dxc_get_vxc)
  call profiling_in(prof, "dXC_LOCAL")

  ASSERT(present(ex) .eqv. present(ec))
  calc_energy = present(ex)

  xc_get_vxc_counter = xc_get_vxc_counter + 1
  !xprint*, "xc_get_vxc call number ", xc_get_vxc_counter

  !Pointer-shortcut for xcs%functl
  !It helps to remember that for xcs%functl(:,:)
  ! (1,:) => exchange,    (2,:) => correlation
  ! (:,1) => unpolarized, (:,2) => polarized
  if(ispin == UNPOLARIZED) then
    functl => xcs%functl(:, 1)
  else
    functl => xcs%functl(:, 2)
  end if

  ! is there anything to do ?
  families = XC_FAMILY_LDA + XC_FAMILY_GGA + XC_FAMILY_HYB_GGA + XC_FAMILY_MGGA
  if(iand(xcs%family, families) == 0) then
    POP_SUB(dxc_get_vxc)
    call profiling_out(prof)
    return
  endif

  n_block = 1000

  ! initialize a couple of handy variables
  gga  = iand(xcs%family, XC_FAMILY_GGA + XC_FAMILY_HYB_GGA + XC_FAMILY_MGGA).ne.0
  mgga = iand(xcs%family, XC_FAMILY_MGGA).ne.0

  !Read the spin channels
  spin_channels = functl(FUNC_X)%spin_channels
  
  if(xcs%xc_density_correction == LR_X) then
    SAFE_ALLOCATE(vx(1:der%mesh%np))
  end if

  call lda_init()
  if(gga .or. xcs%xc_density_correction == LR_X) call  gga_init()
  if(mgga) call mgga_init()

  ! Get the gradient and the Laplacian of the density and the kinetic-energy density
  ! We do it here instead of doing it in gga_init and mgga_init in order to 
  ! avoid calling the subroutine states_calc_quantities twice
  if((gga .and. (.not. mgga)) .or. xcs%xc_density_correction == LR_X) then
    ! get gradient of the density (this is faster than calling states_calc_quantities)
    do isp = 1, spin_channels 
      call dderivatives_grad(der, dens(:, isp), gdens(:, :, isp)) 
    end do
  else if(mgga) then
    if (xc_get_vxc_counter .le. 2) then 
      ipot_to_pass = M_ONE
    else
      ipot_to_pass = ioniz_pot
    end if
    !print *, "Ioniz potential to pass =", ipot_to_pass
    ! We calculate everything from the wavefunctions to benefit from
    ! the error cancellation between the gradient of the density and
    ! tau.
    !
    ! FIXME: Probably this is wrong for non-local corrections or other
    ! cases when the density is not directly generated by the
    ! orbitals.
    call states_calc_quantities(der, st, gi_kinetic_energy_density = tau, density_gradient = gdens, density_laplacian = ldens)

    ! We have to symmetrize everything as they are calculated from the
    ! wavefunctions.
    if(st%symmetrize_density) then
      SAFE_ALLOCATE(symmtmp(1:der%mesh%np, 1:3))
      call symmetrizer_init(symmetrizer, der%mesh)
      do isp = 1, spin_channels
        call dsymmetrizer_apply(symmetrizer, tau(:, isp), symmtmp(:, 1))
        tau(1:der%mesh%np, isp) = symmtmp(1:der%mesh%np, 1)
        call dsymmetrizer_apply(symmetrizer, ldens(:, isp), symmtmp(:, 1))
        ldens(1:der%mesh%np, isp) = symmtmp(1:der%mesh%np, 1)
        call dsymmetrizer_apply_vector(symmetrizer, gdens(:, :, isp), symmtmp)
        gdens(1:der%mesh%np, 1:3, isp) = symmtmp(1:der%mesh%np, 1:3)
      end do

      call symmetrizer_end(symmetrizer)
      SAFE_DEALLOCATE_A(symmtmp)
    end if

    if(functl(FUNC_X)%id == XC_MGGA_X_TB09 .and. der%mesh%sb%periodic_dim == 3) then
      call calc_tb09_c()
    end if

  end if

  space_loop: do ip = 1, der%mesh%np, n_block

    !Resize the dimension of the last block when the number of the mesh points
    !it is not a perfect divider of the dimension of the blocks.
    if(ip + n_block > der%mesh%np) n_block = der%mesh%np - ip + 1 

    ! make a local copy with the correct memory order for libxc
    ib2 = ip
    do ib = 1, n_block
      l_dens(1:spin_channels, ib) = dens(ib2, 1:spin_channels)
      ib2 = ib2 + 1
    end do

    if(gga) then
      ib2 = ip
      do ib = 1, n_block
        l_sigma(1, ib) = sum(gdens(ib2, 1:der%mesh%sb%dim, 1)*gdens(ib2, 1:der%mesh%sb%dim, 1))
        if(ispin /= UNPOLARIZED) then
          ! memo: please check the following indices
          l_sigma(2, ib) = sum(gdens(ib2, 1:der%mesh%sb%dim, 1)*gdens(ib2, 1:der%mesh%sb%dim, 2)) 
          l_sigma(3, ib) = sum(gdens(ib2, 1:der%mesh%sb%dim, 2)*gdens(ib2, 1:der%mesh%sb%dim, 2))
        end if
        ib2 = ib2 + 1
      end do
    end if

    if(mgga) then
      ib2 = ip
      do ib = 1, n_block
        l_tau  (1:spin_channels, ib) =   tau(ib2, 1:spin_channels)
        l_ldens(1:spin_channels, ib) = ldens(ib2, 1:spin_channels)
        ib2 = ib2 + 1
      end do
    end if

    ! Calculate the potential/gradient density in local reference frame.
    functl_loop: do ixc = FUNC_X, FUNC_C

      if(.not. present(vxc)) then ! get only the xc energy

        if(iand(functl(ixc)%flags, XC_FLAGS_HAVE_EXC).ne.0) then
          select case(functl(ixc)%family)

          case(XC_FAMILY_LDA)
            call XC_F90(lda_exc)(functl(ixc)%conf, n_block, l_dens(1,1), l_zk(1))

          case(XC_FAMILY_GGA)
            call XC_F90(gga_exc)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), l_zk(1))

          case(XC_FAMILY_HYB_GGA)
            message(1) = 'Hyb-GGAs are currently disabled.'
            call messages_fatal(1)

          case(XC_FAMILY_MGGA)
            call XC_F90(mgga_exc)(functl(ixc)%conf, n_block, &
              l_dens(1,1), l_sigma(1,1), l_ldens(1,1), l_tau(1,1), l_zk(1))

          case default
            cycle
          end select

        else ! Do not have an energy functional
          l_zk(:) = M_ZERO
        end if

      else if(calc_energy .and. iand(functl(ixc)%flags, XC_FLAGS_HAVE_EXC).ne.0) then
        ! we get the xc energy and potential
        select case(functl(ixc)%family)
        case(XC_FAMILY_LDA)
          call XC_F90(lda_exc_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_zk(1), l_dedd(1,1))

        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          call XC_F90(gga_exc_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), &
            l_zk(1), l_dedd(1,1), l_vsigma(1,1))

        case(XC_FAMILY_MGGA)
          call XC_F90(mgga_exc_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), l_ldens(1,1), l_tau(1,1), &
            l_zk(1), l_dedd(1,1), l_vsigma(1,1), l_dedldens(1,1), l_dedtau(1,1))

        case default
          cycle
        end select

      else ! we just get the potential
        l_zk(:) = M_ZERO

        select case(functl(ixc)%family)
        case(XC_FAMILY_LDA)
          call XC_F90(lda_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_dedd(1,1))

        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          l_vsigma = M_ZERO

          if(functl(ixc)%id == XC_GGA_X_LB) then
            call mesh_r(der%mesh, ip, rr)
            call XC_F90(gga_lb_modified)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), &
              rr, l_dedd(1,1))
          else
            call XC_F90(gga_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), &
              l_dedd(1,1), l_vsigma(1,1))
          end if

        case(XC_FAMILY_MGGA)
          call XC_F90(mgga_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), l_ldens(1,1), l_tau(1,1), &
            l_dedd(1,1), l_vsigma(1,1), l_dedldens(1,1), l_dedtau(1,1))

        case default
          cycle
        end select

      end if

      if(calc_energy) then
        ib2 = ip
        if(functl(ixc)%type == XC_EXCHANGE) then
          do ib = 1, n_block
            ex_per_vol(ib2) = ex_per_vol(ib2) + sum(l_dens(1:spin_channels, ib)) * l_zk(ib)
            ib2 = ib2 + 1
          end do
        else
          do ib = 1, n_block
            ec_per_vol(ib2) = ec_per_vol(ib2) + sum(l_dens(1:spin_channels, ib)) * l_zk(ib)
            ib2 = ib2 + 1
          end do
        end if
      end if

      ! store results
      if(present(vxc)) then

        ib2 = ip
        do ib = 1, n_block
          dedd(ib2, 1:spin_channels) = dedd(ib2, 1:spin_channels) + l_dedd(1:spin_channels, ib)
          ib2 = ib2 + 1
        end do

        ! calculate the spin unpolarized exchange potential for the long range correction
        if(xcs%xc_density_correction == LR_X .and. &
          (functl(ixc)%type == XC_EXCHANGE .or. functl(ixc)%type == XC_EXCHANGE_CORRELATION)) then

          SAFE_ALLOCATE(unp_dens(1:n_block))
          SAFE_ALLOCATE(unp_dedd(1:n_block))

          do ib = 1, n_block
            unp_dens(ib) = sum(l_dens(1:spin_channels, ib))
          end do

          select case(functl(ixc)%family)
          case(XC_FAMILY_LDA)
            call XC_F90(lda_vxc)(xcs%functl(ixc, 1)%conf, n_block, unp_dens(1), unp_dedd(1))

          case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
            l_vsigma = M_ZERO

            call messages_not_implemented('XC density correction for GGA/mGGA')

            if(functl(ixc)%id == XC_GGA_X_LB) then
              call mesh_r(der%mesh, ip, rr)
              call XC_F90(gga_lb_modified)(xcs%functl(ixc, 1)%conf, n_block, unp_dens(1), l_sigma(1,1), &
                rr, unp_dedd(1))
            else
              call XC_F90(gga_vxc)(xcs%functl(ixc, 1)%conf, n_block, unp_dens(1), l_sigma(1,1), &
                unp_dedd(1), l_vsigma(1,1))
            end if
          end select

          ib2 = ip
          do ib = 1, n_block
            vx(ib2) = unp_dedd(ib)
            ib2 = ib2 + 1
          end do

          ! GGA terms are missing here

          SAFE_DEALLOCATE_A(unp_dens)
          SAFE_DEALLOCATE_A(unp_dedd)
        end if

        if((functl(ixc)%family == XC_FAMILY_GGA).or.(functl(ixc)%family == XC_FAMILY_MGGA)) then
          ib2 = ip
          do ib = 1, n_block
            dedgd(ib2,:,1) = dedgd(ib2,:,1) + M_TWO*l_vsigma(1, ib)*gdens(ib2,:,1)
            if(ispin /= UNPOLARIZED) then
              dedgd(ib2,:,1) = dedgd(ib2,:,1) + l_vsigma(2, ib)*gdens(ib2,:,2)
              dedgd(ib2,:,2) = dedgd(ib2,:,2) +  &
                M_TWO*l_vsigma(3, ib)*gdens(ib2,:,2) + l_vsigma(2, ib)*gdens(ib2,:,1)
            end if
            ib2 = ib2 + 1
          end do
        end if

        if(functl(ixc)%family == XC_FAMILY_MGGA) then
          ib2 = ip
          do ib = 1, n_block
            dedldens(ib2, 1:spin_channels) = dedldens(ib2, 1:spin_channels) + l_dedldens(1:spin_channels, ib)
            vtau    (ib2, 1:spin_channels) = vtau    (ib2, 1:spin_channels) + l_dedtau  (1:spin_channels, ib)
            ib2 = ib2 + 1
          end do
        end if
      end if

    end do functl_loop
  end do space_loop

  if(present(deltaxc)) deltaxc = M_ZERO

  if(xcs%xc_density_correction == LR_X) then
    call xc_density_correction_calc(xcs, der, spin_channels, rho, vx, dedd, deltaxc = deltaxc)

    if(calc_energy) then
      ! correct the energy density from Levy-Perdew, note that vx now
      ! contains the correction applied to the xc potential.
      do is = 1, spin_channels
        do ip = 1, der%mesh%np
          ex_per_vol(ip) = ex_per_vol(ip) &
            + vx(ip)*(CNST(3.0)*rho(ip, is) + sum(der%mesh%x(ip, 1:der%mesh%sb%dim)*gdens(ip, 1:der%mesh%sb%dim, is)))
        end do
      end do
    end if
  end if

  ! this has to be done in inverse order
  if(present(vxc)) then
    if(mgga) call mgga_process()
    if( gga) call  gga_process()
    call lda_process()
  end if

  if(calc_energy) then
    ! integrate energies per unit volume
    ex = ex + dmf_integrate(der%mesh, ex_per_vol)
    ec = ec + dmf_integrate(der%mesh, ec_per_vol)
  end if

  ! clean up allocated memory
  call lda_end()
  if(gga .or. xcs%xc_density_correction == LR_X) call  gga_end()
  if(mgga) call mgga_end()

  POP_SUB(dxc_get_vxc)
  call profiling_out(prof)

contains

  ! ---------------------------------------------------------
  ! Takes care of the initialization of the LDA part of the functionals
  !   *) allocates density and dedd, and their local variants
  !   *) calculates the density taking into account nlcc and non-collinear spin
  subroutine lda_init()
    integer :: ii
    FLOAT   :: d(2), dtot, dpol

    PUSH_SUB(xc_get_vxc.lda_init)

    ! allocate some general arrays
    SAFE_ALLOCATE(l_dens(1:spin_channels, 1:n_block))
    SAFE_ALLOCATE(l_zk(1:n_block))

    SAFE_ALLOCATE(dens(1:der%mesh%np_part, 1:spin_channels))
    dens       = M_ZERO

    if(calc_energy) then
      SAFE_ALLOCATE(ex_per_vol(1:der%mesh%np))
      SAFE_ALLOCATE(ec_per_vol(1:der%mesh%np))
      ex_per_vol = M_ZERO
      ec_per_vol = M_ZERO
    end if

    if(present(vxc)) then
      SAFE_ALLOCATE(l_dedd(1:spin_channels, 1:n_block))
      SAFE_ALLOCATE(dedd(1:der%mesh%np_part, 1:spin_channels))
      dedd = M_ZERO
    end if

    do ii = 1, der%mesh%np
      d(1:spin_channels) = rho(ii, 1:spin_channels)

      select case(ispin)
      case(UNPOLARIZED)
        dens(ii, 1) = max(d(1), M_ZERO)
      case(SPIN_POLARIZED)
        dens(ii, 1) = max(d(1), M_ZERO)
        dens(ii, 2) = max(d(2), M_ZERO)
      case(SPINORS)
        dtot = d(1) + d(2)
        dpol = sqrt((d(1) - d(2))**2 + &
          M_FOUR*(rho(ii, 3)**2 + rho(ii, 4)**2))
        dens(ii, 1) = max(M_HALF*(dtot + dpol), M_ZERO)
        dens(ii, 2) = max(M_HALF*(dtot - dpol), M_ZERO)
      end select
    end do

    POP_SUB(xc_get_vxc.lda_init)
  end subroutine lda_init


  ! ---------------------------------------------------------
  ! deallocate variables allocated in lda_init
  subroutine lda_end()
    PUSH_SUB(xc_get_vxc.lda_end)

    SAFE_DEALLOCATE_A(l_dens)
    SAFE_DEALLOCATE_A(l_zk)

    SAFE_DEALLOCATE_A(dens)
    SAFE_DEALLOCATE_A(ex_per_vol)
    SAFE_DEALLOCATE_A(ec_per_vol)
    if(present(vxc)) then
      SAFE_DEALLOCATE_A(l_dedd)
      SAFE_DEALLOCATE_A(dedd)
    end if

    POP_SUB(xc_get_vxc.lda_end)
  end subroutine lda_end


  ! ---------------------------------------------------------
  ! calculates the LDA part of vxc, taking into account non-collinear spin
  subroutine lda_process()
    integer :: ip
    FLOAT :: d(2), dpol, vpol

    PUSH_SUB(xc_get_vxc.lda_process)

    if(ispin == SPINORS) then
      ! rotate back (do not need the rotation matrix for this).
      do ip = 1, der%mesh%np
        d(1:spin_channels) = rho(ip, 1:spin_channels)

        dpol = sqrt((d(1) - d(2))**2 + &
          M_FOUR*(rho(ip, 3)**2 + rho(ip, 4)**2))
        vpol = (dedd(ip, 1) - dedd(ip, 2))*(d(1) - d(2))/(dpol + tiny)

        vxc(ip, 1) = vxc(ip, 1) + M_HALF*(dedd(ip, 1) + dedd(ip, 2) + vpol)
        vxc(ip, 2) = vxc(ip, 2) + M_HALF*(dedd(ip, 1) + dedd(ip, 2) - vpol)
        vxc(ip, 3) = vxc(ip, 3) + (dedd(ip, 1) - dedd(ip, 2))*rho(ip, 3)/(dpol + tiny)
        vxc(ip, 4) = vxc(ip, 4) + (dedd(ip, 1) - dedd(ip, 2))*rho(ip, 4)/(dpol + tiny)
      end do
    elseif(ispin == SPIN_POLARIZED) then
      call lalg_axpy(der%mesh%np, M_ONE, dedd(:, 1), vxc(:, 1))
      call lalg_axpy(der%mesh%np, M_ONE, dedd(:, 2), vxc(:, 2))
    else
      call lalg_axpy(der%mesh%np, M_ONE, dedd(:, 1), vxc(:, 1))
    end if

    POP_SUB(xc_get_vxc.lda_process)
  end subroutine lda_process


  ! ---------------------------------------------------------
  ! initialize GGAs
  !   *) allocates gradient of the density (gdens), dedgd, and its local variants

  subroutine gga_init()
    integer :: ii

    PUSH_SUB(xc_get_vxc.gga_init)

    ii = 1
    if(ispin /= UNPOLARIZED) ii = 3

    ! allocate variables
    SAFE_ALLOCATE(l_sigma(1:ii, 1:n_block))
    SAFE_ALLOCATE(gdens(1:der%mesh%np, 1:3, 1:spin_channels))
    gdens = M_ZERO

    if(present(vxc)) then
      SAFE_ALLOCATE(l_vsigma(1:ii, 1:n_block))
      SAFE_ALLOCATE(dedgd(1:der%mesh%np_part, 1:3, 1:spin_channels))
      dedgd = M_ZERO
    end if

    do ii = 1, 2
      if(functl(ii)%id == XC_GGA_X_LB) then
        call XC_F90(gga_lb_set_par)(functl(ii)%conf, &
          functl(ii)%LB94_modified, functl(ii)%LB94_threshold, ioniz_pot, qtot)
      end if
    end do

    POP_SUB(xc_get_vxc.gga_init)
  end subroutine gga_init


  ! ---------------------------------------------------------
  ! cleans up memory allocated in gga_init
  subroutine gga_end()
    PUSH_SUB(xc_get_vxc.gga_end)

    SAFE_DEALLOCATE_A(l_sigma)
    SAFE_DEALLOCATE_A(gdens)
    if(present(vxc)) then
      SAFE_DEALLOCATE_A(l_vsigma)
      SAFE_DEALLOCATE_A(dedgd)
    end if

    POP_SUB(xc_get_vxc.gga_end)
  end subroutine gga_end


  ! ---------------------------------------------------------
  ! calculates the GGA contribution to vxc
  subroutine gga_process()
    integer :: ip, is
    FLOAT, allocatable :: gf(:,:)

    PUSH_SUB(xc_get_vxc.gga_process)

    ! subtract the divergence of the functional derivative of Exc with respect to
    ! the gradient of the density.
    SAFE_ALLOCATE(gf(1:der%mesh%np, 1:1))
    do is = 1, spin_channels
      call dderivatives_div(der, dedgd(:, :, is), gf(:, 1))
      call lalg_axpy(der%mesh%np, -M_ONE, gf(:,1), dedd(:, is))
    end do
    SAFE_DEALLOCATE_A(gf)

    ! If LB94, we can calculate an approximation to the energy from
    ! Levy-Perdew relation PRA 32, 2010 (1985)
    if(calc_energy .and. functl(1)%id == XC_GGA_X_LB) then
      SAFE_ALLOCATE(gf(1:der%mesh%np, 1:3))

      do is = 1, spin_channels
        call dderivatives_grad(der, dedd(:, is), gf(:,:))
        do ip = 1, der%mesh%np
          ex_per_vol(ip) = ex_per_vol(ip) - dens(ip, is) * sum(der%mesh%x(ip,:)*gf(ip,:))
        end do
      end do

      SAFE_DEALLOCATE_A(gf)
    end if

    POP_SUB(xc_get_vxc.gga_process)
  end subroutine gga_process


  ! ---------------------------------------------------------
  ! initialize meta-GGAs
  !   *) allocate the kinetic-energy density, dedtau, and local variants

  subroutine mgga_init()
    PUSH_SUB(xc_get_vxc.mgga_init)

    ! allocate variables
    SAFE_ALLOCATE( tau(1:der%mesh%np, 1:spin_channels))
    SAFE_ALLOCATE( current (1:der%mesh%np, 1:der%mesh%sb%dim, 1:spin_channels) )
    SAFE_ALLOCATE(ldens(1:der%mesh%np, 1:spin_channels))

    SAFE_ALLOCATE(l_tau  (1:spin_channels, 1:n_block))
    SAFE_ALLOCATE(l_ldens(1:spin_channels, 1:n_block))

    if(present(vxc)) then
      SAFE_ALLOCATE(dedldens(1:der%mesh%np_part, 1:spin_channels))
      dedldens = M_ZERO

      SAFE_ALLOCATE(l_dedtau  (1:spin_channels, 1:n_block))
      SAFE_ALLOCATE(l_dedldens(1:spin_channels, 1:n_block))
    end if

    POP_SUB(xc_get_vxc.mgga_init)
  end subroutine mgga_init


  ! ---------------------------------------------------------

  subroutine calc_tb09_c()
    FLOAT, allocatable :: gnon(:)
    FLOAT gn(MAX_DIM), n, tb09_c
    integer :: ii

    PUSH_SUB(xc_get_vxc.calc_tb09_c)

    SAFE_ALLOCATE(gnon(1:der%mesh%np))

    do ii = 1, der%mesh%np
      if(ispin == UNPOLARIZED) then
        n = dens(ii, 1)
        gn(1:der%mesh%sb%dim) = gdens(ii, 1:der%mesh%sb%dim, 1)
      else
        n = dens(ii, 1) + dens(ii, 2)
        gn(1:der%mesh%sb%dim) = gdens(ii, 1:der%mesh%sb%dim, 1) + gdens(ii, 1:der%mesh%sb%dim, 2)
      end if

      if (n <= CNST(1e-7)) then 
        gnon(ii) = CNST(0.0)
        ! here you will have to print the true gnon(ii) with the correspondent mesh point ii
      else
        gnon(ii) = sqrt(sum((gn(1:der%mesh%sb%dim)/n)**2))
      end if
    end do

    tb09_c =  -CNST(0.012) + CNST(1.023)*sqrt(dmf_integrate(der%mesh, gnon)/der%mesh%sb%rcell_volume)

    write(message(1), '(a,f13.6)') "Info: In the functional TB09 c = ", tb09_c
    call messages_info(1)

    call  XC_F90(mgga_x_tb09_set_par)(functl(1)%conf, tb09_c)

    SAFE_DEALLOCATE_A(gnon)

    POP_SUB(xc_get_vxc.calc_tb09_c)
  end subroutine calc_tb09_c


  ! ---------------------------------------------------------
  ! clean up memory allocates in mgga_init
  subroutine mgga_end()
    PUSH_SUB(xc_get_vxc.mgga_end)

    SAFE_DEALLOCATE_A(tau)
    SAFE_DEALLOCATE_A(current)
    SAFE_DEALLOCATE_A(ldens)

    SAFE_DEALLOCATE_A(l_tau)
    SAFE_DEALLOCATE_A(l_ldens)

    if(present(vxc)) then
      SAFE_DEALLOCATE_A(dedldens)

      SAFE_DEALLOCATE_A(l_dedtau)
      SAFE_DEALLOCATE_A(l_dedldens)
    end if

    POP_SUB(xc_get_vxc.mgga_end)
  end subroutine mgga_end


  ! ---------------------------------------------------------
  ! calculate the mgga contribution to vxc
  subroutine mgga_process()
    integer :: is
    FLOAT, allocatable :: lf(:,:)

    PUSH_SUB(xc_get_vxc.mgga_process)

    ! add the Laplacian of the functional derivative of Exc with respect to
    ! the laplacian of the density.

    SAFE_ALLOCATE(lf(1:der%mesh%np, 1:1))
    do is = 1, spin_channels
      call dderivatives_lapl(der, dedldens(:, is), lf(:, 1))
      call lalg_axpy(der%mesh%np, M_ONE, lf(:, 1), dedd(:, is))
    end do
    SAFE_DEALLOCATE_A(lf)

    POP_SUB(xc_get_vxc.mgga_process)
  end subroutine mgga_process

end subroutine dxc_get_vxc

! -----------------------------------------------------

subroutine xc_density_correction_calc(xcs, der, nspin, density, refvx, vxc, deltaxc)
  type(xc_t),          intent(in)    :: xcs
  type(derivatives_t), intent(in)    :: der
  integer,             intent(in)    :: nspin
  FLOAT,               intent(in)    :: density(:, :)
  FLOAT,               intent(inout) :: refvx(:)
  FLOAT,               intent(inout) :: vxc(:, :)
  FLOAT, optional,     intent(out)   :: deltaxc

  logical :: find_root, done
  integer :: ip, iunit, ierr
  integer, save :: iter = 0
  FLOAT,   save :: ncsave
  character(len=30) :: number
  FLOAT   :: qxc, ncutoff, qxc_old, ncutoff_old, deriv, deriv_old, qxcfin
  FLOAT   :: x1, x2, x3, f1, f2, f3, dd, vol, mindd, maxdd
  FLOAT, allocatable :: nxc(:), lrvxc(:)
  type(profile_t), save :: prof
  FLOAT, parameter :: thres = CNST(1e-6)

  PUSH_SUB('vxc_inc.xc_density_correction_calc')

  call profiling_in(prof, "XC_DENSITY_CORRECTION")

  SAFE_ALLOCATE(nxc(1:der%mesh%np))
  SAFE_ALLOCATE(lrvxc(1:der%mesh%np_part))

  forall(ip = 1:der%mesh%np) lrvxc(ip) = CNST(-1.0)/(CNST(4.0)*M_PI)*refvx(ip)
  call dderivatives_lapl(der, lrvxc, nxc)

  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "rho", der%mesh, density(:, 1), unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "vxcorig", der%mesh, refvx(:), unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "nxc", der%mesh, nxc, unit_one, ierr)

  if(xcs%xcd_optimize_cutoff) then

    x1 = CNST(1.0e-8)
    qxc = get_qxc(der%mesh, nxc, density(:, 1), x1)
    deriv = HUGE(deriv)
    done = .false.

    INCR(iter, 1)
    if(mpi_world%rank == 0) then
      write(number, '(i4)') iter
      iunit = io_open('qxc.'//trim(adjustl(number)), action='write')
    end if
    do
      if(.not. done) then
        ncutoff_old = x1
        qxc_old = qxc
        deriv_old = deriv
      end if

      x1 = x1*CNST(1.01)

      if(mpi_world%rank == 0) then
        write(iunit, *) x1, qxc
      end if

      if(x1 > CNST(1.0)) exit

      qxc = get_qxc(der%mesh, nxc, density(:, 1), x1)

      if(qxc == qxc_old) cycle

      deriv = (qxc - qxc_old)/(x1 - ncutoff_old)

      if(.not. done .and. abs(qxc) >= 1.0_8) then
        find_root = .true.
        done = .true.
        ncutoff = x1
      end if

      if(xcs%xcd_minimum .and. .not. done .and. abs(qxc_old) - abs(qxc) > thres) then
        find_root = .false.
        done = .true.
        ncutoff = x1
        print*, x1, 0
      end if

    end do

    if(mpi_world%rank == 0) call io_close(iunit)

    if(iter > 1) x3 = ncsave

    if(find_root) then
      x1 = ncutoff
      x2 = ncutoff_old
      x3 = ncutoff
      f1 = 1.0_8 + get_qxc(der%mesh, nxc, density(:, 1), x1)
      f2 = 1.0_8 + get_qxc(der%mesh, nxc, density(:, 1), x2)

      do ip = 1, 20
        if(abs(f1 - f2) < 1e-16_8) exit
        x3 = x2 - f2*(x2 - x1)/(f2 - f1)
        f3 = 1.0_8 + get_qxc(der%mesh, nxc, density(:, 1), x3)
        if(abs(f3) < 1e-6_8) exit
        x1 = x2
        f1 = f2
        x2 = x3
        f2 = f3
      end do

      if(x3 <= ncutoff) ncutoff = x3
    end if

    ncsave = x3

  else
    ncutoff = xcs%xcd_ncutoff
  end if

  qxcfin = get_qxc(der%mesh, nxc, density(:, 1), ncutoff)

  do ip = 1, der%mesh%np
    if(density(ip, 1) < ncutoff) then
      nxc(ip) = -nxc(ip)
    else
      nxc(ip) = M_ZERO
    end if
  end do

  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "nxcmod", der%mesh, nxc, unit_one, ierr)

  if(mpi_world%rank == 0) then
    print*, "Iter",    iter, ncutoff, qxcfin
  end if

  call dpoisson_solve(psolver, lrvxc, nxc)

  if(xcs%xcd_normalize .and. abs(qxcfin) > CNST(1e-10)) then
    do ip = 1, der%mesh%np
      lrvxc(ip) = lrvxc(ip)/abs(qxcfin)
    end do
  end if

  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "fulldiffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_AXIS_Y, "./static", "fulldiffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_AXIS_Z, "./static", "fulldiffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_PLANE_X, "./static", "fulldiffvxc.pl", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_PLANE_Y, "./static", "fulldiffvxc.pl", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_PLANE_Z, "./static", "fulldiffvxc.pl", der%mesh, lrvxc, unit_one, ierr)

  forall(ip = 1:der%mesh%np) 
    vxc(ip, 1:nspin) = vxc(ip, 1:nspin) + lrvxc(ip)
    refvx(ip) = lrvxc(ip)
  end forall

  maxdd = -HUGE(maxdd)
  mindd =  HUGE(maxdd)
  vol = M_ZERO
  do ip = 1, der%mesh%np
    if(density(ip, 1) >= ncutoff) then
      vol = vol + der%mesh%volume_element
      maxdd = max(lrvxc(ip), maxdd)
      mindd = min(lrvxc(ip), mindd)
    else
      lrvxc(ip) = M_ZERO
    end if
  end do

  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "diffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_AXIS_Y, "./static", "diffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
  call dio_function_output(C_OUTPUT_HOW_AXIS_Z, "./static", "diffvxc.ax", der%mesh, lrvxc, unit_one, ierr)

  dd = dmf_integrate(der%mesh, lrvxc)/vol

  if(mpi_world%rank == 0) then
    print*, "DD",  -CNST(2.0)*dd, -CNST(2.0)*mindd, -CNST(2.0)*maxdd
  end if

  if(present(deltaxc)) deltaxc = -CNST(2.0)*dd

  call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "fnxc", der%mesh, nxc, unit_one, ierr)

  call profiling_out(prof)

  SAFE_DEALLOCATE_A(lrvxc)
  SAFE_DEALLOCATE_A(nxc)

  POP_SUB('vxc_inc.xc_density_correction_calc')
end subroutine xc_density_correction_calc

! -----------------------------------------------------

FLOAT function get_qxc(mesh, nxc, density, ncutoff)  result(qxc)
  type(mesh_t), intent(in) :: mesh
  FLOAT,        intent(in) :: nxc(:)
  FLOAT,        intent(in) :: density(:)
  FLOAT,        intent(in) :: ncutoff

  integer :: ip
  FLOAT, allocatable :: nxc2(:)

  PUSH_SUB('vxc_inc.get_qxc')

  SAFE_ALLOCATE(nxc2(1:mesh%np))

  do ip = 1, mesh%np
    if(density(ip) < ncutoff) then
      nxc2(ip) = 0.0
    else
      nxc2(ip) = nxc(ip)
    end if
  end do

  qxc = dmf_integrate(mesh, nxc2)

  SAFE_DEALLOCATE_A(nxc2)

  POP_SUB('vxc_inc.get_qxc')
end function get_qxc

subroutine zxc_complex_lda(mesh, rho, vxc, ex, ec, Imrho, Imvxc, Imex, Imec, cmplxscl_th)
  type(mesh_t), intent(in) :: mesh
  FLOAT, intent(in)        :: rho(:, :)
  FLOAT, intent(inout)     :: vxc(:, :)
  FLOAT, intent(inout)     :: ex
  FLOAT, intent(inout)     :: ec
  FLOAT, intent(in)        :: Imrho(:, :)
  FLOAT, intent(inout)     :: Imvxc(:, :)
  FLOAT, intent(inout)     :: Imex
  FLOAT, intent(inout)     :: Imec
  FLOAT, intent(in)        :: cmplxscl_th
  
  CMPLX :: zex, zec, zrho, zvxc, eps_c, last_zvxc
  INTEGER :: i, N
  CMPLX :: rs, rtrs, Q0, Q1, dQ1drs, dedrs, tmpphase, dimphase, vtrial2, vtrial3
  CMPLX, allocatable :: zvxc_arr(:)

  FLOAT :: C0I, C1, CC1, CC2, IF2, gamma, alpha1, beta1, beta2, beta3, beta4, Cx

  ! LDA constants.
  ! Only C0I is used for spin-paired calculations among these five
  C0I = 0.238732414637843
  C1 = -0.45816529328314287
  CC1 = 1.9236610509315362
  CC2 = 2.5648814012420482
  IF2 = 0.58482236226346462
  
  gamma = 0.031091
  alpha1 = 0.21370
  beta1 = 7.5957
  beta2 = 3.5876
  beta3 = 1.6382
  beta4 = 0.49294

  N = size(rho, 1)
  SAFE_ALLOCATE(zvxc_arr(1:N))

  zex = M_z0
  zec = M_z0

  dimphase = exp(-mesh%sb%dim * M_zI * cmplxscl_th)

  !Cx = -3.0 / 4.0 * (3.0 / M_PI)**(1.0 / 3.0)
  Cx = 0.73855876638202234 

  last_zvxc = M_ONE ! entirely arbitrary

  do i=1, N
     zrho = rho(i, 1) + M_zI * Imrho(i, 1)

     ! "simplified", linear exchange potential
     !zex = zex + 0.5 * lda_exchange_prefactor * zrho * zrho * dimphase
     !zvxc = lda_exchange_prefactor * zrho * dimphase !+ 2.0 * 3.1415926535897931 * M_zI

     ! quadratic positive exchange potential
     !zex = zex - lda_exchange_prefactor * (zrho * dimphase)**3.0 / dimphase / 10.
     !zvxc = -3.0 * lda_exchange_prefactor * (zrho * dimphase)**2.0 / 10.

     ! 3d exchange
     zvxc = -Cx * 4.0 / 3.0 * (zrho * dimphase)**(1.0 / 3.0)

     ! Among the three cube roots, choose the one closest to that of
     ! the last iteration.  This choice is quite arbitrary and
     ! probably wrong.  We will correct it later since it only rotates
     ! the potential by a specific phase.
     vtrial2 = zvxc * exp(M_TWO * M_PI * M_zI / M_THREE)
     vtrial3 = zvxc / exp(M_TWO * M_PI * M_zI / M_THREE)
     if (abs(vtrial2 - last_zvxc).lt.abs(zvxc - last_zvxc)) then
        zvxc = vtrial2
     end if
     if (abs(vtrial3 - last_zvxc).lt.abs(zvxc - last_zvxc)) then
        zvxc = vtrial3
     end if
     last_zvxc = zvxc
     
     zex = zex + 3.0 / 4.0 * zvxc * zrho

     ! 2d exchange
     !zex = zex - (4./3.) * sqrt(2./3.1415926535897931) * zrho**(3.0/2.0)
     !zvxc = -1.5 * (4./3.) * sqrt(2./3.1415926535897931) * zrho**(3.0/2.0)

     ! correlation
     !rs = (C0I / zrho)**(1.0 / 3.0)
     !rtrs = sqrt(rs)
     !Q0 = -2.0 * gamma * (1.0 + alpha1 * rs)
     !Q1 = 2.0 * gamma * rtrs * (beta1 + rtrs * (beta2 + rtrs * (beta3 + rtrs * beta4)))
     !eps_c = Q0 * log(1.0 + 1.0 / Q1)
     zec = M_z0 !!!!zec + eps_c * zrho
     !dQ1drs = gamma * (beta1 / rtrs + 2.0 * beta2 + rtrs * (3.0 * beta3 + 4.0 * beta4 * rtrs))
     !dedrs = -2.0 * gamma * alpha1 * eps_c / Q0 - Q0 * dQ1drs / (Q1 * (Q1 + 1.0))
     !zvxc = zvxc + eps_c - rs * dedrs / 3.0
     
     zvxc_arr(i) = zvxc
  end do
  
  tmpphase = exp(M_TWO * M_PI * M_zI / M_THREE)
  do i=1, 2 ! multiply by the phase up to two times
     if (real(zex * tmpphase).lt.real(zex)) then
        zex = zex * tmpphase
        zvxc_arr(:) = zvxc_arr(:) * tmpphase
     end if
  end do

  vxc(:, 1) = real(zvxc_arr)
  Imvxc(:, 1) = aimag(zvxc_arr)

  zex = zex * mesh%volume_element
  zec = zec * mesh%volume_element

  ex = real(zex)
  ec = real(zec)
  Imex = aimag(zex)
  Imec = aimag(zec)

  SAFE_DEALLOCATE_A(zvxc_arr)
  
  print*, 'lda exchange', zex
  print*, 'lda correlation', zec

  
end subroutine zxc_complex_lda

! ----------------------------------------------------------------------------- 
! This is the complex scaled interface for xc functionals.
! It will eventually be merged with the other one dxc_get_vxc after some test
! -----------------------------------------------------------------------------
subroutine zxc_get_vxc(der, xcs, st, rho, ispin, ioniz_pot, qtot, ex, ec, vxc, vtau, Imrho, Imex, Imec, Imvxc, Imvtau, cmplxscl_th)
  type(derivatives_t),  intent(in)    :: der             !< Discretization and the derivative operators and details
  type(xc_t), target,   intent(in)    :: xcs             !< Details about the xc functional used
  type(states_t),       intent(in)    :: st              !< State of the system (wavefunction,eigenvalues...)
  FLOAT,                intent(in)    :: rho(:, :)       !< Electronic density 
  integer,              intent(in)    :: ispin           !< Number of spin channels 
  FLOAT,                intent(in)    :: ioniz_pot
  FLOAT,                intent(in)    :: qtot 
  FLOAT, optional,      intent(inout) :: ex              !< Exchange energy.
  FLOAT, optional,      intent(inout) :: ec              !< Correlation energy.
  FLOAT, optional,      intent(inout) :: vxc(:,:)        !< XC potential
  FLOAT, optional,      intent(inout) :: vtau(:,:)       !< Derivative wrt (two times kinetic energy density)
  FLOAT,                intent(in)    :: Imrho(:, :)     !< cmplxscl: Electronic density 
  FLOAT, optional,      intent(inout) :: Imex            !< cmplxscl: Exchange energy.
  FLOAT, optional,      intent(inout) :: Imec            !< cmplxscl: Correlation energy
  FLOAT, optional,      intent(inout) :: Imvxc(:,:)      !< cmplxscl: XC potential
  FLOAT, optional,      intent(inout) :: Imvtau(:,:)     !< cmplxscl: Derivative wrt (two times kinetic energy density)
  FLOAT,                intent(in)    :: cmplxscl_th     !< complex scaling angle

  
  CMPLX, pointer :: zpot(:), zrho_tot(:)
  CMPLX          :: ztmp
  Integer        :: isp
  type(xc_functl_t), pointer :: functl(:)
  logical         :: calc_energy

  PUSH_SUB(zxc_get_vxc)

  print *, "LDA calc energy exc"
  ASSERT(present(ex) .eqv. present(ec))
  calc_energy = present(ex)

  !Pointer-shortcut for xcs%functl
  !It helps to remember that for xcs%functl(:,:)
  ! (1,:) => exchange,    (2,:) => correlation
  ! (:,1) => unpolarized, (:,2) => polarized
  if(ispin == UNPOLARIZED) then
    functl => xcs%functl(:, 1)
  else
    functl => xcs%functl(:, 2)
  end if


  
  if(functl(1)%id == XC_LDA_XC_CMPLX) then
    
    call zxc_complex_lda(der%mesh, rho, vxc, ex, ec, Imrho, Imvxc, Imex, Imec, cmplxscl_th)

    ! Exact exchange for 2 particles [vxc(r) = 1/2 * vh(r)]
    ! we keep it here for debug purposes
    if(.false.) then
      SAFE_ALLOCATE(zpot(1:size(vxc,1)))
      SAFE_ALLOCATE(zrho_tot(1:size(vxc,1)))

      zrho_tot = M_z0
      do isp = 1, ispin
        zrho_tot(:) = zrho_tot(:)+ rho(:,isp) +M_zI * Imrho(:,isp)
      end do

      call zpoisson_solve(psolver, zpot, zrho_tot, theta = cmplxscl_th)

      zpot = - zpot /CNST(2.0)
      vxc(:,1) = real(zpot(:)) 
      Imvxc(:,1) = aimag(zpot(:))
  
      if(calc_energy) then
        ztmp = M_HALF *zmf_dotp(der%mesh, zrho_tot, zpot, dotu = .true. )
        ex =  real(ztmp)
        Imex = aimag(ztmp)
        ec   = M_ZERO
        Imec = M_ZERO    
      end if
      SAFE_DEALLOCATE_P(zrho_tot)
      SAFE_DEALLOCATE_P(zpot)
    end if
    
  else if(functl(1)%family == XC_FAMILY_NONE) then

    vxc = M_ZERO
    if(calc_energy) then
      ex   = M_ZERO
      Imex = M_ZERO
      ec   = M_ZERO
      Imec = M_ZERO
    end if

  else  
    write(message(1), '(a)') 'The selected XCFunctional will not work with ComplexScaling = yes.'
    write(message(2), '(a)') 'Use XCFunctional = lda_xc_cmplx.'
    call messages_fatal(2)     
  end if

  POP_SUB(zxc_get_vxc)
end subroutine zxc_get_vxc


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

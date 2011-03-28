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

! ---------------------------------------------------------
subroutine xc_get_vxc(der, xcs, st, rho, ispin, ioniz_pot, qtot, ex, ec, vxc, vtau)
  type(derivatives_t),  intent(inout) :: der             !< Discretization and the derivative operators and details
  type(xc_t), target,   intent(in)    :: xcs             !< Details about the xc functional used
  type(states_t),       intent(inout) :: st              !< State of the system (wavefunction,eigenvalues...)
  FLOAT,                intent(in)    :: rho(:, :)       !< Electronic density 
  integer,              intent(in)    :: ispin           !< Number of spin channels 
  FLOAT,                intent(in)    :: ioniz_pot
  FLOAT,                intent(in)    :: qtot 
  FLOAT, optional,      intent(inout) :: ex              !< Exchange energy.
  FLOAT, optional,      intent(inout) :: ec              !< Correlation energy.
  FLOAT, optional,      intent(inout) :: vxc(:,:)        !< XC potential
  FLOAT, optional,      intent(inout) :: vtau(:,:)       !< Derivative wrt (two times kinetic energy density)

  integer :: n_block

  FLOAT, allocatable :: l_zk(:)        ! Local block of the energy functional (with the correct memory order for libxc)
  FLOAT, allocatable :: l_dens(:,:)    ! Local block for the density 
  FLOAT, allocatable :: l_dedd(:,:)    ! Local block of the xchange or correl. potential(with the correct memory order for libxc)
  FLOAT, allocatable :: l_sigma(:,:)   
  FLOAT, allocatable :: l_vsigma(:,:)  
  FLOAT, allocatable :: l_tau(:,:)
  FLOAT, allocatable :: l_ldens(:,:)
  FLOAT, allocatable :: l_dedtau(:,:)
  FLOAT, allocatable :: l_dedldens(:,:)

  FLOAT, allocatable :: dens(:,:)      ! Density
  FLOAT, allocatable :: dedd(:,:)      ! (Functional) Derivative of the xchange or correlation energy with
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

  integer :: ib, ib2, ip, isp, families, ixc, spin_channels
  FLOAT   :: rr
  logical :: gga, mgga
  type(profile_t), save :: prof
  logical :: calc_energy
  type(xc_functl_t), pointer :: functl(:)
  type(symmetrizer_t) :: symmetrizer

  PUSH_SUB(xc_get_vxc)
  call profiling_in(prof, "XC_LOCAL")

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

  ! is there anything to do ?
  families = XC_FAMILY_LDA + XC_FAMILY_GGA + XC_FAMILY_HYB_GGA + XC_FAMILY_MGGA
  if(iand(xcs%family, families) == 0) then
    POP_SUB(xc_get_vxc)
    call profiling_out(prof)
    return
  endif

  n_block = 1000

  ! initialize a couple of handy variables
  gga  = iand(xcs%family, XC_FAMILY_GGA + XC_FAMILY_HYB_GGA + XC_FAMILY_MGGA).ne.0
  mgga = iand(xcs%family, XC_FAMILY_MGGA).ne.0

  !Read the spin channels
  !Index 1 refers to "exchange mode" of type "xc_functl_t". For this purpose using 1 or 2 makes no difference. 
  spin_channels = functl(1)%spin_channels

  call lda_init()
  if( gga) call  gga_init()
  if(mgga) call mgga_init()

  ! Get the gradient and the Laplacian of the density and the kinetic-energy density
  ! We do it here instead of doing it in gga_init and mgga_init in order to 
  ! avoid calling the subroutine states_calc_quantities twice
  if(gga .and. (.not. mgga)) then
    ! get gradient of the density (this is faster than calling states_calc_quantities)
    do isp = 1, spin_channels 
      call dderivatives_grad(der, dens(:, isp), gdens(:, :, isp)) 
    end do
  else if(mgga) then

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

    if(functl(1)%id == XC_MGGA_X_TB09 .and. der%mesh%sb%periodic_dim == 3) then
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
          ! memo: please check the following indexes
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
    functl_loop: do ixc = 1, 2

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

          if(functl(ixc)%id == XC_GGA_C_LB) then
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
  if( gga) call  gga_end()
  if(mgga) call mgga_end()


  POP_SUB(xc_get_vxc)
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

    write(message(1), '(a,f8.6)') "Info: In the functional TB09 c = ", tb09_c
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

end subroutine xc_get_vxc

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

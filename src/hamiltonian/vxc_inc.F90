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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

subroutine xc_get_vxc(der, xcs, st, rho, ispin, ioniz_pot, qtot, vxc, ex, ec, deltaxc, vtau)
  type(derivatives_t),  intent(in)    :: der             !< Discretization and the derivative operators and details
  type(xc_t), target,   intent(in)    :: xcs             !< Details about the xc functional used
  type(states_t),       intent(in)    :: st              !< State of the system (wavefunction,eigenvalues...)
  FLOAT,                intent(in)    :: rho(:, :)       !< Electronic density 
  integer,              intent(in)    :: ispin           !< Number of spin channels 
  FLOAT,                intent(in)    :: ioniz_pot
  FLOAT,                intent(in)    :: qtot 
  FLOAT,                intent(inout) :: vxc(:,:)        !< XC potential
  FLOAT, optional,      intent(inout) :: ex              !< Exchange energy.
  FLOAT, optional,      intent(inout) :: ec              !< Correlation energy.
  FLOAT, optional,      intent(inout) :: deltaxc         !< The XC derivative discontinuity
  FLOAT, optional,      intent(inout) :: vtau(:,:)       !< Derivative wrt (two times kinetic energy density)

  ! Formerly vxc was optional, but I removed this since we always pass vxc, and this simplifies the routine
  ! and avoids some optimization problems. --DAS

  integer, parameter :: N_BLOCK_MAX = 1000
  integer :: n_block

  ! Local blocks (with the correct memory order for libxc):
  !  Input quantities
  FLOAT, allocatable :: l_dens(:,:)     ! Density 
  FLOAT, allocatable :: l_sigma(:,:)    ! Modulus squared of the gradient of the density
  FLOAT, allocatable :: l_ldens(:,:)    ! Laplacian of the density
  FLOAT, allocatable :: l_tau(:,:)      ! Kinetic energy density
  !  Energy
  FLOAT, allocatable :: l_zk(:)
  !  First order (functional) derivatives
  FLOAT, allocatable :: l_dedd(:,:)     ! Derivative of the energy wrt the density
  FLOAT, allocatable :: l_vsigma(:,:)   ! Derivative of the energy wrt sigma
  FLOAT, allocatable :: l_dedldens(:,:) ! Derivative of the energy wrt the laplacian of the density
  FLOAT, allocatable :: l_dedtau(:,:)   ! Derivative of the energy wrt tau

  ! Global arrays
  !  Input quantities
  FLOAT, allocatable :: dens(:,:)      ! Density
  FLOAT, allocatable :: gdens(:,:,:)   ! Gradient of the density
  FLOAT, allocatable :: ldens(:,:)     ! Laplacian of the density
  FLOAT, allocatable :: tau(:,:)       ! Kinetic energy density
  !  Energies
  FLOAT, allocatable :: ex_per_vol(:)  ! Exchange energy per unit volume 
  FLOAT, allocatable :: ec_per_vol(:)  ! Correlation energy per unit volume 
  !  First order (functional) derivatives
  FLOAT, allocatable :: dedd(:,:)      ! Derivative of the exchange or correlation energy wrt the density
  FLOAT, allocatable :: dedgd(:,:,:)   ! Derivative of the exchange or correlation energy wrt the gradient of the density
  FLOAT, allocatable :: dedldens(:,:)  ! Derivative of the exchange or correlation energy wrt the laplacian of the density

  FLOAT, allocatable :: symmtmp(:, :)  ! Temporary vector for the symmetrizer
  FLOAT, allocatable :: vx(:)
  FLOAT, allocatable :: unp_dens(:), unp_dedd(:)

  integer :: ib, ip, isp, families, ixc, spin_channels, is
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

  !Pointer-shortcut for xcs%functional
  !It helps to remember that for xcs%functional(:,:)
  ! (1,:) => exchange,    (2,:) => correlation
  ! (:,1) => unpolarized, (:,2) => polarized
  if(ispin == UNPOLARIZED) then
    functl => xcs%functional(:, 1)
  else
    functl => xcs%functional(:, 2)
  end if

  ! is there anything to do ?
  families = XC_FAMILY_LDA + XC_FAMILY_GGA + XC_FAMILY_HYB_GGA + XC_FAMILY_MGGA + XC_FAMILY_HYB_MGGA
  if(iand(xcs%family, families) == 0) then
    POP_SUB(xc_get_vxc)
    call profiling_out(prof)
    return
  endif

  ! initialize a couple of handy variables
  gga  = iand(xcs%family, XC_FAMILY_GGA + XC_FAMILY_HYB_GGA + XC_FAMILY_MGGA + XC_FAMILY_HYB_MGGA) /= 0
  mgga = iand(xcs%family, XC_FAMILY_MGGA + XC_FAMILY_HYB_MGGA) /= 0

  !Read the spin channels
  spin_channels = functl(FUNC_X)%spin_channels
  
  if(xcs%xc_density_correction == LR_X) then
    SAFE_ALLOCATE(vx(1:der%mesh%np))
  end if

  call lda_init()
  if(gga .or. xcs%xc_density_correction == LR_X) call gga_init()
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
    ! We calculate everything from the wavefunctions to benefit from
    ! the error cancellation between the gradient of the density and
    ! tau.
    !
    ! FIXME: Probably this is wrong for non-local corrections or other
    ! cases when the density is not directly generated by the
    ! orbitals.
    if (xcs%use_gi_ked) then
      call states_calc_quantities(der, st, gi_kinetic_energy_density = tau, density_gradient = gdens, density_laplacian = ldens)
    else
      call states_calc_quantities(der, st, kinetic_energy_density = tau, density_gradient = gdens, density_laplacian = ldens)
    end if

    ! We have to symmetrize everything as they are calculated from the
    ! wavefunctions.
    if(st%symmetrize_density) then
      SAFE_ALLOCATE(symmtmp(1:der%mesh%np, 1:der%mesh%sb%dim))
      call symmetrizer_init(symmetrizer, der%mesh)
      do isp = 1, spin_channels
        call dsymmetrizer_apply(symmetrizer, field = tau(:, isp), symmfield = symmtmp(:, 1))
        tau(1:der%mesh%np, isp) = symmtmp(1:der%mesh%np, 1)
        call dsymmetrizer_apply(symmetrizer, field = ldens(:, isp), symmfield = symmtmp(:, 1))
        ldens(1:der%mesh%np, isp) = symmtmp(1:der%mesh%np, 1)
        call dsymmetrizer_apply(symmetrizer, field_vector = gdens(:, :, isp), symmfield_vector = symmtmp)
        gdens(1:der%mesh%np, 1:der%mesh%sb%dim, isp) = symmtmp(1:der%mesh%np, 1:der%mesh%sb%dim)
      end do

      call symmetrizer_end(symmetrizer)
      SAFE_DEALLOCATE_A(symmtmp)
    end if

    if(functl(FUNC_X)%id == XC_MGGA_X_TB09 .and. der%mesh%sb%periodic_dim == 3) then
      call calc_tb09_c()
    end if

  end if

  call local_allocate()

  space_loop: do ip = 1, der%mesh%np, N_BLOCK_MAX

    call space_loop_init(n_block)

    ! Calculate the potential/gradient density in local reference frame.
    functl_loop: do ixc = FUNC_X, FUNC_C

      if(functl(ixc)%family == XC_FAMILY_NONE) cycle

      if(calc_energy .and. iand(functl(ixc)%flags, XC_FLAGS_HAVE_EXC) /= 0) then
        ! we get the xc energy and potential
        select case(functl(ixc)%family)
        case(XC_FAMILY_LDA)
          call XC_F90(lda_exc_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_zk(1), l_dedd(1,1))

        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          call XC_F90(gga_exc_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), &
            l_zk(1), l_dedd(1,1), l_vsigma(1,1))

        case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
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

        case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
          call XC_F90(mgga_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), l_ldens(1,1), l_tau(1,1), &
            l_dedd(1,1), l_vsigma(1,1), l_dedldens(1,1), l_dedtau(1,1))

        case default
          cycle
        end select

      end if

      if(calc_energy) then
        if(functl(ixc)%type == XC_EXCHANGE) then
          do ib = 1, n_block
            ex_per_vol(ib + ip - 1) = ex_per_vol(ib + ip - 1) + sum(l_dens(1:spin_channels, ib)) * l_zk(ib)
          end do
        else
          do ib = 1, n_block
            ec_per_vol(ib + ip - 1) = ec_per_vol(ib + ip - 1) + sum(l_dens(1:spin_channels, ib)) * l_zk(ib)
          end do
        end if
      end if

      call copy_local_to_global(l_dedd, dedd, n_block, spin_channels, ip)
      
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
          call XC_F90(lda_vxc)(xcs%functional(ixc, 1)%conf, n_block, unp_dens(1), unp_dedd(1))
          
        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
          l_vsigma = M_ZERO
          
          call messages_not_implemented('XC density correction for GGA/mGGA')
          
          if(functl(ixc)%id == XC_GGA_X_LB) then
            call mesh_r(der%mesh, ip, rr)
            call XC_F90(gga_lb_modified)(xcs%functional(ixc, 1)%conf, n_block, unp_dens(1), l_sigma(1,1), &
              rr, unp_dedd(1))
          else
            call XC_F90(gga_vxc)(xcs%functional(ixc, 1)%conf, n_block, unp_dens(1), l_sigma(1,1), &
              unp_dedd(1), l_vsigma(1,1))
          end if
        end select
        
        do ib = 1, n_block
          vx(ib + ip - 1) = unp_dedd(ib)
        end do
        
        ! GGA terms are missing here
        
        ! seems it would be better to allocate and deallocate these arrays outside the space loop.
        SAFE_DEALLOCATE_A(unp_dens)
        SAFE_DEALLOCATE_A(unp_dedd)
      end if

      if((functl(ixc)%family == XC_FAMILY_GGA).or.(functl(ixc)%family == XC_FAMILY_MGGA)) then
        do ib = 1, n_block
          dedgd(ib + ip - 1,:,1) = dedgd(ib + ip - 1,:,1) + M_TWO*l_vsigma(1, ib)*gdens(ib + ip - 1,:,1)
          if(ispin /= UNPOLARIZED) then
            dedgd(ib + ip - 1,:,1) = dedgd(ib + ip - 1,:,1) + l_vsigma(2, ib)*gdens(ib + ip - 1,:,2)
            dedgd(ib + ip - 1,:,2) = dedgd(ib + ip - 1,:,2) +  &
                 M_TWO*l_vsigma(3, ib)*gdens(ib + ip - 1,:,2) + l_vsigma(2, ib)*gdens(ib + ip - 1,:,1)
          end if
        end do
      end if

      if(functl(ixc)%family == XC_FAMILY_MGGA) then
        call copy_local_to_global(l_dedldens, dedldens, n_block, spin_channels, ip)
        call copy_local_to_global(l_dedtau, vtau, n_block, spin_channels, ip)
      end if

    end do functl_loop
  end do space_loop

  call local_deallocate()

  ! Definition of tau in libxc is different, so we need to divide vtau by a factor of two
  if (present(vtau)) vtau = vtau / M_TWO

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
  if(mgga) call mgga_process()
  if( gga) call  gga_process()
  call lda_process()

  if(calc_energy) then
    ! integrate energies per unit volume
    ex = ex + dmf_integrate(der%mesh, ex_per_vol)
    ec = ec + dmf_integrate(der%mesh, ec_per_vol)
  end if

  ! clean up allocated memory
  call lda_end()
  if(gga .or. xcs%xc_density_correction == LR_X) call  gga_end()
  if(mgga) call mgga_end()

  POP_SUB(xc_get_vxc)
  call profiling_out(prof)

contains

  ! ---------------------------------------------------------
  !> make a local copy with the correct memory order for libxc
  subroutine copy_global_to_local(global, local, n_block, spin_channels, ip)
    FLOAT,   intent(in)  :: global(:,:)
    FLOAT,   intent(out) :: local(:,:)
    integer, intent(in)  :: n_block
    integer, intent(in)  :: spin_channels
    integer, intent(in)  :: ip

    integer :: ib

    PUSH_SUB(xc_get_vxc.copy_global_to_local)

    do ib = 1, n_block
      local(1:spin_channels, ib) = global(ib + ip - 1, 1:spin_channels)
    end do

    POP_SUB(xc_get_vxc.copy_global_to_local)
  end subroutine copy_global_to_local

  ! ---------------------------------------------------------
  subroutine copy_local_to_global(local, global, n_block, spin_channels, ip)
    FLOAT,   intent(in)    :: local(:,:)
    FLOAT,   intent(inout) :: global(:,:)
    integer, intent(in)    :: n_block
    integer, intent(in)    :: spin_channels
    integer, intent(in)    :: ip

    integer :: ib

    PUSH_SUB(xc_get_vxc.copy_local_to_global)

    do ib = 1, n_block
      global(ib + ip - 1, 1:spin_channels) = global(ib + ip - 1, 1:spin_channels) + local(1:spin_channels, ib)
    end do

    POP_SUB(xc_get_vxc.copy_local_to_global)
  end subroutine copy_local_to_global

  ! ---------------------------------------------------------
  subroutine space_loop_init(nblock)
    integer, intent(out) :: nblock

    PUSH_SUB(xc_get_vxc.space_loop_init)

    !Resize the dimension of the last block when the number of the mesh points
    !it is not a perfect divisor of the dimension of the blocks.
    nblock = min(der%mesh%np - ip + 1, N_BLOCK_MAX)

    ! make a local copy with the correct memory order for libxc
    call copy_global_to_local(dens, l_dens, nblock, spin_channels, ip)

    if(gga) then
      do ib = 1, nblock
        l_sigma(1, ib) = sum(gdens(ib + ip - 1, 1:der%mesh%sb%dim, 1)**2)
        if(ispin /= UNPOLARIZED) then
          ! memo: please check the following indices
          l_sigma(2, ib) = sum(gdens(ib + ip - 1, 1:der%mesh%sb%dim, 1)*gdens(ib + ip - 1, 1:der%mesh%sb%dim, 2)) 
          l_sigma(3, ib) = sum(gdens(ib + ip - 1, 1:der%mesh%sb%dim, 2)**2)
        end if
      end do
    end if

    if(mgga) then
      call copy_global_to_local(tau, l_tau, nblock, spin_channels, ip)
      ! we adjust for the different definition of tau in libxc
      l_tau(1:spin_channels, 1:nblock) = l_tau(1:spin_channels, 1:nblock) / M_TWO
      call copy_global_to_local(ldens, l_ldens, nblock, spin_channels, ip)
    end if

    POP_SUB(xc_get_vxc.space_loop_init)
  end subroutine space_loop_init

  ! ---------------------------------------------------------
  !> Takes care of the initialization of the LDA part of the functionals
  !!   *) allocates density and dedd, and their local variants
  !!   *) calculates the density taking into account nlcc and non-collinear spin
  subroutine lda_init()
    integer :: ii
    FLOAT   :: d(2), dtot, dpol

    PUSH_SUB(xc_get_vxc.lda_init)

    ! allocate some general arrays

    SAFE_ALLOCATE(dens(1:der%mesh%np_part, 1:spin_channels))
    dens       = M_ZERO

    if(calc_energy) then
      SAFE_ALLOCATE(ex_per_vol(1:der%mesh%np))
      SAFE_ALLOCATE(ec_per_vol(1:der%mesh%np))
      ex_per_vol = M_ZERO
      ec_per_vol = M_ZERO
    end if

    SAFE_ALLOCATE(dedd(1:der%mesh%np_part, 1:spin_channels))
    dedd = M_ZERO
    
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
  !> deallocate variables allocated in lda_init
  subroutine lda_end()
    PUSH_SUB(xc_get_vxc.lda_end)

    SAFE_DEALLOCATE_A(dens)
    SAFE_DEALLOCATE_A(ex_per_vol)
    SAFE_DEALLOCATE_A(ec_per_vol)
    SAFE_DEALLOCATE_A(dedd)

    POP_SUB(xc_get_vxc.lda_end)
  end subroutine lda_end


  ! ---------------------------------------------------------
  !> calculates the LDA part of vxc, taking into account non-collinear spin
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
  !> initialize GGAs
  !!   *) allocates gradient of the density (gdens), dedgd, and its local variants
  subroutine gga_init()
    integer :: ii

    PUSH_SUB(xc_get_vxc.gga_init)

    ! allocate variables
    SAFE_ALLOCATE(gdens(1:der%mesh%np, 1:der%mesh%sb%dim, 1:spin_channels))
    gdens = M_ZERO

    SAFE_ALLOCATE(dedgd(1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:spin_channels))
    dedgd = M_ZERO

    do ii = 1, 2
      if(functl(ii)%id == XC_GGA_X_LB) then
        call XC_F90(gga_lb_set_par)(functl(ii)%conf, &
          functl(ii)%LB94_modified, functl(ii)%LB94_threshold, ioniz_pot, qtot)
      end if
    end do

    POP_SUB(xc_get_vxc.gga_init)
  end subroutine gga_init


  ! ---------------------------------------------------------
  !> cleans up memory allocated in gga_init
  subroutine gga_end()
    PUSH_SUB(xc_get_vxc.gga_end)

    SAFE_DEALLOCATE_A(gdens)
    SAFE_DEALLOCATE_A(dedgd)

    POP_SUB(xc_get_vxc.gga_end)
  end subroutine gga_end


  ! ---------------------------------------------------------
  !> calculates the GGA contribution to vxc
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
      SAFE_ALLOCATE(gf(1:der%mesh%np, 1:der%mesh%sb%dim))

      do is = 1, spin_channels
        call dderivatives_grad(der, dedd(:, is), gf(:,:))
        do ip = 1, der%mesh%np
          ex_per_vol(ip) = ex_per_vol(ip) - dens(ip, is) * sum(der%mesh%x(ip, 1:der%mesh%sb%dim)*gf(ip, 1:der%mesh%sb%dim))
        end do
      end do

      SAFE_DEALLOCATE_A(gf)
    end if

    POP_SUB(xc_get_vxc.gga_process)
  end subroutine gga_process


  ! ---------------------------------------------------------
  !> initialize meta-GGAs
  !!   *) allocate the kinetic-energy density, dedtau, and local variants
  subroutine mgga_init()
    PUSH_SUB(xc_get_vxc.mgga_init)

    ! allocate variables
    SAFE_ALLOCATE( tau(1:der%mesh%np, 1:spin_channels))
    SAFE_ALLOCATE(ldens(1:der%mesh%np, 1:spin_channels))

    SAFE_ALLOCATE(dedldens(1:der%mesh%np_part, 1:spin_channels))
    dedldens = M_ZERO

    POP_SUB(xc_get_vxc.mgga_init)
  end subroutine mgga_init


  ! ---------------------------------------------------------

  subroutine calc_tb09_c()
    FLOAT, allocatable :: gnon(:)
    FLOAT :: gn(MAX_DIM), n, tb09_c
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
        gnon(ii) = M_ZERO
      else
        gnon(ii) = sqrt(sum((gn(1:der%mesh%sb%dim)/n)**2))
      end if
    end do

    tb09_c =  -CNST(0.012) + CNST(1.023)*sqrt(dmf_integrate(der%mesh, gnon)/der%mesh%sb%rcell_volume)

    call XC_F90(mgga_x_tb09_set_par)(functl(1)%conf, tb09_c)

    SAFE_DEALLOCATE_A(gnon)

    POP_SUB(xc_get_vxc.calc_tb09_c)
  end subroutine calc_tb09_c


  ! ---------------------------------------------------------
  !> clean up memory allocated in mgga_init
  subroutine mgga_end()
    PUSH_SUB(xc_get_vxc.mgga_end)

    SAFE_DEALLOCATE_A(tau)
    SAFE_DEALLOCATE_A(ldens)

    SAFE_DEALLOCATE_A(dedldens)

    POP_SUB(xc_get_vxc.mgga_end)
  end subroutine mgga_end

  ! ---------------------------------------------------------
  !> THREADSAFE (no SAFE ALLOCATE or PUSH/POP SUB)
  subroutine local_allocate()
    integer :: ii

    allocate(l_dens(1:spin_channels, 1:N_BLOCK_MAX))
    allocate(l_zk(1:N_BLOCK_MAX))
    allocate(l_dedd(1:spin_channels, 1:N_BLOCK_MAX))

    if(gga .or. xcs%xc_density_correction == LR_X) then
      ii = 1
      if(ispin /= UNPOLARIZED) ii = 3
      
      allocate(l_sigma(1:ii, 1:N_BLOCK_MAX))
      allocate(l_vsigma(1:ii, 1:N_BLOCK_MAX))
    end if

    if(mgga) then
      allocate(l_tau  (1:spin_channels, 1:N_BLOCK_MAX))
      allocate(l_ldens(1:spin_channels, 1:N_BLOCK_MAX))
      allocate(l_dedtau  (1:spin_channels, 1:N_BLOCK_MAX))
      allocate(l_dedldens(1:spin_channels, 1:N_BLOCK_MAX))
    end if

    end subroutine local_allocate

  ! ---------------------------------------------------------
  !> THREADSAFE (no SAFE ALLOCATE or PUSH/POP SUB)
  subroutine local_deallocate()

    deallocate(l_dens)
    deallocate(l_zk)
    deallocate(l_dedd)

    if(gga .or. xcs%xc_density_correction == LR_X) then
      deallocate(l_sigma)
      deallocate(l_vsigma)
    end if

    if(mgga) then
      deallocate(l_tau)
      deallocate(l_ldens)
      deallocate(l_dedtau)
      deallocate(l_dedldens)
    end if

  end subroutine local_deallocate

  ! ---------------------------------------------------------
  !> calculate the mgga contribution to vxc
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

  PUSH_SUB(xc_density_correction_calc)

  call profiling_in(prof, "XC_DENSITY_CORRECTION")

  SAFE_ALLOCATE(nxc(1:der%mesh%np))
  SAFE_ALLOCATE(lrvxc(1:der%mesh%np_part))

  forall(ip = 1:der%mesh%np) lrvxc(ip) = CNST(-1.0)/(CNST(4.0)*M_PI)*refvx(ip)
  call dderivatives_lapl(der, lrvxc, nxc)

  if(in_debug_mode) then
    call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "rho", der%mesh, density(:, 1), unit_one, ierr)
    call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "vxcorig", der%mesh, refvx(:), unit_one, ierr)
    call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "nxc", der%mesh, nxc, unit_one, ierr)
  end if

  if(xcs%xcd_optimize_cutoff) then

    x1 = CNST(1.0e-8)
    qxc = get_qxc(der%mesh, nxc, density(:, 1), x1)
    deriv = HUGE(deriv)
    done = .false.

    INCR(iter, 1)
    if(in_debug_mode) then
      if(mpi_world%rank == 0) then
        write(number, '(i4)') iter
        iunit = io_open('qxc.'//trim(adjustl(number)), action='write')
      end if
    end if
    do
      if(.not. done) then
        ncutoff_old = x1
        qxc_old = qxc
        deriv_old = deriv
      end if

      x1 = x1*CNST(1.01)
      if(in_debug_mode) then
        if(mpi_world%rank == 0) then
          write(iunit, *) x1, qxc
        end if
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

  if(in_debug_mode) then
    call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "nxcmod", der%mesh, nxc, unit_one, ierr)
  
    if(mpi_world%rank == 0) then
      print*, "Iter",    iter, ncutoff, qxcfin
    end if
  end if

  call dpoisson_solve(psolver, lrvxc, nxc)

  if(xcs%xcd_normalize .and. abs(qxcfin) > CNST(1e-10)) then
    do ip = 1, der%mesh%np
      lrvxc(ip) = lrvxc(ip)/abs(qxcfin)
    end do
  end if

  if(in_debug_mode) then
    call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "fulldiffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(C_OUTPUT_HOW_AXIS_Y, "./static", "fulldiffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(C_OUTPUT_HOW_AXIS_Z, "./static", "fulldiffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(C_OUTPUT_HOW_PLANE_X, "./static", "fulldiffvxc.pl", der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(C_OUTPUT_HOW_PLANE_Y, "./static", "fulldiffvxc.pl", der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(C_OUTPUT_HOW_PLANE_Z, "./static", "fulldiffvxc.pl", der%mesh, lrvxc, unit_one, ierr)
  end if

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

  if(in_debug_mode) then
    call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "diffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(C_OUTPUT_HOW_AXIS_Y, "./static", "diffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(C_OUTPUT_HOW_AXIS_Z, "./static", "diffvxc.ax", der%mesh, lrvxc, unit_one, ierr)
  end if
  
  dd = dmf_integrate(der%mesh, lrvxc)/vol

  if(in_debug_mode) then
    if(mpi_world%rank == 0) then
      print*, "DD",  -CNST(2.0)*dd, -CNST(2.0)*mindd, -CNST(2.0)*maxdd
    end if
  end if
  
  if(present(deltaxc)) deltaxc = -CNST(2.0)*dd

  if(in_debug_mode) then
    call dio_function_output(C_OUTPUT_HOW_AXIS_X, "./static", "fnxc", der%mesh, nxc, unit_one, ierr)
  end if
  
  call profiling_out(prof)

  SAFE_DEALLOCATE_A(lrvxc)
  SAFE_DEALLOCATE_A(nxc)

  POP_SUB(xc_density_correction_calc)
end subroutine xc_density_correction_calc

! -----------------------------------------------------

FLOAT function get_qxc(mesh, nxc, density, ncutoff)  result(qxc)
  type(mesh_t), intent(in) :: mesh
  FLOAT,        intent(in) :: nxc(:)
  FLOAT,        intent(in) :: density(:)
  FLOAT,        intent(in) :: ncutoff

  integer :: ip
  FLOAT, allocatable :: nxc2(:)

  PUSH_SUB(get_qxc)

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

  POP_SUB(get_qxc)
end function get_qxc

!------------------------------------------------------------
!
! Complex scaled XC
!
!------------------------------------------------------------

!> Subroutine to stitch discontinuous values of a multiple-valued function
!! together to a single continuous, single-valued function by smoothly
!! joining at the branch cuts.
subroutine stitch(get_branch, functionvalues, startpoint)
  ! Function for getting values of multiple-valued functions.
  ! Each value of the parameter 'branch' corresponds to one such value.
  interface
    CMPLX function get_branch(x, branch)
      implicit none
      CMPLX,   intent(in) :: x
      integer, intent(in) :: branch
    end function get_branch
  end interface

  CMPLX,   intent(inout) :: functionvalues(:, :, :)
  integer, intent(in)    :: startpoint(3)
  
  integer :: i, j, imax, jmax

  imax = size(functionvalues, 1)
  jmax = size(functionvalues, 2)

  call stitchline(get_branch, functionvalues, startpoint, 1)
  do i=1, imax
    call stitchline(get_branch, functionvalues, (/i, startpoint(2), startpoint(3)/), 2)
    do j=1, jmax
      call stitchline(get_branch, functionvalues, (/i, j, startpoint(3)/), 3)
    end do
  end do
end subroutine stitch


!> Like stitch, but stitches along one line only.
subroutine stitchline(get_branch, functionvalues, startpoint, direction, startbranch)
  
  ! Function for getting values of multiple-valued functions.
  ! Each value of the parameter 'branch' corresponds to one such value.
  interface 
    CMPLX function get_branch(x, branch)
      implicit none
      CMPLX,   intent(in) :: x
      integer, intent(in) :: branch
    end function get_branch
  end interface

  CMPLX,             intent(inout) :: functionvalues(:, :, :)
  integer,           intent(in)    :: startpoint(3)
  integer, optional, intent(in)    :: direction
  integer, optional, intent(in)    :: startbranch

  integer :: stitchedpoints, direction_, startbranch_
  
  integer :: currentbranch, npts, i
  integer :: currentlocation(3)
  CMPLX :: prev_value

  PUSH_SUB(stitchline)

  stitchedpoints = 0

  currentlocation = startpoint

  if (present(direction)) then
    direction_ = direction
  else
    direction_ = 1
  end if

  if (present(startbranch)) then
    startbranch_ = startbranch
  else
    startbranch_ = 0
  end if

  npts = size(functionvalues, direction_)

  ! First loop forwards from zero and stitch along the way
  currentbranch = startbranch_
  prev_value = functionvalues(startpoint(1), startpoint(2), startpoint(3))
  do i=startpoint(direction_) + 1, npts
    call stitch_single_point()
  end do
  
  ! Now loop backwards
  currentbranch = startbranch_
  prev_value = functionvalues(startpoint(1), startpoint(2), startpoint(3))
  do i=startpoint(direction_) - 1, 1, -1
    call stitch_single_point()
  end do
  
  POP_SUB(stitchline)
contains

  !recursive 
  subroutine stitch_single_point()
    CMPLX   :: v1, v2, v3, v
    integer :: adj
    
    stitchedpoints = stitchedpoints + 1

    PUSH_SUB(newstitch.stitch_single_point)
    
    currentlocation(direction_) = i

    v1 = get_branch(functionvalues(currentlocation(1), currentlocation(2), currentlocation(3)), currentbranch)
    v2 = get_branch(functionvalues(currentlocation(1), currentlocation(2), currentlocation(3)), currentbranch - 1)
    v3 = get_branch(functionvalues(currentlocation(1), currentlocation(2), currentlocation(3)), currentbranch + 1)
    
    adj = 0
    v = v1
    if (abs(v2 - prev_value) < abs(v - prev_value)) then
      v = v2
      adj = -1
    end if
    if (abs(v3 - prev_value) < abs(v - prev_value)) then
      v = v3
      adj = +1
    end if
    currentbranch = currentbranch + adj
    functionvalues(currentlocation(1), currentlocation(2), currentlocation(3)) = v
    prev_value = v

    POP_SUB(newstitch.stitch_single_point)
  end subroutine stitch_single_point
  
end subroutine stitchline

! For evaluating values of multiple-valued functions when one value,
! e.g. the principal value, is known.  Used to stitch
CMPLX function get_root2_branch(x, branch) result(y)
  CMPLX, intent(in)   :: x
  integer, intent(in) :: branch
  
  y = x * exp(branch * M_zI * M_PI)
end function get_root2_branch

CMPLX function get_root3_branch(x, branch) result(y)
  CMPLX,   intent(in) :: x
  integer, intent(in) :: branch
  
  y = x * exp(branch * M_TWO * M_zI * M_PI / M_THREE)
end function get_root3_branch

CMPLX function get_root6_branch(x, branch) result(y)
  CMPLX,   intent(in) :: x
  integer, intent(in) :: branch
  
  y = x * exp(branch * M_zI * M_PI / M_THREE)
end function get_root6_branch

CMPLX function get_logarithm_branch(x, branch) result(y)
  CMPLX,   intent(in) :: x
  integer, intent(in) :: branch
  
  y = x + branch * M_TWO * M_zI * M_PI
end function get_logarithm_branch


subroutine zxc_complex_lda_gamma(mesh, rootrs, epsc, depsdrs, gamma, alpha1, beta1, beta2, beta3, beta4)
  type(mesh_t), intent(in)  :: mesh
  CMPLX,        intent(in)  :: rootrs(:)
  CMPLX,        intent(out) :: epsc(:)
  CMPLX,        intent(out) :: depsdrs(:)
  FLOAT,        intent(in)  :: gamma
  FLOAT,        intent(in)  :: alpha1
  FLOAT,        intent(in)  :: beta1
  FLOAT,        intent(in)  :: beta2
  FLOAT,        intent(in)  :: beta3
  FLOAT,        intent(in)  :: beta4
  
  CMPLX, allocatable :: Q0(:), Q1(:), dQ1drs(:)

  PUSH_SUB(zxc_complex_lda_gamma)
  
  SAFE_ALLOCATE(Q0(1:mesh%np))
  SAFE_ALLOCATE(Q1(1:mesh%np))
  SAFE_ALLOCATE(dQ1drs(1:mesh%np))

  Q0(:) = -M_TWO * gamma * (M_ONE + alpha1 * rootrs(:)**2)
  Q1(:) = M_TWO * gamma * rootrs(:) * (beta1 + rootrs(:) * (beta2 + rootrs(:) * (beta3 + rootrs(:) * beta4)))
  dQ1drs(:) = gamma * (beta1 / rootrs(:) + M_TWO * beta2 + rootrs(:) * (M_THREE * beta3 + M_FOUR * beta4 * rootrs(:)))

  epsc(:) = log(M_ONE + M_ONE / Q1(:))
  call localstitch(mesh, epsc, get_logarithm_branch)
  epsc(:) = Q0(:) * epsc(:)

  depsdrs(:) = -M_TWO * gamma * alpha1 * epsc(:) / Q0(:) &
    - Q0(:) * dQ1drs(:) / (Q1(:) * (Q1(:) + M_ONE))

  SAFE_DEALLOCATE_A(Q0)
  SAFE_DEALLOCATE_A(Q1)
  SAFE_DEALLOCATE_A(dQ1drs)
  
  POP_SUB(zxc_complex_lda_gamma)
end subroutine zxc_complex_lda_gamma

subroutine localstitch(mesh, array, get_branch)
  type(mesh_t), intent(in)    :: mesh
  CMPLX,        intent(inout) :: array(:)

  interface 
    CMPLX function get_branch(x, branch)
      implicit none
      CMPLX,   intent(in) :: x
      integer, intent(in) :: branch
    end function get_branch
  end interface
  
  type(cube_t)          :: cube
  type(cube_function_t) :: cf
  CMPLX, allocatable :: stitchbuffer(:, :, :)

  PUSH_SUB(localstitch)

  ! XXX this will use mpi_world
  call cube_init(cube, mesh%idx%ll, mesh%sb)
  call cube_function_null(cf)
  call zcube_function_alloc_rs(cube, cf)
  call zmesh_to_cube(mesh, array, cube, cf, .true.)
  
  SAFE_ALLOCATE(stitchbuffer(1:cube%rs_n_global(1), 1:cube%rs_n_global(2), 1:cube%rs_n_global(3)))
  stitchbuffer(:, :, :) = cf%zRS(:, :, :)
  call stitch(get_branch, stitchbuffer, cube%center)
  cf%zRS(:, :, :) = stitchbuffer(:, :, :)
  SAFE_DEALLOCATE_A(stitchbuffer)
  
  call zcube_to_mesh(cube, cf, mesh, array, .true.)
  call zcube_function_free_rs(cube, cf)
  call cube_end(cube)

  POP_SUB(localstitch)
end subroutine 


subroutine zxc_lda_exchange(mesh, zrho, zvx, zex, theta)
  type(mesh_t), intent(in)  :: mesh
  CMPLX,        intent(in)  :: zrho(:)
  CMPLX,        intent(out) :: zvx(:)
  CMPLX,        intent(out) :: zex
  FLOAT,        intent(in)  :: theta
  
  FLOAT, parameter :: exchangeprefactor = CNST(-0.98474502184269641)

  PUSH_SUB(zxc_lda_exchange)
  zvx(:) = exchangeprefactor * (zrho(:) * exp(-mesh%sb%dim * M_zI * theta))**(M_ONE / M_THREE)
  call localstitch(mesh, zvx, get_root3_branch)
  zex = M_THREE / M_FOUR * sum(zvx(:) * zrho(:)) * mesh%volume_element
  POP_SUB(zxc_lda_exchange)
end subroutine zxc_lda_exchange

subroutine zxc_lda_correlation(mesh, zrho, zvc, zec, theta, polarized)
  type(mesh_t), intent(in)  :: mesh
  CMPLX,        intent(in)  :: zrho(:, :)
  CMPLX,        intent(out) :: zvc(:, :)
  CMPLX,        intent(out) :: zec
  FLOAT,        intent(in)  :: theta
  logical,      intent(in)  :: polarized

  CMPLX, allocatable    :: rootrs(:), epsc(:), depsdrs(:), epsc1(:), zrhototal(:), &
    depsdrs1(:), alpha(:), dalphadrs(:), zeta(:), xplus(:), xminus(:)
  
  FLOAT, parameter :: C0I = CNST(0.238732414637843), C1 = CNST(-0.45816529328314287), &
    CC1 = CNST(1.9236610509315362), CC2 = CNST(2.5648814012420482), IF2 = CNST(0.58482236226346462)
  integer          :: gg
  FLOAT            :: ff, f1, xx, zeta3, zeta4, decdrs, decdzeta
  
  PUSH_SUB(zxc_lda_correlation)

  SAFE_ALLOCATE(rootrs(1:mesh%np))
  SAFE_ALLOCATE(zrhototal(1:mesh%np))
  SAFE_ALLOCATE(epsc(1:mesh%np))
  SAFE_ALLOCATE(depsdrs(1:mesh%np))

  zrhototal(:) = sum(zrho, 2)
  rootrs(:) = (M_THREE / (CNST(1e-30) + M_FOUR * M_PI * exp(-mesh%sb%dim * M_zI * theta) * zrhototal(:)))**(M_ONE / M_SIX)
  call localstitch(mesh, rootrs, get_root6_branch)
  
  call zxc_complex_lda_gamma(mesh, rootrs, epsc, depsdrs, &
    CNST(0.031091), CNST(0.21370), CNST(7.5957), CNST(3.5876), CNST(1.6382), CNST(0.49294))

  if(polarized) then
    SAFE_ALLOCATE(zeta(1:mesh%np))
    SAFE_ALLOCATE(epsc1(1:mesh%np))
    SAFE_ALLOCATE(depsdrs1(1:mesh%np))
    SAFE_ALLOCATE(alpha(1:mesh%np))
    SAFE_ALLOCATE(dalphadrs(1:mesh%np))
    SAFE_ALLOCATE(xplus(1:mesh%np))
    SAFE_ALLOCATE(xminus(1:mesh%np))

    zeta(:) = (zrho(:, 1) - zrho(:, 2)) / (CNST(1e-30) + zrhototal(:))

    call zxc_complex_lda_gamma(mesh, rootrs, epsc1, depsdrs1, &
      CNST(0.015545), CNST(0.20548), CNST(14.1189), CNST(6.1977), CNST(3.3662), CNST(0.62517))
    call zxc_complex_lda_gamma(mesh, rootrs, alpha, dalphadrs, &
      CNST(0.016887), CNST(0.11125), CNST(10.357), CNST(3.6231), CNST(0.88026), CNST(0.49671))
    alpha(:) = -alpha(:)
    dalphadrs(:) = -dalphadrs(:)
    
    xplus(:) = (M_ONE + zeta(:))**(M_ONE / M_THREE)
    xminus(:) = (M_ONE - zeta(:))**(M_ONE / M_THREE)
    call localstitch(mesh, xplus, get_root3_branch)
    call localstitch(mesh, xminus, get_root3_branch)
    
    do gg=1, mesh%np
      zeta3 = zeta(gg)**3
      zeta4 = zeta3 * zeta(gg)

      ff = CC1 * ((M_ONE + zeta(gg)) * xplus(gg) + (M_ONE - zeta(gg)) * xminus(gg) - M_TWO)
      f1 = CC2 * (xplus(gg) - xminus(gg))
      xx = M_ONE - zeta4
      decdrs = depsdrs(gg) * (M_ONE - ff * zeta4) + depsdrs(gg) * ff * zeta4 + dalphadrs(gg) * ff * xx * IF2
      
      decdzeta = M_FOUR * zeta3 * ff * (epsc1(gg) - epsc(gg) - alpha(gg) * IF2) + &
        f1 * (zeta4 * epsc1(gg) - zeta4 * epsc(gg) + xx * alpha(gg) * IF2)
      
      zvc(gg, 1) = epsc(gg) - rootrs(gg)**2 * depsdrs(gg) / M_THREE - (zeta(gg) - M_ONE) * decdzeta
      zvc(gg, 2) = epsc(gg) - rootrs(gg)**2 * depsdrs(gg) / M_THREE - (zeta(gg) + M_ONE) * decdzeta
    end do
    zec = sum(epsc(:) * zrhototal(:)) * mesh%volume_element
    
    SAFE_DEALLOCATE_A(epsc1)
    SAFE_DEALLOCATE_A(depsdrs1)
    SAFE_DEALLOCATE_A(alpha)
    SAFE_DEALLOCATE_A(dalphadrs)
    SAFE_DEALLOCATE_A(xplus)
    SAFE_DEALLOCATE_A(xminus)
    SAFE_DEALLOCATE_A(zeta)
  else
    zvc(:, 1) = epsc(:) - rootrs(:)**2 * depsdrs(:) / M_THREE
    zec = sum(epsc(:) * zrho(:, 1)) * mesh%volume_element
  end if

  SAFE_DEALLOCATE_A(epsc)
  SAFE_DEALLOCATE_A(depsdrs)
  SAFE_DEALLOCATE_A(rootrs)
  SAFE_DEALLOCATE_A(zrhototal)

  POP_SUB(zxc_lda_correlation)
end subroutine zxc_lda_correlation

subroutine zxc_complex_lda(mesh, ispin, theta, zrho, zvx, zvc, zex, zec)
  type(mesh_t),    intent(in)    :: mesh
  integer,         intent(in)    :: ispin
  FLOAT,           intent(in)    :: theta
  CMPLX,           intent(in)    :: zrho(:, :)
  CMPLX,           intent(inout) :: zvx(:, :)
  CMPLX,           intent(inout) :: zvc(:, :)
  CMPLX,           intent(inout) :: zex
  CMPLX,           intent(inout) :: zec

  ! LDA correlation parameters
  integer               :: N
  logical               :: calc_energy
  CMPLX                 :: zenergies1(2), zenergies2(2)
  CMPLX                 :: zextmp

  PUSH_SUB(zxc_complex_lda)
  
  if(ispin == UNPOLARIZED) then
    call zxc_lda_exchange(mesh, zrho(:, 1), zvx(:, 1), zex, theta)
  else
    call zxc_lda_exchange(mesh, M_TWO * zrho(:, 1), zvx(:, 1), zex, theta)
    call zxc_lda_exchange(mesh, M_TWO * zrho(:, 2), zvx(:, 2), zextmp, theta)
    zex = M_HALF * (zex + zextmp)
  end if

  call zxc_lda_correlation(mesh, zrho, zvc, zec, theta, ispin == SPIN_POLARIZED)

  call cmplxscl_energy_comm_sum(zex, zec)
  POP_SUB(zxc_complex_lda)
  
end subroutine zxc_complex_lda


subroutine cmplxscl_energy_comm_sum(zex, zec)
  CMPLX, intent(inout) :: zex
  CMPLX, intent(inout) :: zec
  
  CMPLX                :: zenergies1(2), zenergies2(2)

  PUSH_SUB(cmplxscl_energy_comm_sum)

#ifdef HAVE_MPI
  zenergies1(1) = zex
  zenergies1(2) = zec
  call MPI_Allreduce(zenergies1, zenergies2, 2, MPI_CMPLX, MPI_SUM, mpi_world%comm, mpi_err) ! XXX world
  zex = zenergies2(1)
  zec = zenergies2(2)
#endif

  POP_SUB(cmplxscl_energy_comm_sum)
end subroutine cmplxscl_energy_comm_sum

subroutine cmplxscl_add_pbe_div_term(der, mesh, theta, depsdsigma, gradient, zvxc)
  type(derivatives_t), intent(in) :: der
  type(mesh_t),        intent(in) :: mesh
  FLOAT,               intent(in) :: theta
  CMPLX,               intent(in) :: depsdsigma(:)
  CMPLX,               intent(in) :: gradient(:, :)
  CMPLX,            intent(inout) :: zvxc(:, :)
  
  integer :: nn
  CMPLX, allocatable :: divbuffer(:, :), divoutput(:)

  PUSH_SUB(cmplxscl_add_pbe_div_term)
  
  SAFE_ALLOCATE(divbuffer(1:mesh%np_part, 1:mesh%sb%dim))
  divbuffer(mesh%np:mesh%np_part, :) = M_ZERO
  do nn=1, mesh%sb%dim
    divbuffer(1:mesh%np, nn) = depsdsigma(1:mesh%np) * gradient(1:mesh%np, nn)
  end do

  SAFE_ALLOCATE(divoutput(1:mesh%np))
  call zderivatives_div(der, divbuffer, divoutput, .true., .true.)
  SAFE_DEALLOCATE_A(divbuffer)
  
  zvxc(:, 1) = zvxc(:, 1) - M_TWO * divoutput(:) * exp(-M_zI * theta)
  SAFE_DEALLOCATE_A(divoutput)
  
  POP_SUB(cmplxscl_add_pbe_div_term)
end subroutine cmplxscl_add_pbe_div_term

subroutine zxc_complex_lb94(der, mesh, ispin, theta, zrho, zvx, zvc, zex, zec)
  type(derivatives_t),  intent(in)    :: der
  type(mesh_t),         intent(in)    :: mesh
  integer,              intent(in)    :: ispin
  FLOAT,                intent(in)    :: theta
  CMPLX,                intent(in)    :: zrho(:, :)
  CMPLX,                intent(inout) :: zvx(:, :)
  CMPLX,                intent(inout) :: zvc(:, :)
  CMPLX,                intent(inout) :: zex
  CMPLX,                intent(inout) :: zec
  
  CMPLX, allocatable :: gradient(:, :), ncuberoot(:), xbuf(:), x2plus1buf(:), fbuf(:)
  FLOAT, parameter   :: beta = CNST(0.05), xi = CNST(1.2599210498948732) ! 2^(1/3)

  integer :: nn

  PUSH_SUB(zxc_complex_lb94)
  ASSERT(ispin == UNPOLARIZED)
  
  call zxc_complex_lda(mesh, ispin, theta, zrho, zvx, zvc, zex, zec)
  
  SAFE_ALLOCATE(gradient(1:mesh%np, 1:mesh%sb%dim))
  call cmplxscl_get_gradient(der, zrho, theta, gradient)

  SAFE_ALLOCATE(xbuf(1:mesh%np))

  xbuf(:) = sqrt(sum(gradient**2, 2))
  SAFE_DEALLOCATE_A(gradient)
  
  ! Most horrible and awful!  Unstitchable when gradient is 0 at the origin.  But luckily we
  ! will multiply this by some logarithmic expression which may be equally stomach-turning,
  ! possibly cancelling out this twisted perversion of mathematics.
  call localstitch(mesh, xbuf, get_root2_branch) ! Yuck !!!

  SAFE_ALLOCATE(ncuberoot(1:mesh%np))
  ncuberoot(:) = (zrho(:, 1) * exp(-M_zI * theta * mesh%sb%dim))**(M_ONE / M_THREE)
  call localstitch(mesh, ncuberoot, get_root3_branch)
  xbuf(:) = xbuf(:) / ncuberoot(:)**4 ! We could also just do n * ncuberoot of course...

  do nn=1, mesh%np
    if(abs(xbuf(nn)) > CNST(1e4)) then
      xbuf(nn) = xbuf(nn) / abs(xbuf(nn)) ! XXX ?
    end if
  end do
  
  SAFE_ALLOCATE(x2plus1buf(1:mesh%np))
  x2plus1buf(:) = sqrt(xbuf(:)**2 + M_ONE)
  call localstitch(mesh, x2plus1buf, get_root2_branch) ! Less likely to cause trouble
  
  SAFE_ALLOCATE(fbuf(1:mesh%np))
  fbuf(:) = log(xi * (xbuf(:) + x2plus1buf(:)))
  SAFE_DEALLOCATE_A(x2plus1buf)
  call localstitch(mesh, fbuf, get_logarithm_branch)
  fbuf(:) = -beta * xi * xbuf(:)**2 / (M_ONE + M_THREE * beta * xi * xbuf(:) * fbuf(:))
  SAFE_DEALLOCATE_A(xbuf)

  fbuf(:) = fbuf(:) * ncuberoot(:)
  zvx(:, 1) = zvx(:, 1) + fbuf(:)

  SAFE_DEALLOCATE_A(fbuf)
  SAFE_DEALLOCATE_A(ncuberoot)

  ! We would have to sum up the energy over nodes, except with this functional we
  ! do not adjust the energy from the LDA value which was already summed.

  POP_SUB(zxc_complex_lb94)
end subroutine zxc_complex_lb94

subroutine cmplxscl_get_gradient(der, zrho, theta, gradient)
  type(derivatives_t),  intent(in) :: der
  CMPLX,                intent(in) :: zrho(:, :)
  FLOAT,                intent(in) :: theta
  CMPLX,               intent(out) :: gradient(:, :)
  
  CMPLX, allocatable :: zrhotemp(:)
  CMPLX              :: dimphase

  dimphase = exp(M_zI * theta * der%mesh%sb%dim)

  PUSH_SUB(cmplxscl_get_gradient)

  ASSERT(size(zrho, 2) == 1) ! no spin polarization for now
  
  SAFE_ALLOCATE(zrhotemp(1:der%mesh%np_part))
  ! Divide by exp(i dim theta) to get Jacobian right and further by exp(i theta) due to derivative
  zrhotemp(1:der%mesh%np) = zrho(:, 1) * exp(-M_zI * theta * (M_ONE + der%mesh%sb%dim))
  zrhotemp(der%mesh%np + 1:der%mesh%np_part) = M_ZERO
  call zderivatives_grad(der, zrhotemp, gradient)
  SAFE_DEALLOCATE_A(zrhotemp)
  POP_SUB(cmplxscl_get_gradient)
end subroutine cmplxscl_get_gradient


subroutine zxc_complex_pbe(der, mesh, ispin, theta, zrho, zvx, zvc, zex, zec)
  type(derivatives_t),  intent(in)    :: der
  type(mesh_t),         intent(in)    :: mesh
  integer,              intent(in)    :: ispin
  FLOAT,                intent(in)    :: theta
  CMPLX,                intent(in)    :: zrho(:, :)
  CMPLX,                intent(inout) :: zvx(:, :)
  CMPLX,                intent(inout) :: zvc(:, :)
  CMPLX,                intent(inout) :: zex
  CMPLX,                intent(inout) :: zec

  integer               :: nn
  logical               :: calc_energy
  CMPLX                 :: zenergies1(2), zenergies2(2) ! XXX sum/broadcast
  CMPLX                 :: zextmp, cc, dFxds2, dexda2

  ! XXX Some of the float parameters are repeated in LDA.  Define globally or something?
  FLOAT, parameter      :: &
    mu = CNST(0.2195149727645171), kappa = CNST(0.804), C0I = CNST(0.238732414637843), &
    C1 = CNST(-0.45816529328314287), C2 = CNST(0.26053088059892404), C3 = CNST(0.10231023756535741), &
    gamma = CNST(0.0310906908697), CC1 = CNST(1.9236610509315362), CC2 = CNST(2.5648814012420482), &
    IF2 = CNST(0.58482236226346462), beta = CNST(0.06672455060314922)

  CMPLX, allocatable    :: rootrs(:), gradient(:, :), a2(:), epsc(:), decdrs(:), zrhotemp(:), HH(:), depscdsigma(:), &
    depsxdsigma(:)
  CMPLX                 :: rs, ex, dexdrs, s2, xx, yy, Fx, t2, At2, num, denom, AA, dimphase
  CMPLX                 :: tmp, tmp2, dAdrs, dHdA, dHdt2

  PUSH_SUB(zxc_complex_pbe)

  ASSERT(ispin == UNPOLARIZED)

  dimphase = exp(mesh%sb%dim * M_zI * theta)

  SAFE_ALLOCATE(rootrs(1:mesh%np))
  rootrs(:) = (M_THREE / (M_FOUR * M_PI / dimphase * zrho(:, 1)))**(M_ONE / M_SIX)
  call localstitch(mesh, rootrs, get_root6_branch)

  SAFE_ALLOCATE(gradient(1:mesh%np, mesh%sb%dim))

  call cmplxscl_get_gradient(der, zrho, theta, gradient)

  SAFE_ALLOCATE(a2(1:mesh%np))
  a2(:) = sum(gradient**2, 2)
  SAFE_ALLOCATE(epsc(1:mesh%np))
  SAFE_ALLOCATE(decdrs(1:mesh%np))

  call zxc_complex_lda_gamma(mesh, rootrs, epsc, decdrs, gamma, &
    CNST(0.21370), CNST(7.5957), CNST(3.5876), CNST(1.6382), CNST(0.49294))

  SAFE_ALLOCATE(depscdsigma(1:mesh%np_part))
  SAFE_ALLOCATE(depsxdsigma(1:mesh%np_part))
  depscdsigma(mesh%np:mesh%np_part) = M_ZERO
  depsxdsigma(mesh%np:mesh%np_part) = M_ZERO
  
  SAFE_ALLOCATE(HH(1:mesh%np))
  do nn=1, mesh%np
    rs = rootrs(nn)**2
    ex = C1 / rs
    dexdrs = -ex / rs
    cc = (C2 * rs / (zrho(nn, 1) / dimphase))**2
    s2 = a2(nn) * cc
    xx = M_ONE + mu * s2 / kappa
    Fx = M_ONE + kappa - kappa / xx
    
    dFxds2 = mu / xx**2

    dexdrs = dexdrs * Fx + ex * dFxds2 * M_EIGHT * cc * a2(nn) / rs
    dexda2 = ex * dFxds2 * cc
    ex = ex * Fx
    
    zvx(nn, 1) = ex - rs * dexdrs / M_THREE
    zex = zex + ex * zrho(nn, 1)

    depsxdsigma(nn) = zrho(nn, 1) * dexda2 / dimphase

    t2 = C3 * a2(nn) * rs / (zrho(nn, 1) / dimphase)**2
    yy = -epsc(nn) / gamma
    xx = exp(yy)
    if (yy.eq.M_ZERO) then
      AA = M_ZERO
    else
      AA = beta / (gamma * (xx - M_ONE))
    end if
    At2 =  AA * t2
    num = M_ONE + At2
    denom = num + At2**2
    HH(nn) = log(M_ONE + beta * t2 * num / (denom * gamma))
    if(yy.eq.M_ZERO) then
      HH(nn) = M_ZERO
    end if
        
    tmp = gamma * beta / (denom * (beta * t2 * num + gamma * denom))
    tmp2 = AA * AA * xx / beta
    dAdrs = tmp2 * decdrs(nn)
    dHdA = -At2 * t2 * t2 * (M_TWO + At2) * tmp
    dHdt2 = (M_ONE + M_TWO * At2) * tmp
    decdrs(nn) = decdrs(nn) + dHdt2 * M_SEVEN * t2 / rs + dHdA * dAdrs
    
    depscdsigma(nn) = dHdt2 * C3 * rs / (zrho(nn, 1) / dimphase) ! zrho ?
  end do


  call localstitch(mesh, HH, get_logarithm_branch)
  HH(:) = HH(:) * gamma
  
  epsc(:) = epsc(:) + HH(:)

  zvc(:, :) = M_ZERO
  call cmplxscl_add_pbe_div_term(der, mesh, theta, depsxdsigma, gradient, zvx)
  call cmplxscl_add_pbe_div_term(der, mesh, theta, depscdsigma, gradient, zvc)
  SAFE_DEALLOCATE_A(depscdsigma)
  SAFE_DEALLOCATE_A(depsxdsigma)
  SAFE_DEALLOCATE_A(gradient)

  zex = zex * mesh%volume_element
  zec = sum(epsc(:) * zrho(:, 1)) * mesh%volume_element
  zvc(:, 1) = zvc(:, 1) + epsc(:) - rootrs(:)**2 * decdrs(:) / M_THREE

  call cmplxscl_energy_comm_sum(zex, zec)

  SAFE_DEALLOCATE_A(HH)
  SAFE_DEALLOCATE_A(epsc)
  SAFE_DEALLOCATE_A(decdrs)

  SAFE_DEALLOCATE_A(a2)
  SAFE_DEALLOCATE_A(rootrs)
  POP_SUB(zxc_complex_pbe)  
end subroutine zxc_complex_pbe

! ----------------------------------------------------------------------------- 
!> This is the complex scaled interface for xc functionals.
!! It will eventually be merged with the other one dxc_get_vxc after some test
!! -----------------------------------------------------------------------------
subroutine xc_get_vxc_cmplx(der, xcs, ispin, rho, Imrho, vxc, Imvxc, theta, ex, ec, Imex, Imec)
  type(derivatives_t),  intent(in)    :: der             !< Discretization and the derivative operators and details
  type(xc_t), target,   intent(in)    :: xcs             !< Details about the xc functional used
  integer,              intent(in)    :: ispin           !< Number of spin channels 
  FLOAT,                intent(in)    :: rho(:, :)       !< Electronic density 
  FLOAT,                intent(in)    :: Imrho(:, :)     !< cmplxscl: Electronic density 
  FLOAT,                intent(inout) :: vxc(:,:)        !< XC potential
  FLOAT,                intent(inout) :: Imvxc(:,:)      !< cmplxscl: XC potential
  FLOAT,                intent(in)    :: theta           !< complex scaling angle
  FLOAT, optional,      intent(inout) :: ex              !< Exchange 1energy.
  FLOAT, optional,      intent(inout) :: ec              !< Correlation energy.
  FLOAT, optional,      intent(inout) :: Imex            !< cmplxscl: Exchange energy.
  FLOAT, optional,      intent(inout) :: Imec            !< cmplxscl: Correlation energy


  CMPLX, allocatable         :: zrho(:, :), zvx(:, :), zvc(:, :)
  CMPLX                      :: zex, zec
  type(xc_functl_t), pointer :: functl(:)
  logical                    :: calc_energy

  PUSH_SUB(xc_get_vxc_cmplx)

  ASSERT(present(ex) .eqv. present(ec))
  ASSERT(present(ex) .eqv. present(Imex))
  ASSERT(present(ec) .eqv. present(Imec))
  calc_energy = present(ex)
  
  SAFE_ALLOCATE(zrho(1:size(vxc, 1), 1:size(vxc, 2)))
  SAFE_ALLOCATE(zvx(1:size(vxc, 1), 1:size(vxc, 2)))
  SAFE_ALLOCATE(zvc(1:size(vxc, 1), 1:size(vxc, 2)))

  zrho(:, :) = rho(:, :) + M_zI * Imrho(:, :) + CNST(1e-40)
  zvx(:, :) = M_ZERO
  zvc(:, :) = M_ZERO
  zex = M_ZERO
  zec = M_ZERO

  !Pointer-shortcut for xcs%functional
  !It helps to remember that for xcs%functional(:,:)
  ! (1,:) => exchange,    (2,:) => correlation
  ! (:,1) => unpolarized, (:,2) => polarized
  if(ispin == UNPOLARIZED) then
    functl => xcs%functional(:, 1)
  else
    functl => xcs%functional(:, 2)
  end if
  
  if(functl(1)%id == XC_LDA_XC_CMPLX) then
      ! XXX vxc optional parameter?
    call zxc_complex_lda(der%mesh, ispin, theta, zrho, zvx, zvc, zex, zec)
  else if(functl(1)%id == XC_PBE_XC_CMPLX) then
    call zxc_complex_pbe(der, der%mesh, ispin, theta, zrho, zvx, zvc, zex, zec)
  else if(functl(1)%id == XC_LB94_XC_CMPLX) then
    call zxc_complex_lb94(der, der%mesh, ispin, theta, zrho, zvx, zvc, zex, zec)
  else if(functl(1)%id == XC_HALF_HARTREE) then
    ! Exact exchange for 2 particles [vxc(r) = 1/2 * vh(r)]
    ! we keep it here for debug purposes
    call zpoisson_solve(psolver, zvx(:, 1), zrho(:, 1))
    zvx = -M_HALF * zvx  
    if(calc_energy) then
      zex = M_HALF * zmf_dotp(der%mesh, zrho(:, 1), zvx(:, 1), dotu = .true. )
    end if
  else if(functl(1)%family == XC_FAMILY_NONE) then
    ! Everything was initialized to zero, so do nothing.
    ! This block is just here as a reminder.
  else  
    write(message(1), '(a)') 'The selected XCFunctional will not work with ComplexScaling = yes.'
    write(message(2), '(a)') 'Use XCFunctional = lda_xc_cmplx, for example.'
    call messages_fatal(2)     
  end if

  vxc(:, :) = real(zvx(:, :)) + real(zvc(:, :))
  Imvxc(:, :) = aimag(zvx(:, :)) + aimag(zvc(:, :))

  if(calc_energy) then
    ex = real(zex)
    Imex = aimag(zex)
    ec = real(zec)
    Imec = aimag(zec)
  end if
  
  SAFE_DEALLOCATE_A(zrho)
  SAFE_DEALLOCATE_A(zvx)
  SAFE_DEALLOCATE_A(zvc)

  POP_SUB(xc_get_vxc_cmplx)
end subroutine xc_get_vxc_cmplx


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

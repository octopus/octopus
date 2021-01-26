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

subroutine xc_get_vxc(der, xcs, st, psolver, namespace, rho, ispin, vxc, ex, ec, deltaxc, vtau, ex_density, ec_density)
  type(derivatives_t),    intent(in)    :: der             !< Discretization and the derivative operators and details
  type(xc_t), target,     intent(inout) :: xcs             !< Details about the xc functional used
  type(states_elec_t),    intent(in)    :: st              !< State of the system (wavefunction,eigenvalues...)
  type(poisson_t),        intent(in)    :: psolver
  type(namespace_t),      intent(in)    :: namespace
  FLOAT,                  intent(in)    :: rho(:, :)       !< Electronic density 
  integer,                intent(in)    :: ispin           !< Number of spin channels 
  FLOAT, optional,        intent(inout) :: vxc(:,:)        !< XC potential
  FLOAT, optional,        intent(inout) :: ex              !< Exchange energy.
  FLOAT, optional,        intent(inout) :: ec              !< Correlation energy.
  FLOAT, optional,        intent(inout) :: deltaxc         !< The XC derivative discontinuity
  FLOAT, optional,        intent(inout) :: vtau(:,:)       !< Derivative wrt (two times kinetic energy density)
  FLOAT, optional, target, intent(out)  :: ex_density(:)   !< The exchange energy density
  FLOAT, optional, target, intent(out)  :: ec_density(:)   !< The correlation energy density

  integer, parameter :: N_BLOCK_MAX = 10000
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
  FLOAT, pointer :: ex_per_vol(:)  ! Exchange energy per unit volume 
  FLOAT, pointer :: ec_per_vol(:)  ! Correlation energy per unit volume 
  !  First order (functional) derivatives
  FLOAT, allocatable :: dedd(:,:)      ! Derivative of the exchange or correlation energy wrt the density
  FLOAT, allocatable :: dedgd(:,:,:)   ! Derivative of the exchange or correlation energy wrt the gradient of the density
  FLOAT, allocatable :: dedldens(:,:)  ! Derivative of the exchange or correlation energy wrt the laplacian of the density

  FLOAT, allocatable :: vx(:)
  FLOAT, allocatable :: unp_dens(:), unp_dedd(:)

  integer :: ib, ip, isp, families, ixc, spin_channels, is, idir, ipstart, ipend
  FLOAT   :: energy(1:2)
  logical :: gga, mgga, mgga_withexc, libvdwxc
  type(profile_t), save :: prof, prof_libxc
  logical :: calc_energy
  type(xc_functl_t), pointer :: functl(:)
  type(distributed_t) :: distribution
  type(profile_t), save :: prof_gather
  
  PUSH_SUB(xc_get_vxc)
  call profiling_in(prof, TOSTRING(XC_LOCAL))

  nullify(ex_per_vol)
  nullify(ec_per_vol)
  
  ASSERT(present(ex) .eqv. present(ec))
  calc_energy = present(ex) .or. present(ex_density) .or. present(ec_density)

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
  families = XC_FAMILY_LDA + XC_FAMILY_GGA + XC_FAMILY_HYB_GGA + XC_FAMILY_MGGA + XC_FAMILY_HYB_MGGA + XC_FAMILY_LIBVDWXC
  if(bitand(xcs%family, families) == 0) then
    POP_SUB(xc_get_vxc)
    call profiling_out(prof)
    return
  end if

  do ixc = 1, 2
    if(functl(ixc)%family /= XC_FAMILY_NONE .and. bitand(functl(ixc)%family, XC_FAMILY_OEP) == 0) then
      ASSERT(bitand(functl(ixc)%flags, XC_FLAGS_HAVE_VXC) /= 0)
    end if
  end do

  ! initialize a couple of handy variables
  gga  = family_is_gga(xcs%family)
  mgga = family_is_mgga(xcs%family)
  mgga_withexc = family_is_mgga_with_exc(xcs)
  if(mgga) then
    ASSERT(gga)
  end if

  ! libvdwxc counts as a GGA because it needs the density gradient.
  ! However it only calls libxc for LDA correlation and therefore
  ! never initializes l_vsigma in the libxc space/functional loop.
  ! Thus, it must never use l_vsigma!!  libvdwxc adds its own gradient
  ! corrections from the nonlocal part afterwards.
  libvdwxc = bitand(xcs%family, XC_FAMILY_LIBVDWXC) /= 0

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
  ! avoid calling the subroutine states_elec_calc_quantities twice
  if((gga .and. (.not. mgga)) .or. xcs%xc_density_correction == LR_X) then
    ! get gradient of the density (this is faster than calling states_elec_calc_quantities)
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
    if (allocated(st%rho_core)) then
      call messages_not_implemented("MGGA with nonlinear core correction", namespace=namespace)
    end if

    if (xcs%use_gi_ked) then
      call states_elec_calc_quantities(der, st, .true., gi_kinetic_energy_density = tau, &
                                  density_gradient = gdens, density_laplacian = ldens)
    else
      call states_elec_calc_quantities(der, st, .true., kinetic_energy_density = tau, &
 density_gradient = gdens, density_laplacian = ldens)
    end if

    if(functl(FUNC_X)%id == XC_MGGA_X_TB09 .and. der%mesh%sb%periodic_dim == 3) then
      call calc_tb09_c()
    end if

  end if

  if(xcs%functional(FUNC_C,1)%family == XC_FAMILY_HYB_GGA .or. &
         xcs%functional(FUNC_C,1)%family == XC_FAMILY_HYB_MGGA) then

    if (xcs%functional(FUNC_C,1)%id == XC_HYB_GGA_XC_MVORB_HSE06  &
     .or. xcs%functional(FUNC_C,1)%id == XC_HYB_GGA_XC_MVORB_PBEH) then
      call calc_mvorb_alpha()
    end if

  end if



  if(xcs%parallel) then
    call distributed_init(distribution, der%mesh%np, st%st_kpt_mpi_grp%comm)
    ipstart = distribution%start
    ipend = distribution%end
  else
    ipstart = 1
    ipend = der%mesh%np
  end if

  call local_allocate()
 
  space_loop: do ip = ipstart, ipend, N_BLOCK_MAX

    call space_loop_init(ip, ipend, n_block)

    ! Calculate the potential/gradient density in local reference frame.
    functl_loop: do ixc = FUNC_X, FUNC_C

      if(functl(ixc)%family == XC_FAMILY_NONE) cycle

      call profiling_in(prof_libxc, TOSTRING(LIBXC))

      if(calc_energy .and. bitand(functl(ixc)%flags, XC_FLAGS_HAVE_EXC) /= 0) then
        ! we get the xc energy and potential
        select case(functl(ixc)%family)
        case(XC_FAMILY_LDA, XC_FAMILY_LIBVDWXC)
          call xc_f03_lda_exc_vxc(functl(ixc)%conf, int(n_block, XC_SIZE_T), l_dens, l_zk, l_dedd)

        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          call xc_f03_gga_exc_vxc(functl(ixc)%conf, int(n_block, XC_SIZE_T), l_dens, l_sigma, l_zk, l_dedd, l_vsigma)
        case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
          call xc_f03_mgga_exc_vxc(functl(ixc)%conf, int(n_block, XC_SIZE_T), l_dens, l_sigma, l_ldens, l_tau, l_zk, l_dedd, &
            l_vsigma, l_dedldens, l_dedtau)

        case default
          call profiling_out(prof_libxc)
          cycle
        end select

      else ! we just get the potential
        l_zk(:) = M_ZERO

        select case(functl(ixc)%family)
        case(XC_FAMILY_LDA, XC_FAMILY_LIBVDWXC)
          call xc_f03_lda_vxc(functl(ixc)%conf, int(n_block, XC_SIZE_T), l_dens, l_dedd)

        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          call xc_f03_gga_vxc(functl(ixc)%conf, int(n_block, XC_SIZE_T), l_dens, l_sigma, l_dedd, l_vsigma)

        case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
          call xc_f03_mgga_vxc(functl(ixc)%conf, int(n_block, XC_SIZE_T), l_dens, l_sigma, l_ldens, l_tau, l_dedd, &
            l_vsigma, l_dedldens, l_dedtau)

        case default
          call profiling_out(prof_libxc)
          cycle
        end select
        
      end if

      call profiling_out(prof_libxc)

      if(calc_energy) then
        if(functl(ixc)%type == XC_EXCHANGE) then
          do ib = 1, n_block
            ex_per_vol(ib + ip - 1) = ex_per_vol(ib + ip - 1) + sum(l_dens(1:spin_channels, ib))*l_zk(ib)
          end do
        else
          do ib = 1, n_block
            ec_per_vol(ib + ip - 1) = ec_per_vol(ib + ip - 1) + sum(l_dens(1:spin_channels, ib))*l_zk(ib)
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
          call xc_f03_lda_vxc(xcs%functional(ixc, 1)%conf, int(n_block, XC_SIZE_T), unp_dens, unp_dedd)
          
        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
          call messages_not_implemented('XC density correction for GGA/mGGA', namespace=namespace)
          
          call xc_f03_gga_vxc(xcs%functional(ixc, 1)%conf, int(n_block, XC_SIZE_T), unp_dens, l_sigma, unp_dedd, l_vsigma)
        end select
        
        do ib = 1, n_block
          vx(ib + ip - 1) = unp_dedd(ib)
        end do
        
        ! GGA terms are missing here
        
        ! seems it would be better to allocate and deallocate these arrays outside the space loop.
        SAFE_DEALLOCATE_A(unp_dens)
        SAFE_DEALLOCATE_A(unp_dedd)
      end if

      if(family_is_gga(functl(ixc)%family).and.functl(ixc)%family /= XC_FAMILY_LIBVDWXC) then
        do ib = 1, n_block
          dedgd(ib + ip - 1,:,1) = dedgd(ib + ip - 1,:,1) + M_TWO*l_vsigma(1, ib)*gdens(ib + ip - 1,:,1)
          if(ispin /= UNPOLARIZED) then
            dedgd(ib + ip - 1,:,1) = dedgd(ib + ip - 1,:,1) + l_vsigma(2, ib)*gdens(ib + ip - 1,:,2)
            dedgd(ib + ip - 1,:,2) = dedgd(ib + ip - 1,:,2) +  &
                 M_TWO*l_vsigma(3, ib)*gdens(ib + ip - 1,:,2) + l_vsigma(2, ib)*gdens(ib + ip - 1,:,1)
          end if
        end do
      end if

      if(family_is_mgga(functl(ixc)%family)) then
        ! XXXXX does this work correctly when functionals belong to
        ! different families and only one is mgga?
        call copy_local_to_global(l_dedldens, dedldens, n_block, spin_channels, ip)
        if(bitand(functl(ixc)%flags, XC_FLAGS_HAVE_EXC) /= 0 ) &
          call copy_local_to_global(l_dedtau, vtau, n_block, spin_channels, ip)
      end if

    end do functl_loop
  end do space_loop

  call local_deallocate()

  ! calculate the energy, we do the integrals directly so when the data
  ! is fully distributed we do not need to allgather first
  if(present(ex) .or. present(ec)) then

    energy(1:2) = CNST(0.0)
    
    if(der%mesh%use_curvilinear) then
      do ip = ipstart, ipend
        energy(1) = energy(1) + ex_per_vol(ip)*der%mesh%vol_pp(ip)
        energy(2) = energy(2) + ec_per_vol(ip)*der%mesh%vol_pp(ip)
      end do
    else
      do ip = ipstart, ipend
        energy(1) = energy(1) + ex_per_vol(ip)
        energy(2) = energy(2) + ec_per_vol(ip)
      end do
    end if
    
    energy(1:2) = energy(1:2)*der%mesh%volume_element

    if(xcs%parallel) then
      call comm_allreduce(st%dom_st_kpt_mpi_grp%comm, energy)
    else if(der%mesh%parallel_in_domains) then
      call comm_allreduce(der%mesh%mpi_grp%comm, energy)
    end if

    ex = ex + energy(1)
    ec = ec + energy(2)
    
  end if

  if(xcs%parallel) then
    if(distribution%parallel) then
      call profiling_in(prof_gather, TOSTRING(XC_GATHER))
      
      do isp = 1, spin_channels
        call distributed_allgather(distribution, dedd(:, isp))
      end do
      
      if(xcs%xc_density_correction == LR_X) then
        call distributed_allgather(distribution, vx)
      end if
      
      if(gga .or. mgga) then
        do idir = 1, der%mesh%sb%dim
          do isp = 1, spin_channels
            call distributed_allgather(distribution, dedgd(:, idir, isp))
          end do
        end do
      end if
      
      if(mgga) then
        do isp = 1, spin_channels
          call distributed_allgather(distribution, dedldens(:, isp))
          if(mgga_withexc) &
          call distributed_allgather(distribution, vtau(:, isp))
        end do
      end if
      
      call profiling_out(prof_gather)
    end if
  
    call distributed_end(distribution)
  end if

  if(functl(FUNC_C)%family == XC_FAMILY_LIBVDWXC) then
    functl(FUNC_C)%libvdwxc%energy = M_ZERO
    call libvdwxc_calculate(functl(FUNC_C)%libvdwxc, namespace, dens, gdens, dedd, dedgd)
    if(present(ec)) then
      ec = ec + functl(FUNC_C)%libvdwxc%energy
    end if
  end if

  ! Definition of tau in libxc is different, so we need to divide vtau by a factor of two
  if (present(vtau) .and. mgga_withexc) vtau = vtau / M_TWO

  if(present(deltaxc)) deltaxc = M_ZERO

  if(xcs%xc_density_correction == LR_X) then
    call xc_density_correction_calc(xcs, der, psolver, namespace, spin_channels, &
      rho, vx, dedd, deltaxc = deltaxc)

    if(calc_energy) then
      ! correct the energy density from Levy-Perdew, note that vx now
      ! contains the correction applied to the xc potential.
      do is = 1, spin_channels
        do ip = 1, der%mesh%np
          ex_per_vol(ip) = vx(ip)*(CNST(3.0)*rho(ip, is) + sum(der%mesh%x(ip, 1:der%mesh%sb%dim)*gdens(ip, 1:der%mesh%sb%dim, is)))
        end do
      end do
    end if

    if(present(ex)) ex = ex + dmf_integrate(der%mesh, ex_per_vol)
  end if

  ! this has to be done in inverse order
  if(mgga) call mgga_process()
  if( gga) call  gga_process()
  call lda_process()

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

    integer :: ib, is

    PUSH_SUB(xc_get_vxc.copy_global_to_local)

    do is = 1, spin_channels
      !$omp parallel do
      do ib = 1, n_block
        local(is, ib) = global(ib + ip - 1, is)
      end do
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

    integer :: ib, is

    PUSH_SUB(xc_get_vxc.copy_local_to_global)

    do is = 1, spin_channels
      !$omp parallel do
      do ib = 1, n_block
        global(ib + ip - 1, is) = global(ib + ip - 1, is) + local(is, ib)
      end do
    end do

    POP_SUB(xc_get_vxc.copy_local_to_global)
  end subroutine copy_local_to_global

  ! ---------------------------------------------------------
  subroutine space_loop_init(ip, np, nblock)
    integer, intent(in)  :: ip
    integer, intent(in)  :: np
    integer, intent(out) :: nblock

    PUSH_SUB(xc_get_vxc.space_loop_init)

    !Resize the dimension of the last block when the number of the mesh points
    !it is not a perfect divisor of the dimension of the blocks.
    nblock = min(np - ip + 1, N_BLOCK_MAX)

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
      if(present(ex_density)) then
        ex_per_vol => ex_density
      else
        SAFE_ALLOCATE(ex_per_vol(1:der%mesh%np))
      end if
      if(present(ec_density)) then
        ec_per_vol => ec_density
      else
        SAFE_ALLOCATE(ec_per_vol(1:der%mesh%np))
      end if
      ex_per_vol = M_ZERO
      ec_per_vol = M_ZERO
    end if

    SAFE_ALLOCATE(dedd(1:der%mesh%np_part, 1:spin_channels))

    select case(ispin)
    case(UNPOLARIZED)
      !$omp parallel do
      do ii = 1, der%mesh%np
        dedd(ii, 1) = CNST(0.0)
        dens(ii, 1) = max(rho(ii, 1), M_ZERO)
      end do

    case(SPIN_POLARIZED)
      !$omp parallel do
      do ii = 1, der%mesh%np
        dedd(ii, 1:2) = CNST(0.0)
        dens(ii, 1) = max(rho(ii, 1), M_ZERO)
        dens(ii, 2) = max(rho(ii, 2), M_ZERO)
      end do

    case(SPINORS)
      do ii = 1, der%mesh%np
        dedd(ii, 1:spin_channels) = CNST(0.0)
        d(1:spin_channels) = rho(ii, 1:spin_channels)
        dtot = d(1) + d(2)
        dpol = sqrt((d(1) - d(2))**2 + &
          M_FOUR*(rho(ii, 3)**2 + rho(ii, 4)**2))
        dens(ii, 1) = max(M_HALF*(dtot + dpol), M_ZERO)
        dens(ii, 2) = max(M_HALF*(dtot - dpol), M_ZERO)
      end do

    end select

    POP_SUB(xc_get_vxc.lda_init)
  end subroutine lda_init

  ! ---------------------------------------------------------
  !> deallocate variables allocated in lda_init
  subroutine lda_end()
    PUSH_SUB(xc_get_vxc.lda_end)

    SAFE_DEALLOCATE_A(dens)
    if(.not. present(ex_density)) then
      SAFE_DEALLOCATE_P(ex_per_vol)
    end if
    if(.not. present(ec_density)) then
      SAFE_DEALLOCATE_P(ec_per_vol)
    end if
    SAFE_DEALLOCATE_A(dedd)

    POP_SUB(xc_get_vxc.lda_end)
  end subroutine lda_end


  ! ---------------------------------------------------------
  !> calculates the LDA part of vxc, taking into account non-collinear spin
  subroutine lda_process()
    integer :: ip
    FLOAT :: d(2), dpol, vpol

    PUSH_SUB(xc_get_vxc.lda_process)

    if(present(vxc)) then
    
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

    end if
      
    POP_SUB(xc_get_vxc.lda_process)
  end subroutine lda_process


  ! ---------------------------------------------------------
  !> initialize GGAs
  !!   *) allocates gradient of the density (gdens), dedgd, and its local variants
  subroutine gga_init()
    PUSH_SUB(xc_get_vxc.gga_init)

    ! allocate variables
    SAFE_ALLOCATE(gdens(1:der%mesh%np, 1:der%mesh%sb%dim, 1:spin_channels))
    gdens = M_ZERO

    SAFE_ALLOCATE(dedgd(1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:spin_channels))
    dedgd = M_ZERO

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
    integer :: is
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

  subroutine calc_mvorb_alpha()
    FLOAT, allocatable :: gnon(:)
    FLOAT :: tb09_c, alpha
    FLOAT :: gn(MAX_DIM), n
    integer :: ii
    FLOAT :: parameters(3)

    PUSH_SUB(xc_get_vxc.calc_mvorb_alpha)

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
        gnon(ii) = sqrt(gnon(ii))
      end if
    end do

    tb09_c =  dmf_integrate(der%mesh, gnon)/der%mesh%sb%rcell_volume

    SAFE_DEALLOCATE_A(gnon)

    select case(functl(FUNC_C)%id)
    case(XC_HYB_GGA_XC_MVORB_HSE06)
      alpha = CNST(0.121983)+CNST(0.130711)*tb09_c**4

      if(alpha > 1) then
        write(message(1), '(a,f6.3,a)') 'MVORB mixing parameter bigger than one (' , alpha ,').'
        call messages_warning(1, namespace=namespace)
        alpha = CNST(0.25)
      end if


      parameters(1) = alpha
      parameters(2) = xcs%cam_omega
      parameters(3) = xcs%cam_omega
      call xc_f03_func_set_ext_params(functl(FUNC_C)%conf, parameters)
      !The name is confusing. Here alpha is the beta of hybrids in functionals, 
      !but is called alpha in the original paper.
      xcs%cam_beta = alpha

    case(XC_HYB_GGA_XC_MVORB_PBEH)
      alpha = -CNST(1.00778)+CNST(1.10507)*tb09_c

      if(alpha > 1) then
        write(message(1), '(a,f6.3,a)') 'MVORB mixing parameter bigger than one (' , alpha ,').'
        call messages_warning(1, namespace=namespace)
        alpha = CNST(0.25)
      end if
      if(alpha < 0) then
        write(message(1), '(a,f6.3,a)') 'MVORB mixing parameter smaller than zero (' , alpha ,').'
        call messages_warning(1, namespace=namespace)
        alpha = CNST(0.25)
      end if

#if defined HAVE_LIBXC5
      parameters(1) = alpha
      call xc_f03_func_set_ext_params(functl(FUNC_C)%conf, parameters)
#else
      call messages_not_implemented("MVORB with PBE0 requires libxc 5", namespace=namespace)
#endif
      xcs%cam_alpha = alpha
    case default
      call messages_not_implemented("MVORB density-based mixing for functionals other than PBE0 and HSE06", namespace=namespace)
    end select

    POP_SUB(xc_get_vxc.calc_mvorb_alpha)
  end subroutine calc_mvorb_alpha


  ! ---------------------------------------------------------

  subroutine calc_tb09_c()
    FLOAT, allocatable :: gnon(:)
    FLOAT :: gn(MAX_DIM), n, parameters(1)
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

    parameters(1) =  -CNST(0.012) + CNST(1.023)*sqrt(dmf_integrate(der%mesh, gnon)/der%mesh%sb%rcell_volume)

    call xc_f03_func_set_ext_params(functl(1)%conf, parameters)

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

  pure logical function family_is_gga(family)
    integer, intent(in) :: family

    family_is_gga = bitand(family, XC_FAMILY_GGA + XC_FAMILY_HYB_GGA + &
      XC_FAMILY_MGGA + XC_FAMILY_HYB_MGGA + XC_FAMILY_LIBVDWXC) /= 0
  end function  family_is_gga

  pure logical function family_is_mgga(family)
    integer, intent(in) :: family

    family_is_mgga = bitand(family, XC_FAMILY_MGGA + XC_FAMILY_HYB_MGGA) /= 0
  end function family_is_mgga

! -----------------------------------------------------

subroutine xc_density_correction_calc(xcs, der, psolver, namespace, nspin, density, refvx, vxc, deltaxc)
  type(xc_t),          intent(in)    :: xcs
  type(derivatives_t), intent(in)    :: der
  type(poisson_t),     intent(in)    :: psolver
  type(namespace_t),   intent(in)    :: namespace
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
  FLOAT   :: qxc, ncutoff, qxc_old, ncutoff_old, deriv, qxcfin
  FLOAT   :: x1, x2, x3, f1, f2, f3, dd, vol, mindd, maxdd
  FLOAT, allocatable :: nxc(:), lrvxc(:)
  type(profile_t), save :: prof
  FLOAT, parameter :: thres = CNST(1e-6)

  PUSH_SUB(xc_density_correction_calc)

  call profiling_in(prof, TOSTRING(XC_DENSITY_CORRECTION))

  SAFE_ALLOCATE(nxc(1:der%mesh%np))
  SAFE_ALLOCATE(lrvxc(1:der%mesh%np_part))

  do ip = 1, der%mesh%np
    lrvxc(ip) = CNST(-1.0)/(CNST(4.0)*M_PI)*refvx(ip)
  end do
  call dderivatives_lapl(der, lrvxc, nxc)

  if(debug%info) then
    call dio_function_output(OPTION__OUTPUTFORMAT__AXIS_X, "./static", "rho", namespace, &
      der%mesh, density(:, 1), unit_one, ierr)
    call dio_function_output(OPTION__OUTPUTFORMAT__AXIS_X, "./static", "vxcorig", namespace, &
      der%mesh, refvx(:), unit_one, ierr)
    call dio_function_output(OPTION__OUTPUTFORMAT__AXIS_X, "./static", "nxc", namespace, &
      der%mesh, nxc, unit_one, ierr)
  end if

  if(xcs%xcd_optimize_cutoff) then

    x1 = CNST(1.0e-8)
    qxc = get_qxc(der%mesh, nxc, density(:, 1), x1)
    deriv = HUGE(deriv)
    done = .false.

    INCR(iter, 1)
    if(debug%info) then
      if(mpi_world%rank == 0) then
        write(number, '(i4)') iter
        iunit = io_open('qxc.'//trim(adjustl(number)), namespace, action='write')
      end if
    end if
    do
      if(.not. done) then
        ncutoff_old = x1
        qxc_old = qxc
      end if

      x1 = x1*CNST(1.01)
      if(debug%info) then
        if(mpi_world%rank == 0) then
          write(iunit, *) x1, qxc
        end if
      end if
      if(x1 > CNST(1.0)) exit

      qxc = get_qxc(der%mesh, nxc, density(:, 1), x1)

      if(qxc == qxc_old) cycle

      deriv = (qxc - qxc_old)/(x1 - ncutoff_old)

      if(.not. done .and. abs(qxc) >= M_ONE) then
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
      f1 = M_ONE + get_qxc(der%mesh, nxc, density(:, 1), x1)
      f2 = M_ONE + get_qxc(der%mesh, nxc, density(:, 1), x2)

      do ip = 1, 20
        if(abs(f1 - f2) < CNST(1.0e-16)) exit
        x3 = x2 - f2*(x2 - x1)/(f2 - f1)
        f3 = M_ONE + get_qxc(der%mesh, nxc, density(:, 1), x3)
        if(abs(f3) < CNST(1.0e-6)) exit
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

  if(debug%info) then
    call dio_function_output(OPTION__OUTPUTFORMAT__AXIS_X, "./static", "nxcmod", namespace, &
      der%mesh, nxc, unit_one, ierr)
  
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

  if(debug%info) then
    call dio_function_output(OPTION__OUTPUTFORMAT__AXIS_X, "./static", "fulldiffvxc.ax", namespace, &
      der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(OPTION__OUTPUTFORMAT__AXIS_Y, "./static", "fulldiffvxc.ax", namespace, &
      der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(OPTION__OUTPUTFORMAT__AXIS_Z, "./static", "fulldiffvxc.ax", namespace, &
      der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(OPTION__OUTPUTFORMAT__PLANE_X, "./static", "fulldiffvxc.pl", namespace, &
      der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(OPTION__OUTPUTFORMAT__PLANE_Y, "./static", "fulldiffvxc.pl", namespace, &
      der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(OPTION__OUTPUTFORMAT__PLANE_Z, "./static", "fulldiffvxc.pl", namespace, &
      der%mesh, lrvxc, unit_one, ierr)
  end if

  do ip = 1, der%mesh%np
    vxc(ip, 1:nspin) = vxc(ip, 1:nspin) + lrvxc(ip)
    refvx(ip) = lrvxc(ip)
  end do

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

  if(debug%info) then
    call dio_function_output(OPTION__OUTPUTFORMAT__AXIS_X, "./static", "diffvxc.ax", namespace, &
      der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(OPTION__OUTPUTFORMAT__AXIS_Y, "./static", "diffvxc.ax", namespace, &
      der%mesh, lrvxc, unit_one, ierr)
    call dio_function_output(OPTION__OUTPUTFORMAT__AXIS_Z, "./static", "diffvxc.ax", namespace, &
      der%mesh, lrvxc, unit_one, ierr)
  end if
  
  dd = dmf_integrate(der%mesh, lrvxc)/vol

  if(debug%info) then
    if(mpi_world%rank == 0) then
      print*, "DD",  -CNST(2.0)*dd, -CNST(2.0)*mindd, -CNST(2.0)*maxdd
    end if
  end if
  
  if(present(deltaxc)) deltaxc = -CNST(2.0)*dd

  if(debug%info) then
    call dio_function_output(OPTION__OUTPUTFORMAT__AXIS_X, "./static", "fnxc", namespace, &
      der%mesh, nxc, unit_one, ierr)
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
      nxc2(ip) = CNST(0.0)
    else
      nxc2(ip) = nxc(ip)
    end if
  end do

  qxc = dmf_integrate(mesh, nxc2)

  SAFE_DEALLOCATE_A(nxc2)

  POP_SUB(get_qxc)
end function get_qxc


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

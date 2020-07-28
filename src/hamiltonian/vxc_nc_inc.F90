!! Copyright (C) 2020-2021 N. Tancogne-Dejean
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

!> This routines is similar to xc_get_vxc but for noncollinear functionals, which are not implemented in libxc
subroutine xc_get_nc_vxc(der, xcs, st, kpoints, namespace, rho, vxc, ex, ec, vtau, ex_density, ec_density)
  type(derivatives_t),    intent(in)    :: der             !< Discretization and the derivative operators and details
  type(xc_t), target,     intent(inout) :: xcs             !< Details about the xc functional used
  type(states_elec_t),    intent(in)    :: st              !< State of the system (wavefunction,eigenvalues...)
  type(kpoints_t),        intent(in)    :: kpoints
  type(namespace_t),      intent(in)    :: namespace
  FLOAT,                  intent(in)    :: rho(:, :)       !< Electronic density 
  FLOAT,                  intent(inout) :: vxc(:,:)        !< XC potential
  FLOAT, optional,        intent(inout) :: ex              !< Exchange energy.
  FLOAT, optional,        intent(inout) :: ec              !< Correlation energy.
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
  FLOAT, allocatable :: l_zk(:,:)
  !  First order (functional) derivatives
  FLOAT, allocatable :: l_dedd(:,:)
  FLOAT, allocatable :: l_vsigma(:,:)
  FLOAT, allocatable :: l_dedldens(:,:)
  FLOAT, allocatable :: l_dedtau(:,:)

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

  integer :: ib, ip, isp, families, ixc, idir, ipstart, ipend
  FLOAT   :: energy(1:2)
  logical :: mgga, mgga_withexc
  type(profile_t), save :: prof
  logical :: calc_energy
  type(xc_functl_t), pointer :: functl(:)
  type(distributed_t) :: distribution
  type(profile_t), save :: prof_gather
  
  PUSH_SUB(xc_get_nc_vxc)
  call profiling_in(prof, TOSTRING(X(XC_NC_LOCAL)))

  nullify(ex_per_vol)
  nullify(ec_per_vol)
  
  ASSERT(present(ex) .eqv. present(ec))
  calc_energy = present(ex) .or. present(ex_density) .or. present(ec_density)

  functl => xcs%functional(:, 1)

  ! is there anything to do ?
  families = XC_FAMILY_NC_LDA + XC_FAMILY_NC_MGGA
  if(bitand(xcs%family, families) == 0) then
    POP_SUB(xc_get_nc_vxc)
    call profiling_out(prof)
    return
  end if

  do ixc = 1, 2
    if(functl(ixc)%family /= XC_FAMILY_NONE) then
      ASSERT(bitand(functl(ixc)%flags, XC_FLAGS_HAVE_VXC) /= 0)
    end if
  end do

  ! initialize a couple of handy variables
  mgga = family_is_nc_mgga(xcs%family)
  mgga_withexc = family_is_mgga_with_exc(xcs)
  if(mgga_withexc) then
    ASSERT(present(vtau))
  end if

  call nc_lda_init()
  if(mgga) then
    call nc_gga_init()
    call nc_mgga_init()
  end if

  if(mgga) then
    ! We calculate everything from the wavefunctions to benefit from
    ! the error cancellation between the gradient of the density and
    ! tau.
    !
    ! FIXME: Probably this is wrong for non-local corrections or other
    ! cases when the density is not directly generated by the
    ! orbitals.
    if(allocated(st%rho_core)) then
      call messages_not_implemented("MGGA with nonlinear core correction", namespace=namespace)
    end if

    if (xcs%use_gi_ked) then
      call states_elec_calc_quantities(der, st, kpoints, .true., gi_kinetic_energy_density = tau, &
                                  density_gradient = gdens, density_laplacian = ldens)
    else
      call states_elec_calc_quantities(der, st, kpoints, .true., kinetic_energy_density = tau, &
                                  density_gradient = gdens, density_laplacian = ldens)
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

      if(calc_energy .and. bitand(functl(ixc)%flags, XC_FLAGS_HAVE_EXC) /= 0) then
        ! we get the xc energy and potential
        select case(functl(ixc)%family)
        !TODO: Implement
        case(XC_FAMILY_NC_LDA)
          ASSERT(.false.)

        case(XC_FAMILY_NC_MGGA)
          call nc_mgga_exc_vxc(functl(ixc), namespace, n_block, l_dens, l_sigma, l_ldens, l_tau, &
                                l_dedd, l_vsigma, l_dedldens, l_dedtau, l_zk)

        case default
          cycle
        end select

      else ! we just get the potential
        l_zk(:,:) = M_ZERO

        select case(functl(ixc)%family)
        !TODO: Implement
        case(XC_FAMILY_NC_LDA)
          ASSERT(.false.)

        case(XC_FAMILY_NC_MGGA)
          call nc_mgga_exc_vxc(functl(ixc), namespace, n_block, l_dens, l_sigma, l_ldens, l_tau, & 
                                l_dedd, l_vsigma, l_dedldens, l_dedtau)

        case default
          cycle
        end select
        
      end if

      if(calc_energy) then
        if(functl(ixc)%type == XC_EXCHANGE) then
          do ib = 1, n_block
            ex_per_vol(ib + ip - 1) = ex_per_vol(ib + ip - 1) + sum(l_dens(1:2, ib)*l_zk(1:2,ib))
          end do
        else
          do ib = 1, n_block
            ec_per_vol(ib + ip - 1) = ec_per_vol(ib + ip - 1) + sum(l_dens(1:2, ib)*l_zk(1:2,ib))
          end do
        end if
      end if

      call copy_local_to_global(l_dedd, dedd, n_block, ip) 

      if(family_is_mgga(functl(ixc)%family)) then
        do ib = 1, n_block
          dedgd(ib+ip-1,:,1) = dedgd(ib+ip-1,:,1) + M_TWO*l_vsigma(1, ib)*gdens(ib + ip - 1,:,1) &
             + M_TWO*(l_vsigma(3, ib)*gdens(ib + ip - 1,:,3)-l_vsigma(4, ib)*gdens(ib + ip - 1,:,4))
          dedgd(ib+ip-1,:,2) = dedgd(ib+ip-1,:,2) + M_TWO*l_vsigma(2, ib)*gdens(ib + ip - 1,:,2) &
             + M_TWO*(l_vsigma(3, ib)*gdens(ib + ip - 1,:,3)-l_vsigma(4, ib)*gdens(ib + ip - 1,:,4))
          dedgd(ib+ip-1,:,3) = dedgd(ib+ip-1,:,3) + (l_vsigma(1, ib)+l_vsigma(2, ib))*gdens(ib + ip - 1,:,3) &
                                   + (gdens(ib + ip - 1,:,1)+gdens(ib + ip - 1,:,2))*l_vsigma(3, ib)
          dedgd(ib+ip-1,:,4) = dedgd(ib+ip-1,:,4) + (l_vsigma(1, ib)+l_vsigma(2, ib))*gdens(ib + ip - 1,:,4) &
                                   - (gdens(ib + ip - 1,:,1)+gdens(ib + ip - 1,:,2))*l_vsigma(4, ib)
        end do
      end if

      if(family_is_mgga(functl(ixc)%family)) then
        ! XXXXX does this work correctly when functionals belong to
        ! different families and only one is mgga?
        call copy_local_to_global(l_dedldens, dedldens, n_block, ip)
        if(bitand(functl(ixc)%flags, XC_FLAGS_HAVE_EXC) /= 0 ) &
          call copy_local_to_global(l_dedtau, vtau, n_block, ip)
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
      call comm_allreduce(st%dom_st_kpt_mpi_grp, energy)
    else if(der%mesh%parallel_in_domains) then
      call comm_allreduce(der%mesh%mpi_grp, energy)
    end if

    ex = ex + energy(1)
    ec = ec + energy(2)
  end if

  if(xcs%parallel) then
    if(distribution%parallel) then
      call profiling_in(prof_gather, TOSTRING(X(XC_GATHER)))

      do isp = 1, 4
        call distributed_allgather(distribution, dedd(:, isp))
      end do

      if(mgga) then
        do idir = 1, der%mesh%sb%dim
          do isp = 1, 4
            call distributed_allgather(distribution, dedgd(:, idir, isp))
          end do
        end do
      end if

      if(mgga) then
        do isp = 1, 4
          call distributed_allgather(distribution, dedldens(:, isp))
          if(mgga_withexc) then
            call distributed_allgather(distribution, vtau(:, isp))
          end if
        end do
      end if

      call profiling_out(prof_gather)
    end if

    call distributed_end(distribution)
  end if

  ! Definition of tau in libxc is different, so we need to divide vtau by a factor of two
!  if (present(vtau) .and. mgga_withexc) vtau = vtau / M_TWO

  ! this has to be done in inverse order
  if(mgga) then
    call nc_mgga_process()
    call nc_gga_process()
  end if
  call nc_lda_process()

  ! clean up allocated memory
  call nc_lda_end()
  if(mgga) then
    call nc_gga_end()
    call nc_mgga_end()
  end if

  POP_SUB(xc_get_nc_vxc)
  call profiling_out(prof)

contains
  !  ---------------------------------------------------------
  !> make a local copy with the correct memory order for libxc
  subroutine copy_global_to_local(global, local, n_block, ip)
    FLOAT,   intent(in)  :: global(:,:)
    FLOAT,   intent(out) :: local(:,:)
    integer, intent(in)  :: n_block
    integer, intent(in)  :: ip

    integer :: ib, is

    PUSH_SUB(xc_get_nc_vxc.copy_global_to_local)

    do is = 1, 4
      !$omp parallel do
      do ib = 1, n_block
        local(is, ib) = global(ib + ip - 1, is)
      end do
    end do

    POP_SUB(xc_get_nc_vxc.copy_global_to_local)
  end subroutine copy_global_to_local

  ! ---------------------------------------------------------
  subroutine copy_local_to_global(local, global, n_block, ip)
    FLOAT,   intent(in)    :: local(:,:)
    FLOAT,   intent(inout) :: global(:,:)
    integer, intent(in)    :: n_block
    integer, intent(in)    :: ip

    integer :: ib, is

    PUSH_SUB(xc_get_nc_vxc.copy_local_to_global)

    do is = 1, 4
      !$omp parallel do
      do ib = 1, n_block
        global(ib + ip - 1, is) = global(ib + ip - 1, is) + local(is, ib)
      end do
    end do

    POP_SUB(xc_get_nc_vxc.copy_local_to_global)
  end subroutine copy_local_to_global

  ! ---------------------------------------------------------
  subroutine space_loop_init(ip, np, nblock)
    integer, intent(in)  :: ip
    integer, intent(in)  :: np
    integer, intent(out) :: nblock

    CMPLX :: tmp(MAX_DIM), tmp_sum

    PUSH_SUB(xc_get_nc_vxc.space_loop_init)

    !Resize the dimension of the last block when the number of the mesh points
    !it is not a perfect divisor of the dimension of the blocks.
    nblock = min(np - ip + 1, N_BLOCK_MAX)

    ! make a local copy with the correct memory order for libxc
    call copy_global_to_local(dens, l_dens, nblock, ip)

    if(mgga) then
      do ib = 1, nblock
        l_sigma(1, ib) = sum(gdens(ib + ip - 1, 1:der%mesh%sb%dim, 1)**2 &
            + gdens(ib + ip - 1, 1:der%mesh%sb%dim, 3)**2 + gdens(ib + ip - 1, 1:der%mesh%sb%dim, 4)**2)
        l_sigma(2, ib) = sum(gdens(ib + ip - 1, 1:der%mesh%sb%dim, 2)**2 &
            + gdens(ib + ip - 1, 1:der%mesh%sb%dim, 3)**2 + gdens(ib + ip - 1, 1:der%mesh%sb%dim, 4)**2)

        tmp =cmplx(gdens(ib + ip - 1, 1:der%mesh%sb%dim, 3),gdens(ib + ip - 1, 1:der%mesh%sb%dim, 4))
        tmp_sum = sum( tmp * (gdens(ib + ip - 1, 1:der%mesh%sb%dim, 1) + gdens(ib + ip - 1, 1:der%mesh%sb%dim, 2)))
        l_sigma(3, ib) = TOFLOAT(tmp_sum)
        l_sigma(4, ib) = aimag(tmp_sum)
      end do
      call copy_global_to_local(tau, l_tau, nblock, ip)
      ! we adjust for the different definition of tau in libxc
    !  l_tau(1:4, 1:nblock) = l_tau(1:4, 1:nblock) / M_TWO
      call copy_global_to_local(ldens, l_ldens, nblock, ip)

    end if

    POP_SUB(xc_get_nc_vxc.space_loop_init)
  end subroutine space_loop_init

  ! ---------------------------------------------------------
  !> Takes care of the initialization of the LDA part of the functionals
  !!   *) allocates density and dedd, and their local variants
  !!   *) calculates the density taking into account nlcc and non-collinear spin
  subroutine nc_lda_init()
    integer :: ii
    FLOAT   :: dtot, dpol

    PUSH_SUB(xc_get_nc_vxc.nc_lda_init)

    ! allocate some general arrays
    SAFE_ALLOCATE(dens(1:der%mesh%np_part, 1:4))
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

    do ii = 1, der%mesh%np
      vxc(ii, 1:4) = M_ZERO
    end do

    SAFE_ALLOCATE(dedd(1:der%mesh%np_part, 1:4))
    dedd(1:der%mesh%np, 1:4) = M_ZERO
    !$omp parallel do
    do ii = 1, der%mesh%np
      dens(ii, 1) = max(rho(ii, 1), M_ZERO)
      dens(ii, 2) = max(rho(ii, 2), M_ZERO)
      dens(ii, 3) = rho(ii, 3)
      dens(ii, 4) = rho(ii, 4)
    end do

    POP_SUB(xc_get_nc_vxc.nc_lda_init)
  end subroutine nc_lda_init

  ! ---------------------------------------------------------
  !> deallocate variables allocated in nc_lda_init
  subroutine nc_lda_end()
    PUSH_SUB(xc_get_nc_vxc.nc_lda_end)

    SAFE_DEALLOCATE_A(dens)
    if(.not. present(ex_density)) then
      SAFE_DEALLOCATE_P(ex_per_vol)
    end if
    if(.not. present(ec_density)) then
      SAFE_DEALLOCATE_P(ec_per_vol)
    end if
    SAFE_DEALLOCATE_A(dedd)

    POP_SUB(xc_get_nc_vxc.nc_lda_end)
  end subroutine nc_lda_end

  ! ---------------------------------------------------------
  !> calculates the LDA part of vxc
  subroutine nc_lda_process()
    integer :: is

    PUSH_SUB(xc_get_nc_vxc.nc_lda_process)

    do is = 1, 4
      call lalg_axpy(der%mesh%np, M_ONE, dedd(:,is), vxc(:,is))
    end do

    POP_SUB(xc_get_nc_vxc.nc_lda_process)
  end subroutine nc_lda_process

  ! ---------------------------------------------------------
  !> initialize GGAs
  !!   *) allocates gradient of the density (gdens), dedgd, and its local variants
  subroutine nc_gga_init()
    PUSH_SUB(xc_get_nc_vxc.nc_gga_init)

    ! allocate variables
    SAFE_ALLOCATE(gdens(1:der%mesh%np, 1:der%mesh%sb%dim, 1:4))
    gdens = M_ZERO

    SAFE_ALLOCATE(dedgd(1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:4))
    dedgd = M_ZERO

    POP_SUB(xc_get_nc_vxc.nc_gga_init)
  end subroutine nc_gga_init


  ! ---------------------------------------------------------
  !> cleans up memory allocated in gga_init
  subroutine nc_gga_end()
    PUSH_SUB(xc_get_nc_vxc.nc_gga_end)

    SAFE_DEALLOCATE_A(gdens)
    SAFE_DEALLOCATE_A(dedgd)

    POP_SUB(xc_get_nc_vxc.nc_gga_end)
  end subroutine nc_gga_end

  ! ---------------------------------------------------------
  !> calculates the GGA contribution to vxc
  subroutine nc_gga_process()
    integer :: is
    FLOAT, allocatable :: gf(:)

    PUSH_SUB(xc_get_nc_vxc.nc_gga_process)

    ! subtract the divergence of the functional derivative of Exc with respect to
    ! the gradient of the density.
    SAFE_ALLOCATE(gf(1:der%mesh%np))
    do is = 1, 4
      call dderivatives_div(der, dedgd(:, :, is), gf(:))
      call lalg_axpy(der%mesh%np, -M_ONE, gf, dedd(:, is))
    end do
    SAFE_DEALLOCATE_A(gf)

    POP_SUB(xc_get_nc_vxc.nc_gga_process)
  end subroutine nc_gga_process

  ! ---------------------------------------------------------
  !> initialize meta-GGAs
  !!   *) allocate the kinetic-energy density, dedtau, and local variants
  subroutine nc_mgga_init()
    PUSH_SUB(xc_get_nc_vxc.nc_mgga_init)

    ! allocate variables
    SAFE_ALLOCATE( tau(1:der%mesh%np, 1:4))
    SAFE_ALLOCATE(ldens(1:der%mesh%np, 1:4))

    SAFE_ALLOCATE(dedldens(1:der%mesh%np_part, 1:4))
    dedldens = M_ZERO

    POP_SUB(xc_get_nc_vxc.nc_mgga_init)
  end subroutine nc_mgga_init

  ! ---------------------------------------------------------
  !> clean up memory allocated in nc_mgga_init
  subroutine nc_mgga_end()
    PUSH_SUB(xc_get_nc_vxc.nc_mgga_end)

    SAFE_DEALLOCATE_A(tau)
    SAFE_DEALLOCATE_A(ldens)
    SAFE_DEALLOCATE_A(dedldens)

    POP_SUB(xc_get_nc_vxc.nc_mgga_end)
  end subroutine nc_mgga_end

  ! ---------------------------------------------------------
  !> calculate the mgga contribution to vxc
  subroutine nc_mgga_process()
    integer :: is
    FLOAT, allocatable :: lf(:)

    PUSH_SUB(xc_get_nc_vxc.nc_mgga_process)

    ! add the Laplacian of the functional derivative of Exc with respect to
    ! the laplacian of the density.

    SAFE_ALLOCATE(lf(1:der%mesh%np))
    do is = 1, 4
      call dderivatives_lapl(der, dedldens(:, is), lf)
      call lalg_axpy(der%mesh%np, M_ONE, lf, dedd(:, is))
    end do
    SAFE_DEALLOCATE_A(lf)

    POP_SUB(xc_get_nc_vxc.nc_mgga_process)
  end subroutine nc_mgga_process


  ! ---------------------------------------------------------
  !> THREADSAFE (no SAFE ALLOCATE or PUSH/POP SUB)
  subroutine local_allocate()

    allocate(l_dens(1:4, 1:N_BLOCK_MAX))
    allocate(l_zk(1:2, 1:N_BLOCK_MAX))

    if(mgga) then
      allocate(l_sigma   (1:4, 1:N_BLOCK_MAX))
      allocate(l_tau     (1:4, 1:N_BLOCK_MAX))
      allocate(l_ldens   (1:4, 1:N_BLOCK_MAX))
      allocate(l_dedd    (1:4, 1:N_BLOCK_MAX))
      allocate(l_vsigma  (1:4, 1:N_BLOCK_MAX))
      allocate(l_dedldens(1:4, 1:N_BLOCK_MAX))
      allocate(l_dedtau  (1:4, 1:N_BLOCK_MAX))
    end if

  end subroutine local_allocate

  ! ---------------------------------------------------------
  !> THREADSAFE (no SAFE ALLOCATE or PUSH/POP SUB)
  subroutine local_deallocate()

    deallocate(l_dens)
    deallocate(l_zk)
    deallocate(l_dedd)

    if(mgga) then
      deallocate(l_sigma)
      deallocate(l_tau)
      deallocate(l_ldens)
      deallocate(l_vsigma)
      deallocate(l_dedldens)
      deallocate(l_dedtau)
    end if

  end subroutine local_deallocate

end subroutine xc_get_nc_vxc


  pure logical function family_is_nc_mgga(family)
    integer, intent(in) :: family

    family_is_nc_mgga = bitand(family, XC_FAMILY_NC_MGGA) /= 0
  end function family_is_nc_mgga

  ! -----------------------------------------------------
  subroutine  nc_mgga_exc_vxc(functl, namespace, n_block, l_dens, l_sigma, l_ldens, l_tau, &
                                l_dedd, l_vsigma, l_deddldens, l_dedtau, l_zk)
    type(xc_functl_t), intent(in)    :: functl
    type(namespace_t), intent(in)    :: namespace
    integer,           intent(in)    :: n_block
    FLOAT,             intent(in)    :: l_dens(:,:)     ! Density 
    FLOAT,             intent(in)    :: l_sigma(:,:)    ! Modulus squared of the gradient of the density
    FLOAT,             intent(in)    :: l_ldens(:,:)    ! Laplacian of the density
    FLOAT,             intent(in)    :: l_tau(:,:)      ! Kinetic energy density
    FLOAT,             intent(out)   :: l_dedd(:,:)     ! Derivative of the energy versus l_dens
    FLOAT,             intent(out)   :: l_vsigma(:,:)   ! Derivative of the energy versus l_sigma
    FLOAT,             intent(out)   :: l_deddldens(:,:)! Derivative of the energy versus the l_dens
    FLOAT,             intent(out)   :: l_dedtau(:,:)   ! Derivative of the energy versus l_tau
    FLOAT, optional,   intent(out)   :: l_zk(:,:)       ! Energy density

    integer :: ib

    PUSH_SUB(nc_mgga_exc_vxc)

    if(present(l_zk)) then 
      ASSERT(bitand(functl%flags, XC_FLAGS_HAVE_EXC) /= 0) 
    end if

    select case(functl%id)
    case(XC_MGGA_X_NC_BR)

      do ib = 1, n_block

        if(present(l_zk)) then
          call nc_br_vxc_exc(l_dens(:, ib), l_sigma(:, ib), l_ldens(:, ib), l_tau(:, ib), &
                  l_dedd(:,ib), l_vsigma(:,ib), l_deddldens(:,ib), l_dedtau(:,ib), l_zk(:,ib)) 
        else
          call nc_br_vxc_exc(l_dens(:, ib), l_sigma(:, ib), l_ldens(:, ib), l_tau(:, ib), &
                  l_dedd(:,ib), l_vsigma(:,ib), l_deddldens(:,ib), l_dedtau(:,ib))
        end if
      end do
    case default
      message(1) = "Unsupported noncollinear functional"
      call messages_fatal(1, namespace=namespace)
    end select

    POP_SUB(nc_mgga_exc_vxc)
  end subroutine

  !Computes the local curvature of the exchange-hole and get the corresponding values of x and b
  !This allows to compute the local Coulomb potential and the energy
  !The exchange potential is finally constructed from the potential
  subroutine nc_br_vxc_exc(l_dens, l_sigma, l_ldens, l_tau, l_dedd, l_vsigma, l_dedldens, l_dedtau, l_zk)
    FLOAT,             intent(in)    :: l_dens(:)     ! Density 
    FLOAT,             intent(in)    :: l_sigma(:)    ! Modulus squared of the gradient of the density
    FLOAT,             intent(in)    :: l_ldens(:)    ! Laplacian of the density
    FLOAT,             intent(in)    :: l_tau(:)      ! Kinetic energy density
    FLOAT,             intent(out)   :: l_dedd(:)     ! Derivative of the energy versus l_dens
    FLOAT,             intent(out)   :: l_vsigma(:)   ! Derivative of the energy versus l_sigma
    FLOAT,             intent(out)   :: l_dedldens(:) ! Derivative of the energy versus the l_dens
    FLOAT,             intent(out)   :: l_dedtau(:)   ! Derivative of the energy versus l_tau
    FLOAT, optional,   intent(out)   :: l_zk(:)       ! Energy density

    integer :: is
    FLOAT :: l_ontop, l_curv, xx_BR, l_dens_tol(2)
    FLOAT :: cnst, U_BR, P_BR, dUdx
    FLOAT :: dtop_dn(3), dcurv_dn(3), dcurv_d_dn(3)
    CMPLX :: offdiag
    FLOAT, parameter :: gamma = CNST(0.8)
    FLOAT, parameter :: tol_Q = CNST(1e-20)
    FLOAT, parameter :: tol_den = M_EPSILON
    FLOAT, parameter :: tol_x = CNST(1e-7)
    !This tolerance here seems crucial. I checked using Matlab that 
    !the error using the Taylor expansion is lower than 1e-8 compared to the
    !full formula below this threshold. However, for x below this threshold, the
    !full formula become unstable and the error grows. NTD
    FLOAT, parameter :: tol_dUdx = CNST(1e-4)

    l_dedd(1:4)     = M_ZERO
    l_vsigma(1:4)   = M_ZERO
    l_dedldens(1:4) = M_ZERO
    l_dedtau(1:4)     = M_ZERO
    
    do is = 1, 2
      l_dens_tol(is) = l_dens(is) + M_EPSILON

      !We compute the value of the on-top effective exchange hole
      l_ontop = l_dens(is) + (l_dens(3)**2 + l_dens(4)**2)/l_dens_tol(is)
      !We also compute of the curvature and we normalize it properly
      l_curv = l_ldens(is) - M_TWO*gamma*(l_tau(is) - M_FOURTH*l_sigma(is)/(l_dens_tol(is)))
      !Off-diagonal term of the curvature
      offdiag = l_dens(3)*(l_ldens(3)-M_TWO*gamma*l_tau(3)) + l_dens(4)*(l_ldens(4)-M_TWO*gamma*l_tau(4))
      l_curv = (l_curv + offdiag / (l_dens_tol(is))) / CNST(6.0)
      !Avoids division by zero
      if(abs(l_ontop) < tol_Q) l_ontop = sign(tol_Q, l_ontop)
      if(abs(l_curv) < tol_Q) l_curv = sign(tol_Q, l_curv)


      !We get the value of x for up and down components
      xx_BR = nc_br_get_x(l_ontop, l_curv)
      if(xx_BR < M_ZERO) then
        message(1) = "Newton-Raphson produced a negative x value"
        call messages_fatal(1)
      end if

      cnst = M_TWO * (M_PI*l_ontop)**M_THIRD

      !We construct the enengy density 
      if (abs(xx_BR) > tol_x) then
        U_BR = -( M_ONE - exp(-xx_BR)*(M_ONE + M_HALF*xx_BR) ) & 
                   / xx_BR * exp(xx_BR*M_THIRD) * cnst
      else
        !Taylor expansion at x=0
        U_BR = -cnst * (M_HALF + xx_BR/CNST(6.0) - xx_BR**2/CNST(18.0))
      end if

      if(present(l_zk)) then
        l_zk(is) = U_BR*M_HALF
      end if

      ! Minus sign cancels with the next equation
      if(abs(xx_BR) > tol_dUdx) then
        dUdx = (xx_BR**2 + M_TWO*xx_BR + exp(xx_BR)*(xx_BR-M_THREE) + M_THREE) &
                 * exp(-xx_BR*M_TWOTHIRD) / (M_THREE * xx_BR**2)
      else
        !Taylor expansion at x=0
        dUdx = M_ONE/CNST(6.0) - xx_BR/CNST(9.0)
      end if
      dUdx = dUdx * cnst

      P_BR = dUdx * (xx_BR-M_TWO)*xx_BR/(xx_BR**2-M_TWO*xx_BR+M_THREE) &
            * M_THREE * M_HALF * l_dens(is) / l_ontop

      ! The off-diagonal terms in the following are defined as the derivative w.r.t. 
      ! the down-up term (which gives the up-down potential)

      ! Derivative of the on-top density wrt the density
      dtop_dn(1) = M_ONE - (l_dens(3)**2 + l_dens(4)**2)/(l_dens_tol(is)**2)
      dtop_dn(2) = l_dens(3)/(l_dens_tol(is))
      dtop_dn(3) = l_dens(4)/(l_dens_tol(is))

      ! Derivative of the curvature wrt the density
      dcurv_dn(1) =-M_ONE/(l_dens_tol(is)**2) * ( M_HALF*gamma*l_sigma(is) &
               +  (l_dens(3)*(l_ldens(3) - gamma*M_TWO*l_tau(3))) &
               +  (l_dens(4)*(l_ldens(4) - gamma*M_TWO*l_tau(4))))
      dcurv_dn(2) = M_HALF/(l_dens_tol(is)) * ( l_ldens(3) - M_TWO*gamma*l_tau(3) )
      dcurv_dn(3) = M_HALF/(l_dens_tol(is)) * ( l_ldens(4) - M_TWO*gamma*l_tau(4) )
      dcurv_dn(1:3) = dcurv_dn(1:3) / CNST(6.0)

      ! Derivative of the energy wrt the density
      l_dedd(is) = (M_ONE + M_THIRD*l_dens(is)/l_ontop*dtop_dn(1))*U_BR &
                  + P_BR*(M_FIVE*M_THIRD*dtop_dn(1) - l_ontop/l_curv*dcurv_dn(1))

      l_dedd(3) = l_dedd(3) + M_THIRD * l_dens(3) * U_BR/l_ontop
      l_dedd(4) = l_dedd(4) + M_THIRD * l_dens(4) * U_BR/l_ontop
      l_dedd(3) = l_dedd(3) + P_BR*( M_FIVE*M_THIRD*dtop_dn(2) - l_ontop/l_curv * dcurv_dn(2)) 
      l_dedd(4) = l_dedd(4) + P_BR*( M_FIVE*M_THIRD*dtop_dn(3) - l_ontop/l_curv * dcurv_dn(3))

      !Here we recompute Pup and Pdn to avoid diving by l_curv, which can be very small
      !We also divide by the density, to avoid doing it just after
      P_BR = dUdx * xx_BR**2/(xx_BR**2-M_TWO*xx_BR+M_THREE) * exp(-xx_BR*M_TWOTHIRD) &
                * M_THREE * M_HALF * M_THREE / (cnst**2)

      ! Derivative of the energy wrt sigma
      l_vsigma(is) =-M_HALF * P_BR / l_ontop

      ! Derivative of the energy wrt the laplacian of the density
      l_dedldens(is) = -P_BR * l_dens(is) / l_ontop
      l_dedldens(3) = l_dedldens(3) - M_HALF * P_BR * l_dens(3)/ l_ontop
      l_dedldens(4) = l_dedldens(4) - M_HALF * P_BR * l_dens(4)/ l_ontop

    end do !is
  
    l_dedd(1:4) = l_dedd(1:4)*M_HALF
    l_vsigma(1:4) = l_vsigma(1:4) * gamma / CNST(12.0)
    l_dedldens(1:4) = l_dedldens(1:4) / CNST(12.0)
    
    ! Derivative of the energy wrt the kinetic energy
    l_dedtau(1:4) = -M_TWO * gamma * l_dedldens(1:4)

  end subroutine nc_br_vxc_exc

  ! Computes the coefficient x from the local density and the local curvature 
  ! of the exchange hole
  function nc_br_get_x(ldens, lcurv) result(br_x)
    FLOAT,             intent(in)    :: ldens      ! On-top exchange hole 
    FLOAT,             intent(in)    :: lcurv      ! Exchange-hole curvature
    FLOAT                            :: br_x 

    FLOAT :: rhs
    FLOAT, parameter :: tol = CNST(1e-13)
    FLOAT, parameter :: tol_Q = CNST(5e-15)

    !Reduced curvature
    rhs = lcurv / (ldens * ldens ** M_TWOTHIRD)
    if(abs(rhs) < tol_Q) rhs = sign(tol_Q, rhs)

    rhs = M_TWOTHIRD * M_PI ** M_TWOTHIRD / rhs

    if(abs(rhs) < tol) rhs = sign(tol, rhs) 
    br_x = nc_br_rtsafe(rhs, tol)

  end function nc_br_get_x

  !This is inspired by the safe Newton-Raphson method from numerical recipies
  !This function returns the value of x that fulfill the equation
  ! x exp(-2/3*x)/(x-2) = rhs 
  FLOAT function nc_br_rtsafe(rhs, tol) result(br_x)
    FLOAT, intent(in) :: rhs
    FLOAT, intent(in) :: tol

    integer, parameter :: maxit = 500
    FLOAT :: ff, dff, emx, xl, xh, dx, dxold, temp
    integer :: it

    !Select the correct branch for the function
    if(rhs < M_ZERO) then
      xh = M_ZERO
      xl = M_TWO - tol
    else
      xh = M_TWO + tol
      xl = M_TWO + M_ONE/rhs
    end if

    dxold = abs(xh-xl)
    dx = dxold
    br_x = M_HALF*(xl+xh)
   
    !Calculate the function and its derivative
    emx = exp(-M_TWOTHIRD*br_x)
    ff = br_x * emx / rhs - (br_x-M_TWO)
    dff = emx * ( M_ONE - M_TWOTHIRD * br_x) / rhs - M_ONE

    do it = 1, maxit
      if(((br_x-xh)*dff-ff) * ((br_x-xl)*dff-ff) > M_ZERO & ! Bisect if Newton out of range,
          .or. abs(M_TWO*ff) > abs(dxold*dff)) then            ! or not decreasing fast enough
        dxold = dx
        dx = M_HALF*(xh-xl)
        br_x = xl+dx
        if(xl == br_x) return !Change in root is negligible
      else ! Newton step acceptable. Take it.
        dxold = dx
        dx = ff/dff
        temp = br_x
        br_x = max(br_x - dx, M_ZERO)
        if(temp == br_x) return
      end if
      if(abs(dx) < tol) return ! Convergence criterion
      
      !Calculate the function and its derivative
      emx = exp(-M_TWOTHIRD*br_x)
      ff = br_x * emx / rhs - (br_x-M_TWO)
      dff = emx * ( M_ONE - M_TWOTHIRD * br_x) / rhs - M_ONE
      if(ff < M_ZERO) then ! Maintain the bracket on the root
        xl = br_x
      else
        xh = br_x
      end if
    end do

    message(1) = "Newton-Raphson did not converged"
    call messages_fatal(1)

  end function nc_br_rtsafe

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

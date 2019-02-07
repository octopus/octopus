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

subroutine xc_get_vxc(der, xcs, st, rho, ispin, ioniz_pot, qtot, vxc, ex, ec, deltaxc, vtau, ex_density, ec_density)
  type(derivatives_t),  intent(in)    :: der             !< Discretization and the derivative operators and details
  type(xc_t), target,   intent(in)    :: xcs             !< Details about the xc functional used
  type(states_t),       intent(in)    :: st              !< State of the system (wavefunction,eigenvalues...)
  FLOAT,                intent(in)    :: rho(:, :)       !< Electronic density 
  integer,              intent(in)    :: ispin           !< Number of spin channels 
  FLOAT,                intent(in)    :: ioniz_pot
  FLOAT,                intent(in)    :: qtot 
  FLOAT, optional,      intent(inout) :: vxc(:,:)        !< XC potential
  FLOAT, optional,      intent(inout) :: ex              !< Exchange energy.
  FLOAT, optional,      intent(inout) :: ec              !< Correlation energy.
  FLOAT, optional,      intent(inout) :: deltaxc         !< The XC derivative discontinuity
  FLOAT, optional,      intent(inout) :: vtau(:,:)       !< Derivative wrt (two times kinetic energy density)
  FLOAT, optional, target, intent(out)   :: ex_density(:)   !< The exchange energy density
  FLOAT, optional, target, intent(out)   :: ec_density(:)   !< The correlation energy density

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
  ! Global arrays
  !  Input quantities
  FLOAT, allocatable :: dens(:,:)      ! Density
  FLOAT, allocatable :: gdens(:,:,:)   ! Gradient of the density
  !  Energies
  FLOAT, pointer :: ex_per_vol(:)  ! Exchange energy per unit volume 
  FLOAT, pointer :: ec_per_vol(:)  ! Correlation energy per unit volume 
  !  First order (functional) derivatives
  FLOAT, allocatable :: dedd(:,:)      ! Derivative of the exchange or correlation energy wrt the density
  FLOAT, allocatable :: dedgd(:,:,:)   ! Derivative of the exchange or correlation energy wrt the gradient of the density

  FLOAT, allocatable :: vx(:)
  FLOAT, allocatable :: unp_dens(:), unp_dedd(:)

  integer :: ib, ip, isp, families, ixc, spin_channels, is, idir, ipstart, ipend
  FLOAT   :: rr, energy(1:2)
  logical :: gga
  type(profile_t), save :: prof, prof_libxc
  logical :: calc_energy
  type(xc_functl_t), pointer :: functl(:)
  type(distributed_t) :: distribution
  type(profile_t), save :: prof_gather
  
  PUSH_SUB(xc_get_vxc)
  call profiling_in(prof, "XC_LOCAL")

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
  families = XC_FAMILY_LDA + XC_FAMILY_GGA + XC_FAMILY_HYB_GGA + XC_FAMILY_MGGA + XC_FAMILY_HYB_MGGA
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

  !Read the spin channels
  spin_channels = functl(FUNC_X)%spin_channels

  if(xcs%xc_density_correction == LR_X) then
    SAFE_ALLOCATE(vx(1:der%mesh%np))
  end if

  call lda_init()
  if(gga .or. xcs%xc_density_correction == LR_X) call gga_init()

  ! Get the gradient and the Laplacian of the density and the kinetic-energy density
  if(gga) then
    do isp = 1, spin_channels 
      call dderivatives_grad(der, dens(:, isp), gdens(:, :, isp)) 
    end do
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

      call profiling_in(prof_libxc, "LIBXC")

      if(calc_energy .and. bitand(functl(ixc)%flags, XC_FLAGS_HAVE_EXC) /= 0) then
        ! we get the xc energy and potential
        select case(functl(ixc)%family)
        case(XC_FAMILY_LDA)
          call XC_F90(lda_exc_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_zk(1), l_dedd(1,1))

        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          call XC_F90(gga_exc_vxc)(functl(ixc)%conf, n_block, l_dens(1,1), l_sigma(1,1), &
            l_zk(1), l_dedd(1,1), l_vsigma(1,1))

        case default
          call profiling_out(prof_libxc)
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

      if(family_is_gga(functl(ixc)%family)) then
        do ib = 1, n_block
          dedgd(ib + ip - 1,:,1) = dedgd(ib + ip - 1,:,1) + M_TWO*l_vsigma(1, ib)*gdens(ib + ip - 1,:,1)
          if(ispin /= UNPOLARIZED) then
            dedgd(ib + ip - 1,:,1) = dedgd(ib + ip - 1,:,1) + l_vsigma(2, ib)*gdens(ib + ip - 1,:,2)
            dedgd(ib + ip - 1,:,2) = dedgd(ib + ip - 1,:,2) +  &
                 M_TWO*l_vsigma(3, ib)*gdens(ib + ip - 1,:,2) + l_vsigma(2, ib)*gdens(ib + ip - 1,:,1)
          end if
        end do
      end if

    end do functl_loop
  end do space_loop

  call local_deallocate()

  ! calculate the energy, we do the integrals directly so when the data
  ! is fully distributed we do not need to allgather first
  if(present(ex) .or. present(ec)) then

    energy(1:2) = CNST(0.0)
    
    do ip = ipstart, ipend
      energy(1) = energy(1) + ex_per_vol(ip)
      energy(2) = energy(2) + ec_per_vol(ip)
    end do
    
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
      call profiling_in(prof_gather, "XC_GATHER")
      
      do isp = 1, spin_channels
        call distributed_allgather(distribution, dedd(:, isp))
      end do
      
      if(xcs%xc_density_correction == LR_X) then
        call distributed_allgather(distribution, vx)
      end if
      
      if(gga) then
        do idir = 1, der%mesh%sb%dim
          do isp = 1, spin_channels
            call distributed_allgather(distribution, dedgd(:, idir, isp))
          end do
        end do
      end if
            
      call profiling_out(prof_gather)
    end if
  
    call distributed_end(distribution)
  end if

  if(present(deltaxc)) deltaxc = M_ZERO

  ! this has to be done in inverse order
  if( gga) call  gga_process()
  call lda_process()

  ! clean up allocated memory
  call lda_end()
  if(gga .or. xcs%xc_density_correction == LR_X) call  gga_end()

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
    integer :: ii
#ifdef HAVE_LIBXC4
    FLOAT :: parameters(4)
#endif

    PUSH_SUB(xc_get_vxc.gga_init)

    ! allocate variables
    SAFE_ALLOCATE(gdens(1:der%mesh%np, 1:der%mesh%sb%dim, 1:spin_channels))
    gdens = M_ZERO

    SAFE_ALLOCATE(dedgd(1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:spin_channels))
    dedgd = M_ZERO

    do ii = 1, 2
      if(functl(ii)%id == XC_GGA_X_LB) then
#ifdef HAVE_LIBXC4
        parameters(1) = functl(ii)%LB94_modified
        parameters(2) = functl(ii)%LB94_threshold
        parameters(3) = ioniz_pot
        parameters(4) = qtot
        call XC_F90(func_set_ext_params)(functl(ii)%conf, parameters(1))        
#else
        call XC_F90(gga_lb_set_par)(functl(ii)%conf, &
          functl(ii)%LB94_modified, functl(ii)%LB94_threshold, ioniz_pot, qtot)
#endif
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

  end subroutine local_deallocate

end subroutine xc_get_vxc

  pure logical function family_is_gga(family)
    integer, intent(in) :: family

    family_is_gga = bitand(family, XC_FAMILY_GGA + XC_FAMILY_HYB_GGA + XC_FAMILY_MGGA + XC_FAMILY_HYB_MGGA) /= 0
  end function  family_is_gga

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

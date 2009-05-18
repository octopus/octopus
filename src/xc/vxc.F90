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
subroutine xc_get_vxc(gr, xcs, st, rho, ispin, ex, ec, ip, qtot, vxc, vtau)
  type(grid_t),       intent(inout) :: gr
  type(xc_t), target, intent(in)    :: xcs
  type(states_t),     intent(inout) :: st
  FLOAT,              intent(in)    :: rho(:, :)
  integer,            intent(in)    :: ispin
  FLOAT,              intent(inout) :: ex, ec
  FLOAT,              intent(in)    :: ip, qtot
  FLOAT, optional,    intent(inout) :: vxc(:,:)
  FLOAT, optional,    intent(inout) :: vtau(:,:)

  FLOAT :: l_dens(MAX_SPIN), l_dedd(MAX_SPIN), l_sigma(3), l_vsigma(3), l_tau(MAX_SPIN), l_dedtau(MAX_SPIN)
  FLOAT, allocatable :: dens(:,:), dedd(:,:), ex_per_vol(:), ec_per_vol(:)
  FLOAT, allocatable :: gdens(:,:,:), dedgd(:,:,:)
  FLOAT, allocatable :: tau(:,:)

  integer :: jj, ixc, spin_channels
  FLOAT   :: e, r
  logical :: gga, mgga

  type(xc_functl_t), pointer :: functl(:)

  call push_sub('vxc.xc_get_vxc')

  if(ispin == UNPOLARIZED) then
    functl => xcs%functl(:, 1)
  else
    functl => xcs%functl(:, 2)
  end if

  ! is there anything to do ?
  jj = XC_FAMILY_LDA + XC_FAMILY_GGA + XC_FAMILY_HYB_GGA + XC_FAMILY_MGGA
  if(iand(xcs%family, jj) == 0) then
    call pop_sub()
    return
  end if

  call profiling_in(C_PROFILING_XC_LOCAL)

  ! initialize a couple of handy variables
  gga           = iand(xcs%family, XC_FAMILY_GGA + XC_FAMILY_HYB_GGA).ne.0
  mgga          = iand(xcs%family, XC_FAMILY_MGGA).ne.0

  ! This is a bit ugly (why functl(1) and not functl(2)?, but for the moment it works.
  spin_channels = functl(1)%spin_channels

  call lda_init()
  if(gga.or.mgga) call  gga_init()
  if(       mgga) call mgga_init()

  space_loop: do jj = 1, gr%mesh%np

    ! make a local copy with the correct memory order
    l_dens(1:spin_channels) = dens(jj, 1:spin_channels)
    if(gga.or.mgga) then
      l_sigma(1) = sum(gdens(jj, 1:gr%mesh%sb%dim, 1)*gdens(jj, 1:gr%mesh%sb%dim, 1))
      if(ispin /= UNPOLARIZED) then
        l_sigma(2) = sum(gdens(jj, 1:gr%mesh%sb%dim, 1)*gdens(jj, 1:gr%mesh%sb%dim, 2))
        l_sigma(3) = sum(gdens(jj, 1:gr%mesh%sb%dim, 2)*gdens(jj, 1:gr%mesh%sb%dim, 2))
      end if
    end if
    if(mgga) l_tau(1:spin_channels) = tau(jj, 1:spin_channels)

    ! Calculate the potential/gradient density in local reference frame.
    functl_loop: do ixc = 1, 2

      if(.not.present(vxc)) then ! get only the xc energy
        select case(functl(ixc)%family)
        case(XC_FAMILY_LDA)
          call XC_F90(lda_exc)(functl(ixc)%conf, l_dens(1), e)
        case(XC_FAMILY_GGA)
          call XC_F90(gga_exc)(functl(ixc)%conf, l_dens(1), l_sigma(1), e)
        case(XC_FAMILY_HYB_GGA)
          message(1) = 'Hyb-GGAs are currently disabled.'
          call write_fatal(1)
          !call XC_F90(hyb_gga_exc)(functl(ixc)%conf, l_dens(1), l_sigma(1), e)
        case(XC_FAMILY_MGGA)
          call XC_F90(mgga_exc)(functl(ixc)%conf, l_dens(1), l_sigma(1), l_tau(1), e)

        case default
          cycle
        end select

      else ! we also get the xc potential
        select case(functl(ixc)%family)
        case(XC_FAMILY_LDA)
          call XC_F90(lda_vxc)(functl(ixc)%conf, l_dens(1), e, l_dedd(1))
          
        case(XC_FAMILY_GGA)
          if(functl(ixc)%id == XC_GGA_XC_LB) then
            call mesh_r(gr%mesh, jj, r)
            call XC_F90(gga_lb_modified)(functl(ixc)%conf, l_dens(1), l_sigma(1), &
              r, l_dedd(1))
            e = M_ZERO; l_vsigma = M_ZERO
          else
            call XC_F90(gga_vxc)(functl(ixc)%conf, l_dens(1), l_sigma(1), &
              e, l_dedd(1), l_vsigma(1))
          end if
          
        case(XC_FAMILY_HYB_GGA)
          call XC_F90(hyb_gga_vxc)(functl(ixc)%conf, l_dens(1), l_sigma(1), &
            e, l_dedd(1), l_vsigma(1))

        case(XC_FAMILY_MGGA)
          call XC_F90(mgga_vxc)(functl(ixc)%conf, l_dens(1), l_sigma(1), l_tau(1), &
            e, l_dedd(1), l_vsigma(1), l_dedtau(1))

        case default
          cycle
        end select

      end if

      if(functl(ixc)%type == XC_EXCHANGE) then
        ex_per_vol(jj) = ex_per_vol(jj) + sum(l_dens(1:spin_channels)) * e
      else
        ec_per_vol(jj) = ec_per_vol(jj) + sum(l_dens(1:spin_channels)) * e
      end if

      ! store results
      if(present(vxc)) then
        dedd(jj,1:spin_channels) = dedd(jj,1:spin_channels) + l_dedd(1:spin_channels)

        if((functl(ixc)%family == XC_FAMILY_GGA).or.(functl(ixc)%family == XC_FAMILY_MGGA)) then
          dedgd(jj,:,1) = dedgd(jj,:,1) + M_TWO*l_vsigma(1)*gdens(jj,:,1)
          if(ispin /= UNPOLARIZED) then
            dedgd(jj,:,1) = dedgd(jj,:,1) + l_vsigma(2)*gdens(jj,:,2)
            dedgd(jj,:,2) = dedgd(jj,:,2) + M_TWO*l_vsigma(3)*gdens(jj,:,2) + l_vsigma(2)*gdens(jj,:,1)
          end if
        end if

        if(functl(ixc)%family == XC_FAMILY_MGGA) then
          vtau(jj, 1:spin_channels) = vtau(jj, 1:spin_channels) + l_dedtau(1:spin_channels)
        end if
      end if

    end do functl_loop
  end do space_loop

  ! this has to be done in inverse order
  if(present(vxc)) then
    if(       mgga) call mgga_process()
    if(gga.or.mgga) call  gga_process()
    call  lda_process()
  end if

  ! integrate eneries per unit volume
  ex = dmf_integrate(gr%mesh, ex_per_vol)
  ec = dmf_integrate(gr%mesh, ec_per_vol)

  ! clean up allocated memory
  call  lda_end()
  if(gga.or.mgga) call  gga_end()
  if(       mgga) call mgga_end()

  call pop_sub()
  call profiling_out(C_PROFILING_XC_LOCAL)

contains

  ! ---------------------------------------------------------
  ! Takes care of the initialization of the LDA part of the functionals
  !   *) allocates density and dedd, and their local variants
  !   *) calculates the density taking into account nlcc and non-collinear spin
  subroutine lda_init()
    integer :: ii
    FLOAT   :: d(2), dtot, dpol

    call push_sub('vxc.xc_get_vxc.lda_init')

    ! allocate some general arrays
    SAFE_ALLOCATE(dens(1:gr%mesh%np_part, 1:spin_channels))
    SAFE_ALLOCATE(ex_per_vol(1:gr%mesh%np))
    SAFE_ALLOCATE(ec_per_vol(1:gr%mesh%np))

    dens       = M_ZERO
    ex_per_vol = M_ZERO
    ec_per_vol = M_ZERO

    if(present(vxc)) then
      SAFE_ALLOCATE(dedd(1:gr%mesh%np_part, 1:spin_channels))
      dedd       = M_ZERO
    end if

    do ii = 1, gr%mesh%np
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
    
    call pop_sub()
  end subroutine lda_init


  ! ---------------------------------------------------------
  ! deallocate variables allocated in lda_init
  subroutine lda_end()
    SAFE_DEALLOCATE_A(dens)
    SAFE_DEALLOCATE_A(ex_per_vol)
    SAFE_DEALLOCATE_A(ec_per_vol)
    if(present(vxc)) then
      SAFE_DEALLOCATE_A(dedd)
    end if
  end subroutine lda_end


  ! ---------------------------------------------------------
  ! calculates the LDA part of vxc, taking into account non-collinear spin
  subroutine lda_process()
    integer :: i
    FLOAT :: d(2), dpol, vpol

    call push_sub('vxc.xc_get_vxc.lda_process')

    if(ispin == SPINORS) then
      ! rotate back (do not need the rotation matrix for this).
      do i = 1, gr%mesh%np
        d(1:spin_channels) = rho(i, 1:spin_channels)

        dpol = sqrt((d(1) - d(2))**2 + &
          M_FOUR*(rho(i, 3)**2 + rho(i, 4)**2))
        vpol = (dedd(i, 1) - dedd(i, 2))*(d(1) - d(2))/(dpol + tiny)

        vxc(i, 1) = vxc(i, 1) + M_HALF*(dedd(i, 1) + dedd(i, 2) + vpol)
        vxc(i, 2) = vxc(i, 2) + M_HALF*(dedd(i, 1) + dedd(i, 2) - vpol)
        vxc(i, 3) = vxc(i, 3) + (dedd(i, 1) - dedd(i, 2))*rho(i, 3)/(dpol + tiny)
        vxc(i, 4) = vxc(i, 4) + (dedd(i, 1) - dedd(i, 2))*rho(i, 4)/(dpol + tiny)
      end do
    elseif(ispin == SPIN_POLARIZED) then
      call lalg_axpy(gr%mesh%np, M_ONE, dedd(:, 1), vxc(:, 1))
      call lalg_axpy(gr%mesh%np, M_ONE, dedd(:, 2), vxc(:, 2))
    else
      call lalg_axpy(gr%mesh%np, M_ONE, dedd(:, 1), vxc(:, 1))
    end if

    call pop_sub()
  end subroutine lda_process


  ! ---------------------------------------------------------
  ! initialize GGAs
  !   *) allocates gradient of the density (gdens), dedgd, and its local variants
  !   *) calculates the gradient of the density
  subroutine gga_init()
    integer :: ii

    call push_sub('vxc.xc_get_vxc.gga_init')

    ! allocate variables
    SAFE_ALLOCATE(gdens(1:gr%mesh%np, 1:3, 1:spin_channels))
    gdens = M_ZERO

    if(present(vxc)) then
      SAFE_ALLOCATE(dedgd(1:gr%mesh%np_part, 1:3, 1:spin_channels))
      dedgd = M_ZERO
    end if

    ! get gradient of the density
    do ii = 1, spin_channels
      call dderivatives_grad(gr%der, dens(:, ii), gdens(:, :, ii))
    end do

    do ii = 1, 2
      if(functl(ii)%id == XC_GGA_XC_LB) then
        call XC_F90(gga_lb_set_par)(functl(ii)%conf, &
          functl(ii)%LB94_modified, functl(ii)%LB94_threshold, ip, qtot)
      end if
    end do

    call pop_sub()
  end subroutine gga_init


  ! ---------------------------------------------------------
  ! cleans up memory allocated in gga_init
  subroutine gga_end()
    SAFE_DEALLOCATE_A(gdens)
    if(present(vxc)) then
      SAFE_DEALLOCATE_A(dedgd)
    end if
  end subroutine gga_end


  ! ---------------------------------------------------------
  ! calculates the GGA contribution to vxc
  subroutine gga_process()
    integer :: i, is
    FLOAT, allocatable :: gf(:,:)

    call push_sub('vxc.xc_get_vxc.gga_process')

    ! subtract the divergence of the functional derivative of Exc with respect to
    ! the gradient of the density.
    SAFE_ALLOCATE(gf(1:gr%mesh%np, 1:1))
    do is = 1, spin_channels
      call dderivatives_div(gr%der, dedgd(:, :, is), gf(:, 1))
      call lalg_axpy(gr%mesh%np, -M_ONE, gf(:,1), dedd(:, is))
    end do
    SAFE_DEALLOCATE_A(gf)

    ! If LB94, we can calculate an approximation to the energy from
    ! Levy-Perdew relation PRA 32, 2010 (1985)
    if(functl(1)%id == XC_GGA_XC_LB) then
      SAFE_ALLOCATE(gf(1:gr%mesh%np, 1:3))

      do is = 1, spin_channels
        call dderivatives_grad(gr%der, dedd(:, is), gf(:,:))
        do i = 1, gr%mesh%np
          ex_per_vol(i) = ex_per_vol(i) - dens(i, is) * sum(gr%mesh%x(i,:)*gf(i,:))
        end do
      end do

      SAFE_DEALLOCATE_A(gf)
    end if

    call pop_sub()
  end subroutine gga_process


  ! ---------------------------------------------------------
  ! initialize meta-GGAs
  !   *) allocate the kinetic energy density, dedtau, and local variants
  !   *) calculates tau either from a GEA or from the orbitals
  subroutine mgga_init()
    SAFE_ALLOCATE(tau(1:gr%mesh%np, 1:spin_channels))

    ! calculate tau
    call states_calc_tau_jp_gn(gr, st, tau=tau)

  end subroutine mgga_init


  ! ---------------------------------------------------------
  ! clean up memory allocates in mgga_init
  subroutine mgga_end()
    SAFE_DEALLOCATE_A(tau)
  end subroutine mgga_end


  ! ---------------------------------------------------------
  ! calculate the mgga contribution to vxc
  subroutine mgga_process()
  end subroutine mgga_process

end subroutine xc_get_vxc

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

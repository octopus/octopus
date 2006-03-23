!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
subroutine xc_get_vxc(gr, xcs, rho, ispin, vxc, ex, ec, ip, qtot)
  type(grid_t),       intent(inout) :: gr
  type(xc_t), target, intent(in)    :: xcs
  FLOAT,              intent(in)    :: rho(:, :)
  integer,            intent(in)    :: ispin
  FLOAT,              intent(inout) :: vxc(:,:), ex, ec
  FLOAT,              intent(in)    :: ip, qtot

  FLOAT, allocatable :: dens(:,:), dedd(:,:), l_dens(:), l_dedd(:), ex_per_vol(:), ec_per_vol(:)
  FLOAT, allocatable :: gdens(:,:,:), dedgd(:,:,:), l_gdens(:,:), l_dedgd(:,:)
  FLOAT, allocatable :: tau(:,:), dedtau(:,:), l_tau(:), l_dedtau(:)

  integer :: i, ixc, spin_channels
  FLOAT   :: e, r
  logical :: gga, mgga

  type(xc_functl_t), pointer :: functl(:)

  if(ispin == UNPOLARIZED) then
    functl => xcs%functl(:, 1)
  else
    functl => xcs%functl(:, 2)
  end if

  ! is there anything to do ?
  if(iand(xcs%family, XC_FAMILY_LDA + XC_FAMILY_GGA + XC_FAMILY_MGGA) == 0) return

  ! really start
  call push_sub('xc_vxc.xc_get_vxc')

  ! initialize a couple of handy variables
  gga           = iand(xcs%family, XC_FAMILY_GGA).ne.0
  mgga          = iand(xcs%family, XC_FAMILY_MGGA).ne.0

  ! This is a bit ugly (why functl(1) and not functl(2)?, but for the moment it works.
  spin_channels = functl(1)%spin_channels

  call  lda_init()
  if(gga.or.mgga) call  gga_init()
  if(       mgga) call mgga_init()

  space_loop: do i = 1, NP

    ! make a local copy with the correct memory order
    l_dens (:)   = dens (i, :)
    if( gga.or.mgga) l_gdens(:,:) = gdens(i, :,:)
    if(        mgga) l_tau  (:)   = tau  (i, :)

    ! Calculate the potential/gradient density in local reference frame.
    functl_loop: do ixc = 1, 2

      select case(functl(ixc)%family)
      case(XC_FAMILY_LDA)
        call xc_lda(functl(ixc)%conf, l_dens(1), e, l_dedd(1))

      case(XC_FAMILY_GGA)
        if(functl(ixc)%id == XC_GGA_XC_LB) then
          call mesh_r(gr%m, i, r)
          call xc_gga_lb(functl(ixc)%conf, l_dens(1), l_gdens(1,1), &
            r, ip, qtot, l_dedd(1))

          e       = M_ZERO
          l_dedgd = M_ZERO
        else
          call xc_gga(functl(ixc)%conf, l_dens(1), l_gdens(1,1), &
            e, l_dedd(1), l_dedgd(1,1))
        end if

      case(XC_FAMILY_MGGA)
        call xc_mgga(functl(ixc)%conf, l_dens(1), l_gdens(1,1), l_tau(1), &
          e, l_dedd(1), l_dedgd(1,1), l_dedtau(1))

      case default
        cycle
      end select

      if(functl(ixc)%id==XC_LDA_X.or.functl(ixc)%id==XC_GGA_X_PBE.or.&
        functl(ixc)%id==XC_MGGA_X_TPSS) then
        ex_per_vol(i) = ex_per_vol(i) + sum(l_dens(:)) * e
      else
        ec_per_vol(i) = ec_per_vol(i) + sum(l_dens(:)) * e
      end if

      ! store results
      dedd(i,:) = dedd(i,:) + l_dedd(:)

      if(functl(ixc)%family==XC_FAMILY_GGA) then
        dedgd(i,:,:) = dedgd(i,:,:) + l_dedgd(:,:)
      end if

      if(functl(ixc)%family==XC_FAMILY_MGGA) then
        dedgd (i,:,:) = dedgd (i,:,:) + l_dedgd(:,:)
        dedtau(i,:)   = dedtau(i,:)   + l_dedtau(:)
      end if

    end do functl_loop
  end do space_loop

  ! this has to be done in inverse order
  if(       mgga) call mgga_process()
  if(gga.or.mgga) call  gga_process()
  call  lda_process()

  ! integrate eneries per unit volume
  ex = dmf_integrate(gr%m, ex_per_vol)
  ec = dmf_integrate(gr%m, ec_per_vol)

  ! clean up allocated memory
  call  lda_end()
  if(gga.or.mgga) call  gga_end()
  if(       mgga) call mgga_end()

  call pop_sub()

contains

  ! ---------------------------------------------------------
  ! Takes care of the initialization of the LDA part of the functionals
  !   *) allocates density and dedd, and their local variants
  !   *) calculates the density taking into account nlcc and non-collinear spin
  subroutine lda_init()
    integer :: i
    FLOAT   :: d(2), f, dtot, dpol

    ! allocate some general arrays
    ALLOCATE(dens(NP_PART, spin_channels), NP_PART*spin_channels)
    ALLOCATE(dedd(NP_PART, spin_channels), NP_PART*spin_channels)
    ALLOCATE(ex_per_vol(NP), NP)
    ALLOCATE(ec_per_vol(NP), NP)
    ALLOCATE(l_dens(spin_channels), spin_channels)
    ALLOCATE(l_dedd(spin_channels), spin_channels)
    dens       = M_ZERO
    dedd       = M_ZERO
    ex_per_vol = M_ZERO
    ec_per_vol = M_ZERO

    ! get the density
    f = M_ONE/real(spin_channels, PRECISION)
    do i = 1, NP
      d(1:spin_channels) = rho(i, 1:spin_channels)

      select case(ispin)
      case(UNPOLARIZED)
        dens(i, 1) = max(d(1), M_ZERO)
      case(SPIN_POLARIZED)
        dens(i, 1) = max(d(1), M_ZERO)
        dens(i, 2) = max(d(2), M_ZERO)
      case(SPINORS)
        dtot = d(1) + d(2)
        dpol = sqrt((d(1) - d(2))**2 + &
          M_FOUR*(rho(i, 3)**2 + rho(i, 4)**2))
        dens(i, 1) = max(M_HALF*(dtot + dpol), M_ZERO)
        dens(i, 2) = max(M_HALF*(dtot - dpol), M_ZERO)
      end select
    end do

  end subroutine lda_init


  ! ---------------------------------------------------------
  ! deallocates variables allocated in lda_init
  subroutine lda_end()
    deallocate(dens, dedd, ex_per_vol, ec_per_vol, l_dens, l_dedd)
  end subroutine lda_end


  ! ---------------------------------------------------------
  ! calculates the LDA part of vxc, taking into account non-collinear spin
  subroutine lda_process()
    integer :: i
    FLOAT :: d(2), f, dtot, dpol, vpol

    f = M_ONE/real(spin_channels, PRECISION)
    if(ispin == SPINORS) then
      ! rotate back (do not need the rotation matrix for this).
      do i = 1, NP
        d(1:spin_channels) = rho(i, 1:spin_channels)

        dtot = d(1) + d(2)
        dpol = sqrt((d(1) - d(2))**2 + &
          M_FOUR*(rho(i, 3)**2 + rho(i, 4)**2))
        vpol = (dedd(i, 1) - dedd(i, 2))*(d(1) - d(2))/(dpol + tiny)

        vxc(i, 1) = vxc(i, 1) + M_HALF*(dedd(i, 1) + dedd(i, 2) + vpol)
        vxc(i, 2) = vxc(i, 2) + M_HALF*(dedd(i, 1) + dedd(i, 2) - vpol)
        vxc(i, 3) = vxc(i, 3) + (dedd(i, 1) - dedd(i, 2))*rho(i, 3)/(dpol + tiny)
        vxc(i, 4) = vxc(i, 4) + (dedd(i, 1) - dedd(i, 2))*rho(i, 4)/(dpol + tiny)
      end do
    elseif(ispin == SPIN_POLARIZED) then
      vxc(1:NP, 1:2) = vxc(1:NP, 1:2) + dedd(1:NP, 1:2)
    else
      vxc(1:NP, 1) = vxc(1:NP, 1) + dedd(1:NP, 1)
    end if

  end subroutine lda_process


  ! ---------------------------------------------------------
  ! initialize GGAs
  !   *) allocates gradient of the density (gdens), dedgd, and its local variants
  !   *) calculates the gradient of the density
  subroutine gga_init()
    ! allocate variables
    ALLOCATE(gdens(NP,      3, spin_channels), NP     *3*spin_channels)
    ALLOCATE(dedgd(NP_PART, 3, spin_channels), NP_PART*3*spin_channels)
    ALLOCATE(l_gdens(3, spin_channels), 3*spin_channels)
    ALLOCATE(l_dedgd(3, spin_channels), 3*spin_channels)
    gdens = M_ZERO
    dedgd = M_ZERO

    ! get gradient of the density
    do i = 1, spin_channels
      call df_gradient(gr%sb, gr%f_der, dens(:,i), gdens(:,:,i))
    end do
  end subroutine gga_init


  ! ---------------------------------------------------------
  ! cleans up memory allocated in gga_init
  subroutine gga_end()
    deallocate(gdens, dedgd, l_gdens, l_dedgd)
  end subroutine gga_end


  ! ---------------------------------------------------------
  ! calculates the GGA contribution to vxc
  subroutine gga_process()
    integer :: i, is
    FLOAT, allocatable :: gf(:,:)

    ! subtract the divergence of the functional derivative of Exc with respect to
    ! the gradient of the density.
    ALLOCATE(gf(NP, 1), NP*1)
    do is = 1, spin_channels
      call df_divergence(gr%sb, gr%f_der, dedgd(:,:,is), gf(:,1))
      call lalg_axpy(NP, -M_ONE, gf(:,1), dedd(:, is))
    end do
    deallocate(gf)

    ! If LB94, we can calculate an approximation to the energy from
    ! Levy-Perdew relation PRA 32, 2010 (1985)
    if(functl(1)%id == XC_GGA_XC_LB) then
      ALLOCATE(gf(NP, 3), NP*3)

      do is = 1, spin_channels
        call df_gradient(gr%sb, gr%f_der, dedd(:, is), gf(:,:))
        do i = 1, NP
          ex_per_vol(i) = ex_per_vol(i) - dens(i, is) * sum(gr%m%x(i,:)*gf(i,:))
        end do
      end do

      deallocate(gf)
    end if

  end subroutine gga_process


  ! ---------------------------------------------------------
  ! initialize meta-GGAs
  !   *) allocate the kinetic energy density, dedtau, and local variants
  !   *) calculates tau either from a GEA or from the orbitals
  subroutine mgga_init()
    integer :: i, is
    FLOAT   :: f, d
    FLOAT, allocatable :: n2dens(:)

    ALLOCATE(tau   (NP,      spin_channels), NP     *spin_channels)
    ALLOCATE(dedtau(NP_PART, spin_channels), NP_PART*spin_channels)
    ALLOCATE(l_tau   (spin_channels), spin_channels)
    ALLOCATE(l_dedtau(spin_channels), spin_channels)
    tau    = M_ZERO
    dedtau = M_ZERO

    ASSERT(xcs%mGGA_implementation==1.or.xcs%mGGA_implementation==2)

    ! calculate tau
    select case(xcs%mGGA_implementation)
    case (1)  ! GEA implementation
      ALLOCATE(n2dens(NP), NP)
      f = CNST(3.0)/CNST(10.0) * (M_SIX*M_PI*M_PI)**M_TWOTHIRD

      do is = 1, spin_channels
        call df_laplacian(gr%sb, gr%f_der, dens(:,is), n2dens(:))

        do i = 1, NP
          d          = max(dens(i, is), CNST(1e-14))
          tau(i, is) = f * d**(M_FIVE/M_THREE) + &
            sum(gdens(i, :, is)**2)/(CNST(72.0)*d) + &
            n2dens(i)/M_SIX
        end do
      end do

      deallocate(n2dens)

    case(2) ! OEP implementation
      message(1) = 'Error: mgga_init: not yet implemented'
      call write_fatal(1)
    end select

  end subroutine mgga_init


  ! ---------------------------------------------------------
  ! clean up memory allocates in mgga_init
  subroutine mgga_end()
    deallocate(tau, dedtau, l_tau, l_dedtau)
  end subroutine mgga_end


  ! ---------------------------------------------------------
  ! calculate the mgga contribution to vxc
  subroutine mgga_process()
    integer :: i, is
    FLOAT   :: d, f
    FLOAT, allocatable :: gf(:)

    ASSERT(xcs%mGGA_implementation==1.or.xcs%mGGA_implementation==2)

    ! calculate tau
    select case(xcs%mGGA_implementation)
    case (1) ! GEA implementation
      f = CNST(3.0)/CNST(10.0) * (M_SIX*M_PI*M_PI)**M_TWOTHIRD
      ALLOCATE(gf(NP), NP)

      do is = 1, spin_channels
        call df_laplacian(gr%sb, gr%f_der, dedtau(:,is), gf(:))

        do i = 1, NP
          d = max(dens(i, is), CNST(1e-14))
          dedd(i, is) = dedd(i, is) + dedtau(i, is) * &
            (f*d**M_TWOTHIRD - sum(gdens(i, :, is)**2)/(CNST(72.0)*d*d))

          ! add the laplacian of the functional derivative of Exc with respect to tau
          dedd(i, is) = dedd(i, is) + dedtau(i, is) * &
            gf(i)/M_SIX

          dedgd(i, :, is) = dedgd(i, :, is) + dedtau(i, is) * &
            M_TWO*gdens(i, :, is)/(CNST(72.0)*d)
        end do
      end do

      deallocate(gf)

    case (2) ! OEP implementation
      message(1) = 'Error: mgga_process: not yet implemented'
      call write_fatal(1)
    end select
  end subroutine mgga_process

end subroutine xc_get_vxc

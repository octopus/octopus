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

subroutine xc_get_vxc(xcs, m, f_der, rho, ispin, vxc, ex, ec, ip, qtot)
  type(xc_type), target, intent(in)    :: xcs
  type(mesh_type),       intent(in)    :: m
  type(f_der_type),      intent(inout) :: f_der
  FLOAT,                 intent(in)    :: rho(:, :)
  integer,               intent(in)    :: ispin
  FLOAT,                 intent(inout) :: vxc(:,:), ex, ec
  FLOAT,                 intent(in)    :: ip, qtot
  
  FLOAT, allocatable :: dens(:,:), dedd(:,:), l_dens(:), l_dedd(:)
  FLOAT, allocatable :: gdens(:,:,:), dedgd(:,:,:), l_gdens(:,:), l_dedgd(:,:)
  FLOAT, allocatable :: tau(:,:), dedtau(:,:), l_tau(:), l_dedtau(:)

  integer :: i, is, ixc, spin_channels
  FLOAT   :: e, dpol, dtot, vpol, r
  logical :: gga, mgga

  type(xc_functl_type), pointer :: functl(:)

  if(ispin == UNPOLARIZED) then
    functl => xcs%functl(:, 1)
  else
    functl => xcs%functl(:, 2)
  end if

  ! is there anything to do ?
  if(.not.( &
     any(functl(:)%family==XC_FAMILY_LDA).or. &
     any(functl(:)%family==XC_FAMILY_GGA).or. &
     any(functl(:)%family==XC_FAMILY_MGGA) )) return
  
  ! really start
  call push_sub('xc_get_vxc')

  ! initialize a couple of handy variables
  gga           = any(functl(:)%family == XC_FAMILY_GGA)
  mgga          = any(functl(:)%family == XC_FAMILY_MGGA)
  ! This is a bit ugly (why functl(1) and not functl(2)?, but for the moment it works.
  spin_channels = functl(1)%spin_channels

                  call  lda_init()
  if(gga.or.mgga) call  gga_init()
  if(       mgga) call mgga_init()

  space_loop: do i = 1, m%np

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
          call mesh_r(m, i, r)
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
        ex = ex + sum(l_dens(:)) * e * m%vol_pp(i)
      else
        ec = ec + sum(l_dens(:)) * e * m%vol_pp(i)
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


  ! clean up allocated memory
                  call  lda_end()
  if(gga.or.mgga) call  gga_end()
  if(       mgga) call mgga_end()

  call pop_sub()

contains
  ! Takes care of the initialization of the LDA part of the functionals
  !   *) allocates dens(ity) and dedd, and their local variants
  !   *) calculates the density taking into account nlcc and non-collinear spin
  subroutine lda_init()
    integer :: i
    FLOAT   :: d(spin_channels), f, dtot, dpol
    
    ! allocate some general arrays
    allocate(dens(m%np, spin_channels), dedd(m%np, spin_channels))
    allocate(l_dens(spin_channels), l_dedd(spin_channels))
    dedd = M_ZERO

    ! get the density
    f = M_ONE/real(spin_channels, PRECISION)
    do i = 1, m%np
      d(:) = rho(i, :)
      
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


  ! deallocates variables allocated in lda_init
  subroutine lda_end()
    deallocate(dens, dedd, l_dens, l_dedd)
  end subroutine lda_end

  
  ! calculates the LDA part of vxc, taking into account non-collinear spin
  subroutine lda_process()
    integer :: i
    FLOAT :: d(spin_channels), f, dtot, dpol, vpol

    f = M_ONE/real(spin_channels, PRECISION)
    if(ispin == SPINORS) then
      ! rotate back (do not need the rotation matrix for this).
      do i = 1, m%np
        d(:) = rho(i, :)
        
        dtot = d(1) + d(2)
        dpol = sqrt((d(1) - d(2))**2 + &
           M_FOUR*(rho(i, 3)**2 + rho(i, 4)**2))
        vpol = (dedd(i, 1) - dedd(i, 2))*(d(1) - d(2))/(dpol + tiny)
        
        vxc(i, 1) = vxc(i, 1) + M_HALF*(dedd(i, 1) + dedd(i, 2) + vpol)
        vxc(i, 2) = vxc(i, 2) + M_HALF*(dedd(i, 1) + dedd(i, 2) - vpol)
        vxc(i, 3) = vxc(i, 3) + (dedd(i, 1) - dedd(i, 2))*rho(i, 3)/(dpol + tiny)
        vxc(i, 4) = vxc(i, 4) + (dedd(i, 1) - dedd(i, 2))*rho(i, 4)/(dpol + tiny)
      end do
    else
      vxc = vxc + dedd
    end if
    
  end subroutine lda_process


  ! initialize GGAs
  !   *) allocates gradient of the density (gdens), dedgd, and its local variants
  !   *) calculates the gradient of the density
  subroutine gga_init()
    ! allocate variables
    allocate(gdens(m%np, 3, spin_channels), dedgd(m%np, 3, spin_channels))
    allocate(l_gdens(3, spin_channels), l_dedgd(3, spin_channels))
    gdens = M_ZERO
    dedgd = M_ZERO

    ! get gradient of the density
    do i = 1, spin_channels
      call df_gradient(f_der, dens(:,i), gdens(:,:,i))
    end do
  end subroutine gga_init


  ! cleans up memory allocated in gga_init
  subroutine gga_end()
    deallocate(gdens, dedgd, l_gdens, l_dedgd)
  end subroutine gga_end


  ! calculates the GGA contribution to vxc
  subroutine gga_process()
    integer :: i, is
    FLOAT, allocatable :: gf(:,:)

    ! subtract the divergence of the functional derivative of Exc with respect to
    ! the gradient of the density.
    allocate(gf(m%np, 1))
    do is = 1, spin_channels
      call df_divergence(f_der, dedgd(:,:,is), gf(:,1))
      call lalg_axpy(m%np, -M_ONE, gf(:,1), dedd(:, is))
    end do
    deallocate(gf)
    
    ! If LB94, we can calculate an approximation to the energy from 
    ! Levy-Perdew relation PRA 32, 2010 (1985)
    if(functl(1)%id == XC_GGA_XC_LB) then
      allocate(gf(m%np, 3))

      do is = 1, spin_channels
        call df_gradient(f_der, dedd(:, is), gf(:,:))
        do i = 1, m%np
          ex = ex - dens(i, is) * sum(m%x(i,:)*gf(i,:)) * m%vol_pp(i)
        end do
      end do

      deallocate(gf)
    end if

  end subroutine gga_process
  

  ! initialize meta-GGAs
  !   *) allocate the kinetic energy density, dedtau, and local variants
  !   *) calculates tau either from a GEA or from the orbitals
  subroutine mgga_init()
    integer :: i, is
    FLOAT   :: f, d
    FLOAT, allocatable :: n2dens(:)

    allocate(tau(m%np, spin_channels), dedtau(m%np, spin_channels))
    allocate(l_tau(spin_channels), l_dedtau(spin_channels))
    tau    = M_ZERO
    dedtau = M_ZERO

    ASSERT(xcs%mGGA_implementation==1.or.xcs%mGGA_implementation==2)

    ! calculate tau
    select case(xcs%mGGA_implementation)
    case (1)  ! GEA implementation
      allocate(n2dens(m%np))
      f = CNST(3.0)/CNST(10.0) * (M_SIX*M_PI*M_PI)**M_TWOTHIRD

      do is = 1, spin_channels
        call df_laplacian(f_der, dens(:,is), n2dens(:))

        do i = 1, m%np
          d          = max(dens(i, is), CNST(1e-14))
          tau(i, is) = f * d**(M_FIVE/M_THREE) + &
             sum(gdens(i, :, is)**2)/(CNST(72.0)*d) + &
             n2dens(i)/M_SIX
        end do
      end do
      
      deallocate(n2dens)

    case(2) ! OEP implementation
      stop 'not yet implemented'
    end select

  end subroutine mgga_init


  ! clean up memory allocates in mgga_init
  subroutine mgga_end()
    deallocate(tau, dedtau, l_tau, l_dedtau)
  end subroutine mgga_end


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
        allocate(gf(m%np))

        do is = 1, spin_channels
          call df_laplacian(f_der, dedtau(:,is), gf(:))

          do i = 1, m%np
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
        stop 'not yet implemented'
      end select
  end subroutine mgga_process

end subroutine xc_get_vxc


subroutine xc_get_vxc_and_axc(xcs, m, f_der, rho, j, ispin, vxc, axc, ex, ec, exc_j, ip, qtot)
  type(xc_type), target, intent(in)    :: xcs
  type(mesh_type),       intent(in)    :: m
  type(f_der_type),      intent(inout) :: f_der
  FLOAT,                 intent(in)    :: rho(:, :), j(:,:,:)
  integer,               intent(in)    :: ispin
  FLOAT,                 intent(inout) :: vxc(:,:), axc(:,:,:)
  FLOAT,                 intent(inout) :: ex, ec, exc_j
  FLOAT,                 intent(in)    :: ip, qtot

  integer :: spin_channels, family, i, is, id
  FLOAT   :: e
  FLOAT, allocatable :: v(:,:,:), f(:,:,:), dedd(:,:), dedv(:,:,:), tmp(:,:)
  FLOAT, allocatable :: l_dens(:), l_v(:,:), l_dedd(:), l_dedv(:,:)

  call push_sub('xc_get_vxc_and_axc')

  !xc energy and potential in the absence of external magnetic fields
  call xc_get_vxc(xcs, m, f_der, rho, ispin, vxc, ex, ec, ip, qtot)

  spin_channels = xcs%j_functl%spin_channels
  family = xcs%j_functl%family

  ! is there anything else to do ?
  if(family /= XC_FAMILY_LCA) then
    call pop_sub()
    return
  end if

  !allocate memory
  allocate(v(m%np, conf%dim, spin_channels), f(m%np, conf%dim, spin_channels))
  allocate(dedd(m%np, spin_channels), dedv(m%np, conf%dim, spin_channels))

  !Compute j/rho and the vorticity
  do is = 1, spin_channels
    do id = 1, conf%dim
      f(:, id, is)  = j(:, id, is)/rho(:, is)
    end do
    call df_curl(f_der, f(:,:,is), v(:,:,is))
  end do

  !
  allocate(l_dens(spin_channels), l_v(conf%dim, spin_channels))
  allocate(l_dedd(spin_channels), l_dedv(conf%dim, spin_channels))
  space_loop: do i = 1, m%np
    ! make a local copy with the correct memory order
    l_dens (:) = rho(i, :)
    l_v(:,:)   = v(i, :,:)
    
    ! Calculate the potential density in local reference frame.
    select case(family)
    case(XC_FAMILY_LCA)
      call xc_lca(xcs%j_functl%conf, l_dens(1), l_v(1,1), &
                  e, l_dedd(1), l_dedv(1,1))
    end select

    exc_j = exc_j + sum(l_dens(:)) * e * m%vol_pp(i)
    
    ! store results
    dedd(i,:) = dedd(i,:) + l_dedd
    dedv(i,:,:) = dedv(i,:,:) + l_dedv

  end do space_loop
  deallocate(l_dens, l_v, l_dedd, l_dedv)

  !
  vxc = vxc + dedd
  allocate(tmp(m%np, conf%dim))
  do is = 1, spin_channels
    call df_curl(f_der, dedv(:,:,is), tmp)
    do id = 1, conf%dim
      axc(:, id, is) = axc(:, id, is) - tmp(:, id)/rho(:, is)
      vxc(:, is) = vxc(:, is) - axc(:, id, is)*f(:, id, is)
    end do
  end do
  deallocate(tmp)

  !deallocate memory
  deallocate(f, v, dedd, dedv)

  call pop_sub()
end subroutine xc_get_vxc_and_axc

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

subroutine xc_get_vxc(nf, functl, nlcc, m, f_der, st, vxc, ex, ec, ip, qtot)
  integer,              intent(in)    :: nf        ! how many functionals in array
  type(xc_functl_type), intent(in)    :: functl(:) ! functl(nf)
  logical,              intent(in)    :: nlcc
  type(mesh_type),      intent(in)    :: m
  type(f_der_type),     intent(inout) :: f_der
  type(states_type),    intent(in)    :: st
  FLOAT,                intent(inout) :: vxc(:,:), ex, ec
  FLOAT,                intent(in)    :: ip, qtot
  
  FLOAT, allocatable :: dens(:,:), dedd(:,:), l_dens(:), l_dedd(:)
  FLOAT, allocatable :: gdens(:,:,:), dedgd(:,:,:), l_gdens(:,:), l_dedgd(:,:)


  integer :: i, is, ixc, ispin, spin_channels
  FLOAT   :: e, dpol, dtot, vpol, r
  logical :: gga

  ! is there anything to do
  if(.not.(any(functl(:)%family==XC_FAMILY_LDA).or. &
     any(functl(:)%family==XC_FAMILY_GGA))) return
  
  ! really start
  call push_sub('xc_gga')

  ! initialize a couple of handy variables
  gga           = any(functl(:)%family == XC_FAMILY_GGA)
  spin_channels = functl(1)%spin_channels
  ispin         = functl(1)%ispin
  
  allocate(dens(m%np, spin_channels), dedd(m%np, spin_channels))
  allocate(l_dens(spin_channels), l_dedd(spin_channels))
  dedd = M_ZERO

  call get_dens()

  if(gga) then
    allocate(gdens(m%np, 3, spin_channels), dedgd(m%np, 3, spin_channels))
    allocate(l_gdens(3, spin_channels), l_dedgd(3, spin_channels))
    gdens = M_ZERO
    dedgd = M_ZERO

    do i = 1, spin_channels
      call df_gradient(f_der, dens(:,i), gdens(:,:,i))
    end do
  end if

  space_loop: do i = 1, m%np

    ! make a local copy with the correct memory order
    l_dens (:) = dens (i, :)
    if(gga) l_gdens(:,:) = gdens(i, :,:)
    
    ! Calculate the potential/gradient density in local reference frame.
    functl_loop: do ixc = 1, nf

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
      case default
        cycle
      end select

      if(functl(ixc)%id==XC_LDA_X.or.functl(ixc)%id==XC_GGA_X_PBE) then
        ex = ex + sum(l_dens(:)) * e * m%vol_pp(i)
      else
        ec = ec + sum(l_dens(:)) * e * m%vol_pp(i)
      end if

      dedd(i,:) = dedd(i,:) + l_dedd(:)
      if(functl(ixc)%family==XC_FAMILY_GGA) then
        dedgd(i,:,:) = dedgd(i,:,:) + l_dedgd(:,:)
      end if
      
    end do functl_loop

  end do space_loop

  ! We now add substract the divergence of the functional derivative of Exc with respect to
  ! the gradient of the density.
  !   Note: gdens is used as a temporary array
  if(gga) then
    do is = 1, spin_channels
      call df_divergence(f_der, dedgd(:,:,is), gdens(:,1,is))
      call lalg_axpy(m%np, -M_ONE, gdens(:,1,is), dedd(:, is))
    end do
  end if
      
  ! If LB94, we can calculate an approximation to the energy from 
  ! Levy-Perdew relation (PRA 32, 2010 (1985))
  if(functl(1)%id == XC_GGA_XC_LB) then
    do is = 1, spin_channels
      call df_gradient(f_der, dedd(:, is), gdens(:,:,is))
      do i = 1, m%np
        ex = ex - dens(i, is) * sum(m%x(i,:)*gdens(i,:,is)) * m%vol_pp(i)
      end do
    end do
  end if

  ! now rotate to get back vxc
  call get_vxc()

  deallocate(dens, dedd, l_dens, l_dedd)
  if(gga) deallocate(gdens, dedgd, l_gdens, l_dedgd)
  call pop_sub()

contains
  subroutine get_dens()
    integer :: i
    FLOAT   :: d(spin_channels), f, dtot, dpol
    
    f = M_ONE/real(spin_channels, PRECISION)
    do i = 1, m%np
      d(:) = st%rho(i, :)
      
      ! add the non-linear core corrections
      if(nlcc) d(:) = d(:) + f*st%rho_core(i)
      
      select case(ispin)
      case(UNPOLARIZED)
        dens(i, 1) = max(d(1), M_ZERO)
      case(SPIN_POLARIZED)
        dens(i, 1) = max(d(1), M_ZERO)
        dens(i, 2) = max(d(2), M_ZERO)
      case(SPINORS)
        dtot = d(1) + d(2)
        dpol = sqrt((d(1) - d(2))**2 + &
           M_FOUR*(st%rho(i, 3)**2 + st%rho(i, 4)**2))
        dens(i, 1) = max(M_HALF*(dtot + dpol), M_ZERO)
        dens(i, 2) = max(M_HALF*(dtot - dpol), M_ZERO)
      end select
    end do
    
  end subroutine get_dens


  ! rotate back (do not need the rotation matrix for this).
  subroutine get_vxc()
    integer :: i
    FLOAT :: d(spin_channels), f, dtot, dpol, vpol

    f = M_ONE/real(spin_channels, PRECISION)
    if(ispin == SPINORS) then
      do i = 1, m%np
        d(:) = st%rho(i, :)
      
        ! add the non-linear core corrections
        if(nlcc) d(:) = d(:) + f*st%rho_core(i)
        
        dtot = d(1) + d(2)
        dpol = sqrt((d(1) - d(2))**2 + &
           M_FOUR*(st%rho(i, 3)**2 + st%rho(i, 4)**2))
        vpol = (dedd(i, 1) - dedd(i, 2))*(d(1) - d(2))/(dpol + tiny)
        
        vxc(i, 1) = vxc(i, 1) + M_HALF*(dedd(i, 1) + dedd(i, 2) + vpol)
        vxc(i, 2) = vxc(i, 2) + M_HALF*(dedd(i, 1) + dedd(i, 2) - vpol)
        vxc(i, 3) = vxc(i, 3) + (dedd(i, 1) - dedd(i, 2))*st%rho(i, 3)/(dpol + tiny)
        vxc(i, 4) = vxc(i, 4) + (dedd(i, 1) - dedd(i, 2))*st%rho(i, 4)/(dpol + tiny)
      end do
    else
      vxc = vxc + dedd
    end if

  end subroutine get_vxc

end subroutine xc_get_vxc

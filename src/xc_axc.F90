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

subroutine xc_get_vxc_and_axc(xcs, m, f_der, rho, j, ispin, vxc, axc, ex, ec, exc_j, ip, qtot)
  type(xc_type), target, intent(in)    :: xcs
  type(mesh_type),       intent(in)    :: m
  type(f_der_type),      intent(inout) :: f_der
  FLOAT,                 intent(in)    :: rho(:, :), j(:,:,:)
  integer,               intent(in)    :: ispin
  FLOAT,                 intent(inout) :: vxc(:,:), axc(:,:,:)
  FLOAT,                 intent(inout) :: ex, ec, exc_j
  FLOAT,                 intent(in)    :: ip, qtot

  integer :: spin_channels, i, is, id
  FLOAT   :: e
  FLOAT, allocatable :: v(:,:,:), f(:,:,:), dedd(:,:), dedv(:,:,:), tmp(:,:)
  FLOAT, allocatable :: l_dens(:), l_v(:,:), l_dedd(:), l_dedv(:,:)

  ! is there anything to do?
  if(iand(xcs%family, XC_FAMILY_LCA) == 0) return 

  call push_sub('xc_get_vxc_and_axc')

  !xc energy and potential in the absence of external magnetic fields
  call xc_get_vxc(xcs, m, f_der, rho, ispin, vxc, ex, ec, ip, qtot)

  spin_channels = xcs%j_functl%spin_channels

  !allocate memory
  allocate(v(m%np, conf%dim, spin_channels), f(m%np, conf%dim, spin_channels))
  allocate(dedd(m%np, spin_channels), dedv(m%np, conf%dim, spin_channels))

  !Compute j/rho and the vorticity
  do is = 1, spin_channels
    do id = 1, conf%dim
      f(:, id, is) = j(:, id, is)/rho(:, is)
    end do
    call df_curl(f_der, f(:,:,is), v(:,:,is))
  end do

  allocate(l_dens(spin_channels), l_v(conf%dim, spin_channels))
  allocate(l_dedd(spin_channels), l_dedv(conf%dim, spin_channels))
  space_loop: do i = 1, m%np
    ! make a local copy with the correct memory order
    l_dens (:) = rho(i, :)
    l_v(:,:)   = v(i, :,:)
    
    ! Calculate the potential density in local reference frame.
    select case(xcs%j_functl%family)
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

  ! add contributions to vxc and axc
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

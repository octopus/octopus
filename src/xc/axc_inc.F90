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
subroutine xc_get_vxc_and_axc(gr, xcs, st, rho, j, ispin, vxc, axc, ex, ec, exc_j, ip, qtot)
  type(grid_t),       intent(inout) :: gr
  type(xc_t),         intent(in)    :: xcs
  type(states_t),     intent(inout) :: st
  FLOAT,              intent(in)    :: rho(:, :), j(:,:,:)
  integer,            intent(in)    :: ispin
  FLOAT,              intent(inout) :: vxc(:,:), axc(:,:,:)
  FLOAT,              intent(inout) :: ex, ec, exc_j
  FLOAT,              intent(in)    :: ip, qtot

  integer :: spin_channels, i, is, id
  FLOAT   :: e
  FLOAT, allocatable :: v(:,:,:), f(:,:,:), dedd(:,:), dedv(:,:,:), tmp(:,:)
  FLOAT, allocatable :: l_dens(:), l_v(:,:), l_dedd(:), l_dedv(:,:)

  call push_sub('axc_inc.xc_get_vxc_and_axc')

  !xc energy and potential in the absence of external magnetic fields
  call xc_get_vxc(gr, xcs, st, rho, ispin, ex, ec, ip, qtot, vxc=vxc)

  !do we have a current-dependent xc?
  if(iand(xcs%family, XC_FAMILY_LCA) == 0) then
    call pop_sub(); return
  end if

  spin_channels = xcs%j_functl%spin_channels

  !allocate memory
  SAFE_ALLOCATE(v   (1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:spin_channels))
  SAFE_ALLOCATE(f   (1:gr%mesh%np_part, 1:gr%mesh%sb%dim, 1:spin_channels))
  SAFE_ALLOCATE(dedd(1:gr%mesh%np, 1:spin_channels))
  SAFE_ALLOCATE(dedv(1:gr%mesh%np_part, 1:gr%mesh%sb%dim, 1:spin_channels))
  dedd = M_ZERO; dedv = M_ZERO

  !Compute j/rho and the vorticity
  do is = 1, spin_channels
    do id = 1, gr%mesh%sb%dim
      f(1:gr%mesh%np, id, is) = j(1:gr%mesh%np, id, is)/rho(1:gr%mesh%np, is)
    end do
    call dderivatives_curl(gr%der, f(:,:,is), v(:,:,is))
  end do

  SAFE_ALLOCATE(l_dens(1:spin_channels))
  SAFE_ALLOCATE(l_v   (1:gr%mesh%sb%dim, 1:spin_channels))
  SAFE_ALLOCATE(l_dedd(1:spin_channels))
  SAFE_ALLOCATE(l_dedv(1:gr%mesh%sb%dim, 1:spin_channels))
  l_dedd = M_ZERO; l_dedv = M_ZERO
  space_loop: do i = 1, gr%mesh%np
    ! make a local copy with the correct memory order
    l_dens (:) = rho(i, :)
    l_v(:,:)   = v(i, :,:)

    ! Calculate the potential density in local reference frame.
    select case(xcs%j_functl%family)
    case(XC_FAMILY_LCA)
      call XC_F90(lca)(xcs%j_functl%conf, l_dens(1), l_v(1,1), &
        e, l_dedd(1), l_dedv(1,1))
    end select

    exc_j = exc_j + sum(l_dens(:)) * e * gr%mesh%vol_pp(i)

    ! store results
    dedd(i,:) = dedd(i,:) + l_dedd
    dedv(i,:,:) = dedv(i,:,:) + l_dedv

  end do space_loop
  SAFE_DEALLOCATE_A(l_dens)
  SAFE_DEALLOCATE_A(l_v)
  SAFE_DEALLOCATE_A(l_dedd)
  SAFE_DEALLOCATE_A(l_dedv)

  ! add contributions to vxc and axc
  vxc = vxc + dedd
  SAFE_ALLOCATE(tmp(1:gr%mesh%np, 1:gr%mesh%sb%dim))
  do is = 1, spin_channels
    call dderivatives_curl(gr%der, dedv(:, :, is), tmp)

    do id = 1, gr%mesh%sb%dim
      axc(1:gr%mesh%np, id, is) = axc(1:gr%mesh%np, id, is) - tmp(1:gr%mesh%np, id)/rho(1:gr%mesh%np, is)
      vxc(1:gr%mesh%np, is) = vxc(1:gr%mesh%np, is) - axc(1:gr%mesh%np, id, is)*f(1:gr%mesh%np, id, is)
    end do
  end do
  SAFE_DEALLOCATE_A(tmp)

  !deallocate memory
  SAFE_DEALLOCATE_A(f)
  SAFE_DEALLOCATE_A(v)
  SAFE_DEALLOCATE_A(dedd)
  SAFE_DEALLOCATE_A(dedv)

  call pop_sub()
end subroutine xc_get_vxc_and_axc

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

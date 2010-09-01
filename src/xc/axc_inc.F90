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
subroutine xc_get_vxc_and_axc(der, xcs, st, rho, current_density, ispin, vxc, axc, ex, ec, exc_j, ioniz_pot, qtot)
  type(derivatives_t), intent(inout) :: der
  type(xc_t),          intent(in)    :: xcs
  type(states_t),      intent(inout) :: st
  FLOAT,               intent(in)    :: rho(:, :), current_density(:,:,:)
  integer,             intent(in)    :: ispin
  FLOAT,               intent(inout) :: vxc(:,:), axc(:,:,:)
  FLOAT,               intent(inout) :: ex, ec, exc_j
  FLOAT,               intent(in)    :: ioniz_pot, qtot

  integer :: spin_channels, ip, is, id
  FLOAT   :: ee
  FLOAT, allocatable :: vorticity(:,:,:), ff(:,:,:), dedd(:,:), dedv(:,:,:), tmp(:,:)
  FLOAT, allocatable :: l_dens(:), l_vorticity(:,:), l_dedd(:), l_dedv(:,:)

  PUSH_SUB(xc_get_vxc_and_axc)

  !XC energy and potential in the absence of external magnetic fields
  call xc_get_vxc(der, xcs, st, rho, ispin, ex, ec, ioniz_pot, qtot, vxc=vxc)

  !do we have a current-dependent XC?
  if(iand(xcs%family, XC_FAMILY_LCA) == 0) then
    POP_SUB(xc_get_vxc_and_axc)
    return
  end if

  spin_channels = xcs%j_functl%spin_channels

  !allocate memory
  SAFE_ALLOCATE(vorticity(1:der%mesh%np, 1:der%mesh%sb%dim, 1:spin_channels))
  SAFE_ALLOCATE(ff  (1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:spin_channels))
  SAFE_ALLOCATE(dedd(1:der%mesh%np, 1:spin_channels))
  SAFE_ALLOCATE(dedv(1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:spin_channels))
  dedd = M_ZERO
  dedv = M_ZERO

  !Compute j/rho and the vorticity
  do is = 1, spin_channels
    do id = 1, der%mesh%sb%dim
      ff(1:der%mesh%np, id, is) = current_density(1:der%mesh%np, id, is)/rho(1:der%mesh%np, is)
    end do
    call dderivatives_curl(der, ff(:,:,is), vorticity(:,:,is))
  end do

  SAFE_ALLOCATE(l_dens(1:spin_channels))
  SAFE_ALLOCATE(l_vorticity(1:der%mesh%sb%dim, 1:spin_channels))
  SAFE_ALLOCATE(l_dedd(1:spin_channels))
  SAFE_ALLOCATE(l_dedv(1:der%mesh%sb%dim, 1:spin_channels))
  l_dedd = M_ZERO
  l_dedv = M_ZERO

  space_loop: do ip = 1, der%mesh%np
    ! make a local copy with the correct memory order
    l_dens (:) = rho(ip, :)
    l_vorticity(:,:) = vorticity(ip, :,:)

    ! Calculate the potential density in local reference frame.
    select case(xcs%j_functl%family)
    case(XC_FAMILY_LCA)
      call XC_F90(lca)(xcs%j_functl%conf, l_dens(1), l_vorticity(1,1), &
        ee, l_dedd(1), l_dedv(1,1))
    end select

    exc_j = exc_j + sum(l_dens(:)) * ee * der%mesh%vol_pp(ip)

    ! store results
    dedd(ip,:) = dedd(ip,:) + l_dedd
    dedv(ip,:,:) = dedv(ip,:,:) + l_dedv

  end do space_loop
  SAFE_DEALLOCATE_A(l_dens)
  SAFE_DEALLOCATE_A(l_vorticity)
  SAFE_DEALLOCATE_A(l_dedd)
  SAFE_DEALLOCATE_A(l_dedv)

  ! add contributions to vxc and axc
  vxc = vxc + dedd
  SAFE_ALLOCATE(tmp(1:der%mesh%np, 1:der%mesh%sb%dim))
  do is = 1, spin_channels
    call dderivatives_curl(der, dedv(:, :, is), tmp)

    do id = 1, der%mesh%sb%dim
      axc(1:der%mesh%np, id, is) = axc(1:der%mesh%np, id, is) - tmp(1:der%mesh%np, id)/rho(1:der%mesh%np, is)
      vxc(1:der%mesh%np, is) = vxc(1:der%mesh%np, is) - axc(1:der%mesh%np, id, is)*ff(1:der%mesh%np, id, is)
    end do
  end do
  SAFE_DEALLOCATE_A(tmp)

  !deallocate memory
  SAFE_DEALLOCATE_A(ff)
  SAFE_DEALLOCATE_A(vorticity)
  SAFE_DEALLOCATE_A(dedd)
  SAFE_DEALLOCATE_A(dedv)

  POP_SUB(xc_get_vxc_and_axc)
end subroutine xc_get_vxc_and_axc

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

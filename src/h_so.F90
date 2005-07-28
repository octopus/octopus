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

subroutine zso (h, gr, psi, hpsi, ik)
  type(hamiltonian_type), intent(in) :: h
  type(grid_type) :: gr
  R_TYPE, intent(in) :: psi(:, :)
  R_TYPE, intent(inout) :: Hpsi(:, :)
  integer, intent(in) :: ik

  integer :: ivnl
  CMPLX, allocatable :: tpsi(:, :)

  call push_sub('zso')

!!$!!!!I will keep this commented for the moment...
!!$  ASSERT(h%d%dim == 2)
!!$
!!$  allocate(lpsi(gr%m%np, 3, h%d%dim), tpsi(gr%m%np, h%d%dim), lspsi(gr%m%np, h%d%dim), &
!!$           gpsi(gr%m%np, 3, h%d%dim))
!!$  do idim = 1, 2
!!$     call zf_gradient(gr%sb, gr%f_der, psi(:, idim), gpsi(:, :, idim))
!!$     gpsi = -M_zI*gpsi
!!$     do i = 1, gr%m%np
!!$        lpsi(i, 1, idim) = (gr%m%x(i, 2)*gpsi(i, 3, idim) - gr%m%x(i, 3)*gpsi(i, 2, idim))
!!$        lpsi(i, 2, idim) = (gr%m%x(i, 3)*gpsi(i, 1, idim) - gr%m%x(i, 1)*gpsi(i, 3, idim))
!!$        lpsi(i, 3, idim) = (gr%m%x(i, 1)*gpsi(i, 2, idim) - gr%m%x(i, 2)*gpsi(i, 1, idim))
!!$     enddo
!!$  enddo
!!$
!!$  lspsi(:, 1) = M_HALF*(lpsi(:, 3, 1) + lpsi(:, 1, 2) - M_zI*lpsi(:, 2, 2))
!!$  lspsi(:, 2) = M_HALF*(lpsi(:, 1, 1) + M_zI*lpsi(:, 2, 1) - lpsi(:, 3, 2))
!!$
!!$  do idim = 1, h%d%dim
!!$     do ivnl = 1, h%ep%nvnl
!!$        call zproject(gr%m, h%ep%so(ivnl), lspsi(:, idim), hpsi(:, idim), &
!!$                      periodic = simul_box_is_periodic(gr%sb), ik = ik)
!!$     enddo
!!$  enddo
!!$
!!$  deallocate(lpsi, tpsi)
!!$!!!!ENDOFNEW

  ASSERT(h%d%dim == 2)

  allocate(tpsi(gr%m%np, 3))

  do ivnl = 1, h%ep%nvnl
        tpsi = M_z0
        call zproject(gr%m, h%ep%lso(3, ivnl), psi(:, 1), tpsi(:, 1), simul_box_is_periodic(gr%sb), ik)
        tpsi(:, 1) = - M_zI*tpsi(:, 1)
        call zproject(gr%m, h%ep%lso(1, ivnl), psi(:, 2), tpsi(:, 2), simul_box_is_periodic(gr%sb), ik)
        tpsi(:, 2) = - M_zI*tpsi(:, 2)
        call zproject(gr%m, h%ep%lso(2, ivnl), psi(:, 2), tpsi(:, 3), simul_box_is_periodic(gr%sb), ik)
        tpsi(:, 3) = - M_zI*tpsi(:, 3)
        hpsi(:, 1) = hpsi(:, 1) + M_HALF * (tpsi(:, 1) + tpsi(:, 2) + M_zI * tpsi(:, 3))

        tpsi = M_z0
        call zproject(gr%m, h%ep%lso(3, ivnl), psi(:, 2), tpsi(:, 1), simul_box_is_periodic(gr%sb), ik)
        tpsi(:, 1) = - M_zI*tpsi(:, 1)
        call zproject(gr%m, h%ep%lso(1, ivnl), psi(:, 1), tpsi(:, 2), simul_box_is_periodic(gr%sb), ik)
        tpsi(:, 2) = - M_zI*tpsi(:, 2)
        call zproject(gr%m, h%ep%lso(2, ivnl), psi(:, 1), tpsi(:, 3), simul_box_is_periodic(gr%sb), ik)
        tpsi(:, 3) = - M_zI*tpsi(:, 3)
        hpsi(:, 2) = hpsi(:, 2) + M_HALF * (-tpsi(:, 1) + tpsi(:, 2) - M_zI * tpsi(:, 3))
  enddo

  call pop_sub()
end subroutine zso

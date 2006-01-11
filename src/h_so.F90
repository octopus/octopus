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

  CMPLX, allocatable :: tpsi(:, :)

  call push_sub('h_so.zso')

  ASSERT(h%d%dim == 2)

  ALLOCATE(tpsi(NP, 3), NP*3)

  tpsi = M_z0
  call zproject(gr%m, h%ep%lso(3, :), h%ep%nvnl, psi(:, 1), tpsi(:, 1), simul_box_is_periodic(gr%sb), ik)
  call zproject(gr%m, h%ep%lso(1, :), h%ep%nvnl, psi(:, 2), tpsi(:, 2), simul_box_is_periodic(gr%sb), ik)
  call zproject(gr%m, h%ep%lso(2, :), h%ep%nvnl, psi(:, 2), tpsi(:, 3), simul_box_is_periodic(gr%sb), ik)
  hpsi(1:NP, 1) = hpsi(1:NP, 1) - M_zI * M_HALF * (tpsi(:, 1) + tpsi(:, 2) + M_zI * tpsi(:, 3))

  tpsi = M_z0
  call zproject(gr%m, h%ep%lso(3, :), h%ep%nvnl, psi(:, 2), tpsi(:, 1), simul_box_is_periodic(gr%sb), ik)
  call zproject(gr%m, h%ep%lso(1, :), h%ep%nvnl, psi(:, 1), tpsi(:, 2), simul_box_is_periodic(gr%sb), ik)
  call zproject(gr%m, h%ep%lso(2, :), h%ep%nvnl, psi(:, 1), tpsi(:, 3), simul_box_is_periodic(gr%sb), ik)
  hpsi(1:NP, 2) = hpsi(1:NP, 2) - M_zI * M_HALF * (-tpsi(:, 1) + tpsi(:, 2) - M_zI * tpsi(:, 3))

  deallocate(tpsi)

  call pop_sub()
end subroutine zso

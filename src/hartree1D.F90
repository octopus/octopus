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

subroutine hartree1D_solve(h, m, pot, dist)
  type(hartree_type), intent(inout) :: h
  type(mesh_type), intent(IN) :: m
  real(r8), dimension(:), intent(inout) :: pot
  real(r8), dimension(:, :), intent(IN) :: dist

  integer  :: i, j
  real(r8) :: x, y

  sub_name = 'hartree1D_solve'; call push_sub()

  do i=1, m%np
     pot(i) = 0.0_r8
     call mesh_x(m, i, x)
     do j=1, m%np
        call mesh_x(m, j, y)
        pot(i) = pot(i) + sum(dist(j, :))/sqrt(1.0_r8 + (x-y)**2)
     enddo
     pot(i) = pot(i)*m%vol_pp
  enddo

  call pop_sub(); return
end subroutine hartree1D_solve

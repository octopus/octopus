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

subroutine hartree2D_solve(h, m, pot, dist)
  type(hartree_type), intent(inout) :: h
  type(mesh_type), intent(IN) :: m
  real(r8), dimension(:), intent(inout) :: pot
  real(r8), dimension(:, :), intent(IN) :: dist

  integer  :: i, ip, j, jp
  real(r8) :: x(2), y(2)

  call push_sub('hartree2D_solve')

  pot = M_ZERO
  do i = 1, m%np
    call mesh_xyz(m, i, x)
    do j = 1, m%np
      if(i == j) then
        pot(i) = pot(i) + M_TWO*sqrt(M_PI)*sum(dist(i, :))/m%h(1)
      else
        call mesh_xyz(m, j, y)
        pot(i) = pot(i) + sum(dist(j, :))/sqrt(sum((x-y)**2))
      end if
    end do
    pot(i) = pot(i)*m%vol_pp
  end do

  call pop_sub()
end subroutine hartree2D_solve

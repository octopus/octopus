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

subroutine poisson2D_solve(m, pot, rho)
  type(mesh_type), intent(IN) :: m
  FLOAT, intent(out) :: pot(m%np)
  FLOAT, intent(in)  :: rho(m%np)

  integer  :: i, ip, j, jp
  FLOAT :: x(2), y(2)

  ASSERT(poisson_solver == -2)
  
  call push_sub('poisson2D_solve')

  pot = M_ZERO
  do i = 1, m%np
    call mesh_xyz(m, i, x)
    do j = 1, m%np
      if(i == j) then
        pot(i) = pot(i) + M_TWO*sqrt(M_PI)*rho(i)/m%h(1)
      else
        call mesh_xyz(m, j, y)
        pot(i) = pot(i) + rho(j)/sqrt(sum((x-y)**2))
      end if
    end do
    pot(i) = pot(i)*m%vol_pp
  end do

  call pop_sub()
end subroutine poisson2D_solve

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

module hartree
use global
use liboct
use mesh

implicit none

private

type hartree_type
  integer :: solver
end type hartree_type

public :: hartree_type, hartree_init, hartree_solve, hartree_end

contains

subroutine hartree_init(h, m)
  type(hartree_type), intent(inout) :: h
  type(mesh_type), intent(inout) :: m

  message(1) = 'Info: Using direct integration method to solve poisson equation'
  call write_info(1)
    
  return
end subroutine hartree_init

subroutine hartree_end(h)
  type(hartree_type), intent(inout) :: h
end subroutine hartree_end

subroutine hartree_solve(h, m, pot, dist)
  type(hartree_type), intent(inout) :: h
  type(mesh_type), intent(IN) :: m
  real(r8), dimension(:), intent(inout) :: pot
  real(r8), dimension(:), intent(IN) :: dist

  integer  :: i, j
  real(r8) :: x, y

  sub_name = 'hartree_solve'; call push_sub()

  do i=1, m%np
     pot(i) = 0.0_r8
     call mesh_x(m, i, x)
     do j=1, m%np
        call mesh_x(m, j, y)
        pot(i) = pot(i) + dist(j)/sqrt(1.0_r8 + (x-y)**2)
     enddo
     pot(i) = pot(i)*m%vol_pp
  enddo

  call pop_sub(); return
end subroutine hartree_solve

end module hartree

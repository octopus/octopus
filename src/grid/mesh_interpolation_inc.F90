!! Copyright (C) 2014 X. Andrade
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

R_TYPE function X(mesh_interpolation_evaluate)(this, values, position) result(interpolated_value)
  type(mesh_interpolation_t), target, intent(in)    :: this
  R_TYPE,                             intent(in)    :: values(:)
  FLOAT,                              intent(in)    :: position(:)

  type(mesh_t), pointer :: mesh
  integer :: nmax(1:MAX_DIM), nmin(1:MAX_DIM)
  FLOAT :: xd(1:MAX_DIM), posrel(1:MAX_DIM)

  PUSH_SUB(X(mesh_interpolation_evaluate))

  mesh => this%mesh

  ASSERT(mesh%sb%dim == 3)

  posrel(1:mesh%sb%dim) = (position(1:mesh%sb%dim) - mesh%sb%box_offset(1:mesh%sb%dim))/mesh%spacing(1:mesh%sb%dim)

  nmin(1:mesh%sb%dim) = floor(posrel(1:mesh%sb%dim))
  nmax(1:mesh%sb%dim) = nmin(1:mesh%sb%dim) + 1

  xd(1:mesh%sb%dim) = posrel(1:mesh%sb%dim) - nmin(1:mesh%sb%dim)

  ASSERT(all(xd(1:mesh%sb%dim) >= CNST(0.0)))
  ASSERT(all(xd(1:mesh%sb%dim) <= CNST(1.0)))

  select case(mesh%sb%dim)
  case(3)
    call interpolation_3d()
  end select
  
  POP_SUB(X(interpolation_point_evaluate))

  contains 
    
    subroutine interpolation_3d()

      R_TYPE :: c00, c10, c01, c11, c0, c1

      ! trilinear interpolation : http://en.wikipedia.org/wiki/Trilinear_interpolation
      
      c00 = values(mesh%idx%lxyz_inv(nmin(1), nmin(2), nmin(3)))*(CNST(1.0) - xd(1)) + &
        values(mesh%idx%lxyz_inv(nmax(1), nmin(2), nmin(3)))*xd(1)
      c10 = values(mesh%idx%lxyz_inv(nmin(1), nmax(2), nmin(3)))*(CNST(1.0) - xd(1)) + &
        values(mesh%idx%lxyz_inv(nmax(1), nmax(2), nmin(3)))*xd(1)
      c01 = values(mesh%idx%lxyz_inv(nmin(1), nmin(2), nmax(3)))*(CNST(1.0) - xd(1)) + &
        values(mesh%idx%lxyz_inv(nmax(1), nmin(2), nmax(3)))*xd(1)
      c11 = values(mesh%idx%lxyz_inv(nmin(1), nmax(2), nmax(3)))*(CNST(1.0) - xd(1)) + &
        values(mesh%idx%lxyz_inv(nmax(1), nmax(2), nmax(3)))*xd(1)

      c0 = c00*(1 - xd(2)) + c10*xd(2)
      c1 = c01*(1 - xd(2)) + c11*xd(2)
      
      interpolated_value = c0*(1 - xd(3)) - c1*xd(3)

    end subroutine interpolation_3d

end function X(mesh_interpolation_evaluate)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

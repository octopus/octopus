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

R_TYPE function X(interpolator_interpolate)(this, values, position) result(interpolated_value)
  type(interpolator_t), intent(inout) :: this
  R_TYPE,               intent(in)    :: values(:)
  FLOAT,                intent(in)    :: position(:)
  
  type(qshep_t) :: interp

  ASSERT(this%ndim == 3)

  ASSERT(ubound(values, dim = 1) >= this%npoints)
  ASSERT(ubound(position, dim = 1) >= this%ndim)

  call qshep_init(interp, this%npoints, values, this%points(:, 1), this%points(:, 2), this%points(:, 3))
  interpolated_value = qshep_interpolate(interp, values, position)
  call qshep_end(interp)

end function X(interpolator_interpolate)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

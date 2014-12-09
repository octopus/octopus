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

R_TYPE function X(interpolation_point_evaluate)(this, values) result(evaluated_value)
  type(interpolation_point_t), intent(inout) :: this
  R_TYPE,                      intent(in)    :: values(:)

  type(qshep_t) :: interp

  PUSH_SUB(X(interpolation_point_evaluate))

  ASSERT(this%interpolation%ndim == 3)

  ASSERT(ubound(values, dim = 1) >= this%interpolation%npoints)

  call qshep_init(interp, this%interpolation%npoints, values, &
    this%interpolation%points(:, 1), this%interpolation%points(:, 2), this%interpolation%points(:, 3))

  evaluated_value = qshep_interpolate(interp, values, this%position)

  call qshep_end(interp)

  POP_SUB(X(interpolation_point_evaluate))

end function X(interpolation_point_evaluate)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

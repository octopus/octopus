!! Copyright (C) 2021 M. Oliveira
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

#include "global.h"

module box_shape_oct_m
  use box_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::         &
    box_shape_t,    &
    box_shape_init, &
    box_shape_end

  !> Base class for more specialized boxes that are defined by a shape and have
  !! a center.
  type, abstract, extends(box_t) :: box_shape_t
    private
    FLOAT, allocatable, public :: center(:) !< where is the box centered
  end type box_shape_t

  FLOAT, parameter, public :: BOX_BOUNDARY_DELTA = CNST(1e-12)
  
contains

  !--------------------------------------------------------------
  subroutine box_shape_init(this, dim, center)
    class(box_shape_t), intent(inout) :: this
    integer,            intent(in)    :: dim
    FLOAT,              intent(in)    :: center(dim)

    this%dim = dim
    SAFE_ALLOCATE(this%center(1:dim))
    this%center(1:dim) = center(1:dim)

  end subroutine box_shape_init

  !--------------------------------------------------------------
  subroutine box_shape_end(this)
    class(box_shape_t), intent(inout) :: this

    SAFE_DEALLOCATE_A(this%center)

  end subroutine box_shape_end

end module box_shape_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module box_parallelepiped_oct_m
  use box_shape_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  public :: box_parallelepiped_t

  !> Class implementing a parallelepiped box. Currently this is restricted to a
  !! rectangular cuboid (all the faces must be rectangular) and the vectors
  !! generating the parallelepiped must be along the Cartesian axes.
  type, extends(box_shape_t) :: box_parallelepiped_t
    private
    FLOAT, allocatable :: half_length(:) !< half the length of the parallelepiped in each direction.
  contains
    procedure :: contains_points => box_parallelepiped_contains_points
    final     :: box_parallelepiped_finalize
  end type box_parallelepiped_t

  interface box_parallelepiped_t
    procedure box_parallelepiped_constructor
  end interface box_parallelepiped_t

contains

  !--------------------------------------------------------------
  function box_parallelepiped_constructor(dim, center, length) result(box)
    integer,            intent(in) :: dim
    FLOAT,              intent(in) :: center(dim)
    FLOAT,              intent(in) :: length(dim) !< length of the parallelepiped along each Cartesian direction
    class(box_parallelepiped_t), pointer :: box

    PUSH_SUB(box_parallelepiped_constructor)

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    call box_shape_init(box, dim, center)
    SAFE_ALLOCATE(box%half_length(1:dim))
    box%half_length = M_HALF*length

    POP_SUB(box_parallelepiped_constructor)
  end function box_parallelepiped_constructor

  !--------------------------------------------------------------
  subroutine box_parallelepiped_finalize(this)
    type(box_parallelepiped_t), intent(inout) :: this

    PUSH_SUB(box_parallelepiped_finalize)

    call box_shape_end(this)
    SAFE_DEALLOCATE_A(this%half_length)

    POP_SUB(box_parallelepiped_finalize)
  end subroutine box_parallelepiped_finalize

  !--------------------------------------------------------------
  function box_parallelepiped_contains_points(this, nn, xx) result(contained)
    class(box_parallelepiped_t), intent(in)  :: this
    integer,                     intent(in)  :: nn
    FLOAT,                       intent(in)  :: xx(:,:)
    logical :: contained(1:nn)

    integer :: ip
    FLOAT :: ulimit(this%dim), llimit(this%dim)

    do ip = 1, nn
      llimit = this%center - this%half_length - BOX_BOUNDARY_DELTA
      ulimit = this%center + this%half_length + BOX_BOUNDARY_DELTA

      contained(ip) = all(xx(ip, 1:this%dim) >= llimit .and. xx(ip, 1:this%dim) <= ulimit) .neqv. this%is_inside_out()
    end do

  end function box_parallelepiped_contains_points

end module box_parallelepiped_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

!! Copyright (C) 2013 X. Andrade
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

#include "global.h"

module polyhedron_m
  use global_m
  use messages_m
  use io_m
  use profiling_m

  implicit none

  private

  public ::                         &
    polyhedron_t,                   &
    polyhedron_init,                &
    polyhedron_end,                 &
    polyhedron_add_point,           &
    polyhedron_add_triangle

  type polyhedron_t
    FLOAT,   pointer :: points(:, :)
    integer, pointer :: point_indices(:)
    integer, pointer :: triangles(:, :)
    integer          :: nalloc_points
    integer          :: nalloc_triangles
    integer          :: npoints
    integer          :: ntriangles
  end type polyhedron_t

contains

  subroutine polyhedron_init(this)
    type(polyhedron_t), intent(out) :: this
    
    PUSH_SUB(polyhedron_init)

    this%nalloc_points = 1000
    this%nalloc_triangles = 1000
    SAFE_ALLOCATE(this%points(1:3, 1:this%nalloc_points))
    SAFE_ALLOCATE(this%point_indices(1:this%nalloc_points))
    SAFE_ALLOCATE(this%triangles(1:3, this%nalloc_triangles))
    this%npoints = 0
    this%ntriangles = 0

    POP_SUB(polyhedron_init)
  end subroutine polyhedron_init

  !-------------------------------------------------------

  subroutine polyhedron_end(this)
    type(polyhedron_t), intent(inout) :: this
    
    PUSH_SUB(polyhedron_end)

    SAFE_DEALLOCATE_P(this%points)
    SAFE_DEALLOCATE_P(this%point_indices)
    SAFE_DEALLOCATE_P(this%triangles)

    POP_SUB(polyhedron_end)
  end subroutine polyhedron_end
  
  !-------------------------------------------------------

  subroutine polyhedron_add_point(this, index, pos)
    type(polyhedron_t), intent(inout) :: this
    integer,            intent(in)    :: index
    FLOAT,              intent(in)    :: pos(:)

    FLOAT, pointer :: new_points(:, :)
    integer, pointer :: new_point_indices(:)

    if(this%npoints == this%nalloc_points) then
      SAFE_ALLOCATE(new_points(1:3, this%nalloc_points*2))
      SAFE_ALLOCATE(new_point_indices(this%nalloc_points*2))
      new_points(1:3, 1:this%nalloc_points) = this%points(1:3, 1:this%nalloc_points)
      new_point_indices(1:this%nalloc_points) = this%point_indices(1:this%nalloc_points)
      this%nalloc_points = this%nalloc_points*2
      SAFE_DEALLOCATE_P(this%points)
      SAFE_DEALLOCATE_P(this%point_indices)
      this%points => new_points
      this%point_indices => new_point_indices
    end if

    this%npoints = this%npoints + 1
    this%points(1:3, this%npoints) = pos(1:3)
    this%point_indices(this%npoints) = index

  end subroutine polyhedron_add_point

  !------------------------------------------------------

  subroutine polyhedron_add_triangle(this, triangle)
    type(polyhedron_t), intent(inout) :: this
    integer,            intent(in)    :: triangle(:)
    
    integer, pointer :: new_triangles(:, :)

    if(this%ntriangles == this%nalloc_triangles) then
      SAFE_ALLOCATE(new_triangles(1:3, this%nalloc_triangles*2))
      new_triangles(1:3, 1:this%nalloc_triangles) = this%triangles(1:3, 1:this%nalloc_triangles)
      this%nalloc_triangles = this%nalloc_triangles*2
      SAFE_DEALLOCATE_P(this%triangles)
      this%triangles => new_triangles
    end if

    this%ntriangles = this%ntriangles + 1
    this%triangles(1:3, this%ntriangles) = triangle(1:3)

  end subroutine polyhedron_add_triangle

  !------------------------------------------------------

end module polyhedron_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

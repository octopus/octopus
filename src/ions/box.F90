!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2013 M. Oliveira
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

module box_m
  use global_m
  use io_m
  use loct_m
  use messages_m
  use mpi_m
  use profiling_m

  implicit none

  private
  public ::         &
    box_t,          &
    box_create,     &
    box_end,        &
    box_inside,     &
    box_inside_vec, &
    box_copy,       &
    box_inside_bader, &
    box_nearest_point


  integer, parameter, public :: &
    BOX_SPHERE         = 1,         &
    BOX_CYLINDER       = 2,         &
    BOX_PARALLELEPIPED = 3,         &
    BOX_BADER          = 4

  type box_t
    private

    integer :: dim
    integer :: shape

    FLOAT :: center(MAX_DIM)  !< where is the box centered

    FLOAT :: rsize          !< the radius of the sphere or of the cylinder
    FLOAT :: xsize          !< the length of the cylinder in the x-direction
    FLOAT :: lsize(MAX_DIM) !< half of the length of the parallelepiped in each direction.
  end type box_t

contains

  !--------------------------------------------------------------
  subroutine box_create(box, shape, dim, sizes, center)
    type(box_t), intent(out) :: box
    integer,     intent(in)  :: shape
    integer,     intent(in)  :: dim
    FLOAT,       intent(in)  :: sizes(MAX_DIM)
    FLOAT,       intent(in)  :: center(dim)

    PUSH_SUB(box_create)

    box%shape = shape
    box%dim = dim
    box%center(1:dim) = center(1:dim)

    select case (shape)
    case (BOX_SPHERE, BOX_BADER)
      box%rsize = sizes(1)

    case (BOX_CYLINDER)
      if (dim == 2) then
        message(1) = "Cannot create a cylinder in 2D. Use sphere if you want a circle."
        call messages_fatal(1)
      end if
      box%rsize = sizes(1)
      box%xsize = sizes(2)

    case (BOX_PARALLELEPIPED)
      box%lsize(1:dim) = sizes(1:dim)

    case default
      message(1) = "Unknown box shape in box_create."
      call messages_fatal(1)

    end select
    
    POP_SUB(box_create)
  end subroutine box_create

  !--------------------------------------------------------------
  subroutine box_end(box)
    type(box_t), intent(inout) :: box

    PUSH_SUB(box_end)

    box%shape = 0
    box%dim = 0

    POP_SUB(box_end)
  end subroutine box_end

  !--------------------------------------------------------------
  !> Checks if a point is inside the box.
  logical function box_inside(box, point) result(inside)
    type(box_t), intent(in) :: box
    FLOAT,       intent(in) :: point(:)

    FLOAT :: xx(1, 1:MAX_DIM)
    logical :: inside2(1)

    ! no push_sub because this function is called very frequently

    xx(1, 1:box%dim) = point(1:box%dim)

    call box_inside_vec(box, 1, xx, inside2)
    inside = inside2(1)

  end function box_inside

  !--------------------------------------------------------------
  !> Checks if a vector of points are inside the box.
  subroutine box_inside_vec(box, npoints, points, inside)
    type(box_t),  intent(in)  :: box
    integer,      intent(in)  :: npoints
    FLOAT,        intent(in)  :: points(:, :)
    logical,      intent(out) :: inside(:)

    integer :: ip
    real(8), parameter :: DELTA = CNST(1e-12)
    real(8) :: llimit(MAX_DIM), ulimit(MAX_DIM)
    FLOAT :: rr
    FLOAT, allocatable :: xx(:, :)

    ! no push_sub because this function is called very frequently

    SAFE_ALLOCATE(xx(1:box%dim, 1:npoints))
    forall(ip = 1:npoints)
      xx(1:box%dim, ip) = points(ip, 1:box%dim) - box%center(1:box%dim)
    end forall

    select case(box%shape)
    case(BOX_SPHERE)
      forall(ip = 1:npoints)
        inside(ip) = sum(xx(1:box%dim, ip)**2) <= (box%rsize + DELTA)**2
      end forall

    case(BOX_CYLINDER)
      do ip = 1, npoints
        rr = sqrt(sum(xx(2:box%dim, ip)**2))
        inside(ip) = (rr <= box%rsize + DELTA .and. abs(xx(1, ip)) <= box%xsize + DELTA)
      end do

    case(BOX_PARALLELEPIPED) 
      llimit(1:box%dim) = -box%lsize(1:box%dim) - DELTA
      ulimit(1:box%dim) =  box%lsize(1:box%dim) + DELTA

      forall(ip = 1:npoints)
        inside(ip) = all(xx(1:box%dim, ip) >= llimit(1:box%dim) .and. xx(1:box%dim, ip) <= ulimit(1:box%dim))
      end forall

    case(BOX_BADER) 
      message(1) = 'Bader Volumes are not yet implemented'
      call messages_fatal(1)
      !TODO: call box_inside_bader(box, npoints, points, inside, ff, mesh)

    end select

    SAFE_DEALLOCATE_A(xx)

  end subroutine box_inside_vec

  ! --------------------------------------------------------------
  recursive subroutine box_copy(boxout, boxin)
    type(box_t), intent(out) :: boxout
    type(box_t), intent(in)  :: boxin

    PUSH_SUB(box_copy)

    boxout%shape               = boxin%shape
    boxout%dim                 = boxin%dim
    boxout%center(1:boxin%dim) = boxin%center(1:boxin%dim)
    boxout%rsize               = boxin%rsize
    boxout%xsize               = boxin%xsize
    boxout%lsize(1:boxin%dim) = boxin%lsize(1:boxin%dim)

    POP_SUB(box_copy)
  end subroutine box_copy

  !--------------------------------------------------------------
  !> Checks if a vector of points are inside the Bader volume defined by box%center.
  ! TODO: Added the skeleton of a new routine that will check if a point is inside a Bader Volume. 
  !--------------------------------------------------------------
  subroutine box_inside_bader(box, npoints, points, inside)
    type(box_t),  intent(in)  :: box
    integer,      intent(in)  :: npoints
    FLOAT,        intent(in)  :: points(:, :)
    logical,      intent(out) :: inside(:)

    integer :: ip, is, ix, iy, iz, rankmin
    real(8), parameter :: DELTA = CNST(1e-12)
    FLOAT, allocatable :: x(:), xx(:)
    FLOAT :: dmin
    logical :: rhomax
    

  end subroutine box_inside_bader

  !---------------------------------------------------------------------
  !> Returns the index of the point which is nearest to a given vector
  !! position pos. Variable dmin will hold, on exit, the distance between
  !! pos and this nearest mesh point. rankmin will be zero, if the mesh is
  !! not partitioned, and the rank of the processor which holds the point
  !! ind if the mesh is partitioned.
  !! This routine is a copy of the mesh_nearest_point in grid/mesh.F90. 
  !! Here the routine does not need a type(mesh_t) variable.
  ! ----------------------------------------------------------------------
  integer function box_nearest_point(npoints, points, pos, dim, dmin, rankmin) result(ind)
    integer,      intent(in)  :: npoints
    FLOAT,        intent(in)  :: points(:, :)
    FLOAT,        intent(in)  :: pos(MAX_DIM)
    integer,      intent(in)  :: dim
    FLOAT,        intent(out) :: dmin
    integer,      intent(out) :: rankmin
    
    FLOAT :: dd
    integer :: imin, ip
#if defined(HAVE_MPI)
    FLOAT :: min_loc_in(2), min_loc_out(2)
#endif
    
    PUSH_SUB(box_nearest_point)
    
    !find the point of the grid that is closer to the atom
    imin = 0
    dmin = M_ZERO
    do ip = 1, npoints
      dd = sum((pos(1:dim) - points(ip, 1:dim))**2)
      if((dd < dmin) .or. (ip == 1)) then 
        imin = ip
        dmin = dd
      end if
    end do
    
    ind = imin
    POP_SUB(box_nearest_point)
  end function box_nearest_point
  
end module box_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

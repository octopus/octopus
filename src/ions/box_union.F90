!! Copyright (C) 2014 M. Oliveira
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

module box_union_m
  use global_m
  use io_m
  use loct_m
  use messages_m
  use mpi_m
  use profiling_m
  use box_m

  implicit none

  private
  public ::               &
    box_union_t,          &
    box_union_init,       &
    box_union_end,        &
    box_union_inside,     &
    box_union_inside_vec

  type box_union_t
    private
    
    !! TODO: make this a linked list, so that boxes can be added and removed efficiently on-the-fly
    integer :: n_boxes
    type(box_t), pointer :: boxes(:)
  end type box_union_t

contains

  !--------------------------------------------------------------
  subroutine box_union_init(union, n_boxes, boxes)
    type(box_union_t), intent(out) :: union
    integer,           intent(in)  :: n_boxes
    type(box_t),       intent(in)  :: boxes(:)

    integer :: ibox

    PUSH_SUB(box_union_init)

    union%n_boxes = n_boxes
    SAFE_ALLOCATE(union%boxes(n_boxes))

    do  ibox = 1,n_boxes
      call box_copy(union%boxes(ibox), boxes(ibox))
    end do
    
    POP_SUB(box_union_init)
  end subroutine box_union_init

  !--------------------------------------------------------------
  subroutine box_union_end(union)
    type(box_union_t), intent(inout) :: union

    integer :: ibox

    PUSH_SUB(box_union_end)

    do ibox = 1, union%n_boxes
      call box_end(union%boxes(ibox))
    end do
    SAFE_DEALLOCATE_P(union%boxes)

    union%n_boxes = 0

    POP_SUB(box_union_end)
  end subroutine box_union_end

  !--------------------------------------------------------------
  !> Checks if a vector of points are inside the box.
  subroutine box_union_inside_vec(union, npoints, points, inside)
    type(box_union_t),  intent(in)  :: union
    integer,            intent(in)  :: npoints
    FLOAT,              intent(in)  :: points(:, :)
    logical,            intent(out) :: inside(:)

    integer :: ibox
    logical, allocatable :: inside2(:)

    ! no push_sub because this function is called very frequently
    SAFE_ALLOCATE(inside2(npoints))

    inside = .false.
    do ibox = 1, union%n_boxes
      call box_inside_vec(union%boxes(ibox), npoints, points, inside2)
      inside = inside .or. inside2
    end do

    SAFE_DEALLOCATE_A(inside2)

  end subroutine box_union_inside_vec
  
  !--------------------------------------------------------------
  !> Checks if a point are inside the union box.
  logical function box_union_inside(union, point) result(inside)
    type(box_union_t),  intent(in)  :: union
    FLOAT,              intent(in)  :: point(:)

    integer :: ibox
    logical :: inside2

    ! no push_sub because this function is called very frequently

    inside = .false.
    do ibox = 1, union%n_boxes
      if(box_inside(union%boxes(ibox), point)) inside = .true.
    end do

  end function box_union_inside

end module box_union_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

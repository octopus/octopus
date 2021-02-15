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

module box_union_oct_m
  use box_oct_m
  use multibox_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::           &
    box_union_t

  !> Class implementing a box that is an union other boxes.
  type, extends(multibox_t) :: box_union_t
    private
  contains
    procedure :: contains_points => box_union_contains_points
    procedure :: write_info => box_union_write_info
    procedure :: write_short_info => box_union_write_short_info
    final     :: box_union_finalize
  end type box_union_t

  interface box_union_t
    procedure box_union_constructor
  end interface box_union_t

contains

  !--------------------------------------------------------------
  function box_union_constructor(dim) result(box)
    integer, intent(in) :: dim
    class(box_union_t), pointer :: box

    PUSH_SUB(box_union_constructor)

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    box%dim = dim

    POP_SUB(box_union_constructor)
  end function box_union_constructor

  !--------------------------------------------------------------
  subroutine box_union_finalize(this)
    type(box_union_t), intent(inout) :: this

    PUSH_SUB(box_union_finalize)

    call multibox_end(this)

    POP_SUB(box_union_finalize)
  end subroutine box_union_finalize

  !--------------------------------------------------------------
  recursive function box_union_contains_points(this, nn, xx) result(contained)
    class(box_union_t), intent(in) :: this
    integer,            intent(in) :: nn
    FLOAT,              intent(in) :: xx(:,:)
    logical :: contained(nn)

    integer :: ip
    FLOAT :: point(1:this%dim)
    type(box_iterator_t) :: iter
    class(box_t), pointer :: box

    ! A point must be inside at least one box to be considered inside an union of boxes
    do ip = 1, nn
      point(1:this%dim) = xx(ip, 1:this%dim)
      contained(ip) = .false.

      call iter%start(this%list)
      do while (iter%has_next())
        box => iter%get_next()
        contained(ip) = box%contains_point(point)
        if (contained(ip)) exit
      end do

      contained(ip) = contained(ip) .neqv. this%is_inside_out()
    end do

  end function box_union_contains_points

  !--------------------------------------------------------------
  subroutine box_union_write_info(this, iunit)
    class(box_union_t), intent(in) :: this
    integer,            intent(in) :: iunit

    PUSH_SUB(box_union_write_info)

    ! Todo: need to decide how best to display the information of the boxes that make the union

    POP_SUB(box_union_write_info)
  end subroutine box_union_write_info

  !--------------------------------------------------------------
  subroutine box_union_write_short_info(this, iunit)
    class(box_union_t), intent(in) :: this
    integer,            intent(in) :: iunit

    PUSH_SUB(box_union_write_short_info)

    ! Todo: need to decide how best to display the information of the boxes that make the union

    POP_SUB(box_union_write_short_info)
  end subroutine box_union_write_short_info

end module box_union_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

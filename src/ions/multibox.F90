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

module multibox_oct_m
  use box_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::      &
    multibox_t,  &
    multibox_end

  !> Abstract class for boxes that are made up of a list of boxes.
  type, abstract, extends(box_t) :: multibox_t
    private
    type(box_list_t), public :: list !< list containing the boxes that make up this multibox
  contains
    procedure :: add_box => multibox_add_box
  end type multibox_t

contains

  !--------------------------------------------------------------
  subroutine multibox_end(this)
    class(multibox_t), intent(inout) :: this

    type(box_iterator_t) :: iter
    class(box_t), pointer :: box

    PUSH_SUB(multibox_end)

    call iter%start(this%list)
    do while (iter%has_next())
      box => iter%get_next()
      SAFE_DEALLOCATE_P(box)
    end do

    POP_SUB(multibox_end)
  end subroutine multibox_end

  !--------------------------------------------------------------
  subroutine multibox_add_box(this, new_box)
    class(multibox_t), intent(inout) :: this
    class(box_t),      intent(in)    :: new_box

    PUSH_SUB(multibox_add_box)

    ASSERT(this%dim == new_box%dim)

    call this%list%add(new_box)

    POP_SUB(multibox_add_box)
  end subroutine multibox_add_box

end module multibox_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

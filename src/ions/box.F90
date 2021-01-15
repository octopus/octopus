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

module box_oct_m
  use global_oct_m
  use linked_list_oct_m

  implicit none

  private
  public ::        &
    box_t,         &
    box_list_t,    &
    box_iterator_t

  !> The purpose of a box is to tell if something is inside or outside of it.
  !! To do that it provides a function that tells if a given list of points are
  !! inside or outside the box. Furthermore, a box might be turned inside out,
  !! i.e., in that case what is usually considered inside becomes outside and
  !! vice-versa.
  type, abstract :: box_t
    private
    integer, public :: dim                    !< dimensions of the space the box lives in
    logical :: inside_out = .false.           !< if the box is inside out or not
  contains
    procedure(box_contains_points), deferred :: contains_points
    procedure, non_overridable :: contains_point => box_contains_point
    procedure, non_overridable :: is_inside_out => box_is_inside_out
    procedure, non_overridable :: turn_inside_out => box_turn_inside_out
  end type box_t

  abstract interface
    !> Given a list of points, this function should return an array indicating
    !! for each point if it is inside the box or not.
    recursive function box_contains_points(this, nn, xx) result(contained)
      import :: box_t
      class(box_t), intent(in) :: this
      integer,      intent(in) :: nn      !< number of points to check
      FLOAT,        intent(in) :: xx(:,:) !< points to check. The sizes are
                                          !! (1:,1:this%dim), so that it is
                                          !! possible to pass an array with more
                                          !! points than the ones we are
                                          !! checking.
      logical :: contained(1:nn)
    end function box_contains_points
  end interface

  !> These classes extends the list and list iterator to create a box list.
  type, extends(linked_list_t) :: box_list_t
    private
  contains
    procedure :: add => box_list_add_node
  end type box_list_t

  type, extends(linked_list_iterator_t) :: box_iterator_t
    private
  contains
    procedure :: get_next => box_iterator_get_next
  end type box_iterator_t

contains

  !!--------------------------------------------------------------
  !> Turn a box inside out.
  subroutine box_turn_inside_out(this)
    class(box_t), intent(inout) :: this

    this%inside_out = .not. this%inside_out

  end subroutine box_turn_inside_out

  !!--------------------------------------------------------------
  !> Is the box inside out?
  logical function box_is_inside_out(this)
    class(box_t), intent(in) :: this

    box_is_inside_out = this%inside_out

  end function box_is_inside_out

  !!---------------------------------------------------------------
  !> Convenience function to check if a single point is inside the box when that
  !! point is passed as a rank-one array.
  recursive logical function box_contains_point(this, xx) result(contained)
    class(box_t),         intent(in) :: this
    FLOAT,        target, intent(in) :: xx(1:this%dim)

    FLOAT, pointer :: xx_ptr(:,:)
    logical :: points_contained(1)

    xx_ptr(1:1, 1:this%dim) => xx(1:this%dim)
    points_contained = this%contains_points(1, xx_ptr)
    contained = points_contained(1)

  end function box_contains_point

  ! ---------------------------------------------------------
  subroutine box_list_add_node(this, box)
    class(box_list_t)    :: this
    class(box_t), target :: box

    select type (box)
    class is (box_t)
      call this%add_ptr(box)
    class default
      ASSERT(.false.)
    end select

  end subroutine box_list_add_node

  ! ---------------------------------------------------------
  function box_iterator_get_next(this) result(box)
    class(box_iterator_t), intent(inout) :: this
    class(box_t),          pointer       :: box

    select type (ptr => this%get_next_ptr())
    class is (box_t)
      box => ptr
    class default
      ASSERT(.false.)
    end select

  end function box_iterator_get_next

end module box_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

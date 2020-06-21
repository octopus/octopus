!! Copyright (C) 2019 M. Oliveira
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

module list_node_oct_m
  use global_oct_m
  use messages_oct_m
  implicit none

  private
  public :: list_node_t, &
            list_node

  type :: list_node_t
    private
    logical :: clone
    class(*),          pointer :: value => null()
    type(list_node_t), pointer :: next_node => null()
  contains
    procedure :: get
    procedure :: next
    procedure :: set_next
    procedure :: is_equal
    final :: finalize
  end type list_node_t

  interface list_node
    procedure constructor
  end interface list_node

contains

  ! ---------------------------------------------------------
  function constructor(value, next, clone)
    class(*),           target     :: value
    class(list_node_t), pointer    :: next
    logical,            intent(in) :: clone
    class(list_node_t), pointer    :: constructor

    ! No safe_allocate macro here, as its counterpart in linked_list.F90
    ! causes an internal compiler error with GCC 6.4.0
    allocate(constructor)
    constructor%next_node => next
    constructor%clone = clone
    if (constructor%clone) then
      allocate(constructor%value, source=value)
    else
      constructor%value => value
    end if

  end function constructor

  ! ---------------------------------------------------------
  function next(this)
    class(list_node_t), intent(in) :: this
    class(list_node_t), pointer    :: next

    next => this%next_node

  end function next

  ! ---------------------------------------------------------
  subroutine set_next(this, next_node)
    class(list_node_t), intent(inout) :: this
    class(list_node_t), pointer       :: next_node

    this%next_node => next_node

  end subroutine set_next

  ! ---------------------------------------------------------
  function get(this)
    class(list_node_t), intent(in) :: this
    class(*),           pointer :: get

    get => this%value

  end function get

  ! ---------------------------------------------------------
  logical function is_equal(this, value)
    class(list_node_t), intent(in) :: this
    class(*),           target     :: value

    ! First try to match the two types and compare the values.
    ! Note that the list of types taken into account might not be exhaustive.
    is_equal = .false.
    select type (ptr => this%value)
    type is (integer)
      select type (value)
      type is (integer)
        is_equal = value == ptr
      end select
    type is (FLOAT)
      select type (value)
      type is (FLOAT)
        is_equal = value == ptr
      end select
    type is (complex)
      select type (value)
      type is (complex)
        is_equal = value == ptr
      end select
    type is (character(len=*))
      select type (value)
      type is (character(len=*))
        is_equal = value == ptr
      end select
    type is (logical)
      select type (value)
      type is (logical)
        is_equal = value .eqv. ptr
      end select
    end select

    ! If we were not able to match the types, then we check if the two values
    ! point to the same target.
    if (.not. is_equal) then
      is_equal = associated(this%value, value)
    end if

  end function is_equal

  subroutine finalize(this)
    type(list_node_t), intent(inout) :: this

    if (associated(this%next_node)) then
      nullify(this%next_node)
    end if
    if (associated(this%value)) then
      if (this%clone) then
        deallocate(this%value)
      else
        nullify(this%value)
      end if
    end if

  end subroutine finalize
  
end module list_node_oct_m

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
  public :: list_node_t

  type :: list_node_t
    private
    logical :: clone
    class(*),          pointer :: value => null()
    type(list_node_t), pointer :: next_node => null()
  contains
    procedure :: get => list_node_get
    procedure :: next => list_node_next
    procedure :: set_next => list_node_set_next
    procedure :: is_equal => list_node_is_equal
    procedure :: copy => list_node_copy
    final :: list_node_finalize
  end type list_node_t

  interface list_node_t
    procedure list_node_constructor
  end interface list_node_t

contains

  ! ---------------------------------------------------------
  function list_node_constructor(value, next, clone) result(constructor)
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

  end function list_node_constructor

  ! ---------------------------------------------------------
  function list_node_copy(this, next)
    class(list_node_t), target  :: this
    class(list_node_t), pointer :: next
    class(list_node_t), pointer :: list_node_copy

    list_node_copy => list_node_constructor(this%value, next, this%clone)

  end function list_node_copy

  ! ---------------------------------------------------------
  function list_node_next(this) result(next)
    class(list_node_t), intent(in) :: this
    class(list_node_t), pointer    :: next

    next => this%next_node

  end function list_node_next

  ! ---------------------------------------------------------
  subroutine list_node_set_next(this, next_node)
    class(list_node_t), intent(inout) :: this
    class(list_node_t), pointer       :: next_node

    this%next_node => next_node

  end subroutine list_node_set_next

  ! ---------------------------------------------------------
  function list_node_get(this) result(get)
    class(list_node_t), intent(in) :: this
    class(*),           pointer :: get

    get => this%value

  end function list_node_get

  ! ---------------------------------------------------------
  logical function list_node_is_equal(this, value) result(is_equal)
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

  end function list_node_is_equal

  subroutine list_node_finalize(this)
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

  end subroutine list_node_finalize
  
end module list_node_oct_m

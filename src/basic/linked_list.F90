!! Copyright (C) 2019-2020 M. Oliveira
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

module linked_list_oct_m
  use global_oct_m
  use list_node_oct_m
  use messages_oct_m
  implicit none

  private
  public :: linked_list_t,          &
            linked_list_iterator_t, &
            list_t,                 &
            list_iterator_t,        &
            integer_list_t,         &
            integer_iterator_t

  !---------------------------------------------------------------------------
  !> The following class implements a linked list of unlimited polymorphic
  !> values. This allows the storage of any type of data. Iterating over the
  !> list is done using the associated iterator. These two classes are not meant
  !> to used as is, but rather to be extended and by providing an add method to the
  !> list and a get_next method to the iterator.
  type :: linked_list_t
    private
    class(list_node_t), pointer :: first_node => null()
    class(list_node_t), pointer :: last_node => null()
  contains
    procedure :: add_node => linked_list_add_node
    procedure :: add_ptr  => linked_list_add_node_ptr
    procedure :: add_copy => linked_list_add_node_copy
    procedure :: delete => linked_list_delete_node
    procedure :: has => linked_list_has
    procedure :: copy => linked_list_copy
    generic   :: assignment(=) => copy
    final     :: linked_list_finalize
  end type linked_list_t

  type :: linked_list_iterator_t
    private
    class(list_node_t), pointer :: next_node => null()
  contains
    procedure :: start    => linked_list_iterator_start
    procedure :: has_next => linked_list_iterator_has_next
    procedure :: get_next_ptr => linked_list_iterator_get_next_ptr
  end type linked_list_iterator_t

  !---------------------------------------------------------------------------
  !> The following class implements a linked list of unlimited polymorphic
  !> values. Iterating over the list is done using the associated iterator and
  !> there are two ways of doing so. The first is by using a "do while"
  !> construct:
  !>
  !>  call iter%start(list)
  !>  do while (iter%has_next())
  !>     value => iter%get_next()
  !>     ...
  !>  end do
  !>
  !> The second method is with a simple "do":
  !>
  !>  call iter%start(list)
  !>  do
  !>     if (.not. iter%has_next()) exit
  !>     value => iter%get_next()
  !>     ...
  !>  end do
  !>
  type, extends(linked_list_t) :: list_t
    private
  contains
    procedure :: add => list_add_node
  end type list_t

  type, extends(linked_list_iterator_t) :: list_iterator_t
  contains
    procedure :: get_next => list_iterator_get_next
  end type list_iterator_t

  !---------------------------------------------------------------------------
  !> The following class implements a linked list of integer values. Note that
  !> the get method returns an integer, not a pointer.
  type, extends(linked_list_t) :: integer_list_t
    private
  contains
    procedure :: add => integer_list_add_node
  end type integer_list_t

  type, extends(linked_list_iterator_t) :: integer_iterator_t
  contains
    procedure :: get_next => integer_iterator_get_next
  end type integer_iterator_t

contains

  ! Linked list
  ! ---------------------------------------------------------
  subroutine linked_list_add_node(this, value, clone)
    class(linked_list_t), intent(inout) :: this
    class(*),             target        :: value
    logical,              intent(in)    :: clone

    class(list_node_t), pointer :: new_node

    if (.not. associated(this%first_node)) then
      this%first_node => list_node_t(value, this%first_node, clone)
      this%last_node => this%first_node
    else
      new_node => list_node_t(value, this%last_node%next(), clone)
      call this%last_node%set_next(new_node)
      this%last_node => new_node
    end if

  end subroutine linked_list_add_node

  ! ---------------------------------------------------------
  subroutine linked_list_add_node_ptr(this, value)
    class(linked_list_t), intent(inout) :: this
    class(*),             target        :: value

    call this%add_node(value, clone=.false.)

  end subroutine linked_list_add_node_ptr

  ! ---------------------------------------------------------
  subroutine linked_list_add_node_copy(this, value)
    class(linked_list_t), intent(inout) :: this
    class(*),             target        :: value

    call this%add_node(value, clone=.true.)

  end subroutine linked_list_add_node_copy

  ! ---------------------------------------------------------
  subroutine linked_list_delete_node(this, value)
    class(linked_list_t), intent(inout) :: this
    class(*),             target        :: value

    class(list_node_t), pointer :: previous, current, next

    previous => null()
    current => null()
    next => this%first_node
    do while (associated(next))
      previous => current
      current => next
      next => next%next()
      if (current%is_equal(value)) then
        if (associated(next) .and. .not. associated(previous)) then
          ! First node
          this%first_node => next
        else if (.not. associated(next) .and. associated(previous)) then
          ! Last node
          call previous%set_next(null())
          this%last_node => previous
        else if (.not. associated(next) .and. .not. associated(previous)) then
          ! List only has one node
          nullify(this%first_node)
          nullify(this%last_node)
        else
          ! Neither the first nor the last node
          call previous%set_next(next)
        end if
        deallocate(current)
        exit
      end if
    end do

  end subroutine linked_list_delete_node

  ! ---------------------------------------------------------
  subroutine linked_list_finalize(this)
    type(linked_list_t), intent(inout) :: this

    class(list_node_t), pointer :: current, next

    current => this%first_node
    do while (associated(current))
      next => current%next()
      deallocate(current)
      current => next
    end do
    nullify(this%first_node)
    nullify(this%last_node)

  end subroutine linked_list_finalize

  ! ---------------------------------------------------------
  subroutine linked_list_copy(lhs, rhs)
    class(linked_list_t), intent(out) :: lhs
    class(linked_list_t), intent(in)  :: rhs

    class(list_node_t), pointer :: current, new_node

    current => rhs%first_node
    do while (associated(current))
      if (.not. associated(lhs%first_node)) then
        lhs%first_node => current%copy(lhs%first_node)
        lhs%last_node => lhs%first_node
      else
        new_node => current%copy(lhs%last_node%next())
        call lhs%last_node%set_next(new_node)
        lhs%last_node => new_node
      end if
      current => current%next()
    end do

  end subroutine linked_list_copy

  ! ---------------------------------------------------------
  logical function linked_list_has(this, value)
    class(linked_list_t),        intent(inout) :: this
    class(*),                    target        :: value

    class(list_node_t), pointer :: current

    current => this%first_node
    linked_list_has = .false.
    do while (associated(current) .and. .not. linked_list_has)
      linked_list_has = current%is_equal(value)
      current => current%next()
    end do

  end function linked_list_has

  ! ---------------------------------------------------------
  subroutine linked_list_iterator_start(this, list)
    class(linked_list_iterator_t),         intent(inout) :: this
    class(linked_list_t),          target, intent(in)    :: list

    this%next_node => list%first_node

  end subroutine linked_list_iterator_start

  ! ---------------------------------------------------------
  logical function linked_list_iterator_has_next(this)
    class(linked_list_iterator_t), intent(in) :: this

    linked_list_iterator_has_next = associated(this%next_node)

  end function linked_list_iterator_has_next

  ! ---------------------------------------------------------
  function linked_list_iterator_get_next_ptr(this) result(value)
    class(linked_list_iterator_t), intent(inout) :: this
    class(*),                      pointer       :: value

    value => this%next_node%get()
    this%next_node => this%next_node%next()

  end function linked_list_iterator_get_next_ptr


  ! Unlimited polymorphic list

  ! ---------------------------------------------------------
  subroutine list_add_node(this, value)
    class(list_t), intent(inout) :: this
    class(*),      target        :: value

    call this%add_ptr(value)

  end subroutine list_add_node

  ! ---------------------------------------------------------
  function list_iterator_get_next(this) result(value)
    class(list_iterator_t), intent(inout) :: this
    class(*),               pointer       :: value

    value => this%get_next_ptr()

  end function list_iterator_get_next

  ! Integer list

  ! ---------------------------------------------------------
  subroutine integer_list_add_node(this, value)
    class(integer_list_t), intent(inout) :: this
    integer,               target        :: value

    call this%add_copy(value)

  end subroutine integer_list_add_node

  ! ---------------------------------------------------------
  function integer_iterator_get_next(this) result(value)
    class(integer_iterator_t), intent(inout) :: this
    integer                                  :: value

    select type (ptr => this%get_next_ptr())
    type is (integer)
      value = ptr
    class default
      ASSERT(.false.)
    end select

  end function integer_iterator_get_next

end module linked_list_oct_m

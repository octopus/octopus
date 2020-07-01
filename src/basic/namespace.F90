!! Copyright (C) 2019 M. Oliveira, S. Ohlmann
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

module namespace_oct_m
  use global_oct_m
  use mpi_oct_m
  implicit none

  private
  public :: namespace_t, &
            global_namespace

  integer, parameter, public :: MAX_NAMESPACE_LEN = 128

  type :: namespace_t
    private
    character(len=MAX_NAMESPACE_LEN) :: name = ""
    type(namespace_t), pointer, public :: parent => NULL()
  contains
    procedure :: get => namespace_get
    procedure :: len => namespace_len
    generic   :: operator(==) => equal
    procedure :: equal => namespace_equal
    procedure :: contains => namespace_contains
  end type namespace_t

  interface namespace_t
    procedure namespace_constructor
  end interface namespace_t

  type(namespace_t) :: global_namespace

contains

  !---------------------------------------------------------
  !> Create namespace from name. If parent is present, the new namespace will be
  !! a child of it. It is also possible to create a namespace with several
  !! ancestors on the fly by providing a dot separated list of names. A
  !! different delimiter can be specified with the delimiter optional argument.
  recursive type(namespace_t) function namespace_constructor(name, parent, delimiter)
    character(len=*),                    intent(in) :: name
    type(namespace_t), optional, target, intent(in) :: parent
    character(len=1),  optional,         intent(in) :: delimiter

    integer :: total_len, parent_len, n_start
    character(len=1) :: delimiter_

    if (present(delimiter)) then
      delimiter_ = delimiter
    else
      delimiter_ = '.'
    end if

    ! Calculate total length of namespace, including the parent
    total_len = len_trim(name)
    if (present(parent)) then
      parent_len = parent%len()
      if (parent_len > 0) then
        total_len = total_len + parent_len + 1
      end if
    end if

    ! If total length is too large, stop and explain the reason
    if (total_len > MAX_NAMESPACE_LEN) then
      write(stderr,'(a)') '*** Fatal Error (description follows)'
      write(stderr,'(a)') 'Trying to create the following namespace:'
      if (present(parent)) then
        if (parent%len() > 0) then
          write(stderr,'(a)') trim(parent%get()) // delimiter_ // name
        end if
      else
        write(stderr,'(a)') name
      end if
      write(stderr,'(a,i4,a)') 'but namespaces are limited to ', MAX_NAMESPACE_LEN, ' characters'
#ifdef HAVE_MPI
      if(mpi_world%comm /= -1) call MPI_Abort(mpi_world%comm, 999, mpi_err)
#endif
      stop
    end if

    ! Find last delimiter in name
    n_start = index(name, delimiter_, back=.true.)

    if (n_start == 0) then
      ! There are no implicit parents in the name

      ! We do not allow the creation of empty namespaces, as that might lead to ambiguous paths
      ASSERT(len_trim(name) > 0)

      namespace_constructor%name = name
      if (present(parent)) then
        namespace_constructor%parent => parent
      else
        nullify(namespace_constructor%parent)
      end if

    else
      ! There are implicit parents in the name, so we create them recursively
      namespace_constructor%name = name(n_start+1:len(name))
      allocate(namespace_constructor%parent)
      namespace_constructor%parent = namespace_t(name(1:n_start-1), parent)
    end if

  end function namespace_constructor

  ! ---------------------------------------------------------
  recursive function namespace_get(this, delimiter) result(name)
    class(namespace_t),           intent(in) :: this
    character(len=1),   optional, intent(in) :: delimiter
    character(len=MAX_NAMESPACE_LEN) :: name

    character(len=1) :: delimiter_

    if (present(delimiter)) then
      delimiter_ = delimiter
    else
      delimiter_ = '.'
    end if

    name = this%name
    if (associated(this%parent)) then
      if (len_trim(name) > 0 .and. this%parent%len() > 0) then
        name = trim(this%parent%get(delimiter_)) // delimiter_ // trim(name)
      end if
    end if

  end function namespace_get

  ! ---------------------------------------------------------
  pure recursive function namespace_len(this)
    class(namespace_t), intent(in) :: this
    integer :: namespace_len

    integer :: parent_len

    namespace_len = len_trim(this%name)
    if (associated(this%parent)) then
      parent_len = this%parent%len()
      if (parent_len > 0) then
        namespace_len = namespace_len + parent_len + 1
      end if
    end if

  end function namespace_len

  ! ---------------------------------------------------------
  elemental logical function namespace_equal(lhs, rhs)
    class(namespace_t), intent(in) :: lhs
    class(namespace_t), intent(in) :: rhs

    namespace_equal = lhs%name == rhs%name

  end function namespace_equal

  ! ---------------------------------------------------------
  logical function namespace_contains(this, namespace)
    class(namespace_t), target, intent(in) :: this
    type(namespace_t),  target, intent(in) :: namespace

    logical :: found_ancestor
    type(namespace_t),  pointer :: this_ancestor, namespace_ancestor

    ! Find if there is a common ancestor
    found_ancestor = .false.
    namespace_ancestor => namespace
    do while (associated(namespace_ancestor))
      found_ancestor = namespace_ancestor == this
      if (found_ancestor) exit
      namespace_ancestor => namespace_ancestor%parent
    end do

    if (found_ancestor) then
      ! Check if the remaining ancestors are also equal
      namespace_contains = .true.
      this_ancestor => this
      do while (associated(namespace_ancestor%parent) .and. associated(this_ancestor%parent))
        namespace_contains = namespace_ancestor%parent == this_ancestor%parent
        if (.not. namespace_contains) exit
        this_ancestor => this_ancestor%parent
        namespace_ancestor => namespace_ancestor%parent
      end do
    else
      ! We did not find a common ancestor
      namespace_contains = .false.
    end if

  end function namespace_contains

end module namespace_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

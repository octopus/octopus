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

  integer, parameter :: MAX_NAMESPACE_LEN = 128

  type :: namespace_t
    private
    character(len=MAX_NAMESPACE_LEN) :: name = ""
    type(namespace_t), pointer :: parent => NULL()
  contains
    procedure :: get => namespace_get
    procedure :: len => namespace_len
  end type namespace_t

  interface namespace_t
    procedure namespace_init
  end interface namespace_t

  type(namespace_t) :: global_namespace

contains

  ! ---------------------------------------------------------
  type(namespace_t) function namespace_init(name, parent)
    character(len=*),                    intent(in) :: name
    type(namespace_t), optional, target, intent(in) :: parent

    integer :: total_len

    ! Calculate total length of namespace, including the parent
    total_len = len_trim(name)
    if (present(parent)) then
      if (parent%len() > 0) then
        total_len = total_len + parent%len() + 1
      end if
    end if

    ! If total length is too large, stop and explain the reason
    if (total_len > MAX_NAMESPACE_LEN) then
      write(stderr,'(a)') '*** Fatal Error (description follows)'
      write(stderr,'(a)') 'Trying to create the following namespace:'
      if (present(parent)) then
        if (parent%len() > 0) then
          write(stderr,'(a)') trim(parent%get()) // "." // name
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

    ! Now initialize the type
    namespace_init%name = name
    if (present(parent)) then
      namespace_init%parent => parent
    else
      nullify(namespace_init%parent)
    end if

  end function namespace_init

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

    name = ""
    if (associated(this%parent)) then
      if (this%parent%len() > 0) then
        name = trim(this%parent%get()) // delimiter_
      end if
    end if
    name = trim(name)//this%name

  end function namespace_get

  ! ---------------------------------------------------------
  pure recursive function namespace_len(this)
    class(namespace_t), intent(in) :: this
    integer :: namespace_len

    namespace_len = len_trim(this%name)
    if (associated(this%parent)) then
      if (this%parent%len() > 0) then
        namespace_len = namespace_len + this%parent%len() + 1
      end if
    end if

  end function namespace_len

end module namespace_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

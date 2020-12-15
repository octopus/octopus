!! Copyright (C)  2020 M. Oliveira
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

module algorithm_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                  &
    algorithmic_operation_t, &
    algorithm_t,             &
    algorithm_iterator_t

  integer, parameter, public :: ALGO_LABEL_LEN = 50

  type :: algorithmic_operation_t
    character(len=ALGO_LABEL_LEN) :: id !< Operation identifier. We use a string instead of an integer to minimize the chance of having duplicated identifiers.
    character(len=ALGO_LABEL_LEN) :: label !< Label describing what the code is doing when performing this operation.
  end type algorithmic_operation_t

  !> An algorithm is a list of algorithmic operations
  type, extends(linked_list_t) :: algorithm_t
    private
  contains
    procedure :: add_operation => algorithm_add_operation
  end type algorithm_t

  !> Iterator to loop over algorithmic operations of an algorithm
  type, extends(linked_list_iterator_t) :: algorithm_iterator_t
    private
  contains
    procedure :: get_next => algorithm_iterator_get_next
  end type algorithm_iterator_t

contains

  ! ---------------------------------------------------------
  subroutine algorithm_add_operation(this, operation)
    class(algorithm_t),             intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    PUSH_SUB(algorithm_add_operation)

    select type (operation)
    class is (algorithmic_operation_t)
      call this%add_copy(operation)
    class default
      ASSERT(.false.)
    end select

    POP_SUB(algorithm_add_operation)
  end subroutine algorithm_add_operation

  ! ---------------------------------------------------------
  function algorithm_iterator_get_next(this) result(operation)
    class(algorithm_iterator_t),  intent(inout) :: this
    type(algorithmic_operation_t)               :: operation

    PUSH_SUB(algorithm_iterator_get_next)

    select type (ptr => this%get_next_ptr())
    class is (algorithmic_operation_t)
      operation = ptr
    class default
      ASSERT(.false.)
    end select

    POP_SUB(algorithm_iterator_get_next)
  end function algorithm_iterator_get_next

end module algorithm_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

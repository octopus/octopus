!! Copyright (C)  2019 N. Tancogne-Dejean
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

module propagator_abst_oct_m
  use global_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                            &
    propagator_abst_t

  type, abstract :: propagator_abst_t
    private

    type(linked_list_t), public :: list

    FLOAT, public :: internal_time
    FLOAT, public :: dt
    logical :: step_done

  contains
    !Below are the list of operations that needs to be implemented
    procedure :: get_td_operation => propagator_get_tdop
    procedure :: step_is_done => propagator_step_is_done
    procedure :: rewind => propagator_rewind
    procedure :: finished => propagator_finished
    procedure(propagator_do_td_op), deferred :: do_td_op
  end type propagator_abst_t

  abstract interface
    subroutine propagator_do_td_op(this)
      import propagator_abst_t
      class(propagator_abst_t), intent(inout) :: this
    end subroutine propagator_do_td_op
  end interface

contains

  subroutine propagator_rewind(this)
    class(propagator_abst_t), intent(inout) :: this

    PUSH_SUB(propagator_rewind)

    call this%list%rewind()
    this%step_done = .false.

    POP_SUB(propagator_rewind)
  end subroutine propagator_rewind

  subroutine propagator_finished(this)
    class(propagator_abst_t), intent(inout) :: this

    PUSH_SUB(propagator_finished)

    this%step_done = .true.

    POP_SUB(propagator_finished)
  end subroutine propagator_finished

  integer function propagator_get_tdop(this) result(tdop)
    class(propagator_abst_t), intent(in) :: this

    class(*), pointer :: current_ops

    PUSH_SUB(propagator_get_tdop)

    current_ops => this%list%current()

    select type(current_ops)
    type is (integer)
      tdop = current_ops
    class default
      message(1) = "Corrupted list."
      call messages_fatal(1)
    end select

    POP_SUB(propagator_get_tdop)
  end function propagator_get_tdop

  logical pure function propagator_step_is_done(this) result(step_is_done)
    class(propagator_abst_t), intent(in) :: this

    step_is_done = this%step_done

  end function propagator_step_is_done

end module propagator_abst_oct_m


!!o, Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

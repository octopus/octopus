!! Copyright (C) 2020 Heiko Appel, M. Oliveira
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

module clock_oct_m
  use global_oct_m
  use loct_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  integer, parameter :: CLOCK_TICK = 1

  public ::                &
    clock_t,               &
    CLOCK_TICK

  type clock_t
    private
    integer :: tick        !< internal clock counter which is incremented by one when the clock is advanced
    FLOAT   :: time_step   !< physical simulation time increment which corresponds to a single clock tick
    FLOAT   :: time_       !< physical simulation time
  contains
    procedure :: print => clock_print                 !< print internal state of the clock
    procedure :: print_str => clock_print_str         !< print internal state of the clock to a string
    procedure :: set_time => clock_set_time           !< set the clock only to the time of a given input clock
    procedure :: copy => clock_copy                   !< set the clock to the state of a given input clock
    procedure :: get_tick => clock_get_tick           !< get value of internal clock counter
    procedure :: time => clock_time                   !< get the current physical simulation time of the clock
    procedure :: reset => clock_reset                 !< set the internal clock counter back to zero
    procedure :: clock_is_equal
    generic   :: operator(.eq.) => clock_is_equal
    procedure :: clock_is_different
    generic   :: operator(/=) => clock_is_different
    procedure :: clock_is_earlier
    generic   :: operator(.lt.) => clock_is_earlier
    procedure :: clock_is_later
    generic   :: operator(.gt.) => clock_is_later
    procedure :: clock_is_equal_or_earlier
    generic   :: operator(.le.) => clock_is_equal_or_earlier
    procedure :: clock_is_equal_or_later
    generic   :: operator(.ge.) => clock_is_equal_or_later
    procedure :: clock_copy
    generic   :: assignment(=) => clock_copy
    procedure :: clock_add_tick
    generic   :: operator(+) => clock_add_tick
    procedure :: clock_subtract_tick
    generic   :: operator(-) => clock_subtract_tick
  end type clock_t

  interface clock_t
    module procedure clock_init
  end interface clock_t

contains

  ! ---------------------------------------------------------
  !> Initialize the clock with a given label and associated physical time step.
  !! The internal clock counter starts at zero or if the optional argument initial_tick is given
  !! at the value of initial_tick.
  type(clock_t) function clock_init(time_step, initial_tick) result(this)
    FLOAT,   optional, intent(in) :: time_step
    integer, optional, intent(in) :: initial_tick

    PUSH_SUB(clock_init)

    this%tick = optional_default(initial_tick, 0)
    this%time_step = optional_default(time_step, M_ZERO)
    if (this%time_step <= M_ZERO) then
      this%time_ = M_ZERO
    else
      this%time_ = this%tick*this%time_step
    end if

    POP_SUB(clock_init)
  end function clock_init

  ! ---------------------------------------------------------
  function clock_print_str(this) result(clock_string)
    class(clock_t), intent(in)    :: this
    character(len=65)             :: clock_string

    PUSH_SUB(clock_print_str)

    write(clock_string,'(A7,F16.6,A,I8.8,A)') &
        '[Clock:',                   &
        this%time_,                  &
        '|',                         &
        this%tick,                   &
        ']'

    POP_SUB(clock_print_str)
  end function clock_print_str

  ! ---------------------------------------------------------
  subroutine clock_print(this)
    class(clock_t), intent(in) :: this

    PUSH_SUB(clock_print)

    message(1) = this%print_str()
    call messages_info(1)

    POP_SUB(clock_print)
  end subroutine clock_print

  ! ---------------------------------------------------------
  subroutine clock_set_time(this, new)
    class(clock_t), intent(inout) :: this
    class(clock_t), intent(in)    :: new

    logical :: commensurable
    integer :: this_granularity, new_granularity

    PUSH_SUB(clock_set_time)

    if (this%time_step > M_ZERO .and. new%time_step > M_ZERO) then
      if (this%time_step >= new%time_step) then
        commensurable = ceiling(this%time_step/new%time_step) == floor(this%time_step/new%time_step)
        this_granularity = ceiling(this%time_step/new%time_step)
        new_granularity = 1
      else
        commensurable = ceiling(new%time_step/this%time_step) == floor(new%time_step/this%time_step)
        this_granularity = 1
        new_granularity = ceiling(new%time_step/this%time_step)
      end if

      if (.not. commensurable) then
        message(1) = 'Cannot set clock new time, as it is not commensurable with clock time-step.'
        call messages_fatal(1)
      end if

      this%tick = (new%tick * new_granularity) / this_granularity
      this%time_ = this%tick*this%time_step
    else
      this%time_ = new%time_
    end if

    POP_SUB(clock_set_time)
  end subroutine clock_set_time

  ! ---------------------------------------------------------
  subroutine clock_copy(this, clock_in)
    class(clock_t), intent(in)    :: clock_in
    class(clock_t), intent(inout) :: this

    PUSH_SUB(clock_copy)

    this%tick = clock_in%tick
    this%time_step = clock_in%time_step
    this%time_ = clock_in%time_

    POP_SUB(clock_copy)
  end subroutine clock_copy

  ! ---------------------------------------------------------
  type(clock_t) function clock_add_tick(clock, tick) result(new_clock)
    class(clock_t), intent(in) :: clock
    integer,        intent(in) :: tick

    PUSH_SUB(clock_add_tick)

    new_clock = clock
    new_clock%tick = new_clock%tick + tick
    new_clock%time_ = new_clock%tick*new_clock%time_step

    POP_SUB(clock_add_tick)
  end function clock_add_tick

  ! ---------------------------------------------------------
  type(clock_t) function clock_subtract_tick(clock, tick) result(new_clock)
    class(clock_t), intent(in) :: clock
    integer,        intent(in) :: tick

    PUSH_SUB(clock_subtract_tick)

    new_clock = clock
    new_clock%tick = new_clock%tick - tick
    new_clock%time_ = new_clock%tick*new_clock%time_step

    POP_SUB(clock_subtract_tick)
  end function clock_subtract_tick

  ! ---------------------------------------------------------
  integer function clock_get_tick(this) result(tick)
    class(clock_t), intent(in) :: this

    PUSH_SUB(clock_get_tick)

    tick = this%tick

    POP_SUB(clock_get_tick)
  end function clock_get_tick

  ! ---------------------------------------------------------
  FLOAT function clock_time(this)
    class(clock_t), intent(in) :: this

    PUSH_SUB(clock_time)

    clock_time = this%time_

    POP_SUB(clock_time)
  end function clock_time

  ! ---------------------------------------------------------
  subroutine clock_reset(this)
    class(clock_t), intent(inout) :: this

    PUSH_SUB(clock_reset)

    this%tick = 0
    this%time_ = M_ZERO

    POP_SUB(clock_reset)
  end subroutine clock_reset

  ! ---------------------------------------------------------
  logical function clock_is_earlier(clock_a, clock_b) result(is_earlier)
    class(clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(clock_is_earlier)

    is_earlier = clock_a%time_ < clock_b%time_

    POP_SUB(clock_is_earlier)
  end function clock_is_earlier

  ! ---------------------------------------------------------
  logical function clock_is_later(clock_a, clock_b) result(is_later)
    class(clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(clock_is_later)

    is_later = clock_a%time_ > clock_b%time_

    POP_SUB(clock_is_later)
  end function clock_is_later

  ! ---------------------------------------------------------
  logical function clock_is_equal_or_earlier(clock_a, clock_b) result(is_earlier)
    class(clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(clock_is_equal_or_earlier)

    is_earlier = clock_a%time_ <= clock_b%time_

    POP_SUB(clock_is_equal_or_earlier)
  end function clock_is_equal_or_earlier

  ! ---------------------------------------------------------
  logical function clock_is_equal_or_later(clock_a, clock_b) result(is_later)
    class(clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(clock_is_equal_or_later)

    is_later = clock_a%time_ >= clock_b%time_

    POP_SUB(clock_is_equal_or_later)
  end function clock_is_equal_or_later

  ! ---------------------------------------------------------
  logical function clock_is_equal(clock_a, clock_b) result(are_equal)
    class(clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(clock_is_equal)

    are_equal = clock_a%time_ == clock_b%time_

    POP_SUB(clock_is_equal)
  end function clock_is_equal

  ! ---------------------------------------------------------
  logical function clock_is_different(clock_a, clock_b) result(are_diff)
    class(clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(clock_is_different)

    are_diff = clock_a%time_ /= clock_b%time_

    POP_SUB(clock_is_different)
  end function clock_is_different

end module clock_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

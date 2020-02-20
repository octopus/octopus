!! Copyright (C) 2020 Heiko Appel
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
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                &
    clock_t,    &
    operator(.eq.),        &
    operator(.lt.),        &
    operator(.gt.),        &
    operator(.le.),        &
    operator(.ge.)

  type clock_t
    private
    integer :: clock_tick
    integer :: granularity
    FLOAT   :: time_step
    type(namespace_t) :: namespace

  contains
    procedure :: print => clock_print
    procedure :: set => clock_set
    procedure :: get_tick => clock_get_tick
    procedure :: get_sim_time => clock_get_sim_time
    procedure :: increment => clock_increment
    procedure :: decrement => clock_decrement
    procedure :: reset => clock_reset
    procedure :: is_later_with_step => clock_is_later_with_step
  end type clock_t

  interface clock_t
    module procedure clock_init
  end interface clock_t

  interface assignment(=)
    procedure clock_set
  end interface

  interface operator(.eq.)
    procedure clock_is_equal
  end interface

  interface operator(.lt.)
    procedure clock_is_earlier
  end interface

  interface operator(.gt.)
    procedure clock_is_later
  end interface

  interface operator(.le.)
    procedure clock_is_equal_or_earlier
  end interface

  interface operator(.ge.)
    procedure clock_is_equal_or_later
  end interface


contains

  ! ---------------------------------------------------------
  type(clock_t) function clock_init(namespace, time_step, smallest_algo_dt, initial_tick) result(this)
    type(namespace_t), intent(in) :: namespace
    FLOAT,             intent(in) :: time_step, smallest_algo_dt
    integer, optional             :: initial_tick

    PUSH_SUB(clock_init)

    ! this needs to be adapted later on for a more sophisticated handling of the clock namespace
    this%namespace = namespace
    this%clock_tick = optional_default(initial_tick, 0)

    if (ceiling(time_step/smallest_algo_dt) == floor(time_step/smallest_algo_dt)) then
      this%granularity = ceiling(time_step/smallest_algo_dt)
    else
      message(1) = 'Timesteps of the clocks are not comensurate.'
      message(2) = 'Please adapt the time steps of your subsystems to make them compatible.'
      call messages_fatal(2)
    endif

    this%time_step = time_step

    POP_SUB(clock_init)
  end function clock_init

  ! ---------------------------------------------------------
  subroutine clock_print(this)
    class(clock_t), intent(in) :: this
    
    PUSH_SUB(clock_print)

    write(message(1),'(A7,A16,A,I8.8,A,I8.8,A,I8.8,A,F8.6,A)') &
        '[Clock:',                         &
        trim(this%namespace%get()),        &
        '|',                               &
        this%clock_tick,                   &
        '|',                               &
        this%granularity,                  &
        '|',                               &
        this%clock_tick*this%granularity,  &
        '|',                               &
        this%time_step,                    &
        ']'
    call messages_info(1)

    POP_SUB(clock_print)
  end subroutine clock_print

  ! ---------------------------------------------------------
  subroutine clock_set(this, clock_in)
    class(clock_t), intent(in) :: clock_in
    class(clock_t), intent(inout) :: this
    
    PUSH_SUB(clock_set)

    this%clock_tick = clock_in%clock_tick
    this%granularity = clock_in%granularity
    this%time_step = clock_in%time_step

    POP_SUB(clock_set)
  end subroutine clock_set

  ! ---------------------------------------------------------
  integer function clock_get_tick(this) result(current_global_tick)
    class(clock_t), intent(in) :: this

    PUSH_SUB(clock_get_tick)

    current_global_tick = this%clock_tick * this%granularity

    POP_SUB(clock_get_tick)
  end function clock_get_tick

  ! ---------------------------------------------------------
  FLOAT function clock_get_sim_time(this) result(current_time)
    class(clock_t), intent(in) :: this

    PUSH_SUB(clock_get_sim_time)

    current_time = this%clock_tick * this%time_step

    POP_SUB(clock_get_sim_time)
  end function clock_get_sim_time

  ! ---------------------------------------------------------
  subroutine clock_increment(this, steps)
    class(clock_t), intent(inout) :: this
    integer, optional,         intent(in)    :: steps

    PUSH_SUB(clock_increment)

    this%clock_tick = this%clock_tick + optional_default(steps, 1)

    POP_SUB(clock_increment)
  end subroutine clock_increment

  ! ---------------------------------------------------------
  subroutine clock_decrement(this, steps)
    class(clock_t), intent(inout) :: this
    integer, optional,         intent(in)    :: steps

    PUSH_SUB(clock_decrement)

    this%clock_tick = this%clock_tick - optional_default(steps, 1)

    POP_SUB(clock_decrement)
  end subroutine clock_decrement

  ! ---------------------------------------------------------
  subroutine clock_reset(this)
    class(clock_t), intent(inout) :: this

    PUSH_SUB(clock_reset)

    this%clock_tick = 0

    POP_SUB(clock_reset)
  end subroutine clock_reset

  ! ---------------------------------------------------------
  logical function clock_is_earlier(clock_a, clock_b) result(is_earlier)
    type(clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(clock_is_earlier)

    is_earlier = clock_a%get_tick() < clock_b%get_tick()

    POP_SUB(clock_is_earlier)
  end function clock_is_earlier

  ! ---------------------------------------------------------
  logical function clock_is_later(clock_a, clock_b) result(is_later)
    class(clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(clock_is_later)

    is_later = clock_a%get_tick() > clock_b%get_tick()

    POP_SUB(clock_is_later)
  end function clock_is_later

  ! ---------------------------------------------------------
  logical function clock_is_equal_or_earlier(clock_a, clock_b) result(is_earlier)
    type(clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(clock_is_equal_or_earlier)

    is_earlier = clock_a%get_tick() <= clock_b%get_tick()

    POP_SUB(clock_is_equal_or_earlier)
  end function clock_is_equal_or_earlier

  ! ---------------------------------------------------------
  logical function clock_is_equal_or_later(clock_a, clock_b) result(is_later)
    class(clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(clock_is_equal_or_later)

    is_later = clock_a%get_tick() >= clock_b%get_tick()

    POP_SUB(clock_is_equal_or_later)
  end function clock_is_equal_or_later

  ! ---------------------------------------------------------
  logical function clock_is_equal(clock_a, clock_b) result(are_equal)
    class(clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(clock_is_equal)

    are_equal = clock_a%get_tick() == clock_b%get_tick()

    POP_SUB(clock_is_equal)
  end function clock_is_equal

  ! ---------------------------------------------------------
  logical function clock_is_later_with_step(clock_a, clock_b) result(is_later_with_step)
    class(clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(clock_is_later_with_step)

    is_later_with_step = (clock_a%get_tick() + clock_a%granularity) > clock_b%get_tick()

    POP_SUB(clock_is_later_with_step)
  end function clock_is_later_with_step

end module clock_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

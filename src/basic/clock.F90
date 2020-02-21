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
    clock_t,               &
    operator(.eq.),        &
    operator(.lt.),        &
    operator(.gt.),        &
    operator(.le.),        &
    operator(.ge.)

  type clock_t
    private
    integer :: clock_tick  !< internal clock counter which is incremented by one when the clock is advanced
    integer :: granularity !< one clock tick corresponds to an advancement of 'granularity' steps on the finest scale in the system
    FLOAT   :: time_step   !< physical simulation time increment which corresponds to a single clock tick
    type(namespace_t) :: namespace !< namespace for labeling the clock

  contains
    procedure :: print => clock_print                 !< print internal state of the clock
    procedure :: print_str => clock_print_str         !< print internal state of the clock to a string
    procedure :: print_message => clock_print_message !< print internal state of the clock together with a given message
    procedure :: set_time => clock_set_time           !< set the clock only to the time of a given input clock
    procedure :: copy => clock_copy                   !< set the clock and the namespace to the state of a given input clock
    procedure :: get_tick => clock_get_tick           !< get value of internal clock counter
    procedure :: get_sim_time => clock_get_sim_time   !< get the current physical simulation time of the clock
    procedure :: increment => clock_increment         !< increment the internal clock counter by one or several steps
    procedure :: decrement => clock_decrement         !< decrement the internal clock counter by one or several steps
    procedure :: reset => clock_reset                 !< set the internal clock counter back to zero
    procedure :: is_later_with_step => clock_is_later_with_step !< compare two clocks and indicate if the first clock plus
                                                                !! one step is later than the second clock
  end type clock_t

  interface clock_t
    module procedure clock_init
  end interface clock_t

  interface assignment(=)
    procedure clock_copy
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
  !> Initialize the clock with the given namespace and associated physical time step.
  !! The internal clock counter starts at zero or if the optional argument initial_tick is given
  !! at the value of initial_tick.
  !! When several clocks are used simultaneously, then smallest_dt is the smallest time step of
  !! all the clocks, or equivalently the time step of the fastest ticking clock. If only a single
  !! clock is used or all clocks run with the same time step, then smallest_dt is identical to
  !! time_step. 
  !! All clocks are compared on the finest scale. The granularity tells how many steps on the finest 
  !! scale an increase of the clock tick corresponds to. In other words, it tells in which 'units'
  !! an increment of the internal clock counter is counting the time on the finest scale.
  !! So far the implementation assumes that all clocks in the system have commensurable time steps.
  type(clock_t) function clock_init(namespace, time_step, smallest_dt, initial_tick) result(this)
    type(namespace_t), intent(in) :: namespace
    FLOAT,             intent(in) :: time_step, smallest_dt
    integer, optional             :: initial_tick

    PUSH_SUB(clock_init)

    ! this needs to be adapted later on for a more sophisticated handling of the clock namespace
    this%namespace = namespace
    this%clock_tick = optional_default(initial_tick, 0)

    if (ceiling(time_step/smallest_dt) == floor(time_step/smallest_dt)) then
      this%granularity = ceiling(time_step/smallest_dt)
    else
      message(1) = 'Timesteps of the clocks are not commensurable.'
      message(2) = 'Please adapt the time steps of your subsystems to make them compatible.'
      call messages_fatal(2)
    endif

    this%time_step = time_step

    POP_SUB(clock_init)
  end function clock_init

  ! ---------------------------------------------------------
  function clock_print_str(this) result(clock_string)
    class(clock_t), intent(in) :: this
    character(len=64)          :: clock_string

    PUSH_SUB(clock_print_str)

    write(clock_string,'(A7,A11,A,F16.6,A,I8.8,A,I8.8,A,I8.8,A)') &
        '[Clock:',                         &
        trim(this%namespace%get()),        &
        '|',                               &
        this%time_step*this%clock_tick,    &
        '|',                               &
        this%clock_tick,                   &
        '|',                               &
        this%granularity,                  &
        '|',                               &
        this%clock_tick*this%granularity,  &
        ']'

    POP_SUB(clock_print_str)
  end function clock_print_str

  ! ---------------------------------------------------------
  function clock_print_message(this, message) result(clock_message)
    class(clock_t),   intent(in) :: this
    character(len=*), intent(in) :: message
    character(len=256)           :: clock_message

    PUSH_SUB(clock_print_message)

    clock_message = trim(this%print_str()) // trim(message)

    POP_SUB(clock_print_message)
  end function clock_print_message

  ! ---------------------------------------------------------
  subroutine clock_print(this)
    class(clock_t), intent(in) :: this

    PUSH_SUB(clock_print)

    message(1) = this%print_str()
    call messages_info(1)

    POP_SUB(clock_print)
  end subroutine clock_print

  ! ---------------------------------------------------------
  subroutine clock_set_time(this, clock_in)
    class(clock_t), intent(in)    :: clock_in
    class(clock_t), intent(inout) :: this

    PUSH_SUB(clock_set_time)

    this%clock_tick = clock_in%clock_tick
    this%granularity = clock_in%granularity
    this%time_step = clock_in%time_step

    POP_SUB(clock_set_time)
  end subroutine clock_set_time

  ! ---------------------------------------------------------
  subroutine clock_copy(this, clock_in)
    class(clock_t), intent(in)    :: clock_in
    class(clock_t), intent(inout) :: this

    PUSH_SUB(clock_copy)

    call this%set_time(clock_in)
    this%namespace = clock_in%namespace

    POP_SUB(clock_copy)
  end subroutine clock_copy

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
    integer, optional, intent(in) :: steps

    PUSH_SUB(clock_increment)

    this%clock_tick = this%clock_tick + optional_default(steps, 1)

    POP_SUB(clock_increment)
  end subroutine clock_increment

  ! ---------------------------------------------------------
  subroutine clock_decrement(this, steps)
    class(clock_t), intent(inout) :: this
    integer, optional, intent(in) :: steps

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
  !> This function returns true if taking clock_a and adding one
  !! time step is later in time than the current time of clock_b
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

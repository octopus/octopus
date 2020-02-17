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

module simulation_clock_oct_m
  use global_oct_m
  use loct_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                &
    simulation_clock_t,    &
    operator(.eq.),        &
    operator(.lt.),        &
    operator(.gt.),        &
    operator(.le.),        &
    operator(.ge.)

  type simulation_clock_t
    integer :: clock_tick
    integer :: granularity
    FLOAT   :: time_step
    integer :: sec_since_epoch
    integer :: usec_since_epoch
    type(namespace_t) :: namespace


  contains
    procedure :: print => simulation_clock_print
    procedure :: set => simulation_clock_set
    procedure :: get_tick => simulation_clock_get_tick
    procedure :: get_sim_time => simulation_clock_get_sim_time
    procedure :: increment => simulation_clock_increment
    procedure :: decrement => simulation_clock_decrement
    procedure :: reset => simulation_clock_reset
    procedure :: is_later_with_step => simulation_clock_is_later_with_step
  end type simulation_clock_t

  interface simulation_clock_t
    module procedure simulation_clock_init
  end interface simulation_clock_t

  interface assignment(=)
    procedure simulation_clock_set
  end interface

  interface operator(.eq.)
    procedure simulation_clock_is_equal
  end interface

  interface operator(.lt.)
    procedure simulation_clock_is_earlier
  end interface

  interface operator(.gt.)
    procedure simulation_clock_is_later
  end interface

  interface operator(.le.)
    procedure simulation_clock_is_equal_or_earlier
  end interface

  interface operator(.ge.)
    procedure simulation_clock_is_equal_or_later
  end interface


contains

  ! ---------------------------------------------------------
  type(simulation_clock_t) function simulation_clock_init(namespace, time_step, smallest_algo_dt) result(this)
    type(namespace_t), intent(in)  :: namespace
    FLOAT, intent(in)              :: time_step, smallest_algo_dt

    integer :: epoch_sec, epoch_usec

    PUSH_SUB(simulation_clock_init)

    this%namespace = namespace
    this%clock_tick = 0
    this%granularity = ceiling(time_step/smallest_algo_dt)
    this%time_step = time_step

    call loct_gettimeofday(epoch_sec, epoch_usec)
    this%sec_since_epoch = epoch_sec
    this%usec_since_epoch = epoch_usec

    POP_SUB(simulation_clock_init)
  end function simulation_clock_init

  ! ---------------------------------------------------------
  subroutine simulation_clock_print(this)
    class(simulation_clock_t), intent(inout) :: this
    
    PUSH_SUB(simulation_clock_print)

    write(message(1),'(A7,A16,A,I8.8,A,I8.8,A,I8.8,A,F8.6,A,I10.10,A,I6.6,A)') &
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
        '|',                               &
        this%sec_since_epoch,              &
        '|',                               &
        this%usec_since_epoch,             &
        ']'
    call messages_info(1)

    POP_SUB(simulation_clock_print)
  end subroutine simulation_clock_print

  ! ---------------------------------------------------------
  subroutine simulation_clock_set(this, clock_in)
    class(simulation_clock_t), intent(in) :: clock_in
    class(simulation_clock_t), intent(inout) :: this
    
    PUSH_SUB(simulation_clock_set)

    this%clock_tick = clock_in%clock_tick
    this%granularity = clock_in%granularity
    this%time_step = clock_in%time_step
    this%sec_since_epoch = clock_in%sec_since_epoch
    this%usec_since_epoch = clock_in%usec_since_epoch

    POP_SUB(simulation_clock_set)
  end subroutine simulation_clock_set

  ! ---------------------------------------------------------
  integer function simulation_clock_get_tick(this) result(current_global_tick)
    class(simulation_clock_t), intent(in) :: this

    PUSH_SUB(simulation_clock_get_tick)

    current_global_tick = this%clock_tick*this%granularity

    POP_SUB(simulation_clock_get_tick)
  end function simulation_clock_get_tick

  ! ---------------------------------------------------------
  FLOAT function simulation_clock_get_sim_time(this) result(current_time)
    class(simulation_clock_t), intent(in) :: this

    PUSH_SUB(simulation_clock_get_sim_time)

    current_time = this%clock_tick*this%time_step

    POP_SUB(simulation_clock_get_sim_time)
  end function simulation_clock_get_sim_time

  ! ---------------------------------------------------------
  subroutine simulation_clock_increment(this)
    class(simulation_clock_t), intent(inout) :: this

    PUSH_SUB(simulation_clock_increment)

    this%clock_tick = this%clock_tick + 1

    POP_SUB(simulation_clock_increment)
  end subroutine simulation_clock_increment

  ! ---------------------------------------------------------
  subroutine simulation_clock_decrement(this)
    class(simulation_clock_t), intent(inout) :: this

    PUSH_SUB(simulation_clock_decrement)

    this%clock_tick = this%clock_tick - 1

    POP_SUB(simulation_clock_decrement)
  end subroutine simulation_clock_decrement

  ! ---------------------------------------------------------
  subroutine simulation_clock_reset(this)
    class(simulation_clock_t), intent(inout) :: this

    PUSH_SUB(simulation_clock_reset)

    this%clock_tick = 0

    POP_SUB(simulation_clock_reset)
  end subroutine simulation_clock_reset

  ! ---------------------------------------------------------
  logical function simulation_clock_is_earlier(clock_a, clock_b) result(is_earlier)
    type(simulation_clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(simulation_clock_is_earlier)

    if(clock_a%get_tick() < clock_b%get_tick()) then
        is_earlier = .true.
    else
        is_earlier = .false.
    end if

    POP_SUB(simulation_clock_is_earlier)
  end function simulation_clock_is_earlier

  ! ---------------------------------------------------------
  logical function simulation_clock_is_later(clock_a, clock_b) result(is_earlier)
    class(simulation_clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(simulation_clock_is_later)

    is_earlier = simulation_clock_is_earlier(clock_b, clock_a)

    POP_SUB(simulation_clock_is_later)
  end function simulation_clock_is_later

  ! ---------------------------------------------------------
  logical function simulation_clock_is_equal_or_earlier(clock_a, clock_b) result(is_earlier)
    type(simulation_clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(simulation_clock_is_earlier)

    if(clock_a%get_tick() <= clock_b%get_tick()) then
        is_earlier = .true.
    else
        is_earlier = .false.
    end if

    POP_SUB(simulation_clock_is_earlier)
  end function simulation_clock_is_equal_or_earlier

  ! ---------------------------------------------------------
  logical function simulation_clock_is_equal_or_later(clock_a, clock_b) result(is_earlier)
    class(simulation_clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(simulation_clock_is_later)

    is_earlier = simulation_clock_is_earlier(clock_b, clock_a)

    POP_SUB(simulation_clock_is_later)
  end function simulation_clock_is_equal_or_later

  ! ---------------------------------------------------------
  logical function simulation_clock_is_equal(clock_a, clock_b) result(are_equal)
    class(simulation_clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(simulation_clock_is_equal)

    if(clock_a%get_tick() == clock_b%get_tick()) then
        are_equal = .true.
    else
        are_equal = .false.
    end if

    POP_SUB(simulation_clock_is_equal)
  end function simulation_clock_is_equal

  ! ---------------------------------------------------------
  logical function simulation_clock_is_later_with_step(clock_a, clock_b) result(is_later_with_step)
    class(simulation_clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(simulation_clock_is_later_with_step)

    if(clock_a%get_tick() + clock_a%granularity > clock_b%get_tick()) then
        is_later_with_step = .true.
    else
        is_later_with_step = .false.
    end if

    POP_SUB(simulation_clock_is_later_with_step)
  end function simulation_clock_is_later_with_step

end module simulation_clock_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                &
    simulation_clock_t

  type simulation_clock_t
    integer :: clock_tick
    integer :: granularity
    FLOAT   :: time_step
    integer :: time_since_epoch

  contains
    procedure :: get_tick => simulation_clock_get_tick
    procedure :: get_sim_time => simulation_clock_get_sim_time
    procedure :: increment => simulation_clock_increment
    procedure :: decrement => simulation_clock_decrement
    procedure :: reset => simulation_clock_reset
    procedure :: set => simulation_clock_set
    procedure :: is_earlier => simulation_clock_is_earlier
    procedure :: is_later => simulation_clock_is_later
    procedure :: is_equal => simulation_clock_is_equal
    procedure :: is_later_with_step => simulation_clock_is_later_with_step
  end type simulation_clock_t

  interface simulation_clock_t
    module procedure simulation_clock_init
  end interface simulation_clock_t

  interface operator(==)
    procedure simulation_clock_is_equal
  end interface

  interface assignment(=)
    procedure simulation_clock_set
  end interface

  interface operator(<)
    procedure simulation_clock_is_earlier
  end interface

  interface operator(>)
    procedure simulation_clock_is_later
  end interface

contains

  ! ---------------------------------------------------------
  type(simulation_clock_t) function simulation_clock_init(time_step, smallest_algo_dt) result(this)
    FLOAT, intent(in)   :: time_step, smallest_algo_dt

    PUSH_SUB(simulation_clock_init)

    this%clock_tick = 0
    this%granularity = ceiling(time_step/smallest_algo_dt)
    this%time_step = time_step
    this%time_since_epoch = 0 ! Fixme: get time since epoch here from OS

    POP_SUB(simulation_clock_init)
  end function simulation_clock_init

  ! ---------------------------------------------------------
  subroutine simulation_clock_set(clock_out, clock_in)
    class(simulation_clock_t), intent(in) :: clock_in
    class(simulation_clock_t), intent(out) :: clock_out
    
    PUSH_SUB(simulation_clock_set)

    clock_out%clock_tick = clock_in%clock_tick
    clock_out%granularity = clock_in%granularity
    clock_out%time_step = clock_in%time_step
    clock_out%time_since_epoch = clock_in%time_since_epoch

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
    class(simulation_clock_t), intent(in) :: clock_a, clock_b

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

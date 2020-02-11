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
    FLOAT   :: time_step
    integer :: time_since_epoch

  contains
    procedure :: get_tick => simulation_clock_get_tick
    procedure :: get_sim_time => simulation_clock_get_sim_time
    procedure :: increment => simulation_clock_increment
    procedure :: reset => simulation_clock_reset
    procedure :: is_smaller => simulation_clock_is_smaller
    procedure :: is_larger => simulation_clock_is_larger
    procedure :: is_equal => simulation_clock_is_equal

  end type simulation_clock_t

  interface simulation_clock_t
    module procedure simulation_clock_init
  end interface simulation_clock_t

contains

  ! ---------------------------------------------------------
  type(simulation_clock_t) function simulation_clock_init(time_step) result(this)
    FLOAT, intent(in)   :: time_step

    PUSH_SUB(simulation_clock_init)

    this%clock_tick = 0
    this%time_step = time_step
    this%time_since_epoch = 0 ! Fixme: get time since epoch here from OS

    POP_SUB(simulation_clock_init)
  end function simulation_clock_init

  ! ---------------------------------------------------------
  integer function simulation_clock_get_tick(this) result(current_tick)
    class(simulation_clock_t), intent(in) :: this

    PUSH_SUB(simulation_clock_get_tick)

    current_tick = this%clock_tick

    POP_SUB(simulation_clock_get_tick)
  end function simulation_clock_get_tick

  ! ---------------------------------------------------------
  FLOAT function simulation_clock_get_sim_time(this) result(current_time)
    class(simulation_clock_t), intent(in) :: this

    PUSH_SUB(simulation_clock_get_tick)

    current_time = this%clock_tick*this%time_step

    POP_SUB(simulation_clock_get_tick)
  end function simulation_clock_get_sim_time

  ! ---------------------------------------------------------
  subroutine simulation_clock_increment(this)
    class(simulation_clock_t), intent(inout) :: this

    PUSH_SUB(simulation_clock_update)

    this%clock_tick = this%clock_tick + 1

    POP_SUB(simulation_clock_update)
  end subroutine simulation_clock_increment

  ! ---------------------------------------------------------
  subroutine simulation_clock_reset(this)
    class(simulation_clock_t), intent(inout) :: this

    PUSH_SUB(simulation_clock_reset)

    this%clock_tick = 0

    POP_SUB(simulation_clock_reset)
  end subroutine simulation_clock_reset

  ! ---------------------------------------------------------
  logical function simulation_clock_is_smaller(clock_a, clock_b) result(is_smaller)
    class(simulation_clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(simulation_clock_is_smaller)

    if(clock_a%get_tick() < clock_b%get_tick()) then
        is_smaller = .true.
    else
        is_smaller = .false.
    end if

    POP_SUB(simulation_clock_is_smaller)
  end function simulation_clock_is_smaller

  ! ---------------------------------------------------------
  logical function simulation_clock_is_larger(clock_a, clock_b) result(is_smaller)
    class(simulation_clock_t), intent(in) :: clock_a, clock_b

    PUSH_SUB(simulation_clock_is_larger)

    is_smaller = simulation_clock_is_smaller(clock_b, clock_a)

    POP_SUB(simulation_clock_is_larger)
  end function simulation_clock_is_larger

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

end module simulation_clock_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

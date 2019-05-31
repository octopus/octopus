!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

!> This module provices a simple timer class which can be used to trigger the writing of a restart
!! file before the requested CPU time is up.
!!
!! It allows to take into account the time required for one iteration and optionally a time margin
!! for completing the restart dump process.

#include "global.h"

module walltimer_oct_m

  use global_oct_m
  use loct_oct_m
  use messages_oct_m
  use parser_oct_m
  
  implicit none
  
  private
  
  FLOAT :: start_time       !< time the timer was started
  FLOAT :: last_tap         !< time of the last call to tap()
  FLOAT :: iteration_time   !< time difference of two calls to tap()
  FLOAT :: margin           !< additional time margin for writing the restart file
  FLOAT :: duration         !< time when the alarm should trigger
  
  logical :: active 
  logical :: auto_tap                  !< if .true., tap() is automatically called in every wakeUp() call.

  public :: walltimer_init, walltimer_end, walltimer_tap, walltimer_alarm
  
contains

  !> initialize the timer
  
  subroutine walltimer_init(auto)

    logical, optional, intent(IN) :: auto   !< automatically call walltimer_tap in walltimer_alarm() if .true.

    FLOAT  :: alarm_time, write_time

    PUSH_SUB(walltimer_init)

    start_time = M_ZERO
    last_tap = M_ZERO
    iteration_time = M_ZERO
    margin = M_ZERO

    active = .false.
    
    auto_tap = .true.
    if(present(auto)) auto_tap = auto


    ! The following have to be moved to the right place, after the names for the variables have been confirmed:

    !%Variable Walltime
    !%Type float
    !%Default 0
    !%Section Execution::IO
    !%Description
    !% Time in minutes before which the restart file will be written. This is to make sure that at least one restart
    !% file can be written before the code might be killed to to exceeding the given CPU time.
    !% If a finite time (in minutes) is specified, the code will write the restart file when the next
    !% iteration (plus the RestartWriteTime) would exceed the given time.
    !% A value less than 1 second (1/60 minutes) will disable the timer.
    !%End0.0
    call parse_Variable('Walltime', M_ZERO, alarm_time)
    call setAlarm(alarm_time*60.d0)
    
    !%Variable RestartWriteTime
    !%Type float
    !%Default 5
    !%Section Execution::IO
    !%Description
    !% The RestartWriteTime (in minutes) will be subtracted from the WallTime to allow time for writing the restart file.
    !% In huge calculations, this value should be increased.
    !%End
    call parse_Variable('RestartWriteTime', CNST(5.0), write_time)
    if(write_time > alarm_time/4) write_time = alarm_time/4
    call setMargin(write_time*60.d0)
    
    call start()
    
    POP_SUB(walltimer_init)

  end subroutine walltimer_init

  
  !> empty destructor
  
  subroutine walltimer_end()

    PUSH_SUB(walltimer_end)

    active = .false.

    POP_SUB(walltimer_end)

  end subroutine walltimer_end


  

  !> set alarm interval in seconds

  subroutine setAlarm(time)

    double precision :: time

    PUSH_SUB(setAlarm)

    duration = time

    POP_SUB(setAlarm)        

  end subroutine setAlarm


  !> set safty margin in seconds

  subroutine setMargin(time)

    double precision :: time

    PUSH_SUB(setMargin)

    margin = time

    POP_SUB(setMargin)

  end subroutine setMargin


  !> start the timer (save starting time)
  
  subroutine start()

    PUSH_SUB(start)

    start_time = loct_clock()
    last_tap = start_time
    if(duration > 1.d0) active = .true.
    
    POP_SUB(start)

  end subroutine start



  !> measure time of on itertion
  
  subroutine walltimer_tap()

    double precision :: now

    PUSH_SUB(walltimer_tap)

    now = loct_clock()
    
    iteration_time = now - last_tap
    last_tap = now
    
    POP_SUB(walltimer_tap)

  end subroutine walltimer_tap


  !> indicate whether time is up

  logical function walltimer_alarm(print)

    logical, optional :: print
    double precision :: now

    PUSH_SUB(walltimer_alarm)

    now = loct_clock()
    
    if(present(print)) then
      if(print) write(*,*) "Walltimer. elapsed time = ", now - start_time, " (", duration, ") ", active 
    end if
    
    if(auto_tap) call walltimer_tap()
    
    if( active .and. (now > start_time + duration - iteration_time - margin) ) then
      walltimer_alarm = .true.
    else
      walltimer_alarm = .false.
    end if
    
    POP_SUB(walltimer_alarm)

  end function walltimer_alarm


  !> get time since start (in seconds)
 
  double precision function walltimer_getElapsedTime()

    PUSH_SUB(walltimer_getElapsedTime)

    walltimer_getElapsedTime = loct_clock() - start_time

    POP_SUB(walltimer_getElapsedTime)

  end function walltimer_getElapsedTime


  !> get iteration time

  double precision function walltimer_getIterationTime()

    PUSH_SUB(walltimer_getIterationTime)

    walltimer_getIterationTime = iteration_time

    POP_SUB(walltimer_getIterationTime)

  end function walltimer_getIterationTime


end module walltimer_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
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

module walltimer_oct_m

  use loct_oct_m
  use parser_oct_m
  
  implicit none
  
  private
  
  double precision :: start_time       !< time the timer was started
  double precision :: last_tap         !< time of the last call to tap()
  double precision :: iteration_time   !< time difference of two calls to tap()
  double precision :: margin           !< additional time margin for writing the restart file
  double precision :: duration         !< time when the alarm should trigger
  
  logical :: active 
  logical :: auto_tap                  !< if .true., tap() is automatically called in every wakeUp() call.

  public :: walltimer_init, walltimer_end, walltimer_tap, walltimer_alarm
  
contains

  !> initialize the timer
  
  subroutine walltimer_init(auto)

    logical, optional, intent(IN) :: auto   !< automatically call walltimer_tap in walltimer_alarm() if .true.

    real  :: alarm_time, write_time

    start_time = 0.d0
    last_tap = 0.d0
    iteration_time = 0.d0
    margin = 0.d0

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
    !%End
    call parse_Variable('Walltime', 0.0, alarm_time)
    call setAlarm(alarm_time*60.d0)
    
    !%Variable RestartWriteTime
    !%Type float
    !%Default 5
    !%Section Execution::IO
    !%Description
    !% The RestartWriteTime (in minutes) will be subtracted from the WallTime to allow time for writing the restart file.
    !% In huge calculations, this value should be increased.
    !%End
    call parse_Variable('RestartWriteTime', 5.0, write_time)
    if(write_time > alarm_time/4) write_time = alarm_time/4
    call setMargin(write_time*60.d0)
    
    call start()
    
    
  end subroutine walltimer_init

  
  !> empty destructor
  
  subroutine walltimer_end()

    active = .false.

  end subroutine walltimer_end


  

  !> set alarm interval in seconds

  subroutine setAlarm(time)

    double precision :: time

    duration = time
        
  end subroutine setAlarm


  !> set safty margin in seconds

  subroutine setMargin(time)

    double precision :: time

    margin = time

  end subroutine setMargin


  !> start the timer (save starting time)
  
  subroutine start()

    
    start_time = loct_clock()
    last_tap = start_time
    if(duration > 1.d0) active = .true.
    
    
  end subroutine start



  !> measure time of on itertion
  
  subroutine walltimer_tap()

    double precision :: now


    now = loct_clock()
    
    iteration_time = now - last_tap
    last_tap = now
    
  end subroutine walltimer_tap


  !> indicate whether time is up

  logical function walltimer_alarm(print)

    logical, optional :: print
    double precision :: now

    now = loct_clock()
    
    if(present(print)) then
       if(print) write(*,*) "Timer::WakeUp called. Time = ", now - start_time, " (", duration, ") ", active 
    endif
    
    if(auto_tap) call walltimer_tap()
    
    if( active .and. (now > start_time + duration - iteration_time - margin) ) then
       walltimer_alarm = .true.
    else
       walltimer_alarm = .false.
    endif
    
  end function walltimer_alarm


  !> get time since start (in seconds)
 
  double precision function walltimer_getElapsedTime()

    walltimer_getElapsedTime = loct_clock() - start_time

  end function walltimer_getElapsedTime


  !> get iteration time

  double precision function walltimer_getIterationTime()

    walltimer_getIterationTime = iteration_time

  end function walltimer_getIterationTime


end module walltimer_oct_m

!! Copyright (C) 2021 M. Lueders
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

module multisystem_debug_oct_m

  use algorithm_oct_m
  use clock_oct_m
  use global_oct_m
  use interaction_partner_oct_m
  use io_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use profiling_oct_m
  
  implicit none

  private

  public &
    multisystem_debug_init,            &
    multisystem_debug_end,             &
    multisystem_debug_start_log,       &
    multisystem_debug_stop_log,        &
    multisystem_debug_write_event_in,  &
    multisystem_debug_write_event_out, &
    event_info_t,                      &
    event_function_call_t,             &
    event_handle_t

  integer, parameter, public :: MAX_INFO_LEN = 128


  !-------------------------------------------------------------------

  type, abstract :: event_info_t
  private
  contains
    procedure(event_info_get_info),  deferred :: get_info  
  end type event_info_t

  abstract interface

    function event_info_get_info(this) result(res)
      import event_info_t
      import MAX_INFO_LEN
      class(event_info_t), intent(in) :: this
      character(len=MAX_INFO_LEN)     :: res
    end function  event_info_get_info

  end interface

  !-------------------------------------------------------------------

  type, extends(event_info_t) :: event_function_call_t
    character(len=MAX_INFO_LEN) :: function_name
    character(len=ALGO_LABEL_LEN) :: op_label
  contains
    procedure :: get_info => event_function_call_get_info
  end type event_function_call_t

  interface event_function_call_t
    procedure :: event_function_call_constructor
  end interface event_function_call_t

  !-------------------------------------------------------------------

  type :: event_handle_t
    character, allocatable :: name(:)
    integer, public :: enter_ID
  end type event_handle_t

  interface event_handle_t
    procedure :: event_handle_constructor
  end interface event_handle_t

  !-------------------------------------------------------------------

  type(mpi_grp_t) :: mpi_grp
  integer iunit
  integer event_ID
  logical log_active

contains


  function event_handle_constructor( id, name ) result(handle)
    integer, intent(in)                   :: id
    character(len=MAX_INFO_LEN), optional :: name
    type(event_handle_t)                  :: handle

    PUSH_SUB(event_handle_constructor)
    handle%enter_ID = id
    POP_SUB(event_handle_constructor)

  end function event_handle_constructor
  !-------------------------------------------------------------------

  function event_function_call_constructor(name, op) result(event)
    character(*) name
    type(algorithmic_operation_t), optional :: op
    type(event_function_call_t), pointer :: event

    PUSH_SUB(event_function_call_constructor)

    SAFE_ALLOCATE(event)
    event%function_name = name

    if(present(op)) then
      event%op_label = op%label
    else
      event%op_label = "NULL"
    endif

    POP_SUB(event_function_call_constructor)

  end function event_function_call_constructor

  
  function event_function_call_get_info(this) result(info)
    class(event_function_call_t), intent(in) :: this
    character(len=MAX_INFO_LEN)  :: info

    info = "function: " // trim(this%function_name) 
    if(this%op_label /= "NULL") then
      info =  "function: " // trim(this%function_name) // "|operation: " // trim(this%op_label)
    endif

  end function event_function_call_get_info

  !-------------------------------------------------------------------

  subroutine multisystem_debug_init(filename, namespace, group)
    character(*),      intent(in)      :: filename
    type(namespace_t), intent(in)      :: namespace
    type(mpi_grp_t),   intent(in)      :: group

    PUSH_SUB(multisystem_debug_init)

    mpi_grp = group

    event_ID = 0
    if(mpi_grp%rank == 0) then
      iunit = io_open(filename, namespace, action="write", status="unknown" )
    end if
    log_active = .false.

    POP_SUB(multisystem_debug_init) 

  end subroutine multisystem_debug_init

  subroutine multisystem_debug_end()

    PUSH_SUB(multisystem_debug_end)
    if(mpi_grp%rank == 0) then
      call io_close(iunit)
    end if
    POP_SUB(multisystem_debug_end)

  end subroutine multisystem_debug_end

  subroutine multisystem_debug_start_log()

    PUSH_SUB(multisystem_debug_start_log)
    log_active = .true.
    POP_SUB(multisystem_debug_start_log)
    
  end subroutine multisystem_debug_start_log

  subroutine multisystem_debug_stop_log()

    PUSH_SUB(multisystem_debug_stop_log)
    log_active = .false.
    POP_SUB(multisystem_debug_stop_log)

  end subroutine multisystem_debug_stop_log


  function multisystem_debug_write_event_in(system, event, extra,  system_clock, prop_clock, & 
                                            interaction_clock, partner_clock, requested_clock) result(handle)
    class(interaction_partner_t), intent(in), optional  :: system
    class(event_info_t), intent(in)                     :: event
    character(*), optional                              :: extra
    type(clock_t), intent(in), optional                 :: system_clock
    type(clock_t), intent(in), optional                 :: prop_clock
    type(clock_t), intent(in), optional                 :: interaction_clock
    type(clock_t), intent(in), optional                 :: partner_clock
    type(clock_t), intent(in), optional                 :: requested_clock
    type(event_handle_t)         :: handle

    character(len = MAX_NAMESPACE_LEN) ::  system_name

    PUSH_SUB(multisystem_debug_write_event_in)

    if( (mpi_grp%rank == 0) .and. log_active) then

      if(present(system)) then
        system_name = '.'//trim(system%namespace%get())
        if (system_name == '.') system_name = ''
      else
        system_name = 'KEEP'
      endif

      handle = event_handle_t(event_ID)

      write(iunit, '("IN  step: ",I10," | system: ",a,"| ",a)' , advance='no' ) event_ID, trim(system_name), trim(event%get_info())

      if( present(extra)) then
        write(iunit, '(" | ",a)' , advance='no')  trim(extra)
      end if

      if (present(system_clock)) then
        write(iunit, '(" | system_clock:", E15.5)' , advance='no')  system_clock%time()
      end if

      if (present(prop_clock)) then
        write(iunit, '(" | prop_clock:", E15.5)' , advance='no')  prop_clock%time()
      end if

      if (present(interaction_clock)) then
        write(iunit, '(" | interaction_clock:", E15.5)' , advance='no')  interaction_clock%time()
      end if

      if (present(partner_clock)) then
        write(iunit, '(" | partner_clock:", E15.5)' , advance='no')  partner_clock%time()
      end if

      if (present(requested_clock)) then
        write(iunit, '(" | requested_clock:", E15.5)' , advance='no')  requested_clock%time()
      end if

      write(iunit, '()' , advance='yes')

      event_ID = event_ID + 1

    endif

    POP_SUB(multisystem_debug_write_event_in)

  end function multisystem_debug_write_event_in

  subroutine multisystem_debug_write_event_out(handle, extra, update, system_clock, prop_clock, &
                                               interaction_clock, partner_clock, requested_clock) 
    class(event_handle_t), intent(in)    :: handle
    character(*), optional               :: extra
    logical, optional                    :: update
    type(clock_t), intent(in), optional  :: system_clock
    type(clock_t), intent(in), optional  :: prop_clock
    type(clock_t), intent(in), optional  :: interaction_clock
    type(clock_t), intent(in), optional  :: partner_clock
    type(clock_t), intent(in), optional  :: requested_clock

    character(17)                        :: update_string

    PUSH_SUB(multisystem_debug_write_event_out)

    if(log_active) then

      if(present(update)) then
        if(update) then
          update_string = " | updated: true"
        else
          update_string = " | updated: false"
        endif
      else
        update_string = ""
      endif

      write(iunit, '("OUT step: ",I10," | closes: ",I10)', advance='no')  &
        event_ID, handle%enter_ID

      if(present(update)) then
        if(update) then
            write(iunit, '(" | updated: true")', advance='no')
          else
            write(iunit, '(" | updated: false")', advance='no')
          endif
      end if
  
      if( present(extra)) then
        write(iunit, '(" | ",a)' , advance='no')  trim(extra)
      end if

      if (present(system_clock)) then
        write(iunit, '(" | system_clock:", E15.5)' , advance='no')  system_clock%time()
      end if

      if (present(prop_clock)) then
        write(iunit, '(" | prop_clock:", E15.5)' , advance='no')  prop_clock%time()
      end if

      if (present(interaction_clock)) then
        write(iunit, '(" | interaction_clock:", E15.5)' , advance='no')  interaction_clock%time()
      end if

      if (present(partner_clock)) then
        write(iunit, '(" | partner_clock:", E15.5)' , advance='no')  partner_clock%time()
      end if

      if (present(requested_clock)) then
        write(iunit, '(" | requested_clock:", E15.5)' , advance='no')  requested_clock%time()
      end if

      write(iunit, '()' , advance='yes')

      event_ID = event_ID + 1

    endif 


    POP_SUB(multisystem_debug_write_event_out)

  end subroutine multisystem_debug_write_event_out


end module multisystem_debug_oct_m
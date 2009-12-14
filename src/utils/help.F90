!! Copyright (C) 2007 Xavier Andrade
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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: help.F90 $

#include "global.h"

program oct_help
  
  use command_line_m
  use global_m
  use messages_m
  use varinfo_m

  implicit none

  integer :: argc, ierr
  character(len=32) :: mode
  character(len=100) :: varname

  integer, parameter :: help_stdout = 6, help_stderr = 0

  call global_init()
  call command_line_version()

  if( .not. command_line_is_available() ) then 
    write(help_stderr, '(a)') "Your fortran compiler doesn't support command line arguments,"
    write(help_stderr, '(a)') "the oct-help command is not available."
    stop
  end if

  argc = command_argument_count()

  if( argc == 0 ) then
    write(help_stderr, '(a)') 'Usage: oct-help { show | search | list }'
    stop
  end if

  call get_command_argument(1, mode)
  
  select case(mode)

  case("show")
    
    if (argc >= 2) then
      call get_command_argument(2, varname)
      call varinfo_print(help_stdout, trim(varname), ierr)

      if (ierr /= 0) then
        write(help_stderr, '(a)') "Error: Variable "//trim(varname)//" not found."
      end if

    else
      write(help_stderr, '(a)') "Gives information about a variable"
      write(help_stderr, '(a)') "Usage: oct-help show variable_name"
    end if

  case("search")

    if (argc >= 2) then
      call get_command_argument(2, varname)
      call varinfo_search(help_stdout, trim(varname), ierr)
    else
      write(help_stderr, '(a)') 'Searches for variable names that contain a certain string'
      write(help_stderr, '(a)') 'Usage: oct-help search string'
    end if

  case("list")
    call varinfo_search(help_stdout, "", ierr)

  case default
    write(help_stderr, '(a,a)') 'Unknown command: ', mode
    write(help_stderr, '(a)') 'Usage: oct-help { show | search | list}'

  end select
    
  call global_end()
    
end program oct_help

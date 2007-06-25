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

  call global_init()

  if( .not. command_line_is_available() ) then 
    message(1) = "Your fortran compiler doesn't support command line arguments,"
    message(2) = "the oct-help command is not available."
    call write_fatal(2)
  end if

  argc = command_argument_count()

  if( argc == 0 ) then
    message(1) = "Usage: oct-test { show | search } variable_name" 
    call write_info(1)
  end if

  call get_command_argument(1, mode)
  
  select case(mode)

  case("show")
    
    if (argc >= 2) then
      call get_command_argument(2, varname)
      call varinfo_print(stdout, trim(varname), ierr)

      if (ierr /= 0) then
        message(1) = "Error: Variable "//trim(varname)//" not found."
        call write_info(1)
      end if

    else
      message(1) = "Gives information about a variable"
      message(2) = "Usage: oct-help show variable_name"
      call write_info(2)
    end if

  case("search")

    if (argc >= 2) then
      call get_command_argument(2, varname)
      call varinfo_search(stdout, trim(varname), ierr)
    else
      message(1) = "Searches for variable names that contains a certain string"
      message(2) = "Usage: oct-help search string"
      call write_info(2)
    end if

  end select
    
  call global_end()
    
end program oct_help

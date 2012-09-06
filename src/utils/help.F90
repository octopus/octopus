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

  integer :: ierr
  character(len=32) :: mode
  character(len=100) :: varname

  integer, parameter :: help_stdout = 6, help_stderr = 0

  call global_init()

  call getopt_init(ierr)
  if(ierr.ne.0) then
    write(stderr, '(a)') "Your Fortran compiler doesn't support command-line arguments;"
    write(stderr, '(a)') "the oct-help command is not available."
    stop
  end if

  mode    = " "
  varname = " "
  call getopt_help(mode, varname)
  call getopt_end()

  select case(mode)
  case("print")
    call varinfo_print(help_stdout, trim(varname), ierr)
    if (ierr /= 0) then
      write(help_stderr, '(a)') "Error: Variable "//trim(varname)//" not found."
    end if

  case("search")
    call varinfo_search(help_stdout, trim(varname), ierr)

  case("list")
    call varinfo_search(help_stdout, "", ierr)
  end select
    
  call global_end()
    
end program oct_help

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: cross_section.F90 2781 2007-03-23 10:58:32Z lorenzen $

#include "global.h"

module getopt_m
#ifdef FC_COMMAND_LINE_MODULE
  use FC_COMMAND_LINE_MODULE
#endif

  implicit none

#ifdef FC_COMMAND_LINE_INCLUDE
include FC_COMMAND_LINE_INCLUDE
#endif

  private
  public :: getopt_init,                &
            getopt_oscillator_strength


  interface 
    subroutine set_number_clarg(argc)
      integer :: argc
    end subroutine set_number_clarg
    subroutine set_clarg(i, argstring)
      integer :: i
      character(len=*) :: argstring
    end subroutine set_clarg
    subroutine getopt_oscillator_strength(omega, filename, searchinterval)
      real(8) :: omega
      character(len=*) :: filename
      real(8) :: searchinterval
    end subroutine getopt_oscillator_strength
  end interface


  contains


  subroutine getopt_init(ierr)
    integer, intent(out) :: ierr
    integer :: argc, i
    character(len=100), allocatable :: argstring(:)
#ifdef FC_COMMAND_LINE_ARGUMENTS
    argc = command_argument_count()
    allocate(argstring(0:argc))
    call set_number_clarg(argc)
    do i = 0, argc
      call get_command_argument(i, argstring(i))
      call set_clarg(i, argstring(i))
    end do  
    deallocate(argstring)
    ierr = 0
#else
    ierr = -1
#endif
  end subroutine getopt_init



end module getopt_m



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

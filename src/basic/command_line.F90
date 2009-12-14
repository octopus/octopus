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

module command_line_m

#ifdef FC_COMMAND_LINE_MODULE
  use FC_COMMAND_LINE_MODULE
#endif

  use global_m
  use messages_m
  use mpi_m
  use profiling_m

  implicit none

#ifdef FC_COMMAND_LINE_INCLUDE
  include FC_COMMAND_LINE_INCLUDE
#endif

  private
  public :: getopt_init,                &
            getopt_oscillator_strength, &
            command_line_is_available,  &
            command_line_version

#if FC_COMMAND_LINE_ARGUMENTS != 2003
  public :: command_argument_count, get_command_argument
#endif

! If Fortran 2003 interface to command line arguments is not
! available, define it using an interface over Fortran 77 API
#if FC_COMMAND_LINE_ARGUMENTS == 77 && ! defined(FC_COMMAND_LINE_INTRINSIC)

  interface command_argument_count
#ifdef FC_COMMAND_LINE_IMPLICIT
     integer function iargc() 
     end function iargc
#else
     module procedure iargc
#endif
  end interface

  interface get_command_argument
#ifdef FC_COMMAND_LINE_IMPLICIT
     subroutine getarg(c, a)
       integer :: c
       character(len=*) :: a
     end subroutine getarg
#else
     module procedure getarg
#endif
  end interface
  
#endif /* FC_COMMAND_LINE_ARGUMENTS == 77 */

  interface 
    subroutine set_number_clarg(argc)
      integer :: argc
    end subroutine set_number_clarg
    subroutine set_clarg(i, argstring)
      integer :: i
      character(len=*) :: argstring
    end subroutine set_clarg
    subroutine getopt_oscillator_strength(mode, omega, searchinterval, &
                                          order, nresonances, nfrequencies, time, &
                                          l, m, damping, file)
      integer :: mode
      real(8) :: omega
      real(8) :: searchinterval
      integer :: order, nresonances, nfrequencies
      real(8) :: time
      integer :: print_omega_file
      integer :: l, m
      real(8) :: damping
      character(len=*) :: file
    end subroutine getopt_oscillator_strength
    subroutine getopt_harmonic_spectrum(w0, m, pol)
      real(8)          :: w0
      integer          :: m
      character(len=*) :: pol
    end subroutine getopt_harmonic_spectrum
  end interface

  contains

    pure logical function command_line_is_available()
#ifdef FC_COMMAND_LINE_ARGUMENTS
      command_line_is_available = .true.
#else
      command_line_is_available = .false.
#endif
    end function command_line_is_available

  subroutine getopt_init(ierr)
    integer, intent(out) :: ierr
    integer :: argc, i
    character(len=100), allocatable :: argstring(:)
#ifdef FC_COMMAND_LINE_ARGUMENTS
    argc = command_argument_count()
    SAFE_ALLOCATE(argstring(0:argc))
    call set_number_clarg(argc)
    do i = 0, argc
      call get_command_argument(i, argstring(i))
      call set_clarg(i, argstring(i))
    end do  
    SAFE_DEALLOCATE_A(argstring)
    ierr = 0
#else
    ierr = -1
#endif
  end subroutine getopt_init


!if there is no way to access command line, define some dummy
!functions to avoid problems when linking

#if defined(FC_COMMAND_LINE_INTRINSIC) || ! defined(FC_COMMAND_LINE_ARGUMENTS)

  integer function command_argument_count()
#if defined(FC_COMMAND_LINE_INTRINSIC)
    command_argument_count = iargc()
#else
    command_argument_count = 0
#endif
  end function command_argument_count

  subroutine get_command_argument(c, a)
    integer,          intent(in)     :: c
    character(len=*), intent(out)    :: a
#if defined(FC_COMMAND_LINE_INTRINSIC)
    call getarg(c, a)
#endif
  end subroutine get_command_argument

#endif

  subroutine command_line_version()
    integer :: ii
    character(len=100) :: arg
    
    do ii = 1, command_argument_count()
      call get_command_argument(ii, arg)
      if(arg == "--version") then
        write(message(1), '(5a)') 'octopus ', trim(conf%version), ' (svn version ', trim(conf%latest_svn), ')'
        call write_info(1)

        if(flush_messages .and. mpi_grp_is_root(mpi_world)) then
          close(iunit_err)
        end if

#ifdef HAVE_MPI
        call MPI_Finalize(mpi_err)
#endif
        stop
      end if
    end do

  end subroutine command_line_version

end module command_line_m



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

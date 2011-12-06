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

! ---------------------------------------------------------
!
! This module:
!
! (i) provides an interface to the command line, as defined by the Fortran 2003 routines
! "get_command_argument", "command_argument_count", and "get_command". These routines
! may not exist in all compilers (octopus presently does not require Fortran 2003), and
! in that case they are defined here. But it may happen that the compiler may not be able
! to access the command line arguments at all, and in that case the utilities that need them
! are not usable.
!
! (ii) provides an interface to the getopt C library, so that the command line arguments
! can be managed in this "standard" way.
!
! ---------------------------------------------------------
module command_line_m

  ! ---------------------------------------------------------
  ! The compilation depends on several macros defined at configure time (the m4 code that
  ! takes care of testing the compiler in order to build this macros is in m4/fc_command_line_m4):
  !
  ! FC_COMMAND_LINE_ARGUMENTS maybe be "2003", "77", or undefined. In the first case, 
  !   the "get_command_argument", "command_argument_count", and "get_command" routines
  !   are defined. In the second case, they are coded here, by accessing the non-standard
  !   FORTRAN 77 "getarg" and "iargc", present in many compiler extensions. If this macro
  !   is undefined, then the command line arguments cannot be accessed at all.
  ! FC_COMMAND_LINE_MODULE is the module that should be used in order to have access to
  !   the getarg and iargc procedures, in case "FC_COMMAND_LINE_ARGUMENTS = 77".  This 
  !   module depends on the compiler, and in some cases it is none at all.
  ! FC_COMMAND_LINE_INCLUDE is the file to be included in order to have access to the getarg
  !   and iargc procedures, in case "FC_COMMAND_LINE_ARGUMENTS = 77". This file depends on
  !   the compiler, and in some case it is none at all.
  ! FC_COMMAND_LINE_INTRINSIC is defined if "FC_COMMAND_LINE_ARGUMENTS = 77", and the 
  !   procedures iargc and getarg are intrinsi, i.e. there is no need to use any module,
  !   include any file, or declare them.
  ! FC_COMMAND_LINE_IMPLICIT is defined if "FC_COMMAND_LINE_ARGUMENTS = 77", and the procedures
  !   iargc and getarg are implicit, i.e. there is no need to use any module or include any
  !   file, but have to be declared.
  ! ---------------------------------------------------------

#ifdef FC_COMMAND_LINE_MODULE
  use FC_COMMAND_LINE_MODULE
#endif

  implicit none

#ifdef FC_COMMAND_LINE_INCLUDE
  include FC_COMMAND_LINE_INCLUDE
#endif

  private
  public :: getopt_init,                 &
            getopt_end,                  &
            getopt_octopus,              &
            getopt_casida_spectrum,      &
            getopt_center_geom,          &
            getopt_propagation_spectrum, &
            getopt_rotatory_strength,    &
            getopt_vibrational,          &
            getopt_xyz_anim,             &
            getopt_oscillator_strength,  &
            getopt_harmonic_spectrum,    &
            getopt_help
#if FC_COMMAND_LINE_ARGUMENTS != 2003
  public :: command_argument_count,     &
            get_command_argument
#endif


  ! ---------------------------------------------------------
  ! First, the public interfaces.


  ! Each program/utility that needs to use the getopt features should have
  ! an interface here -- the definition of the procedure should be given in the
  ! getopt_f.c file.
  interface
    subroutine getopt_octopus
    end subroutine getopt_octopus

    subroutine getopt_casida_spectrum
    end subroutine getopt_casida_spectrum

    subroutine getopt_center_geom
    end subroutine getopt_center_geom

    subroutine getopt_dielectric_function
    end subroutine getopt_dielectric_function

    subroutine getopt_propagation_spectrum(fname)
      character(len=*) :: fname
    end subroutine getopt_propagation_spectrum

    subroutine getopt_rotatory_strength
    end subroutine getopt_rotatory_strength

    subroutine getopt_vibrational(mode)
      integer :: mode
    end subroutine getopt_vibrational

    subroutine getopt_xyz_anim
    end subroutine getopt_xyz_anim

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

    subroutine getopt_harmonic_spectrum(w0, m, ar, x, y, z, pol)
      real(8)          :: w0
      integer          :: m
      integer          :: ar
      character(len=*) :: pol
      real(8)          :: x
      real(8)          :: y
      real(8)          :: z
    end subroutine getopt_harmonic_spectrum

    subroutine getopt_help(mode, name)
      character(len=*) :: mode
      character(len=*) :: name
    end subroutine getopt_help
    
    subroutine getopt_photoelectron_spectrum(mode,interp)
      integer          :: mode
      integer          :: interp
    end subroutine getopt_photoelectron_spectrum

  end interface

  ! If Fortran 2003 interface to command line arguments is not
  ! available, define it using an interface over Fortran 77 API.
  !
  ! This cannot be done in case the compiler defines the Fortran 77 API through
  ! intrinsic procedures. That case is taken care below.
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


  ! ---------------------------------------------------------
  ! The following interfaces are private to this module, and should
  ! not be called from outside.

  interface 
    subroutine set_number_clarg(argc)
      integer :: argc
    end subroutine set_number_clarg
    subroutine set_clarg(i, argstring)
      integer :: i
      character(len=*) :: argstring
    end subroutine set_clarg
    subroutine clean_clarg()
    end subroutine clean_clarg
  end interface


  contains


  ! ---------------------------------------------------------
  !> Initializes the getopt machinery. Must be called before attempting
  !! to parse the options. On output, ierr is zero if the call is successful,
  !! and -1 if the command line arguments cannot be accessed.
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


  subroutine getopt_end
#ifdef FC_COMMAND_LINE_ARGUMENTS
    call clean_clarg()
#endif
  end subroutine getopt_end

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


end module command_line_m



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

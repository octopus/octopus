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

#include "global.h"

module loct_oct_m

  implicit none

  !> Define which routines can be seen from the outside
  private
  public ::                  &
    loct_clock,              &
    loct_gettimeofday,       &
    loct_nanosleep,          &
    loct_getcwd,             &
    loct_dirname,            &
    loct_basename,           &
    loct_realpath,           &
    loct_sysname,            &
    loct_search_file_lr,     &
    loct_mkdir,              &
    loct_stat,               &
    loct_rm,                 &
    loct_dir_exists,         &
    loct_number_of_lines,    &
    loct_break_C_string,     &
    loct_getenv,             &
    loct_isinstringlist,     &
    loct_progress_bar,       &
    loct_printRecipe,        &
    loct_strerror,           &
    loct_get_memory_usage,   &
    loct_exit_failure

#if defined(HAVE_GDLIB)
  public ::                    &
    loct_gdimage_create_from,  &
    loct_gdimage_sx,           &
    loct_gdimage_sy,           &
    loct_gdimage_get_pixel_rgb, &
    loct_gdimagedestroy
#endif

  ! ---------------------------------------------------------
  !> System information (time, memory, sysname)

  interface loct_strerror
    subroutine oct_strerror(errno, res)
      implicit none
      integer, intent(in) :: errno
      character(len=*), intent(out)  :: res
    end subroutine oct_strerror
  end interface loct_strerror

  interface loct_clock
    function oct_clock()
     implicit none
      real(8) :: oct_clock
    end function oct_clock
  end interface loct_clock

  interface loct_gettimeofday
    subroutine oct_gettimeofday(sec, usec)
      implicit none
      integer, intent(out) :: sec, usec
    end subroutine oct_gettimeofday
  end interface loct_gettimeofday

  interface loct_nanosleep
    subroutine oct_nanosleep(sec, nsec)
      implicit none
      integer, intent(in) :: sec  !< number of seconds
      integer, intent(in) :: nsec !< + number of nanoseconds
    end subroutine oct_nanosleep
  end interface loct_nanosleep

  interface loct_sysname
    subroutine oct_sysname(name)
      implicit none
      character(len=*), intent(out) :: name
    end subroutine oct_sysname
  end interface loct_sysname

  interface loct_getcwd
    subroutine oct_getcwd(name)
      implicit none
      character(len=*), intent(out) :: name
    end subroutine oct_getcwd
  end interface loct_getcwd

  interface loct_realpath
    subroutine oct_realpath(fnam, rnam)
      character(len=*), intent(in)  :: fnam
      character(len=*), intent(out) :: rnam
    end subroutine oct_realpath
  end interface

  interface loct_dirname
    subroutine oct_dirname(fnam, dnam)
      character(len=*), intent(in)  :: fnam
      character(len=*), intent(out) :: dnam
    end subroutine oct_dirname
  end interface

  interface loct_basename
     subroutine oct_basename(fnam, dnam)
       character(len=*), intent(in)  :: fnam
       character(len=*), intent(out) :: dnam
     end subroutine oct_basename
  end interface


  ! ---------------------------------------------------------
  !> File-handling
  interface loct_mkdir
    subroutine oct_mkdir(name)
      implicit none
      character(len=*), intent(in) :: name
    end subroutine oct_mkdir
  end interface loct_mkdir

  interface loct_stat
    subroutine oct_stat(ierr, name, mod_time)
      implicit none
      integer,          intent(out) :: ierr
      character(len=*), intent(in)  :: name
      character(len=*), intent(out) :: mod_time
    end subroutine oct_stat
  end interface loct_stat

  interface loct_rm
    subroutine oct_rm(name)
      implicit none
      character(len=*), intent(in) :: name
    end subroutine oct_rm
  end interface loct_rm

  interface loct_number_of_lines
    integer function oct_number_of_lines(filename)
      implicit none
      character(len=*), intent(in) :: filename
    end function oct_number_of_lines
  end interface loct_number_of_lines

  interface loct_break_C_string
    subroutine oct_break_C_string(str, s, line)
      use iso_c_binding
      implicit none
      type(c_ptr),       intent(in)    :: str
      type(c_ptr),       intent(inout) :: s
      character(len=*),  intent(out)   :: line
    end subroutine oct_break_C_string
  end interface loct_break_C_string

  interface loct_search_file_lr
    subroutine oct_search_file_lr(freq, tag, ierr, dirname)
      implicit none
      REAL_DOUBLE,      intent(inout) :: freq
      integer,          intent(in)    :: tag
      integer,          intent(out)   :: ierr
      character(len=*), intent(in)    :: dirname
    end subroutine oct_search_file_lr
  end interface loct_search_file_lr

  ! ---------------------------------------------------------
  !> Varia
  interface loct_getenv
    subroutine oct_getenv(var, val)
      implicit none
      character(len=*), intent(in)  :: var
      character(len=*), intent(out) :: val
    end subroutine oct_getenv
  end interface loct_getenv

  interface loct_progress_bar
    subroutine oct_progress_bar(a, maxcount)
      implicit none
      integer, intent(in) :: a, maxcount
    end subroutine oct_progress_bar
  end interface loct_progress_bar

  interface loct_printRecipe
    subroutine oct_printRecipe(dir, filename)
      implicit none
      character(len=*), intent(in)  :: dir
      character(len=*), intent(out) :: filename
    end subroutine oct_printRecipe
  end interface loct_printRecipe

  interface loct_exit_failure
    subroutine oct_exit_failure()
      implicit none
    end subroutine oct_exit_failure
  end interface loct_exit_failure

  interface loct_wfs_list
    subroutine oct_wfs_list(str, l)
      implicit none
      character(len=*), intent(in)  :: str
      integer,          intent(out) :: l !< array
    end subroutine oct_wfs_list
  end interface loct_wfs_list

  ! ---------------------------------------------------------
  !> GD library
#if defined(HAVE_GDLIB)
  interface loct_gdimage_create_from
    function oct_gdimage_create_from(filename)
      use iso_c_binding
      implicit none
      type(c_ptr) :: oct_gdimage_create_from
      character(len=*), intent(in) :: filename
    end function oct_gdimage_create_from
  end interface loct_gdimage_create_from

  interface loct_gdimage_sx
    function oct_gdimage_sx(im)
      use iso_c_binding
      implicit none
      integer :: oct_gdimage_sx
      type(c_ptr), intent(in) :: im
    end function oct_gdimage_sx
  end interface loct_gdimage_sx

  interface loct_gdimage_sy
    function oct_gdimage_sy(im)
      use iso_c_binding
      implicit none
      integer :: oct_gdimage_sy
      type(c_ptr), intent(in) :: im
    end function oct_gdimage_sy
  end interface loct_gdimage_sy

  interface loct_gdimage_get_pixel_rgb
    subroutine oct_gdimage_get_pixel_rgb(im, x, y, r, g, b)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in)  :: im
      integer,     intent(in)  :: x, y
      integer,     intent(out) :: r, g, b
    end subroutine oct_gdimage_get_pixel_rgb
  end interface loct_gdimage_get_pixel_rgb

  interface loct_gdimagedestroy
    subroutine oct_gdimagedestroy(im)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: im
    end subroutine oct_gdimagedestroy
  end interface loct_gdimagedestroy
#endif

 interface loct_get_memory_usage
   integer(SIZEOF_VOIDP) function oct_get_memory_usage()
     implicit none
   end function oct_get_memory_usage
 end interface loct_get_memory_usage

contains

  logical function loct_isinstringlist(a, s) result(inlist)
    integer,          intent(in) :: a
    character(len=*), intent(in) :: s

    integer, allocatable :: list(:)

    allocate(list(2**14))

    call loct_wfs_list(s, list(1))
    inlist = .false.
    if (list(a) == 1) inlist = .true.

    deallocate(list)

  end function loct_isinstringlist


  logical function loct_dir_exists(dirname) result(exists)
    character(len=*), intent(in) :: dirname

    interface oct_dir_exists
      integer function oct_dir_exists(dirname)
        implicit none
        character(len=*), intent(in)    :: dirname
      end function oct_dir_exists
    end interface oct_dir_exists

    exists = oct_dir_exists(dirname) /= 0

  end function loct_dir_exists

end module loct_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

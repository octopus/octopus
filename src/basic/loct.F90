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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module loct_m

  implicit none

  ! Define which routines can be seen from the outside
  private
  public ::                  &
    loct_clock,              &
    loct_gettimeofday,       &
    loct_nanosleep,          &
    loct_getcwd,             &
    loct_sysname,            &
    loct_mkdir,              &
    loct_stat,               &
    loct_rm,                 &
    loct_rm_status_files,    &
    loct_number_of_lines,    &
    loct_getenv,             &
    loct_isinstringlist,     &
    loct_progress_bar,       &
    loct_printRecipe,        &
    loct_strerror,           &
    get_memory_usage,        &
    loct_exit_failure,       &
    loct_pointer_copy

#if defined(HAVE_GDLIB)
  public ::                    &
    loct_gdimage_create_from,  &
    loct_gdimage_sx,           &
    loct_gdimage_sy,           &
    loct_gdimage_get_pixel_rgb
#endif

  interface loct_pointer_copy
    module procedure sloct_pointer_copy_1
    module procedure sloct_pointer_copy_2
    module procedure sloct_pointer_copy_3
    module procedure sloct_pointer_copy_4
    module procedure dloct_pointer_copy_1
    module procedure dloct_pointer_copy_2
    module procedure dloct_pointer_copy_3
    module procedure dloct_pointer_copy_4
    module procedure cloct_pointer_copy_1
    module procedure cloct_pointer_copy_2
    module procedure cloct_pointer_copy_3
    module procedure cloct_pointer_copy_4
    module procedure zloct_pointer_copy_1
    module procedure zloct_pointer_copy_2
    module procedure zloct_pointer_copy_3
    module procedure zloct_pointer_copy_4
    module procedure iloct_pointer_copy_1
    module procedure iloct_pointer_copy_2
    module procedure iloct_pointer_copy_3
    module procedure iloct_pointer_copy_4
    module procedure aloct_pointer_copy_1
    module procedure aloct_pointer_copy_2
    module procedure aloct_pointer_copy_3
    module procedure aloct_pointer_copy_4
    module procedure lloct_pointer_copy_1
    module procedure lloct_pointer_copy_2
    module procedure lloct_pointer_copy_3
    module procedure lloct_pointer_copy_4
  end interface loct_pointer_copy


  ! ---------------------------------------------------------
  !> System information (time, memory, sysname)

  interface loct_strerror
    subroutine oct_strerror(errno, res)
      integer, intent(in) :: errno
      character(len=*), intent(out)  :: res
    end subroutine oct_strerror
  end interface

  interface loct_clock
    function oct_clock()
      real(8) :: oct_clock
    end function oct_clock
  end interface

  interface loct_gettimeofday
    subroutine oct_gettimeofday(sec, usec)
      integer, intent(out) :: sec, usec
    end subroutine oct_gettimeofday
  end interface

  interface loct_nanosleep
    subroutine oct_nanosleep(sec, usec)
      integer, intent(in) :: sec, usec
    end subroutine oct_nanosleep
  end interface

  interface loct_sysname
    subroutine oct_sysname(name)
      character(len=*), intent(out) :: name
    end subroutine oct_sysname
  end interface

  interface loct_getcwd
    subroutine oct_getcwd(name)
      character(len=*), intent(out) :: name
    end subroutine oct_getcwd
  end interface


  ! ---------------------------------------------------------
  !> File-handling
  interface loct_mkdir
    subroutine oct_mkdir(name)
      character(len=*), intent(in) :: name
    end subroutine oct_mkdir
  end interface

  interface loct_stat
    subroutine oct_stat(ierr, name)
      integer,          intent(out) :: ierr
      character(len=*), intent(in)  :: name
    end subroutine oct_stat
  end interface

  interface loct_rm
    subroutine oct_rm(name)
      character(len=*), intent(in) :: name
    end subroutine oct_rm
  end interface

  interface loct_number_of_lines
    integer function number_of_lines(filename)
      character(len=*), intent(in) :: filename
    end function number_of_lines
  end interface


  ! ---------------------------------------------------------
  !> Varia
  interface loct_getenv
    subroutine oct_getenv(var, value)
      character(len=*), intent(in)  :: var
      character(len=*), intent(out) :: value
    end subroutine oct_getenv
  end interface

  interface loct_progress_bar
    subroutine oct_progress_bar(a, max)
      integer, intent(in) :: a, max
    end subroutine oct_progress_bar
  end interface

  interface loct_printRecipe
    subroutine oct_printRecipe(dir, filename)
      character(len=*), intent(in)  :: dir
      character(len=*), intent(out) :: filename
    end subroutine oct_printRecipe
  end interface

  interface 
    subroutine loct_exit_failure()
    end subroutine loct_exit_failure
  end interface

  ! ---------------------------------------------------------
  !> GD library
#if defined(HAVE_GDLIB)
  interface loct_gdimage_create_from
    function oct_gdimage_create_from(filename)
      use c_pointer_m
      type(c_ptr) :: oct_gdimage_create_from
      character(len=*), intent(in) :: filename
    end function oct_gdimage_create_from
  end interface

  interface loct_gdimage_sx
    function oct_gdimage_sx(im)
      use c_pointer_m
      integer :: oct_gdimage_sx
      type(c_ptr), intent(in) :: im
    end function oct_gdimage_sx
  end interface

  interface loct_gdimage_sy
    function oct_gdimage_sy(im)
      use c_pointer_m
      integer :: oct_gdimage_sy
      type(c_ptr), intent(in) :: im
    end function oct_gdimage_sy
  end interface

  interface loct_gdimage_get_pixel_rgb
    subroutine oct_gdimage_get_pixel_rgb(im, x, y, r, g, b)
      use c_pointer_m
      type(c_ptr), intent(in)  :: im
      integer,   intent(in)  :: x, y
      integer,   intent(out) :: r, g, b
    end subroutine oct_gdimage_get_pixel_rgb
  end interface
#endif

 interface
   integer(SIZEOF_VOIDP) function get_memory_usage()
   end function get_memory_usage
  end interface

contains

  logical function loct_isinstringlist(a, s) result(inlist)
    integer,          intent(in) :: a
    character(len=*), intent(in) :: s

    integer, allocatable :: list(:)

    allocate(list(2**14))

    call oct_wfs_list(s, list)
    inlist = .false.
    if (list(a).eq.1) inlist = .true.

    deallocate(list)

  end function loct_isinstringlist

  subroutine loct_rm_status_files(current_label)
    character(len=*), intent(in) :: current_label

    call loct_rm('exec/'//trim(current_label)//'oct-status-running')
    call loct_rm('exec/'//trim(current_label)//'oct-status-finished')
    call loct_rm('exec/'//trim(current_label)//'oct-status-aborted')

  end subroutine loct_rm_status_files

#  define TYPE real(4)
#  define SUBNAME(x) s ## x
#  include "loct_inc.F90"
#  undef SUBNAME
#  undef TYPE

#  define TYPE real(8)
#  define SUBNAME(x) d ## x
#  include "loct_inc.F90"
#  undef SUBNAME
#  undef TYPE

#  define TYPE complex(4)
#  define SUBNAME(x) c ## x
#  include "loct_inc.F90"
#  undef SUBNAME
#  undef TYPE

#  define TYPE complex(8)
#  define SUBNAME(x) z ## x
#  include "loct_inc.F90"
#  undef SUBNAME
#  undef TYPE

#  define TYPE  integer
#  define SUBNAME(x) i ## x
#  include "loct_inc.F90"
#  undef SUBNAME
#  undef TYPE

#  define TYPE  character(len=*)
#  define SUBNAME(x) a ## x
#  include "loct_inc.F90"
#  undef SUBNAME
#  undef TYPE

#  define TYPE  logical
#  define SUBNAME(x) l ## x
#  include "loct_inc.F90"
#  undef SUBNAME
#  undef TYPE

end module loct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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

module write_iter_m

  implicit none

  ! Define which routines can be seen from the outside.
  private

  public ::                  &
    write_iter_init,         &
    write_iter_clear,        &
    write_iter_flush,        &
    write_iter_end,          &
    write_iter_start,        &
    write_iter_string,       &
    write_iter_header_start, &
    write_iter_header,       &
    write_iter_nl,           &
    write_iter_double,       &
    write_iter_int

  ! ---------------------------------------------------------
  !> write_iter functions
  interface
    subroutine write_iter_init(out,  iter, factor, file)
      use c_pointer_m
      type(c_ptr)      :: out
      integer          :: iter
      FLOAT            :: factor
      character(len=*) :: file
    end subroutine write_iter_init
    subroutine write_iter_clear(out)
      use c_pointer_m
      type(c_ptr) :: out
    end subroutine write_iter_clear
    subroutine write_iter_flush(out)
      use c_pointer_m
      type(c_ptr) :: out
    end subroutine write_iter_flush
    subroutine write_iter_end(out)
      use c_pointer_m
      type(c_ptr) :: out
    end subroutine write_iter_end
    subroutine write_iter_start(out)
      use c_pointer_m
      type(c_ptr) :: out
    end subroutine write_iter_start
    subroutine write_iter_string(out, string)
      use c_pointer_m
      type(c_ptr)      :: out
      character(len=*) :: string
    end subroutine write_iter_string
    subroutine write_iter_header_start(out)
      use c_pointer_m
      type(c_ptr) :: out
    end subroutine write_iter_header_start
    subroutine write_iter_header(out, string)
      use c_pointer_m
      type(c_ptr)      :: out
      character(len=*) :: string
    end subroutine write_iter_header
    subroutine write_iter_nl(out)
      use c_pointer_m
      type(c_ptr) :: out
    end subroutine write_iter_nl
  end interface

  interface write_iter_double
    subroutine write_iter_double_1(out, d, n)
      use c_pointer_m
      type(c_ptr) :: out
      integer   :: n
      real(8)   :: d
    end subroutine write_iter_double_1
    subroutine write_iter_double_n(out, d, n)
      use c_pointer_m
      type(c_ptr) :: out
      integer   :: n
      real(8)   :: d(n)
    end subroutine write_iter_double_n
    subroutine write_iter_float_1(out, d, n)
      use c_pointer_m
      type(c_ptr) :: out
      integer   :: n
      real(4)   :: d
    end subroutine write_iter_float_1
    subroutine write_iter_float_n(out, d, n)
      use c_pointer_m
      type(c_ptr) :: out
      integer   :: n
      real(4)   :: d(n)
    end subroutine write_iter_float_n
  end interface write_iter_double

  interface write_iter_int
    subroutine write_iter_int_1(out, i, n)
      use c_pointer_m
      type(c_ptr) :: out
      integer   :: n
      integer   :: i
    end subroutine write_iter_int_1
    subroutine write_iter_int_n(out, i, n)
      use c_pointer_m
      type(c_ptr) :: out
      integer   :: n
      integer   :: i(n)
    end subroutine write_iter_int_n
  end interface write_iter_int

end module write_iter_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

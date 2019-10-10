!! Copyright (C) 2019 S. Ohlmann
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

! -----------------------------------------------------------------------
!> This module contains interfaces for routines in allocate_special.c
! -----------------------------------------------------------------------
module allocate_special_oct_m
  use iso_c_binding

  implicit none

  public ! only interfaces in this module

  interface
    function dallocate_special(size) bind(c, name='dallocate_special')
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: dallocate_special
    end function dallocate_special
  end interface

  interface
    function zallocate_special(size) bind(c, name='zallocate_special')
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: zallocate_special
    end function zallocate_special
  end interface

  interface
    function sallocate_special(size) bind(c, name='sallocate_special')
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: sallocate_special
    end function sallocate_special
  end interface

  interface
    function callocate_special(size) bind(c, name='callocate_special')
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: callocate_special
    end function callocate_special
  end interface

  interface
    subroutine deallocate_special(array) bind(c, name='deallocate_special')
      import :: c_ptr
      type(c_ptr), value :: array
    end subroutine deallocate_special
  end interface

end module allocate_special_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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
!> This module contains interfaces for routines in allocate_hardware_aware.c
! -----------------------------------------------------------------------
module allocate_hardware_aware_oct_m
  use iso_c_binding

  implicit none

  public ! only interfaces in this module

  interface
    function dallocate_hardware_aware(size) bind(c)
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: dallocate_hardware_aware
    end function dallocate_hardware_aware
  end interface

  interface
    function zallocate_hardware_aware(size) bind(c)
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: zallocate_hardware_aware
    end function zallocate_hardware_aware
  end interface

  interface
    function sallocate_hardware_aware(size) bind(c)
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: sallocate_hardware_aware
    end function sallocate_hardware_aware
  end interface

  interface
    function callocate_hardware_aware(size) bind(c)
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: callocate_hardware_aware
    end function callocate_hardware_aware
  end interface

  interface
    subroutine deallocate_hardware_aware(array) bind(c)
      import :: c_ptr
      type(c_ptr), value :: array
    end subroutine deallocate_hardware_aware
  end interface

end module allocate_hardware_aware_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

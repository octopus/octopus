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
  use accel_oct_m
  use iso_c_binding

  implicit none

  private

  public :: &
    zallocate_hardware_aware, &
    dallocate_hardware_aware, &
    callocate_hardware_aware, &
    sallocate_hardware_aware, &
    deallocate_hardware_aware

  interface
    function dallocate_aligned(size) bind(c)
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: dallocate_aligned
    end function dallocate_aligned
  end interface

  interface
    function zallocate_aligned(size) bind(c)
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: zallocate_aligned
    end function zallocate_aligned
  end interface

  interface
    function sallocate_aligned(size) bind(c)
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: sallocate_aligned
    end function sallocate_aligned
  end interface

  interface
    function callocate_aligned(size) bind(c)
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: callocate_aligned
    end function callocate_aligned
  end interface

  interface
    subroutine deallocate_aligned(array) bind(c)
      import :: c_ptr
      type(c_ptr), value :: array
    end subroutine deallocate_aligned
  end interface

  interface
    function dallocate_pinned(size) bind(c)
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: dallocate_pinned
    end function dallocate_pinned
  end interface

  interface
    function zallocate_pinned(size) bind(c)
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: zallocate_pinned
    end function zallocate_pinned
  end interface

  interface
    function sallocate_pinned(size) bind(c)
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: sallocate_pinned
    end function sallocate_pinned
  end interface

  interface
    function callocate_pinned(size) bind(c)
      import :: c_ptr, c_int
      integer(c_int), value :: size
      type(c_ptr) :: callocate_pinned
    end function callocate_pinned
  end interface

  interface
    subroutine deallocate_pinned(array) bind(c)
      import :: c_ptr
      type(c_ptr), value :: array
    end subroutine deallocate_pinned
  end interface

contains

  function zallocate_hardware_aware(size)
    integer(c_int) :: size
    type(c_ptr) :: zallocate_hardware_aware

    ! allocate pinned memory for GPU runs, otherwise aligned memory
    if(accel_is_enabled()) then
      zallocate_hardware_aware = zallocate_pinned(size)
    else
      zallocate_hardware_aware = zallocate_aligned(size)
    end if
  end function zallocate_hardware_aware

  function dallocate_hardware_aware(size)
    integer(c_int) :: size
    type(c_ptr) :: dallocate_hardware_aware

    ! allocate pinned memory for GPU runs, otherwise aligned memory
    if(accel_is_enabled()) then
      dallocate_hardware_aware = dallocate_pinned(size)
    else
      dallocate_hardware_aware = dallocate_aligned(size)
    end if
  end function dallocate_hardware_aware

  function callocate_hardware_aware(size)
    integer(c_int) :: size
    type(c_ptr) :: callocate_hardware_aware

    ! allocate pinned memory for GPU runs, otherwise aligned memory
    if(accel_is_enabled()) then
      callocate_hardware_aware = callocate_pinned(size)
    else
      callocate_hardware_aware = callocate_aligned(size)
    end if
  end function callocate_hardware_aware

  function sallocate_hardware_aware(size)
    integer(c_int) :: size
    type(c_ptr) :: sallocate_hardware_aware

    ! allocate pinned memory for GPU runs, otherwise aligned memory
    if(accel_is_enabled()) then
      sallocate_hardware_aware = sallocate_pinned(size)
    else
      sallocate_hardware_aware = sallocate_aligned(size)
    end if
  end function sallocate_hardware_aware

  subroutine deallocate_hardware_aware(array)
    type(c_ptr), value :: array

    ! deallocate pinned memory for GPU runs, otherwise aligned memory
    if(accel_is_enabled()) then
      call deallocate_pinned(array)
    else
      call deallocate_aligned(array)
    end if
  end subroutine deallocate_hardware_aware

end module allocate_hardware_aware_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

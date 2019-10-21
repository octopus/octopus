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

! This is the fortran part of the wrapper around the NVTX 
! (NVIDIA Tools Extension) profiling functions.

#include "global.h"

module nvtx_oct_m

  implicit none

  private
  public ::          &
    nvtx_range_push, &
    nvtx_range_pop

  interface

    subroutine nvtx_range_push(range_name)
      use iso_c_binding
      implicit none
      
      character(len=*), intent(in)    :: range_name
    end subroutine nvtx_range_push

    subroutine nvtx_range_pop()
      use iso_c_binding
      implicit none
    end subroutine nvtx_range_pop
  end interface
  
end module nvtx_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

!! Copyright (C) 2008 X. Andrade
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
!! $Id: global.F90 3856 2008-03-10 14:52:36Z marques $

#include "global.h"
module hardware_m
  
  private
  
  public ::             &
    hardware,           &
    cache_t,            &
    hardware_t,         &
    hardware_init,      &
    hardware_end

  type cache_t
    integer :: size
    integer :: line_size
  end type cache_t

  type hardware_t 
    type(cache_t) :: l1
    type(cache_t) :: l2
  end type hardware_t
  
  type(hardware_t) :: hardware

contains

  subroutine hardware_init

    !for the moment we will use fixed values

    hardware%l1%size = 32*1024

    hardware%l1%line_size = 64

    hardware%l2%size = 4096*1024

    hardware%l2%line_size = 64

  end subroutine hardware_init

  subroutine hardware_end

  end subroutine hardware_end

end module hardware_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

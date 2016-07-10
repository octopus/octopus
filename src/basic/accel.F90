!! Copyright (C) 2016 X. Andrade
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
!! $Id: opencl.F90 15203 2016-03-19 13:15:05Z xavier $

#include "global.h"

module accel_oct_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use types_oct_m
  
  implicit none 

  private

  public ::                 &
    accel_context_t,        &
    accel_device_t,         &
    accel_mem_t,            &
    accel_kernel_t,         &
    accel_t
  
  type accel_context_t
#ifdef HAVE_OPENCL
    type(cl_context) :: cl_context
#else
    integer          :: dummy
#endif
  end type accel_context_t

  type accel_device_t
#ifdef HAVE_OPENCL
    type(cl_device_id) :: cl_device
#else
    integer         :: dummy
#endif
  end type accel_device_t

  type accel_t 
#ifdef HAVE_OPENCL
    type(accel_context_t)  :: context
    type(cl_command_queue) :: command_queue
    type(accel_device_t)   :: device
#endif
    integer                :: max_workgroup_size
    integer                :: local_memory_size
    logical                :: enabled
    logical                :: shared_mem
  end type accel_t

  type accel_mem_t
#ifdef HAVE_OPENCL
    type(cl_mem)           :: mem
#endif
    integer(SIZEOF_SIZE_T) :: size
    type(type_t)           :: type
    integer                :: flags
  end type accel_mem_t

  type accel_kernel_t
#ifdef HAVE_OPENCL
    type(cl_kernel)               :: kernel
#endif
    logical                       :: initialized = .false.
    type(accel_kernel_t), pointer :: next
    integer                       :: arg_count
  end type accel_kernel_t

  type(accel_t), public :: accel

end module accel_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

!! Copyright (C) 2019 X. Andrade
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

module cuda_oct_m

  implicit none

  private
  public ::                             &
    cuda_init,                          &
    cuda_end,                           &
    cuda_module_map_init,               &
    cuda_module_map_end,                &
    cuda_build_program,                 &
    cuda_create_kernel,                 &
    cuda_release_module,                &
    cuda_release_kernel,                &
    cuda_device_max_threads_per_block,  &
    cuda_device_total_memory,           &
    cuda_device_shared_memory,          &
    cuda_mem_alloc,                     &
    cuda_mem_free,                      &
    cuda_alloc_arg_array,               &
    cuda_free_arg_array,                &
    cuda_kernel_set_arg_buffer,         &
    cuda_context_synchronize,           &
    cuda_launch_kernel,                 &
    cuda_device_name,                   &
    cuda_device_capability,             &
    cuda_driver_version

  interface

    subroutine cuda_init(context, device, stream, device_number, rank)
      use iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: context
      type(c_ptr), intent(inout) :: device
      type(c_ptr), intent(inout) :: stream
      integer,     intent(out)   :: device_number
      integer,     intent(out)   :: rank
    end subroutine cuda_init
    
    ! -------------------------------------------------

    subroutine cuda_end(context, device)
      use iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: context
      type(c_ptr), intent(inout) :: device
    end subroutine cuda_end

    ! -------------------------------------------------

    subroutine cuda_module_map_init(module_map)
      use iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: module_map
    end subroutine cuda_module_map_init

    ! -------------------------------------------------

    subroutine cuda_module_map_end(module_map)
      use iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: module_map
    end subroutine cuda_module_map_end
    
    ! -------------------------------------------------

    subroutine cuda_build_program(module_map, modul, device, fname, flags)
      use iso_c_binding
      implicit none
      
      type(c_ptr),      intent(inout) :: module_map
      type(c_ptr),      intent(inout) :: modul
      type(c_ptr),      intent(inout) :: device
      character(len=*), intent(in)    :: fname
      character(len=*), intent(in)    :: flags
    end subroutine cuda_build_program
    
    ! -------------------------------------------------

    subroutine cuda_create_kernel(kernel, modul, kernel_name)
      use iso_c_binding
      implicit none
      
      type(c_ptr),      intent(inout) :: kernel
      type(c_ptr),      intent(inout) :: modul
      character(len=*), intent(in)    :: kernel_name
    end subroutine cuda_create_kernel
    ! -------------------------------------------------

    subroutine cuda_release_module(modul)
      use iso_c_binding
      implicit none
      
      type(c_ptr),      intent(inout) :: modul
    end subroutine cuda_release_module

    ! -------------------------------------------------

    subroutine cuda_release_kernel(kernel)
      use iso_c_binding
      implicit none
      
      type(c_ptr),      intent(inout) :: kernel
    end subroutine cuda_release_kernel
    
    ! -------------------------------------------------

    subroutine cuda_device_max_threads_per_block(device, max_threads)
      use iso_c_binding
      implicit none
      
      type(c_ptr),      intent(inout) :: device
      integer,          intent(out)   :: max_threads
    end subroutine cuda_device_max_threads_per_block

    ! -------------------------------------------------

    subroutine cuda_device_total_memory(device, total_memory)
      use iso_c_binding
      implicit none
      
      type(c_ptr),      intent(inout) :: device
      integer(8),       intent(out)   :: total_memory
    end subroutine cuda_device_total_memory

    ! -------------------------------------------------

    subroutine cuda_device_shared_memory(device, shared_memory)
      use iso_c_binding
      implicit none
      
      type(c_ptr),      intent(inout) :: device
      integer(8),       intent(out)   :: shared_memory
    end subroutine cuda_device_shared_memory
    
    ! -------------------------------------------------
    
    subroutine cuda_mem_alloc(cuda_ptr, size)
      use iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: cuda_ptr
      integer(8),  intent(in)    :: size
    end subroutine cuda_mem_alloc

    ! -------------------------------------------------
    
    subroutine cuda_mem_free(cuda_ptr)
      use iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: cuda_ptr
    end subroutine cuda_mem_free
    
    ! -------------------------------------------------
    
    subroutine cuda_alloc_arg_array(arg_array)
      use iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: arg_array
    end subroutine cuda_alloc_arg_array

    ! -------------------------------------------------
    
    subroutine cuda_free_arg_array(arg_array)
      use iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: arg_array
    end subroutine cuda_free_arg_array

    ! -------------------------------------------------
    
    subroutine cuda_kernel_set_arg_buffer(arg_array, cuda_ptr, arg_index)
      use iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: arg_array
      type(c_ptr), intent(in)    :: cuda_ptr
      integer,     intent(in)    :: arg_index      
    end subroutine cuda_kernel_set_arg_buffer

    ! -------------------------------------------------
    
    subroutine cuda_context_synchronize()
      implicit none
    end subroutine cuda_context_synchronize

    ! -------------------------------------------------
        
    subroutine cuda_launch_kernel(kernel, griddim, blockdim, shared_mem, arg_array)
      use iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: kernel
      integer(8),  intent(in)    :: griddim
      integer(8),  intent(in)    :: blockdim
      integer(8),  intent(in)    :: shared_mem
      type(c_ptr), intent(inout) :: arg_array
    end subroutine cuda_launch_kernel

    ! -------------------------------------------------

    subroutine cuda_device_name(device, name)
      use iso_c_binding
      implicit none
      
      type(c_ptr),      intent(inout) :: device
      character(len=*), intent(inout) :: name
    end subroutine cuda_device_name

    ! -------------------------------------------------

    subroutine cuda_device_capability(device, major, minor)
      use iso_c_binding
      implicit none
      
      type(c_ptr),      intent(inout) :: device
      integer,          intent(out)   :: major
      integer,          intent(out)   :: minor
    end subroutine cuda_device_capability
    
    ! -------------------------------------------------

    subroutine cuda_driver_version(version)
      use iso_c_binding
      implicit none
      
      integer,       intent(out)   :: version
    end subroutine cuda_driver_version

    subroutine cuda_deref(cuda_ptr, cuda_deref_ptr)
      use iso_c_binding
      implicit none

      type(c_ptr), intent(in)  :: cuda_ptr
      type(c_ptr), intent(out) :: cuda_deref_ptr
    end subroutine cuda_deref

  end interface
  
end module cuda_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

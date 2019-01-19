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
  use global_oct_m
  use io_oct_m
  use loct_oct_m
  use messages_oct_m
  use profiling_oct_m
  
  implicit none

  private
  public ::                    &
    cuda_init,                 &
    cuda_end,                  &
    cuda_mem_alloc,            &
    cuda_module_map_init,      &
    cuda_module_map_end,       &
    cuda_build_program,        &
    cuda_create_kernel,        &
    cuda_release_module,       &
    cuda_release_kernel
  
  interface

    subroutine cuda_init(context, device)
      use iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: context
      type(c_ptr), intent(inout) :: device
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
    
    subroutine cuda_mem_alloc(cuda_ptr, size)
      use iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: cuda_ptr
      integer(8),  intent(in)    :: size
    end subroutine cuda_mem_alloc
    
  end interface
  
end module cuda_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

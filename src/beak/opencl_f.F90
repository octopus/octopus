!! Copyright (C) 2010 X. Andrade
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
!! $Id: opencl_f.F90 3587 2007-11-22 16:43:00Z xavier $

#include "global.h"

module opencl_m
  use c_pointer_m
  use global_m

  implicit none 

  private

  public ::                   &
    opencl_t,                 &
    opencl_init,              &
    opencl_end,               &
    opencl_mem_t,             &
    dopencl_create_buffer,    &
    zopencl_create_buffer,    &
    iopencl_create_buffer,    &
    opencl_write_buffer,      &
    opencl_release_buffer

  type opencl_t 
    type(c_ptr) :: env
    type(c_ptr) :: kernel_vpsi
  end type opencl_t

  type opencl_mem_t
    type(c_ptr) :: mem
  end type opencl_mem_t

  type(opencl_t), public :: opencl

  ! this values are copied from OpenCL include CL/cl.h
  integer, parameter, public ::        &
    CL_MEM_READ_WRITE = 1,             &
    CL_MEM_WRITE_ONLY = 2,             &
    CL_MEM_READ_ONLY  = 4

  interface
    subroutine f90_opencl_env_init(this, source_path)
      use c_pointer_m

      type(c_ptr),      intent(out) :: this
      character(len=*), intent(in)  :: source_path
    end subroutine f90_opencl_env_init

    ! ----------------------------------------------------

    subroutine f90_opencl_env_end(this)
      use c_pointer_m

      type(c_ptr), intent(inout) :: this
    end subroutine f90_opencl_env_end

    ! ----------------------------------------------------

    subroutine f90_opencl_kernel_init(this, env, source_file, kernel_base_name)
      use c_pointer_m

      type(c_ptr),      intent(out)   :: this
      type(c_ptr),      intent(inout) :: env
      character(len=*), intent(in)    :: source_file
      character(len=*), intent(in)    :: kernel_base_name
    end subroutine f90_opencl_kernel_init

    ! ----------------------------------------------------

    subroutine f90_opencl_kernel_end(this)
      use c_pointer_m

      type(c_ptr), intent(inout) :: this
    end subroutine f90_opencl_kernel_end

    ! ----------------------------------------------------

    subroutine f90_opencl_create_buffer(this, env, flags, size)
      use c_pointer_m

      type(c_ptr),            intent(inout) :: this
      type(c_ptr),            intent(inout) :: env
      integer,                intent(in)    :: flags
      integer(SIZEOF_SIZE_T), intent(in)    :: size
    end subroutine f90_opencl_create_buffer
    
    ! ----------------------------------------------------

    subroutine f90_opencl_release_buffer(this)
      use c_pointer_m

      type(c_ptr),            intent(inout) :: this
    end subroutine f90_opencl_release_buffer

    ! ----------------------------------------------------

  end interface

  interface opencl_write_buffer
    module procedure iopencl_write_buffer, dopencl_write_buffer, zopencl_write_buffer
  end interface

  contains
    
    subroutine opencl_init(this)
      type(opencl_t), intent(out) :: this
      type(c_ptr) :: buf
      call f90_opencl_env_init(this%env, trim(conf%share)//'/opencl/')      

      ! now initialize the kernels
      call f90_opencl_kernel_init(this%kernel_vpsi, this%env, "vpsi.cl", "vpsi");
    end subroutine opencl_init

    ! ------------------------------------------

    subroutine opencl_end(this)
      type(opencl_t), intent(inout) :: this
      
      call f90_opencl_kernel_end(this%kernel_vpsi);
      call f90_opencl_env_end(this%env)
    end subroutine opencl_end

    ! ------------------------------------------

    subroutine opencl_release_buffer(this)
      type(opencl_mem_t), intent(inout) :: this

      call f90_opencl_release_buffer(this%mem)
    end subroutine opencl_release_buffer

    ! ------------------------------------------

#include "undef.F90"
#include "real.F90"
#include "opencl_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "opencl_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "opencl_inc.F90"

end module opencl_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

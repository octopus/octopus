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
  use types_m
  use profiling_m

  implicit none 

  private

  public ::                   &
    opencl_is_available,      &
    opencl_init,              &
    opencl_end,               &
    opencl_mem_t,             &
    opencl_create_buffer,     &
    opencl_write_buffer,      &
    opencl_read_buffer,       &
    opencl_release_buffer,    &
    opencl_padded_size,       &
    opencl_finish,            &
    opencl_set_kernel_arg,    &
    opencl_kernel_run

  type opencl_t 
    type(c_ptr) :: env
  end type opencl_t

  type opencl_mem_t
    type(c_ptr) :: mem
  end type opencl_mem_t

  type(opencl_t) :: opencl
  type(c_ptr), public :: kernel_vpsi

  ! this values are copied from OpenCL include CL/cl.h
  integer, parameter, public ::        &
    CL_MEM_READ_WRITE = 1,             &
    CL_MEM_WRITE_ONLY = 2,             &
    CL_MEM_READ_ONLY  = 4

  ! this function are defined in opencl_low.c
  interface
    ! ---------------------------------------------------

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

    subroutine f90_opencl_finish(this)
      use c_pointer_m

      type(c_ptr),            intent(inout) :: this
    end subroutine f90_opencl_finish

    ! ----------------------------------------------------

    subroutine f90_opencl_set_kernel_arg_buf(kernel, type, index, buffer)
      use c_pointer_m

      type(c_ptr), intent(inout) :: kernel
      integer,     intent(in)    :: type
      integer,     intent(in)    :: index
      type(c_ptr), intent(in)    :: buffer
    end subroutine f90_opencl_set_kernel_arg_buf
    
    ! ----------------------------------------------------

    subroutine f90_opencl_kernel_run(kernel, env, type, dim, globalsizes, localsizes)
      use c_pointer_m

      type(c_ptr),            intent(inout) :: kernel
      type(c_ptr),            intent(inout) :: env
      integer,                intent(in)    :: type
      integer,                intent(in)    :: dim
      integer(SIZEOF_SIZE_T), intent(in)    :: globalsizes
      integer(SIZEOF_SIZE_T), intent(in)    :: localsizes
    end subroutine f90_opencl_kernel_run

  end interface

  interface opencl_create_buffer
    module procedure opencl_create_buffer_4, opencl_create_buffer_8
  end interface

  interface opencl_write_buffer
    module procedure iopencl_write_buffer, dopencl_write_buffer, zopencl_write_buffer
  end interface

  interface opencl_read_buffer
    module procedure iopencl_read_buffer, dopencl_read_buffer, zopencl_read_buffer
  end interface

  interface opencl_set_kernel_arg
    module procedure                 &
      opencl_set_kernel_arg_buffer,  &
      iopencl_set_kernel_arg_data,   &
      dopencl_set_kernel_arg_data,   &
      zopencl_set_kernel_arg_data
  end interface
  
  type(profile_t), save :: prof_read, prof_write, prof_kernel_run
  
  contains

    pure logical function opencl_is_available() result(available)
    
      available = .true.
    end function opencl_is_available

    ! ------------------------------------------

    subroutine opencl_init()
      call f90_opencl_env_init(opencl%env, trim(conf%share)//'/opencl/')      

      ! now initialize the kernels
      call f90_opencl_kernel_init(kernel_vpsi, opencl%env, "vpsi.cl", "vpsi");
    end subroutine opencl_init

    ! ------------------------------------------

    subroutine opencl_end()
      
      call f90_opencl_kernel_end(kernel_vpsi)
      call f90_opencl_env_end(opencl%env)
    end subroutine opencl_end

    ! ------------------------------------------

    subroutine opencl_create_buffer_4(this, flags, type, size)
      type(opencl_mem_t), intent(inout) :: this
      integer,            intent(in)    :: flags
      integer,            intent(in)    :: type
      integer,            intent(in)    :: size
      
      integer(SIZEOF_SIZE_T) :: fsize
      
      fsize = size*types_get_size(type)
      
      call f90_opencl_create_buffer(this%mem, opencl%env, flags, fsize)
      
    end subroutine opencl_create_buffer_4

    ! ------------------------------------------

    subroutine opencl_create_buffer_8(this, flags, type, size)
      type(opencl_mem_t),     intent(inout) :: this
      integer,                intent(in)    :: flags
      integer,                intent(in)    :: type
      integer(SIZEOF_SIZE_T), intent(in)    :: size
      
      integer(SIZEOF_SIZE_T) :: fsize
      
      fsize = size*types_get_size(type)

      call f90_opencl_create_buffer(this%mem, opencl%env, flags, fsize)
      
    end subroutine opencl_create_buffer_8

    ! ------------------------------------------

    subroutine opencl_release_buffer(this)
      type(opencl_mem_t), intent(inout) :: this

      call f90_opencl_release_buffer(this%mem)
    end subroutine opencl_release_buffer

    ! ------------------------------------------

    integer(SIZEOF_SIZE_T) function opencl_padded_size(nn) result(psize)
      integer,        intent(in) :: nn

      integer :: modnn, bsize

      bsize = 512

      psize = nn
      modnn = mod(nn, bsize)
      if(modnn /= 0) psize = psize + bsize - modnn

    end function opencl_padded_size

    ! ------------------------------------------

    subroutine opencl_finish()
      call f90_opencl_finish(opencl%env)
    end subroutine opencl_finish

    ! ------------------------------------------

    subroutine opencl_set_kernel_arg_buffer(kernel, type, narg, buffer)
      type(c_ptr),        intent(inout) :: kernel
      integer,            intent(in)    :: type
      integer,            intent(in)    :: narg
      type(opencl_mem_t), intent(in)    :: buffer
      
      call f90_opencl_set_kernel_arg_buf(kernel, type, narg, buffer%mem)

    end subroutine opencl_set_kernel_arg_buffer

    ! ------------------------------------------

    subroutine opencl_kernel_run(kernel, type, globalsizes, localsizes)
      type(c_ptr),        intent(inout) :: kernel
      integer,            intent(in)    :: type
      integer,            intent(in)    :: globalsizes(:)
      integer,            intent(in)    :: localsizes(:)
      
      integer :: dim
      integer(SIZEOF_SIZE_T) :: gsizes(1:3)
      integer(SIZEOF_SIZE_T) :: lsizes(1:3)

      call profiling_in(prof_kernel_run, "CL_KERNEL_RUN")

      dim = ubound(globalsizes, dim = 1)
      ASSERT(dim == ubound(localsizes, dim = 1))
      
      gsizes(1:dim) = int(globalsizes(1:dim), SIZEOF_SIZE_T)
      lsizes(1:dim) = int(localsizes(1:dim), SIZEOF_SIZE_T)

      ASSERT(all(mod(gsizes, lsizes) == 0))

      call f90_opencl_kernel_run(kernel, opencl%env, type, dim, gsizes(1), lsizes(1))
      call opencl_finish()
      call profiling_out(prof_kernel_run)

    end subroutine opencl_kernel_run

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

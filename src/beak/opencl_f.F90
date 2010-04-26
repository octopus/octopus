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
    opencl_padded_size

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

    integer(SIZEOF_SIZE_T) function opencl_padded_size(nn, type) result(psize)
      integer,        intent(in) :: nn
      integer,        intent(in) :: type

      integer :: modnn, bsize

      bsize = 1024/types_get_size(type)

      psize = nn
      modnn = mod(nn, bsize)
      if(modnn /= 0) psize = psize + bsize - modnn
      
    end function opencl_padded_size

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

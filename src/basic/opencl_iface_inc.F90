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
  
  
  ! these values are copied from OpenCL include CL/cl.h
  
  integer, parameter  ::  CL_SUCCESS = 0
  integer, parameter  ::  CL_DEVICE_NOT_FOUND = -1
  integer, parameter  ::  CL_DEVICE_NOT_AVAILABLE = -2
  integer, parameter  ::  CL_COMPILER_NOT_AVAILABLE = -3
  integer, parameter  ::  CL_MEM_OBJECT_ALLOCATION_FAIL = -4
  integer, parameter  ::  CL_OUT_OF_RESOURCES = -5
  integer, parameter  ::  CL_OUT_OF_HOST_MEMORY = -6
  integer, parameter  ::  CL_PROFILING_INFO_NOT_AVAILABLE = -7
  integer, parameter  ::  CL_MEM_COPY_OVERLAP = -8
  integer, parameter  ::  CL_IMAGE_FORMAT_MISMATCH = -9
  integer, parameter  ::  CL_IMAGE_FORMAT_NOT_SUPPORTED = -10
  integer, parameter  ::  CL_BUILD_PROGRAM_FAILURE = -11
  integer, parameter  ::  CL_MAP_FAILURE = -12
  integer, parameter  ::  CL_INVALID_VALUE = -30
  integer, parameter  ::  CL_INVALID_DEVICE_TYPE = -31
  integer, parameter  ::  CL_INVALID_PLATFORM = -32
  integer, parameter  ::  CL_INVALID_DEVICE = -33
  integer, parameter  ::  CL_INVALID_CONTEXT = -34
  integer, parameter  ::  CL_INVALID_QUEUE_PROPERTIES = -35
  integer, parameter  ::  CL_INVALID_COMMAND_QUEUE = -36
  integer, parameter  ::  CL_INVALID_HOST_PTR = -37
  integer, parameter  ::  CL_INVALID_MEM_OBJECT = -38
  integer, parameter  ::  CL_INVALID_IMAGE_FORMAT_DESC = -39
  integer, parameter  ::  CL_INVALID_IMAGE_SIZE = -40
  integer, parameter  ::  CL_INVALID_SAMPLER = -41
  integer, parameter  ::  CL_INVALID_BINARY = -42
  integer, parameter  ::  CL_INVALID_BUILD_OPTIONS = -43
  integer, parameter  ::  CL_INVALID_PROGRAM = -44
  integer, parameter  ::  CL_INVALID_PROGRAM_EXECUTABLE = -45
  integer, parameter  ::  CL_INVALID_KERNEL_NAME = -46
  integer, parameter  ::  CL_INVALID_KERNEL_DEFINITION = -47
  integer, parameter  ::  CL_INVALID_KERNEL = -48
  integer, parameter  ::  CL_INVALID_ARG_INDEX = -49
  integer, parameter  ::  CL_INVALID_ARG_VALUE = -50
  integer, parameter  ::  CL_INVALID_ARG_SIZE = -51
  integer, parameter  ::  CL_INVALID_KERNEL_ARGS = -52
  integer, parameter  ::  CL_INVALID_WORK_DIMENSION = -53
  integer, parameter  ::  CL_INVALID_WORK_GROUP_SIZE = -54
  integer, parameter  ::  CL_INVALID_WORK_ITEM_SIZE = -55
  integer, parameter  ::  CL_INVALID_GLOBAL_OFFSET = -56
  integer, parameter  ::  CL_INVALID_EVENT_WAIT_LIST = -57
  integer, parameter  ::  CL_INVALID_EVENT = -58
  integer, parameter  ::  CL_INVALID_OPERATION = -59
  integer, parameter  ::  CL_INVALID_GL_OBJECT = -60
  integer, parameter  ::  CL_INVALID_BUFFER_SIZE = -61
  integer, parameter  ::  CL_INVALID_MIP_LEVEL = -62
  integer, parameter  ::  CL_INVALID_GLOBAL_WORK_SIZE = -63


  integer, parameter, public ::        &
    CL_MEM_READ_WRITE = 1,             &
    CL_MEM_WRITE_ONLY = 2,             &
    CL_MEM_READ_ONLY  = 4

  ! this function are defined in opencl_low.c
  interface
    ! ---------------------------------------------------

    subroutine f90_cl_env_init(env, source_path)
      use c_pointer_m

      implicit none

      type(c_ptr),      intent(out) :: env
      character(len=*), intent(in)  :: source_path
    end subroutine f90_cl_env_init

    ! ----------------------------------------------------

    subroutine f90_cl_env_end(env)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: env
    end subroutine f90_cl_env_end

    ! ----------------------------------------------------

    integer function f90_cl_max_workgroup_size(env)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: env
    end function f90_cl_max_workgroup_size

    ! ----------------------------------------------------

    subroutine f90_cl_build_program(prog, env, source_file)
      use c_pointer_m

      implicit none

      type(c_ptr),      intent(out)   :: prog
      type(c_ptr),      intent(inout) :: env
      character(len=*), intent(in)    :: source_file
    end subroutine f90_cl_build_program

    ! ----------------------------------------------------

    subroutine f90_cl_release_program(prog, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: prog
      integer,     intent(out)   :: ierr
    end subroutine f90_cl_release_program

    ! ----------------------------------------------------

    subroutine f90_cl_create_kernel(kernel, prog, kernel_name)
      use c_pointer_m

      implicit none

      type(c_ptr),      intent(out)   :: kernel
      type(c_ptr),      intent(inout) :: prog
      character(len=*), intent(in)    :: kernel_name
    end subroutine f90_cl_create_kernel

    ! ----------------------------------------------------

    subroutine f90_cl_release_kernel(kernel)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: kernel
    end subroutine f90_cl_release_kernel

    ! ----------------------------------------------------

    integer function f90_cl_kernel_wgroup_size(kernel, env)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: kernel
      type(c_ptr), intent(inout) :: env
    end function f90_cl_kernel_wgroup_size

    ! ----------------------------------------------------

    subroutine f90_cl_create_buffer(this, env, flags, size)
      use c_pointer_m

      implicit none

      type(c_ptr),            intent(inout) :: this
      type(c_ptr),            intent(inout) :: env
      integer,                intent(in)    :: flags
      integer(SIZEOF_SIZE_T), intent(in)    :: size
    end subroutine f90_cl_create_buffer

    ! ----------------------------------------------------

    subroutine f90_cl_release_buffer(this)
      use c_pointer_m

      implicit none

      type(c_ptr),            intent(inout) :: this
    end subroutine f90_cl_release_buffer

    ! ----------------------------------------------------

    subroutine f90_cl_finish(this)
      use c_pointer_m

      implicit none

      type(c_ptr),            intent(inout) :: this
    end subroutine f90_cl_finish

    ! ----------------------------------------------------

    subroutine f90_cl_set_kernel_arg_buf(kernel, index, buffer, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: kernel
      integer,     intent(in)    :: index
      type(c_ptr), intent(in)    :: buffer
      integer,     intent(in)    :: ierr
    end subroutine f90_cl_set_kernel_arg_buf

    ! ----------------------------------------------------

    subroutine f90_cl_set_kernel_arg_local(kernel, index, size, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: kernel
      integer,     intent(in)    :: index
      integer,     intent(in)    :: size
      integer,     intent(in)    :: ierr
    end subroutine f90_cl_set_kernel_arg_local

    ! ----------------------------------------------------

    subroutine f90_cl_kernel_run(kernel, env, dim, globalsizes, localsizes, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr),            intent(inout) :: kernel
      type(c_ptr),            intent(inout) :: env
      integer,                intent(in)    :: dim
      integer(SIZEOF_SIZE_T), intent(in)    :: globalsizes
      integer(SIZEOF_SIZE_T), intent(in)    :: localsizes
      integer,                intent(in)    :: ierr
    end subroutine f90_cl_kernel_run

  end interface

  !! Local Variables:
  !! mode: f90
  !! coding: utf-8
  !! End:

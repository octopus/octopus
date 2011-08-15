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
  
  integer, parameter, public ::  CL_SUCCESS = 0
  integer, parameter, public ::  CL_DEVICE_NOT_FOUND = -1
  integer, parameter, public ::  CL_DEVICE_NOT_AVAILABLE = -2
  integer, parameter, public ::  CL_COMPILER_NOT_AVAILABLE = -3
  integer, parameter, public ::  CL_MEM_OBJECT_ALLOCATION_FAIL = -4
  integer, parameter, public ::  CL_OUT_OF_RESOURCES = -5
  integer, parameter, public ::  CL_OUT_OF_HOST_MEMORY = -6
  integer, parameter, public ::  CL_PROFILING_INFO_NOT_AVAILABLE = -7
  integer, parameter, public ::  CL_MEM_COPY_OVERLAP = -8
  integer, parameter, public ::  CL_IMAGE_FORMAT_MISMATCH = -9
  integer, parameter, public ::  CL_IMAGE_FORMAT_NOT_SUPPORTED = -10
  integer, parameter, public ::  CL_BUILD_PROGRAM_FAILURE = -11
  integer, parameter, public ::  CL_MAP_FAILURE = -12
  integer, parameter, public ::  CL_INVALID_VALUE = -30
  integer, parameter, public ::  CL_INVALID_DEVICE_TYPE = -31
  integer, parameter, public ::  CL_INVALID_PLATFORM = -32
  integer, parameter, public ::  CL_INVALID_DEVICE = -33
  integer, parameter, public ::  CL_INVALID_CONTEXT = -34
  integer, parameter, public ::  CL_INVALID_QUEUE_PROPERTIES = -35
  integer, parameter, public ::  CL_INVALID_COMMAND_QUEUE = -36
  integer, parameter, public ::  CL_INVALID_HOST_PTR = -37
  integer, parameter, public ::  CL_INVALID_MEM_OBJECT = -38
  integer, parameter, public ::  CL_INVALID_IMAGE_FORMAT_DESC = -39
  integer, parameter, public ::  CL_INVALID_IMAGE_SIZE = -40
  integer, parameter, public ::  CL_INVALID_SAMPLER = -41
  integer, parameter, public ::  CL_INVALID_BINARY = -42
  integer, parameter, public ::  CL_INVALID_BUILD_OPTIONS = -43
  integer, parameter, public ::  CL_INVALID_PROGRAM = -44
  integer, parameter, public ::  CL_INVALID_PROGRAM_EXECUTABLE = -45
  integer, parameter, public ::  CL_INVALID_KERNEL_NAME = -46
  integer, parameter, public ::  CL_INVALID_KERNEL_DEFINITION = -47
  integer, parameter, public ::  CL_INVALID_KERNEL = -48
  integer, parameter, public ::  CL_INVALID_ARG_INDEX = -49
  integer, parameter, public ::  CL_INVALID_ARG_VALUE = -50
  integer, parameter, public ::  CL_INVALID_ARG_SIZE = -51
  integer, parameter, public ::  CL_INVALID_KERNEL_ARGS = -52
  integer, parameter, public ::  CL_INVALID_WORK_DIMENSION = -53
  integer, parameter, public ::  CL_INVALID_WORK_GROUP_SIZE = -54
  integer, parameter, public ::  CL_INVALID_WORK_ITEM_SIZE = -55
  integer, parameter, public ::  CL_INVALID_GLOBAL_OFFSET = -56
  integer, parameter, public ::  CL_INVALID_EVENT_WAIT_LIST = -57
  integer, parameter, public ::  CL_INVALID_EVENT = -58
  integer, parameter, public ::  CL_INVALID_OPERATION = -59
  integer, parameter, public ::  CL_INVALID_GL_OBJECT = -60
  integer, parameter, public ::  CL_INVALID_BUFFER_SIZE = -61
  integer, parameter, public ::  CL_INVALID_MIP_LEVEL = -62
  integer, parameter, public ::  CL_INVALID_GLOBAL_WORK_SIZE = -63


  integer, parameter, public ::        &
    CL_MEM_READ_WRITE = 1,             &
    CL_MEM_WRITE_ONLY = 2,             &
    CL_MEM_READ_ONLY  = 4

  ! this function are defined in opencl_low.c
  interface

    ! ---------------------------------------------------

    integer function flGetPlatformIDs(iplatform, platform_id)
      use c_pointer_m

      implicit none
      integer,          intent(in)  :: iplatform
      type(c_ptr),      intent(out) :: platform_id
    end function flGetPlatformIDs

    ! ---------------------------------------------------

    subroutine f90_cl_init_context(platform_id, context)
      use c_pointer_m

      implicit none
      type(c_ptr),      intent(in)   :: platform_id
      type(c_ptr),      intent(out)  :: context
    end subroutine f90_cl_init_context
    
    ! ---------------------------------------------------

    subroutine f90_cl_init_device(idevice, platform, context, device)
      use c_pointer_m

      implicit none

      integer,          intent(in)    :: idevice
      type(c_ptr),      intent(inout) :: platform
      type(c_ptr),      intent(inout) :: context
      type(c_ptr),      intent(out)   :: device
    end subroutine f90_cl_init_device

    ! ----------------------------------------------------

    subroutine flReleaseContext(context)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: context
    end subroutine flReleaseContext

    ! ----------------------------------------------------

    subroutine flCreateCommandQueue(command_queue, context, device, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: command_queue
      type(c_ptr), intent(inout) :: context
      type(c_ptr), intent(inout) :: device
      integer,     intent(out)   :: ierr
      
    end subroutine flCreateCommandQueue

    ! ----------------------------------------------------

    subroutine flReleaseCommandQueue(command_queue, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: command_queue
      integer,     intent(out)   :: ierr
      
    end subroutine flReleaseCommandQueue

    ! ----------------------------------------------------

    integer function f90_cl_max_workgroup_size(device)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: device
    end function f90_cl_max_workgroup_size


    ! ----------------------------------------------------

    integer function f90_cl_device_has_extension(device, extension)
      use c_pointer_m

      implicit none

      type(c_ptr),      intent(inout) :: device
      character(len=*), intent(in)    :: extension
    end function f90_cl_device_has_extension

    ! ----------------------------------------------------

    integer function f90_cl_device_local_mem_size(device)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: device
    end function f90_cl_device_local_mem_size

    ! ----------------------------------------------------

    integer function f90_cl_device_max_constant_buffer_size(device)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: device
    end function f90_cl_device_max_constant_buffer_size

    ! ----------------------------------------------------

    subroutine f90_cl_create_program_from_file(prog, context, source_file)
      use c_pointer_m

      implicit none

      type(c_ptr),      intent(out)   :: prog
      type(c_ptr),      intent(inout) :: context
      character(len=*), intent(in)    :: source_file
    end subroutine f90_cl_create_program_from_file

    ! ----------------------------------------------------

    subroutine f90_cl_build_program(prog, context, device, flags)
      use c_pointer_m

      implicit none

      type(c_ptr),      intent(inout) :: prog
      type(c_ptr),      intent(inout) :: context
      type(c_ptr),      intent(inout) :: device
      character(len=*), intent(in)    :: flags
    end subroutine f90_cl_build_program

    ! ----------------------------------------------------

    subroutine f90_cl_release_program(prog, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: prog
      integer,     intent(out)   :: ierr
    end subroutine f90_cl_release_program

    ! ----------------------------------------------------

    subroutine f90_cl_create_kernel(kernel, prog, kernel_name, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr),      intent(out)   :: kernel
      type(c_ptr),      intent(inout) :: prog
      character(len=*), intent(in)    :: kernel_name
      integer,          intent(out)   :: ierr
    end subroutine f90_cl_create_kernel

    ! ----------------------------------------------------

    subroutine f90_cl_release_kernel(kernel, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: kernel
      integer,     intent(out)   :: ierr
    end subroutine f90_cl_release_kernel

    ! ----------------------------------------------------

    integer function f90_cl_kernel_wgroup_size(kernel, device)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: kernel
      type(c_ptr), intent(inout) :: device
    end function f90_cl_kernel_wgroup_size

    ! ----------------------------------------------------

    subroutine f90_cl_create_buffer(this, context, flags, size, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr),            intent(inout) :: this
      type(c_ptr),            intent(inout) :: context
      integer,                intent(in)    :: flags
      integer(SIZEOF_SIZE_T), intent(in)    :: size
      integer,                intent(out)   :: ierr
    end subroutine f90_cl_create_buffer

    ! ----------------------------------------------------

    subroutine f90_cl_release_buffer(this, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr),            intent(inout) :: this
      integer,                intent(out)   :: ierr
    end subroutine f90_cl_release_buffer

    ! ----------------------------------------------------

    subroutine flFinish(command_queue, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr),            intent(inout) :: command_queue
      integer,                intent(out)   :: ierr
    end subroutine flFinish
    
    ! ----------------------------------------------------

    subroutine f90_cl_set_kernel_arg_buf(kernel, index, buffer, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: kernel
      integer,     intent(in)    :: index
      type(c_ptr), intent(in)    :: buffer
      integer,     intent(out)   :: ierr
    end subroutine f90_cl_set_kernel_arg_buf

    ! ----------------------------------------------------

    subroutine f90_cl_set_kernel_arg_local(kernel, index, size, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: kernel
      integer,     intent(in)    :: index
      integer,     intent(in)    :: size
      integer,     intent(out)   :: ierr
    end subroutine f90_cl_set_kernel_arg_local

    ! ----------------------------------------------------

    subroutine flEnqueueNDRangeKernel(kernel, command_queue, dim, globalsizes, localsizes, ierr)
      use c_pointer_m

      implicit none

      type(c_ptr),            intent(inout) :: kernel
      type(c_ptr),            intent(inout) :: command_queue
      integer,                intent(in)    :: dim
      integer(SIZEOF_SIZE_T), intent(in)    :: globalsizes
      integer(SIZEOF_SIZE_T), intent(in)    :: localsizes
      integer,                intent(out)   :: ierr
    end subroutine flEnqueueNDRangeKernel

  end interface

  !! Local Variables:
  !! mode: f90
  !! coding: utf-8
  !! End:

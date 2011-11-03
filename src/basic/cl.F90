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
  !! $Id: cl.F90 3587 2007-11-22 16:43:00Z xavier $

#include "global.h"
  
module cl_m
  use c_pointer_m

  implicit none 
  
  private

  public ::                          &
    flGetPlatformIDs,                &
    flEnqueueNDRangeKernel,          &
    flGetDeviceInfo,                 &
    flReleaseContext,                &
    flCreateCommandQueue,            &
    flReleaseCommandQueue,           &
    flFinish,                        &
    f90_cl_get_number_of_devices,    &
    f90_cl_init_context,             &
    f90_cl_init_device,              &
    f90_cl_create_program_from_file, &
    f90_cl_build_program,            &
    f90_cl_release_program,          &
    f90_cl_create_kernel,            &
    f90_cl_release_kernel,           &
    f90_cl_kernel_wgroup_size,       &
    f90_cl_create_buffer,            &
    f90_cl_set_kernel_arg_buf,       &
    f90_cl_set_kernel_arg_local

#include "cl_constants_inc.F90"

  interface

    integer function flGetPlatformIDs(iplatform, platform_id)
      use c_pointer_m

      implicit none
      integer,          intent(in)  :: iplatform
      type(c_ptr),      intent(out) :: platform_id
    end function flGetPlatformIDs


    ! ---------------------------------------------------

    integer function f90_cl_get_number_of_devices(platform_id)
      use c_pointer_m

      implicit none
      type(c_ptr),      intent(in)   :: platform_id
    end function f90_cl_get_number_of_devices


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

  ! ---------------------------------------------------

  interface flGetDeviceInfo

    subroutine flgetdeviceinfo_str(device, param_name, param_value)
      use c_pointer_m
      
      implicit none
      type(c_ptr),      intent(in)   :: device
      integer,          intent(in)   :: param_name
      character(len=*), intent(out)  :: param_value
    end subroutine flgetdeviceinfo_str

    subroutine flgetdeviceinfo_int(device, param_name, param_value)
      use c_pointer_m
      
      implicit none
      type(c_ptr),      intent(in)   :: device
      integer,          intent(in)   :: param_name
      integer,          intent(out)  :: param_value
    end subroutine flgetdeviceinfo_int

    subroutine flgetdeviceinfo_int64(device, param_name, param_value)
      use c_pointer_m
      
      implicit none
      type(c_ptr),      intent(in)   :: device
      integer,          intent(in)   :: param_name
      integer(8),       intent(out)  :: param_value
    end subroutine flgetdeviceinfo_int64

    module procedure flgetdeviceinfo_logical

  end interface flGetDeviceInfo
  
  ! ---------------------------------------------------

  contains

    subroutine flgetdeviceinfo_logical(device, param_name, param_value)
      type(c_ptr),      intent(in)   :: device
      integer,          intent(in)   :: param_name
      logical,          intent(out)  :: param_value

      integer(8) :: param_value_64

      call flgetdeviceinfo_int64(device, param_name, param_value_64)

      param_value = param_value_64 /= 0

    end subroutine flgetdeviceinfo_logical


end module cl_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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

#include "config_F90.h"

module cl_types
  implicit none 
  
  type :: cl_platform_id
    private 
    integer, pointer :: p 
  end type cl_platform_id
  
  type :: cl_device_id
    private 
    integer, pointer :: p 
  end type cl_device_id

  type :: cl_context
    private 
    integer, pointer :: p 
  end type cl_context

  type :: cl_command_queue
    private 
    integer, pointer :: p 
  end type cl_command_queue

  type :: cl_mem
    private 
    integer, pointer :: p 
  end type cl_mem

  type :: cl_program
    private 
    integer, pointer :: p 
  end type cl_program

  type :: cl_kernel
    private 
    integer, pointer :: p 
  end type cl_kernel

  type :: cl_event
    private 
    integer, pointer :: p 
  end type cl_event

  type :: cl_sampler
    private 
    integer, pointer :: p 
  end type cl_sampler

end module cl_types
 
module cl_m
  use cl_types

  implicit none 

  private

  ! the datatypes
  public ::                          &
    cl_platform_id,                  &
    cl_device_id,                    &
    cl_context,                      &
    cl_command_queue,                &
    cl_mem,                          &
    cl_program,                      &
    cl_kernel

  ! the functions
  public ::                          &
    flGetPlatformIDs,                &
    flGetPlatformInfo,               &
    flCreateContext,                 &
    flReleaseContext,                &
    flCreateCommandQueue,            &
    flReleaseCommandQueue,           &
    flGetDeviceInfo,                 &
    flGetDeviceIDs,                  &
    flEnqueueNDRangeKernel,          &
    flFinish,                        &
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

    subroutine f90_cl_init_context(platform_id, context)
      use cl_types

      implicit none
      type(cl_platform_id), intent(in)   :: platform_id
      type(cl_context),     intent(out)  :: context
    end subroutine f90_cl_init_context

    ! ---------------------------------------------------

    subroutine f90_cl_init_device(idevice, platform, context, device)
      use cl_types

      implicit none

      integer,              intent(in)    :: idevice
      type(cl_platform_id), intent(inout) :: platform
      type(cl_context),     intent(inout) :: context
      type(cl_device_id),   intent(out)   :: device
    end subroutine f90_cl_init_device

    ! ----------------------------------------------------

    subroutine flReleaseContext(context, status)
      use cl_types

      implicit none

      type(cl_context), intent(inout) :: context
      integer,          intent(out)   :: status
    end subroutine flReleaseContext

    ! ----------------------------------------------------

    subroutine flReleaseCommandQueue(command_queue, status)
      use cl_types

      implicit none

      type(cl_command_queue), intent(inout) :: command_queue
      integer,                intent(out)   :: status

    end subroutine flReleaseCommandQueue


    ! ----------------------------------------------------

    subroutine f90_cl_create_program_from_file(prog, context, source_file)
      use cl_types

      implicit none

      type(cl_program), intent(out)   :: prog
      type(cl_context), intent(inout) :: context
      character(len=*), intent(in)    :: source_file
    end subroutine f90_cl_create_program_from_file

    ! ----------------------------------------------------

    subroutine f90_cl_build_program(prog, context, device, flags)
      use cl_types

      implicit none

      type(cl_program),   intent(inout) :: prog
      type(cl_context),   intent(inout) :: context
      type(cl_device_id), intent(inout) :: device
      character(len=*),   intent(in)    :: flags
    end subroutine f90_cl_build_program

    ! ----------------------------------------------------

    subroutine f90_cl_release_program(prog, status)
      use cl_types

      implicit none

      type(cl_program), intent(inout) :: prog
      integer,          intent(out)   :: status
    end subroutine f90_cl_release_program

    ! ----------------------------------------------------

    subroutine f90_cl_create_kernel(kernel, prog, kernel_name, status)
      use cl_types

      implicit none

      type(cl_kernel),  intent(out)   :: kernel
      type(cl_program), intent(inout) :: prog
      character(len=*), intent(in)    :: kernel_name
      integer,          intent(out)   :: status
    end subroutine f90_cl_create_kernel

    ! ----------------------------------------------------

    subroutine f90_cl_release_kernel(kernel, status)
      use cl_types

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(out)   :: status
    end subroutine f90_cl_release_kernel

    ! ----------------------------------------------------

    integer function f90_cl_kernel_wgroup_size(kernel, device)
      use cl_types

      implicit none

      type(cl_kernel),    intent(inout) :: kernel
      type(cl_device_id), intent(inout) :: device
    end function f90_cl_kernel_wgroup_size

    ! ----------------------------------------------------

    subroutine f90_cl_create_buffer(this, context, flags, size, status)
      use cl_types

      implicit none

      type(cl_mem),           intent(inout) :: this
      type(cl_context),       intent(inout) :: context
      integer,                intent(in)    :: flags
      integer(SIZEOF_SIZE_T), intent(in)    :: size
      integer,                intent(out)   :: status
    end subroutine f90_cl_create_buffer

    ! ----------------------------------------------------

    subroutine f90_cl_release_buffer(this, status)
      use cl_types

      implicit none

      type(cl_mem),           intent(inout) :: this
      integer,                intent(out)   :: status
    end subroutine f90_cl_release_buffer

    ! ----------------------------------------------------

    subroutine flFinish(command_queue, status)
      use cl_types

      implicit none

      type(cl_command_queue), intent(inout) :: command_queue
      integer,                intent(out)   :: status
    end subroutine flFinish

    ! ----------------------------------------------------

    subroutine f90_cl_set_kernel_arg_buf(kernel, index, buffer, status)
      use cl_types

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(in)    :: index
      type(cl_mem),    intent(in)    :: buffer
      integer,         intent(out)   :: status
    end subroutine f90_cl_set_kernel_arg_buf

    ! ----------------------------------------------------

    subroutine f90_cl_set_kernel_arg_local(kernel, index, size, status)
      use cl_types

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(in)    :: index
      integer,         intent(in)    :: size
      integer,         intent(out)   :: status
    end subroutine f90_cl_set_kernel_arg_local

    ! ----------------------------------------------------

    subroutine flEnqueueNDRangeKernel(command_queue, kernel, dim, globalsizes, localsizes, status)
      use cl_types

      implicit none

      type(cl_command_queue), intent(inout) :: command_queue
      type(cl_kernel),        intent(inout) :: kernel
      integer,                intent(in)    :: dim
      integer(SIZEOF_SIZE_T), intent(in)    :: globalsizes
      integer(SIZEOF_SIZE_T), intent(in)    :: localsizes
      integer,                intent(out)   :: status
    end subroutine flEnqueueNDRangeKernel

  end interface

  ! ---------------------------------------------------

  interface flGetPlatformIDs

    subroutine flgetplatformids_num(num_platforms, status)
      use cl_types

      implicit none
      integer,              intent(out)  :: num_platforms
      integer,              intent(out)  :: status
    end subroutine flgetplatformids_num

    module procedure flgetplatformids_list

  end interface flGetPlatformIDs

  ! ---------------------------------------------------

  interface

    subroutine flGetPlatformInfo(platform, param_name, param_value, status)
      use cl_types

      implicit none
      type(cl_platform_id), intent(in)   :: platform
      integer,              intent(in)   :: param_name
      character(len=*),     intent(out)  :: param_value
      integer,              intent(out)  :: status
    end subroutine flGetPlatformInfo

  end interface

  ! ---------------------------------------------------

  interface flGetDeviceIDs

    subroutine flgetdeviceids_num(platform, device_type, num_devices, status)
      use cl_types

      implicit none
      type(cl_platform_id), intent(in)   :: platform
      integer,              intent(in)   :: device_type
      integer,              intent(out)  :: num_devices
      integer,              intent(out)  :: status
    end subroutine flgetdeviceids_num

    module procedure flgetdeviceids_list

  end interface flGetDeviceIDs

  ! ---------------------------------------------------

  interface flGetDeviceInfo

    subroutine flgetdeviceinfo_str(device, param_name, param_value, status)
      use cl_types

      implicit none
      type(cl_device_id), intent(in)   :: device
      integer,            intent(in)   :: param_name
      character(len=*),   intent(out)  :: param_value
      integer,            intent(out)  :: status
    end subroutine flgetdeviceinfo_str

    subroutine flgetdeviceinfo_int(device, param_name, param_value, status)
      use cl_types

      implicit none
      type(cl_device_id), intent(in)   :: device
      integer,            intent(in)   :: param_name
      integer,            intent(out)  :: param_value
      integer,            intent(out)  :: status
    end subroutine flgetdeviceinfo_int

    subroutine flgetdeviceinfo_int64(device, param_name, param_value, status)
      use cl_types

      implicit none
      type(cl_device_id), intent(in)   :: device
      integer,            intent(in)   :: param_name
      integer(8),         intent(out)  :: param_value
      integer,            intent(out)  :: status
    end subroutine flgetdeviceinfo_int64

    module procedure flgetdeviceinfo_logical

  end interface flGetDeviceInfo

  ! ---------------------------------------------------

contains

  subroutine flgetplatformids_list(num_entries, platforms, num_platforms, status)
    integer,              intent(out)  :: num_entries
    type(cl_platform_id), intent(out)  :: platforms(:)
    integer,              intent(out)  :: num_platforms
    integer,              intent(out)  :: status

#ifdef HAVE_OPENCL
    integer                         :: iplatform
    type(cl_platform_id), allocatable :: plat(:)

    interface
      subroutine flgetplatformids_listall(num_entries, platforms, num_platforms, status)
        use cl_types

        implicit none

        integer,              intent(out)  :: num_entries
        type(cl_platform_id), intent(out)  :: platforms
        integer,              intent(out)  :: num_platforms
        integer,              intent(out)  :: status
      end subroutine flgetplatformids_listall

      subroutine flgetplatformids_getplat(allplatforms, iplatform, platform)
        use cl_types

        implicit none

        type(cl_platform_id), intent(in)   :: allplatforms
        integer,              intent(in)   :: iplatform
        type(cl_platform_id), intent(out)  :: platform
      end subroutine flgetplatformids_getplat
    end interface

    ! since our cl_platform_id type might be longer than the C
    ! cl_platform_id type we need to get all the values in an array
    ! and the copy them explicitly to the return array

    allocate(plat(1:num_entries))

    call flgetplatformids_listall(num_entries, plat(1), num_platforms, status)

    do iplatform = 1, num_platforms
      call flgetplatformids_getplat(plat(1), iplatform - 1, platforms(iplatform))
    end do

    deallocate(plat)
#endif
  end subroutine flgetplatformids_list

  ! ----------------------------------------------------------

  subroutine flgetdeviceids_list(platform, device_type, num_entries, devices, num_devices, status)
    type(cl_platform_id), intent(in)   :: platform
    integer,              intent(in)   :: device_type
    integer,              intent(out)  :: num_entries
    type(cl_device_id),   intent(out)  :: devices(:)
    integer,              intent(out)  :: num_devices
    integer,              intent(out)  :: status

#ifdef HAVE_OPENCL
    integer                         :: idevice
    type(cl_device_id), allocatable :: dev(:)

    interface
      subroutine flgetdeviceids_listall(platform, device_type, num_entries, devices, num_devices, status)
        use cl_types

        implicit none

        type(cl_platform_id), intent(in)   :: platform
        integer,              intent(in)   :: device_type
        integer,              intent(out)  :: num_entries
        type(cl_device_id),   intent(out)  :: devices
        integer,              intent(out)  :: num_devices
        integer,              intent(out)  :: status
      end subroutine flgetdeviceids_listall

      subroutine flgetdeviceids_getdev(alldevices, idevice, device)
        use cl_types

        implicit none

        type(cl_device_id),   intent(in)   :: alldevices
        integer,              intent(in)   :: idevice
        type(cl_device_id),   intent(out)  :: device
      end subroutine flgetdeviceids_getdev
    end interface

    ! since our cl_device_id type might be longer than the C
    ! cl_device_id type we need to get all the values in an array
    ! and the copy them explicitly to the return array

    allocate(dev(1:num_entries))

    call flgetdeviceids_listall(platform, device_type, num_entries, dev(1), num_devices, status)

    do idevice = 1, num_devices
      call flgetdeviceids_getdev(dev(1), idevice - 1, devices(idevice))
    end do

    deallocate(dev)
#endif
  end subroutine flgetdeviceids_list

  ! ----------------------------------------------------------

  subroutine flgetdeviceinfo_logical(device, param_name, param_value, status)
    type(cl_device_id), intent(in)   :: device
    integer,            intent(in)   :: param_name
    logical,            intent(out)  :: param_value
    integer,            intent(out)  :: status

    integer(8) :: param_value_64

#ifdef HAVE_OPENCL
    call flgetdeviceinfo_int64(device, param_name, param_value_64, status)
#endif

    param_value = param_value_64 /= 0

  end subroutine flgetdeviceinfo_logical

  ! ---------------------------------------------------------

  type(cl_context) function flCreateContext(platform, num_devices, devices, errcode_ret) result(context)
    type(cl_platform_id), intent(in)   :: platform
    integer,              intent(in)   :: num_devices
    type(cl_device_id),   intent(in)   :: devices(:)
    integer,              intent(out)  :: errcode_ret

    interface
      subroutine flcreatecontext_low(platform, num_devices, devices, errcode_ret, context)
        use cl_types

        implicit none

        type(cl_platform_id), intent(in)   :: platform
        integer,              intent(in)   :: num_devices
        type(cl_device_id),   intent(in)   :: devices
        integer,              intent(out)  :: errcode_ret
        type(cl_context),     intent(out)  :: context
      end subroutine flcreatecontext_low

      subroutine flgetdeviceids_setdev(alldevices, idevice, device)
        use cl_types

        implicit none

        type(cl_device_id),   intent(out)   :: alldevices
        integer,              intent(in)   :: idevice
        type(cl_device_id),   intent(in)  :: device
      end subroutine flgetdeviceids_setdev
    end interface

    integer :: idev
    type(cl_device_id), allocatable :: devs(:)
    
    allocate(devs(1:num_devices))

    do idev = 1, num_devices
#ifdef HAVE_OPENCL
      call flgetdeviceids_setdev(devs(1), idev - 1, devices(idev))
#endif
    end do
    
#ifdef HAVE_OPENCL
    call flcreatecontext_low(platform, num_devices, devs(1), errcode_ret, context)
#endif

    deallocate(devs)

  end function flCreateContext

  type(cl_command_queue) function flCreateCommandQueue(context, device, properties, errcode_ret) result(command_queue)
      type(cl_context),       intent(inout) :: context
      type(cl_device_id),     intent(inout) :: device
      integer,                intent(in)    :: properties
      integer,                intent(out)   :: errcode_ret

    interface
      subroutine flcreatecommandqueue_low(context, device, properties, errcode_ret, command_queue)
        use cl_types
        
        implicit none
        
        type(cl_context),       intent(inout) :: context
        type(cl_device_id),     intent(inout) :: device
        integer,                intent(in)    :: properties
        integer,                intent(out)   :: errcode_ret
        type(cl_command_queue), intent(inout) :: command_queue
      end subroutine flcreatecommandqueue_low
    end interface

#ifdef HAVE_OPENCL
    call flcreatecommandqueue_low(context, device, properties, errcode_ret, command_queue)
#endif

  end function flCreateCommandQueue

end module cl_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

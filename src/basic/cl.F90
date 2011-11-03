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
    flEnqueueNDRangeKernel,          &
    flGetDeviceInfo,                 &
    flGetDeviceIDs,                  &
    flReleaseContext,                &
    flCreateCommandQueue,            &
    flReleaseCommandQueue,           &
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

    integer function flGetPlatformIDs(iplatform, platform_id)
      use cl_types

      implicit none
      integer,              intent(in)  :: iplatform
      type(cl_platform_id), intent(out) :: platform_id
    end function flGetPlatformIDs

    ! ---------------------------------------------------

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

    subroutine flReleaseContext(context)
      use cl_types

      implicit none

      type(cl_context), intent(inout) :: context
    end subroutine flReleaseContext

    ! ----------------------------------------------------

    subroutine flCreateCommandQueue(command_queue, context, device, ierr)
      use cl_types

      implicit none

      type(cl_command_queue), intent(inout) :: command_queue
      type(cl_context),       intent(inout) :: context
      type(cl_device_id),     intent(inout) :: device
      integer,                intent(out)   :: ierr

    end subroutine flCreateCommandQueue

    ! ----------------------------------------------------

    subroutine flReleaseCommandQueue(command_queue, ierr)
      use cl_types

      implicit none

      type(cl_command_queue), intent(inout) :: command_queue
      integer,                intent(out)   :: ierr

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

    subroutine f90_cl_release_program(prog, ierr)
      use cl_types

      implicit none

      type(cl_program), intent(inout) :: prog
      integer,          intent(out)   :: ierr
    end subroutine f90_cl_release_program

    ! ----------------------------------------------------

    subroutine f90_cl_create_kernel(kernel, prog, kernel_name, ierr)
      use cl_types

      implicit none

      type(cl_kernel),  intent(out)   :: kernel
      type(cl_program), intent(inout) :: prog
      character(len=*), intent(in)    :: kernel_name
      integer,          intent(out)   :: ierr
    end subroutine f90_cl_create_kernel

    ! ----------------------------------------------------

    subroutine f90_cl_release_kernel(kernel, ierr)
      use cl_types

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(out)   :: ierr
    end subroutine f90_cl_release_kernel

    ! ----------------------------------------------------

    integer function f90_cl_kernel_wgroup_size(kernel, device)
      use cl_types

      implicit none

      type(cl_kernel),    intent(inout) :: kernel
      type(cl_device_id), intent(inout) :: device
    end function f90_cl_kernel_wgroup_size

    ! ----------------------------------------------------

    subroutine f90_cl_create_buffer(this, context, flags, size, ierr)
      use cl_types

      implicit none

      type(cl_mem),           intent(inout) :: this
      type(cl_context),       intent(inout) :: context
      integer,                intent(in)    :: flags
      integer(SIZEOF_SIZE_T), intent(in)    :: size
      integer,                intent(out)   :: ierr
    end subroutine f90_cl_create_buffer

    ! ----------------------------------------------------

    subroutine f90_cl_release_buffer(this, ierr)
      use cl_types

      implicit none

      type(cl_mem),           intent(inout) :: this
      integer,                intent(out)   :: ierr
    end subroutine f90_cl_release_buffer

    ! ----------------------------------------------------

    subroutine flFinish(command_queue, ierr)
      use cl_types

      implicit none

      type(cl_command_queue), intent(inout) :: command_queue
      integer,                intent(out)   :: ierr
    end subroutine flFinish

    ! ----------------------------------------------------

    subroutine f90_cl_set_kernel_arg_buf(kernel, index, buffer, ierr)
      use cl_types

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(in)    :: index
      type(cl_mem),    intent(in)    :: buffer
      integer,         intent(out)   :: ierr
    end subroutine f90_cl_set_kernel_arg_buf

    ! ----------------------------------------------------

    subroutine f90_cl_set_kernel_arg_local(kernel, index, size, ierr)
      use cl_types

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(in)    :: index
      integer,         intent(in)    :: size
      integer,         intent(out)   :: ierr
    end subroutine f90_cl_set_kernel_arg_local

    ! ----------------------------------------------------

    subroutine flEnqueueNDRangeKernel(kernel, command_queue, dim, globalsizes, localsizes, ierr)
      use cl_types

      implicit none

      type(cl_kernel),        intent(inout) :: kernel
      type(cl_command_queue), intent(inout) :: command_queue
      integer,                intent(in)    :: dim
      integer(SIZEOF_SIZE_T), intent(in)    :: globalsizes
      integer(SIZEOF_SIZE_T), intent(in)    :: localsizes
      integer,                intent(out)   :: ierr
    end subroutine flEnqueueNDRangeKernel

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


end module cl_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

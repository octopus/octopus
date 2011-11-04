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

module cl_device_m
  use cl_types_m

  implicit none 

  private

  public ::                          &
    flGetDeviceInfo,                 &
    flGetDeviceIDs

  ! ---------------------------------------------------

  interface flGetDeviceIDs

    subroutine flgetdeviceids_num(platform, device_type, num_devices, status)
      use cl_types_m

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
      use cl_types_m

      implicit none
      type(cl_device_id), intent(in)   :: device
      integer,            intent(in)   :: param_name
      character(len=*),   intent(out)  :: param_value
      integer,            intent(out)  :: status
    end subroutine flgetdeviceinfo_str

    subroutine flgetdeviceinfo_int(device, param_name, param_value, status)
      use cl_types_m

      implicit none
      type(cl_device_id), intent(in)   :: device
      integer,            intent(in)   :: param_name
      integer,            intent(out)  :: param_value
      integer,            intent(out)  :: status
    end subroutine flgetdeviceinfo_int

    subroutine flgetdeviceinfo_int64(device, param_name, param_value, status)
      use cl_types_m

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
        use cl_types_m

        implicit none

        type(cl_platform_id), intent(in)   :: platform
        integer,              intent(in)   :: device_type
        integer,              intent(out)  :: num_entries
        type(cl_device_id),   intent(out)  :: devices
        integer,              intent(out)  :: num_devices
        integer,              intent(out)  :: status
      end subroutine flgetdeviceids_listall

      subroutine flgetdeviceids_getdev(alldevices, idevice, device)
        use cl_types_m

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

end module cl_device_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

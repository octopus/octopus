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
 
module cl_command_queue_m
  use cl_types_m

  implicit none 

  private

  ! the functions
  public ::                          &
    flCreateCommandQueue,            &
    flReleaseCommandQueue,           &
    flEnqueueNDRangeKernel,          &
    flFinish

  interface

    ! ----------------------------------------------------

    subroutine flReleaseCommandQueue(command_queue, status)
      use cl_types_m

      implicit none

      type(cl_command_queue), intent(inout) :: command_queue
      integer,                intent(out)   :: status

    end subroutine flReleaseCommandQueue

    ! ----------------------------------------------------

    subroutine flFinish(command_queue, status)
      use cl_types_m

      implicit none

      type(cl_command_queue), intent(inout) :: command_queue
      integer,                intent(out)   :: status
    end subroutine flFinish

    ! ----------------------------------------------------

    subroutine flEnqueueNDRangeKernel(command_queue, kernel, dim, globalsizes, localsizes, status)
      use cl_types_m

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

contains

  ! --------------------------------------------------------

  type(cl_command_queue) function flCreateCommandQueue(context, device, properties, errcode_ret) result(command_queue)
      type(cl_context),       intent(inout) :: context
      type(cl_device_id),     intent(inout) :: device
      integer,                intent(in)    :: properties
      integer,                intent(out)   :: errcode_ret

    interface
      subroutine flcreatecommandqueue_low(context, device, properties, errcode_ret, command_queue)
        use cl_types_m
        
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

end module cl_command_queue_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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

module cl_kernel_m
  use cl_types_m

  implicit none 

  private

  public ::                          &
    clCreateKernel,                  &
    clReleaseKernel,                 &
    clKernelWorkGroupInfo,           &
    clSetKernelArgLocal,             &
    clSetKernelArg

  interface

    subroutine clReleaseKernel(kernel, status)
      use cl_types_m

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(out)   :: status
    end subroutine clReleaseKernel

  end interface


  interface clSetKernelArg

    subroutine clsetkernelarg_buf(kernel, arg_index, arg_value, status)
      use cl_types_m

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(in)    :: arg_index
      type(cl_mem),    intent(in)    :: arg_value
      integer,         intent(out)   :: status
    end subroutine clsetkernelarg_buf

    ! ----------------------------------------------------

    subroutine clsetkernelarg_int(kernel, arg_index, arg_value, status)
      use cl_types_m

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(in)    :: arg_index
      integer(4),      intent(in)    :: arg_value
      integer,         intent(out)   :: status
    end subroutine clsetkernelarg_int

    ! ----------------------------------------------------

    subroutine clsetkernelarg_float(kernel, arg_index, arg_value, status)
      use cl_types_m

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(in)    :: arg_index
      real(4),         intent(in)    :: arg_value
      integer,         intent(out)   :: status
    end subroutine clsetkernelarg_float

    ! ----------------------------------------------------

    subroutine clsetkernelarg_double(kernel, arg_index, arg_value, status)
      use cl_types_m

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(in)    :: arg_index
      real(8),         intent(in)    :: arg_value
      integer,         intent(out)   :: status
    end subroutine clsetkernelarg_double

  end interface clSetKernelArg

  ! ----------------------------------------------------

  interface
    subroutine clSetKernelArgLocal(kernel, arg_index, arg_size, status)
      use cl_types_m

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(in)    :: arg_index
      integer(8),      intent(in)    :: arg_size
      integer,         intent(out)   :: status
    end subroutine clSetKernelArgLocal
  end interface

  interface clKernelWorkGroupInfo

    subroutine clkernelworkgroupinfo_int64(kernel, device, param_name, param_value, retcode_err)
      use cl_types_m
      
      implicit none

      type(cl_kernel),    intent(in)    :: kernel
      type(cl_device_id), intent(in)    :: device
      integer,            intent(in)    :: param_name
      integer(8),         intent(out)   :: param_value
      integer,            intent(out)   :: retcode_err
    end subroutine clkernelworkgroupinfo_int64

  end interface clKernelWorkGroupInfo

contains

  type(cl_kernel) function clCreateKernel(program, kernel_name, errcode_ret) result(kernel)
    type(cl_program), intent(inout) :: program
    character(len=*), intent(in)    :: kernel_name
    integer,          intent(out)   :: errcode_ret

    interface
      subroutine clcreatekernel_low(program, kernel_name, errcode_ret, kernel)
        use cl_types_m

        implicit none

        type(cl_program), intent(inout) :: program
        character(len=*), intent(in)    :: kernel_name
        integer,          intent(out)   :: errcode_ret
        type(cl_kernel),  intent(out)   :: kernel
      end subroutine clcreatekernel_low
    end interface

#ifdef HAVE_OPENCL
    call clcreatekernel_low(program, kernel_name, errcode_ret, kernel)
#endif

  end function clCreateKernel

end module cl_kernel_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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
    f90_cl_create_kernel,            &
    clReleaseKernel,                 &
    f90_cl_kernel_wgroup_size,       &
    f90_cl_set_kernel_arg_buf,       &
    f90_cl_set_kernel_arg_local

  interface

    ! ----------------------------------------------------

    subroutine f90_cl_create_kernel(kernel, prog, kernel_name, status)
      use cl_types_m

      implicit none

      type(cl_kernel),  intent(out)   :: kernel
      type(cl_program), intent(inout) :: prog
      character(len=*), intent(in)    :: kernel_name
      integer,          intent(out)   :: status
    end subroutine f90_cl_create_kernel

    ! ----------------------------------------------------

    subroutine clReleaseKernel(kernel, status)
      use cl_types_m
      
      implicit none
      
      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(out)   :: status
    end subroutine clReleaseKernel

    ! ----------------------------------------------------

    integer function f90_cl_kernel_wgroup_size(kernel, device)
      use cl_types_m

      implicit none

      type(cl_kernel),    intent(inout) :: kernel
      type(cl_device_id), intent(inout) :: device
    end function f90_cl_kernel_wgroup_size

    ! ----------------------------------------------------

    subroutine f90_cl_set_kernel_arg_buf(kernel, index, buffer, status)
      use cl_types_m

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(in)    :: index
      type(cl_mem),    intent(in)    :: buffer
      integer,         intent(out)   :: status
    end subroutine f90_cl_set_kernel_arg_buf

    ! ----------------------------------------------------

    subroutine f90_cl_set_kernel_arg_local(kernel, index, size, status)
      use cl_types_m

      implicit none

      type(cl_kernel), intent(inout) :: kernel
      integer,         intent(in)    :: index
      integer,         intent(in)    :: size
      integer,         intent(out)   :: status
    end subroutine f90_cl_set_kernel_arg_local

  end interface

end module cl_kernel_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

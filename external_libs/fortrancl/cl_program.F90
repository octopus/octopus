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
 
module cl_program_m
  use cl_types_m

  implicit none 

  private

  public ::                          &
    f90_cl_create_program_from_file, &
    clBuildProgram,                  &
    clReleaseProgram,                &
    clGetProgramBuildInfo

  interface

    ! ----------------------------------------------------

    subroutine f90_cl_create_program_from_file(prog, context, source_file)
      use cl_types_m

      implicit none

      type(cl_program), intent(out)   :: prog
      type(cl_context), intent(inout) :: context
      character(len=*), intent(in)    :: source_file
    end subroutine f90_cl_create_program_from_file

    ! ----------------------------------------------------

    subroutine clReleaseProgram(prog, status)
      use cl_types_m

      implicit none

      type(cl_program), intent(inout) :: prog
      integer,          intent(out)   :: status
    end subroutine clReleaseProgram

  end interface

  interface clBuildProgram

    subroutine clBuildProgram_nodevices(program, options, retcode_err)
      use cl_types_m

      implicit none

      type(cl_program),   intent(inout) :: program
      character(len=*),   intent(in)    :: options
      integer,            intent(in)    :: retcode_err
    end subroutine clBuildProgram_nodevices

  end interface clBuildProgram

  interface clGetProgramBuildInfo
    
    subroutine clGetProgramBuildInfo_str(program, device, param_name, param_value, retcode_err)
      use cl_types_m

      implicit none

      type(cl_program),   intent(in)    :: program
      type(cl_device_id), intent(in)    :: device
      integer,            intent(in)    :: param_name
      character(len=*),   intent(out)   :: param_value
      integer,            intent(out)   :: retcode_err
    end subroutine clGetProgramBuildInfo_str

  end interface clGetProgramBuildInfo

  ! ----------------------------------------------------

end module cl_program_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

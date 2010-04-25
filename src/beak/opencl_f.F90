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

#include "global.h"

module opencl_m
  use c_pointer_m
  use global_m

  implicit none 

  private

  public ::             &
    opencl_t,           &
    opencl_init,        &
    opencl_end

  type opencl_t 
    type(c_ptr) :: env
    type(c_ptr) :: kernel_vpsi
  end type opencl_t

  type(opencl_t), public :: opencl

  interface
    subroutine f90_opencl_env_init(this, source_path)
      use c_pointer_m

      type(c_ptr),      intent(out) :: this
      character(len=*), intent(in)  :: source_path
    end subroutine f90_opencl_env_init

    subroutine f90_opencl_env_end(this)
      use c_pointer_m

      type(c_ptr), intent(inout) :: this
    end subroutine f90_opencl_env_end

    subroutine f90_opencl_kernel_init(this, env, source_file, kernel_base_name)
      use c_pointer_m

      type(c_ptr),      intent(out)   :: this
      type(c_ptr),      intent(inout) :: env
      character(len=*), intent(in)    :: source_file
      character(len=*), intent(in)    :: kernel_base_name
    end subroutine f90_opencl_kernel_init

    subroutine f90_opencl_kernel_end(this)
      use c_pointer_m

      type(c_ptr), intent(inout) :: this
    end subroutine f90_opencl_kernel_end
  end interface

  contains
    
    subroutine opencl_init(this)
      type(opencl_t), intent(out) :: this

      call f90_opencl_env_init(this%env, trim(conf%share)//'/opencl/')      

      ! now initialize the kernels
      call f90_opencl_kernel_init(this%kernel_vpsi, this%env, "vpsi.cl", "vpsi");
    end subroutine opencl_init

    ! ------------------------------------------
    subroutine opencl_end(this)
      type(opencl_t), intent(inout) :: this
      
      call f90_opencl_kernel_end(this%kernel_vpsi);
      call f90_opencl_env_end(this%env)
    end subroutine opencl_end

end module opencl_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

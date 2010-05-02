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
  use datasets_m
  use global_m
  use messages_m
  use types_m
  use parser_m
  use profiling_m

  implicit none 

  private

  public ::                       &
    opencl_is_enabled,            &
    opencl_init,                  &
    opencl_end,                   &
    opencl_mem_t,                 &
    opencl_create_buffer,         &
    opencl_write_buffer,          &
    opencl_read_buffer,           &
    opencl_release_buffer,        &
    opencl_padded_size,           &
    opencl_finish,                &
    opencl_set_kernel_arg,        &
    opencl_max_workgroup_size,    &
    opencl_kernel_workgroup_size, &
    opencl_kernel_run

  type opencl_t 
    type(c_ptr) :: env
    integer     :: max_workgroup_size
    logical     :: enabled
  end type opencl_t

  type opencl_mem_t
    private
    type(c_ptr)            :: mem
    integer(SIZEOF_SIZE_T) :: size
    integer                :: type
  end type opencl_mem_t

  type(opencl_t) :: opencl

  ! the kernels
  type(c_ptr), public :: dvpsi
  type(c_ptr), public :: zvpsi
  type(c_ptr), public :: zvpsi_spinors
  type(c_ptr), public :: set_zero
  type(c_ptr), public :: dset_zero_part
  type(c_ptr), public :: zset_zero_part
  type(c_ptr), public :: doperate
  type(c_ptr), public :: zoperate

  ! this values are copied from OpenCL include CL/cl.h
  integer, parameter, public ::        &
    CL_MEM_READ_WRITE = 1,             &
    CL_MEM_WRITE_ONLY = 2,             &
    CL_MEM_READ_ONLY  = 4

  ! this function are defined in opencl_low.c
  interface
    ! ---------------------------------------------------

    subroutine f90_opencl_env_init(env, source_path)
      use c_pointer_m

      implicit none

      type(c_ptr),      intent(out) :: env
      character(len=*), intent(in)  :: source_path
    end subroutine f90_opencl_env_init

    ! ----------------------------------------------------

    subroutine f90_opencl_env_end(env)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: env
    end subroutine f90_opencl_env_end

    ! ----------------------------------------------------

    integer function f90_opencl_max_workgroup_size(env)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: env
    end function f90_opencl_max_workgroup_size

    ! ----------------------------------------------------

    subroutine f90_opencl_build_program(prog, env, source_file)
      use c_pointer_m

      implicit none

      type(c_ptr),      intent(out)   :: prog
      type(c_ptr),      intent(inout) :: env
      character(len=*), intent(in)    :: source_file
    end subroutine f90_opencl_build_program

    ! ----------------------------------------------------

    subroutine f90_opencl_release_program(prog)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: prog
    end subroutine f90_opencl_release_program

    ! ----------------------------------------------------

    subroutine f90_opencl_create_kernel(kernel, prog, kernel_name)
      use c_pointer_m

      implicit none

      type(c_ptr),      intent(out)   :: kernel
      type(c_ptr),      intent(inout) :: prog
      character(len=*), intent(in)    :: kernel_name
    end subroutine f90_opencl_create_kernel

    ! ----------------------------------------------------

    subroutine f90_opencl_release_kernel(kernel)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: kernel
    end subroutine f90_opencl_release_kernel

    ! ----------------------------------------------------

    integer function f90_opencl_kernel_wgroup_size(kernel, env)
      use c_pointer_m
      
      implicit none
      
      type(c_ptr), intent(inout) :: kernel
      type(c_ptr), intent(inout) :: env
    end function f90_opencl_kernel_wgroup_size

    ! ----------------------------------------------------

    subroutine f90_opencl_create_buffer(this, env, flags, size)
      use c_pointer_m

      implicit none

      type(c_ptr),            intent(inout) :: this
      type(c_ptr),            intent(inout) :: env
      integer,                intent(in)    :: flags
      integer(SIZEOF_SIZE_T), intent(in)    :: size
    end subroutine f90_opencl_create_buffer
    
    ! ----------------------------------------------------

    subroutine f90_opencl_release_buffer(this)
      use c_pointer_m

      implicit none

      type(c_ptr),            intent(inout) :: this
    end subroutine f90_opencl_release_buffer

    ! ----------------------------------------------------

    subroutine f90_opencl_finish(this)
      use c_pointer_m

      implicit none

      type(c_ptr),            intent(inout) :: this
    end subroutine f90_opencl_finish

    ! ----------------------------------------------------

    subroutine f90_opencl_set_kernel_arg_buf(kernel, index, buffer)
      use c_pointer_m

      implicit none

      type(c_ptr), intent(inout) :: kernel
      integer,     intent(in)    :: index
      type(c_ptr), intent(in)    :: buffer
    end subroutine f90_opencl_set_kernel_arg_buf
    
    ! ----------------------------------------------------

    subroutine f90_opencl_kernel_run(kernel, env, dim, globalsizes, localsizes)
      use c_pointer_m

      implicit none

      type(c_ptr),            intent(inout) :: kernel
      type(c_ptr),            intent(inout) :: env
      integer,                intent(in)    :: dim
      integer(SIZEOF_SIZE_T), intent(in)    :: globalsizes
      integer(SIZEOF_SIZE_T), intent(in)    :: localsizes
    end subroutine f90_opencl_kernel_run

  end interface

  interface opencl_create_buffer
    module procedure opencl_create_buffer_4, opencl_create_buffer_8
  end interface

  interface opencl_write_buffer
    module procedure iopencl_write_buffer, dopencl_write_buffer, zopencl_write_buffer, &
      iopencl_write_buffer_2, dopencl_write_buffer_2, zopencl_write_buffer_2
  end interface

  interface opencl_read_buffer
    module procedure iopencl_read_buffer, dopencl_read_buffer, zopencl_read_buffer
  end interface

  interface opencl_set_kernel_arg
    module procedure                 &
      opencl_set_kernel_arg_buffer,  &
      iopencl_set_kernel_arg_data,   &
      dopencl_set_kernel_arg_data,   &
      zopencl_set_kernel_arg_data
  end interface
  
  type(profile_t), save :: prof_read, prof_write, prof_kernel_run
  
  contains

    pure logical function opencl_is_enabled() result(enabled)
      
      enabled = opencl%enabled
    end function opencl_is_enabled

    ! ------------------------------------------

    subroutine opencl_init()
      type(c_ptr) :: prog
      logical  :: disable

      call push_sub('opencl.opencl_init')
      
      !%Variable DisableOpenCL
      !%Type logical
      !%Default yes
      !%Section Execution::OpenCL
      !%Description
      !% If Octopus was compiled with OpenCL support, it will try to
      !% initialize and use an OpenCL device. By setting this variable
      !% to <tt>yes</tt> you tell Octopus not to use OpenCL.
      !%End
      call parse_logical(datasets_check('DisableOpenCL'), .false., disable)
      opencl%enabled = .not. disable
      
      if(disable) then
        call pop_sub('opencl.opencl_init')
        return
      end if

      call f90_opencl_env_init(opencl%env, trim(conf%share)//'/opencl/')   
      
      opencl%max_workgroup_size = f90_opencl_max_workgroup_size(opencl%env)
      
      ! now initialize the kernels
      call f90_opencl_build_program(prog, opencl%env, "vpsi.cl")
      call f90_opencl_create_kernel(dvpsi, prog, "dvpsi")
      call f90_opencl_create_kernel(zvpsi, prog, "zvpsi")
      call f90_opencl_create_kernel(zvpsi_spinors, prog, "zvpsi_spinors")
      call f90_opencl_release_program(prog)
      
      call f90_opencl_build_program(prog, opencl%env, "set_zero.cl")
      call f90_opencl_create_kernel(set_zero, prog, "set_zero")
      call f90_opencl_create_kernel(dset_zero_part, prog, "dset_zero_part")
      call f90_opencl_create_kernel(zset_zero_part, prog, "zset_zero_part")
      call f90_opencl_release_program(prog)
      
      call f90_opencl_build_program(prog, opencl%env, "operate.cl")
      call f90_opencl_create_kernel(doperate, prog, "doperate")
      call f90_opencl_create_kernel(zoperate, prog, "zoperate")
      call f90_opencl_release_program(prog)

      call pop_sub('opencl.opencl_init')
    end subroutine opencl_init

    ! ------------------------------------------

    subroutine opencl_end()

      call push_sub('opencl.opencl_end')

      if(opencl_is_enabled()) then
        call f90_opencl_release_kernel(dvpsi)
        call f90_opencl_release_kernel(zvpsi)
        call f90_opencl_env_end(opencl%env)
      end if

      call pop_sub('opencl.opencl_end')
    end subroutine opencl_end

    ! ------------------------------------------

    subroutine opencl_create_buffer_4(this, flags, type, size)
      type(opencl_mem_t), intent(inout) :: this
      integer,            intent(in)    :: flags
      integer,            intent(in)    :: type
      integer,            intent(in)    :: size
      
      integer(SIZEOF_SIZE_T) :: fsize

      call push_sub('opencl.opencl_create_buffer_4')

      this%type = type
      this%size = size      
      fsize = size*types_get_size(type)
      
      call f90_opencl_create_buffer(this%mem, opencl%env, flags, fsize)

      call pop_sub('opencl.opencl_create_buffer_4')
    end subroutine opencl_create_buffer_4

    ! ------------------------------------------

    subroutine opencl_create_buffer_8(this, flags, type, size)
      type(opencl_mem_t),     intent(inout) :: this
      integer,                intent(in)    :: flags
      integer,                intent(in)    :: type
      integer(SIZEOF_SIZE_T), intent(in)    :: size
      
      integer(SIZEOF_SIZE_T) :: fsize
      
      call push_sub('opencl.opencl_create_buffer_8')

      this%type = type
      this%size = size

      fsize = size*types_get_size(type)

      call f90_opencl_create_buffer(this%mem, opencl%env, flags, fsize)
      
      call pop_sub('opencl.opencl_create_buffer_8')
    end subroutine opencl_create_buffer_8

    ! ------------------------------------------

    subroutine opencl_release_buffer(this)
      type(opencl_mem_t), intent(inout) :: this

      this%size = 0
      call f90_opencl_release_buffer(this%mem)

    end subroutine opencl_release_buffer

    ! ------------------------------------------

    integer(SIZEOF_SIZE_T) pure function opencl_get_buffer_size(this) result(size)
      type(opencl_mem_t), intent(in) :: this

      size = this%size
    end function opencl_get_buffer_size

    ! -----------------------------------------

    integer pure function opencl_get_buffer_type(this) result(type)
      type(opencl_mem_t), intent(in) :: this

      type = this%type
    end function opencl_get_buffer_type

    ! -----------------------------------------

    integer(SIZEOF_SIZE_T) function opencl_padded_size(nn) result(psize)
      integer,        intent(in) :: nn

      integer :: modnn, bsize

      bsize = opencl_max_workgroup_size()

      psize = nn
      modnn = mod(nn, bsize)
      if(modnn /= 0) psize = psize + bsize - modnn

    end function opencl_padded_size

    ! ------------------------------------------

    subroutine opencl_finish()
      call f90_opencl_finish(opencl%env)
    end subroutine opencl_finish

    ! ------------------------------------------

    subroutine opencl_set_kernel_arg_buffer(kernel, narg, buffer)
      type(c_ptr),        intent(inout) :: kernel
      integer,            intent(in)    :: narg
      type(opencl_mem_t), intent(in)    :: buffer
      
      call push_sub('opencl.opencl_set_kernel_arg_buffer')

      call f90_opencl_set_kernel_arg_buf(kernel, narg, buffer%mem)

      call pop_sub('opencl.opencl_set_kernel_arg_buffer')

    end subroutine opencl_set_kernel_arg_buffer

    ! ------------------------------------------

    subroutine opencl_kernel_run(kernel, globalsizes, localsizes)
      type(c_ptr),        intent(inout) :: kernel
      integer,            intent(in)    :: globalsizes(:)
      integer,            intent(in)    :: localsizes(:)
      
      integer :: dim
      integer(SIZEOF_SIZE_T) :: gsizes(1:3)
      integer(SIZEOF_SIZE_T) :: lsizes(1:3)

      call push_sub('opencl.opencl_kernel_run')
      call profiling_in(prof_kernel_run, "CL_KERNEL_RUN")

      dim = ubound(globalsizes, dim = 1)

      ASSERT(dim == ubound(localsizes, dim = 1))
      ASSERT(all(localsizes <= opencl_max_workgroup_size()))
      ASSERT(all(mod(globalsizes, localsizes) == 0))
     
      gsizes(1:dim) = int(globalsizes(1:dim), SIZEOF_SIZE_T)
      lsizes(1:dim) = int(localsizes(1:dim), SIZEOF_SIZE_T)

      call f90_opencl_kernel_run(kernel, opencl%env, dim, gsizes(1), lsizes(1))
      call opencl_finish()
      call profiling_out(prof_kernel_run)
      call pop_sub('opencl.opencl_kernel_run')
    end subroutine opencl_kernel_run

    ! -----------------------------------------------

    integer pure function opencl_max_workgroup_size() result(max_workgroup_size)
      max_workgroup_size = opencl%max_workgroup_size
    end function opencl_max_workgroup_size

    ! -----------------------------------------------
    integer function opencl_kernel_workgroup_size(kernel) result(workgroup_size)
      type(c_ptr), intent(inout) :: kernel
      
      workgroup_size = f90_opencl_kernel_wgroup_size(kernel, opencl%env)
    end function opencl_kernel_workgroup_size

#include "undef.F90"
#include "real.F90"
#include "opencl_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "opencl_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "opencl_inc.F90"

end module opencl_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

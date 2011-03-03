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
!! $Id: opencl.F90 3587 2007-11-22 16:43:00Z xavier $

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
    opencl_is_enabled

#ifdef HAVE_OPENCL
  public ::                       &
    opencl_t,                     &
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
    opencl_kernel_run,            &
    opencl_build_program,         &
    opencl_release_program,       &
    opencl_release_kernel,        &
    opencl_create_kernel,         &
    f90_cl_device_local_mem_size

  type opencl_t 
    type(c_ptr) :: platform_id
    type(c_ptr) :: context
    type(c_ptr) :: command_queue
    type(c_ptr) :: device
    integer     :: max_workgroup_size
    integer     :: local_memory_size
    logical     :: enabled
  end type opencl_t

  type opencl_mem_t
    private
    type(c_ptr)            :: mem
    integer(SIZEOF_SIZE_T) :: size
    type(type_t)           :: type
  end type opencl_mem_t

  type(opencl_t), public :: opencl

  ! the kernels
  type(c_ptr), public :: kernel_vpsi
  type(c_ptr), public :: kernel_vpsi_spinors
  type(c_ptr), public :: set_zero
  type(c_ptr), public :: set_zero_part
  type(c_ptr), public :: kernel_daxpy
  type(c_ptr), public :: kernel_zaxpy
  type(c_ptr), public :: kernel_copy
  type(c_ptr), public :: kernel_projector_bra
  type(c_ptr), public :: kernel_projector_ket
  type(c_ptr), public :: kernel_projector_ket_copy
  type(c_ptr), public :: dpack
  type(c_ptr), public :: zpack
  type(c_ptr), public :: dunpack
  type(c_ptr), public :: zunpack
  type(c_ptr), public :: kernel_subarray_gather
  type(c_ptr), public :: kernel_density_real
  type(c_ptr), public :: kernel_density_complex
  type(c_ptr), public :: kernel_phase

#include "opencl_iface_inc.F90"

  interface opencl_create_buffer
    module procedure opencl_create_buffer_4, opencl_create_buffer_8
  end interface

  interface opencl_write_buffer
    module procedure iopencl_write_buffer, dopencl_write_buffer, zopencl_write_buffer, &
      iopencl_write_buffer_2, dopencl_write_buffer_2, zopencl_write_buffer_2
  end interface

  interface opencl_read_buffer
    module procedure iopencl_read_buffer, dopencl_read_buffer, zopencl_read_buffer, &
      iopencl_read_buffer_2, dopencl_read_buffer_2, zopencl_read_buffer_2
  end interface

  interface opencl_set_kernel_arg
    module procedure                 &
      opencl_set_kernel_arg_buffer,  &
      iopencl_set_kernel_arg_data,   &
      dopencl_set_kernel_arg_data,   &
      zopencl_set_kernel_arg_data,   &
      opencl_set_kernel_arg_local
  end interface

  type(profile_t), save :: prof_read, prof_write, prof_kernel_run
#endif

  contains

    pure logical function opencl_is_enabled() result(enabled)
#ifdef HAVE_OPENCL
      enabled = opencl%enabled
#else
      enabled = .false.
#endif
    end function opencl_is_enabled

#ifdef HAVE_OPENCL
    ! ------------------------------------------

    subroutine opencl_init()
      type(c_ptr) :: prog
      logical  :: disable
      integer  :: ierr, idevice, iplatform

      PUSH_SUB(opencl_init)
      
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

      !%Variable OpenCLPlatform
      !%Type integer
      !%Default 0
      !%Section Execution::OpenCL
      !%Description
      !% This variable selects the OpenCL platform that Octopus will
      !% use. Platform 0 is used by default.
      !%End
      call parse_integer(datasets_check('OpenCLPlatform'), 0, iplatform)

     !%Variable OpenCLDevice
      !%Type integer
      !%Default 0
      !%Section Execution::OpenCL
      !%Description
      !% This variable selects the OpenCL device that Octopus will
      !% use. Device 0 is used by default.
      !%End
      call parse_integer(datasets_check('OpenCLDevice'), 0, idevice)

      if(disable) then
        POP_SUB(opencl_init)
        return
      end if

      call messages_print_stress(stdout, "OpenCL")

      ierr = flGetPlatformIDs(iplatform, opencl%platform_id)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "GetPlatformIDs")

      call f90_cl_init_context(opencl%platform_id, opencl%context)
      call f90_cl_init_device(idevice, opencl%context, opencl%device)

      call flCreateCommandQueue(opencl%command_queue, opencl%context, opencl%device, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "CreateCommandQueue")
      
      opencl%max_workgroup_size = f90_cl_max_workgroup_size(opencl%device)
      opencl%local_memory_size = f90_cl_device_local_mem_size(opencl%device)

      ! now initialize the kernels
      call opencl_build_program(prog, trim(conf%share)//'/opencl/vpsi.cl')
      call opencl_create_kernel(kernel_vpsi, prog, "vpsi")
      call opencl_create_kernel(kernel_vpsi_spinors, prog, "vpsi_spinors")
      call opencl_release_program(prog)
      
      call opencl_build_program(prog, trim(conf%share)//'/opencl/set_zero.cl')
      call opencl_create_kernel(set_zero, prog, "set_zero")
      call opencl_create_kernel(set_zero_part, prog, "set_zero_part")
      call opencl_release_program(prog)
      
      call opencl_build_program(prog, trim(conf%share)//'/opencl/axpy.cl')
      call opencl_create_kernel(kernel_daxpy, prog, "daxpy")
      call opencl_create_kernel(kernel_zaxpy, prog, "zaxpy")
      call opencl_release_program(prog)

      call opencl_build_program(prog, trim(conf%share)//'/opencl/projector.cl')
      call opencl_create_kernel(kernel_projector_bra, prog, "projector_bra")
      call opencl_create_kernel(kernel_projector_ket, prog, "projector_ket")
      call opencl_create_kernel(kernel_projector_ket_copy, prog, "projector_ket_copy")
      call opencl_release_program(prog)

      call opencl_build_program(prog, trim(conf%share)//'/opencl/pack.cl')
      call opencl_create_kernel(dpack, prog, "dpack")
      call opencl_create_kernel(zpack, prog, "zpack")
      call opencl_create_kernel(dunpack, prog, "dunpack")
      call opencl_create_kernel(zunpack, prog, "zunpack")
      call opencl_release_program(prog)

      call opencl_build_program(prog, trim(conf%share)//'/opencl/copy.cl')
      call opencl_create_kernel(kernel_copy, prog, "copy")
      call opencl_release_program(prog)

      call opencl_build_program(prog, trim(conf%share)//'/opencl/subarray.cl')
      call opencl_create_kernel(kernel_subarray_gather, prog, "subarray_gather")
      call opencl_release_program(prog)

      call opencl_build_program(prog, trim(conf%share)//'/opencl/density.cl')
      call opencl_create_kernel(kernel_density_real, prog, "density_real")
      call opencl_create_kernel(kernel_density_complex, prog, "density_complex")
      call opencl_release_program(prog)

      call opencl_build_program(prog, trim(conf%share)//'/opencl/phase.cl')
      call opencl_create_kernel(kernel_phase, prog, "phase")
      call opencl_release_program(prog)

      call messages_print_stress(stdout)

      POP_SUB(opencl_init)
    end subroutine opencl_init

    ! ------------------------------------------

    subroutine opencl_end()
      integer :: ierr

      PUSH_SUB(opencl_end)

      if(opencl_is_enabled()) then
        call opencl_release_kernel(kernel_vpsi)
        call opencl_release_kernel(kernel_vpsi_spinors)
        call opencl_release_kernel(set_zero)
        call opencl_release_kernel(set_zero_part)
        call opencl_release_kernel(kernel_daxpy)
        call opencl_release_kernel(kernel_zaxpy)
        call opencl_release_kernel(kernel_copy)
        call opencl_release_kernel(kernel_projector_bra)
        call opencl_release_kernel(kernel_projector_ket)
        call opencl_release_kernel(kernel_projector_ket_copy)
        call opencl_release_kernel(dpack)
        call opencl_release_kernel(zpack)
        call opencl_release_kernel(dunpack)
        call opencl_release_kernel(zunpack)
        call opencl_release_kernel(kernel_subarray_gather)
        call opencl_release_kernel(kernel_density_real)
        call opencl_release_kernel(kernel_density_complex)
        call opencl_release_kernel(kernel_phase)

        call flReleaseCommandQueue(opencl%command_queue, ierr)

        if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "ReleaseCommandQueue")
        call flReleaseContext(opencl%context)
      end if

      POP_SUB(opencl_end)
    end subroutine opencl_end

    ! ------------------------------------------

    subroutine opencl_create_buffer_4(this, flags, type, size)
      type(opencl_mem_t), intent(inout) :: this
      integer,            intent(in)    :: flags
      type(type_t),       intent(in)    :: type
      integer,            intent(in)    :: size
      
      integer(SIZEOF_SIZE_T) :: fsize
      integer :: ierr

      PUSH_SUB(opencl_create_buffer_4)

      this%type = type
      this%size = size      
      fsize = size*types_get_size(type)
      
      call f90_cl_create_buffer(this%mem, opencl%context, flags, fsize, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "create_buffer")

      POP_SUB(opencl_create_buffer_4)
    end subroutine opencl_create_buffer_4

    ! ------------------------------------------

    subroutine opencl_create_buffer_8(this, flags, type, size)
      type(opencl_mem_t),     intent(inout) :: this
      integer,                intent(in)    :: flags
      type(type_t),           intent(in)    :: type
      integer(SIZEOF_SIZE_T), intent(in)    :: size
      
      integer(SIZEOF_SIZE_T) :: fsize
      integer :: ierr

      PUSH_SUB(opencl_create_buffer_8)

      this%type = type
      this%size = size

      fsize = size*types_get_size(type)

      call f90_cl_create_buffer(this%mem, opencl%context, flags, fsize, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "create_buffer")

      POP_SUB(opencl_create_buffer_8)
    end subroutine opencl_create_buffer_8

    ! ------------------------------------------

    subroutine opencl_release_buffer(this)
      type(opencl_mem_t), intent(inout) :: this

      integer :: ierr

      this%size = 0
      call f90_cl_release_buffer(this%mem, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "release_buffer")
    end subroutine opencl_release_buffer

    ! ------------------------------------------

    integer(SIZEOF_SIZE_T) pure function opencl_get_buffer_size(this) result(size)
      type(opencl_mem_t), intent(in) :: this

      size = this%size
    end function opencl_get_buffer_size

    ! -----------------------------------------

    type(type_t) pure function opencl_get_buffer_type(this) result(type)
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
      integer :: ierr

      call profiling_in(prof_kernel_run, "CL_KERNEL_RUN")

      call flFinish(opencl%command_queue, ierr)

      call profiling_out(prof_kernel_run)

      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, 'cl_finish') 
    end subroutine opencl_finish

    ! ------------------------------------------

    subroutine opencl_set_kernel_arg_buffer(kernel, narg, buffer)
      type(c_ptr),        intent(inout) :: kernel
      integer,            intent(in)    :: narg
      type(opencl_mem_t), intent(in)    :: buffer
      
      integer :: ierr

      PUSH_SUB(opencl_set_kernel_arg_buffer)

      call f90_cl_set_kernel_arg_buf(kernel, narg, buffer%mem, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "set_kernel_arg_buf")

      POP_SUB(opencl_set_kernel_arg_buffer)

    end subroutine opencl_set_kernel_arg_buffer

    ! ------------------------------------------

    subroutine opencl_set_kernel_arg_local(kernel, narg, type, size)
      type(c_ptr),        intent(inout) :: kernel
      integer,            intent(in)    :: narg
      type(type_t),       intent(in)    :: type
      integer,            intent(in)    :: size

      integer :: ierr
      integer :: size_in_bytes
      PUSH_SUB(opencl_set_kernel_arg_local)

      size_in_bytes = size*types_get_size(type)
      
      if(size_in_bytes > opencl%local_memory_size) then
        write(message(1), '(a,f12.6,a)') "CL Error: requested local memory: ", dble(size_in_bytes)/1024.0, " Kb"
        write(message(2), '(a,f12.6,a)') "          available local memory: ", dble(opencl%local_memory_size)/1024.0, " Kb"
        call write_fatal(2)
      end if

      call f90_cl_set_kernel_arg_local(kernel, narg, size_in_bytes, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "set_kernel_arg_buf")

      POP_SUB(opencl_set_kernel_arg_local)

    end subroutine opencl_set_kernel_arg_local

    ! ------------------------------------------

    subroutine opencl_kernel_run(kernel, globalsizes, localsizes)
      type(c_ptr),        intent(inout) :: kernel
      integer,            intent(in)    :: globalsizes(:)
      integer,            intent(in)    :: localsizes(:)
      
      integer :: dim, ierr
      integer(SIZEOF_SIZE_T) :: gsizes(1:3)
      integer(SIZEOF_SIZE_T) :: lsizes(1:3)

      PUSH_SUB(opencl_kernel_run)
      call profiling_in(prof_kernel_run, "CL_KERNEL_RUN")

      dim = ubound(globalsizes, dim = 1)

      ASSERT(dim == ubound(localsizes, dim = 1))
      ASSERT(all(localsizes <= opencl_max_workgroup_size()))
      ASSERT(all(mod(globalsizes, localsizes) == 0))
     
      gsizes(1:dim) = int(globalsizes(1:dim), SIZEOF_SIZE_T)
      lsizes(1:dim) = int(localsizes(1:dim), SIZEOF_SIZE_T)

      call flEnqueueNDRangeKernel(kernel, opencl%command_queue, dim, gsizes(1), lsizes(1), ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "EnqueueNDRangeKernel")

      call profiling_out(prof_kernel_run)
      POP_SUB(opencl_kernel_run)
    end subroutine opencl_kernel_run

    ! -----------------------------------------------

    integer pure function opencl_max_workgroup_size() result(max_workgroup_size)
      max_workgroup_size = opencl%max_workgroup_size
    end function opencl_max_workgroup_size

    ! -----------------------------------------------
    integer function opencl_kernel_workgroup_size(kernel) result(workgroup_size)
      type(c_ptr), intent(inout) :: kernel
      
      workgroup_size = f90_cl_kernel_wgroup_size(kernel, opencl%device)
    end function opencl_kernel_workgroup_size

    ! -----------------------------------------------

    subroutine opencl_build_program(prog, filename)
      type(c_ptr),      intent(inout) :: prog
      character(len=*), intent(in)    :: filename

      call f90_cl_build_program(prog, opencl%context, opencl%device, filename)
    end subroutine opencl_build_program

    ! -----------------------------------------------

    subroutine opencl_release_program(prog)
      type(c_ptr),      intent(inout) :: prog

      integer :: ierr

      call f90_cl_release_program(prog, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "release_program")
    end subroutine opencl_release_program

    ! -----------------------------------------------

    subroutine opencl_release_kernel(prog)
      type(c_ptr),      intent(inout) :: prog

      integer :: ierr

      call f90_cl_release_kernel(prog, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "release_kernel")
    end subroutine opencl_release_kernel

    ! -----------------------------------------------
    subroutine opencl_create_kernel(kernel, prog, name)
      type(c_ptr),      intent(inout) :: kernel
      type(c_ptr),      intent(inout) :: prog
      character(len=*), intent(in)    :: name

      integer :: ierr

      call f90_cl_create_kernel(kernel, prog, name, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "create_kernel")
    end subroutine opencl_create_kernel

    ! ------------------------------------------------
    
    subroutine opencl_print_error(ierr, name)
      integer,          intent(in) :: ierr
      character(len=*), intent(in) :: name

      character(len=40) :: errcode
    
      select case(ierr)
      case(CL_SUCCESS); errcode = 'CL_SUCCESS '
      case(CL_DEVICE_NOT_FOUND); errcode = 'CL_DEVICE_NOT_FOUND '
      case(CL_DEVICE_NOT_AVAILABLE); errcode = 'CL_DEVICE_NOT_AVAILABLE '
      case(CL_COMPILER_NOT_AVAILABLE); errcode = 'CL_COMPILER_NOT_AVAILABLE '
      case(CL_MEM_OBJECT_ALLOCATION_FAIL); errcode = 'CL_MEM_OBJECT_ALLOCATION_FAILURE '
      case(CL_OUT_OF_RESOURCES); errcode = 'CL_OUT_OF_RESOURCES '
      case(CL_OUT_OF_HOST_MEMORY); errcode = 'CL_OUT_OF_HOST_MEMORY '
      case(CL_PROFILING_INFO_NOT_AVAILABLE); errcode = 'CL_PROFILING_INFO_NOT_AVAILABLE '
      case(CL_MEM_COPY_OVERLAP); errcode = 'CL_MEM_COPY_OVERLAP '
      case(CL_IMAGE_FORMAT_MISMATCH); errcode = 'CL_IMAGE_FORMAT_MISMATCH '
      case(CL_IMAGE_FORMAT_NOT_SUPPORTED); errcode = 'CL_IMAGE_FORMAT_NOT_SUPPORTED '
      case(CL_BUILD_PROGRAM_FAILURE); errcode = 'CL_BUILD_PROGRAM_FAILURE '
      case(CL_MAP_FAILURE); errcode = 'CL_MAP_FAILURE '
      case(CL_INVALID_VALUE); errcode = 'CL_INVALID_VALUE '
      case(CL_INVALID_DEVICE_TYPE); errcode = 'CL_INVALID_DEVICE_TYPE '
      case(CL_INVALID_PLATFORM); errcode = 'CL_INVALID_PLATFORM '
      case(CL_INVALID_DEVICE); errcode = 'CL_INVALID_DEVICE '
      case(CL_INVALID_CONTEXT); errcode = 'CL_INVALID_CONTEXT '
      case(CL_INVALID_QUEUE_PROPERTIES); errcode = 'CL_INVALID_QUEUE_PROPERTIES '
      case(CL_INVALID_COMMAND_QUEUE); errcode = 'CL_INVALID_COMMAND_QUEUE '
      case(CL_INVALID_HOST_PTR); errcode = 'CL_INVALID_HOST_PTR '
      case(CL_INVALID_MEM_OBJECT); errcode = 'CL_INVALID_MEM_OBJECT '
      case(CL_INVALID_IMAGE_FORMAT_DESC); errcode = 'CL_INVALID_IMAGE_FORMAT_DESCRIPTOR '
      case(CL_INVALID_IMAGE_SIZE); errcode = 'CL_INVALID_IMAGE_SIZE '
      case(CL_INVALID_SAMPLER); errcode = 'CL_INVALID_SAMPLER '
      case(CL_INVALID_BINARY); errcode = 'CL_INVALID_BINARY '
      case(CL_INVALID_BUILD_OPTIONS); errcode = 'CL_INVALID_BUILD_OPTIONS '
      case(CL_INVALID_PROGRAM); errcode = 'CL_INVALID_PROGRAM '
      case(CL_INVALID_PROGRAM_EXECUTABLE); errcode = 'CL_INVALID_PROGRAM_EXECUTABLE '
      case(CL_INVALID_KERNEL_NAME); errcode = 'CL_INVALID_KERNEL_NAME '
      case(CL_INVALID_KERNEL_DEFINITION); errcode = 'CL_INVALID_KERNEL_DEFINITION '
      case(CL_INVALID_KERNEL); errcode = 'CL_INVALID_KERNEL '
      case(CL_INVALID_ARG_INDEX); errcode = 'CL_INVALID_ARG_INDEX '
      case(CL_INVALID_ARG_VALUE); errcode = 'CL_INVALID_ARG_VALUE '
      case(CL_INVALID_ARG_SIZE); errcode = 'CL_INVALID_ARG_SIZE '
      case(CL_INVALID_KERNEL_ARGS); errcode = 'CL_INVALID_KERNEL_ARGS '
      case(CL_INVALID_WORK_DIMENSION); errcode = 'CL_INVALID_WORK_DIMENSION '
      case(CL_INVALID_WORK_GROUP_SIZE); errcode = 'CL_INVALID_WORK_GROUP_SIZE '
      case(CL_INVALID_WORK_ITEM_SIZE); errcode = 'CL_INVALID_WORK_ITEM_SIZE '
      case(CL_INVALID_GLOBAL_OFFSET); errcode = 'CL_INVALID_GLOBAL_OFFSET '
      case(CL_INVALID_EVENT_WAIT_LIST); errcode = 'CL_INVALID_EVENT_WAIT_LIST '
      case(CL_INVALID_EVENT); errcode = 'CL_INVALID_EVENT '
      case(CL_INVALID_OPERATION); errcode = 'CL_INVALID_OPERATION '
      case(CL_INVALID_GL_OBJECT); errcode = 'CL_INVALID_GL_OBJECT '
      case(CL_INVALID_BUFFER_SIZE); errcode = 'CL_INVALID_BUFFER_SIZE '
      case(CL_INVALID_MIP_LEVEL); errcode = 'CL_INVALID_MIP_LEVEL '
      case(CL_INVALID_GLOBAL_WORK_SIZE); errcode = 'CL_INVALID_GLOBAL_WORK_SIZE '
      end select

      message(1) = 'Error: OpenCL '//trim(name)//' '//trim(errcode)
      call write_fatal(1)
  
    end subroutine opencl_print_error

#include "undef.F90"
#include "real.F90"
#include "opencl_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "opencl_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "opencl_inc.F90"
#endif
end module opencl_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

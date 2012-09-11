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
#ifdef HAVE_OPENCL
  use cl
#endif
#ifdef HAVE_CLAMDBLAS
  use clAmdBlas
#endif
#ifdef HAVE_CLAMDFFT
  use clAmdFft
#endif
  use datasets_m
  use global_m
  use io_m
  use messages_m
  use mpi_m
  use types_m
  use parser_m
  use profiling_m
  use unit_system_m

  implicit none 

  private

  public ::                       &
    opencl_is_enabled,            &
    opencl_t,                     &
    opencl_init,                  &
    opencl_end,                   &
    opencl_padded_size,           &
    opencl_mem_t

#ifdef HAVE_OPENCL
  public ::                       &
    opencl_create_buffer,         &
    opencl_write_buffer,          &
    opencl_read_buffer,           &
    opencl_release_buffer,        &
    opencl_finish,                &
    opencl_set_kernel_arg,        &
    opencl_max_workgroup_size,    &
    opencl_kernel_workgroup_size, &
    opencl_kernel_run,            &
    opencl_build_program,         &
    opencl_release_program,       &
    opencl_release_kernel,        &
    opencl_create_kernel,         &
    opencl_print_error,           &
    clblas_print_error,           &
    clfft_print_error,            &
    opencl_set_buffer_to_zero
#endif

  type opencl_t 
#ifdef HAVE_OPENCL
    type(cl_platform_id)   :: platform_id
    type(cl_context)       :: context
    type(cl_command_queue) :: command_queue
    type(cl_device_id)     :: device
#endif
    integer                :: max_workgroup_size
    integer                :: local_memory_size
    logical                :: enabled
  end type opencl_t

  type opencl_mem_t
#ifdef HAVE_OPENCL
    type(cl_mem)           :: mem
#endif
    integer(SIZEOF_SIZE_T) :: size
    type(type_t)           :: type
  end type opencl_mem_t

  type(opencl_t), public :: opencl

#ifdef HAVE_OPENCL
  ! the kernels
  type(cl_kernel), public :: kernel_vpsi
  type(cl_kernel), public :: kernel_vpsi_spinors
  type(cl_kernel), public :: set_zero_part
  type(cl_kernel), public :: kernel_daxpy
  type(cl_kernel), public :: kernel_zaxpy
  type(cl_kernel), public :: kernel_copy
  type(cl_kernel), public :: kernel_projector_bra
  type(cl_kernel), public :: kernel_projector_ket
  type(cl_kernel), public :: dpack
  type(cl_kernel), public :: zpack
  type(cl_kernel), public :: dunpack
  type(cl_kernel), public :: zunpack
  type(cl_kernel), public :: kernel_subarray_gather
  type(cl_kernel), public :: kernel_density_real
  type(cl_kernel), public :: kernel_density_complex
  type(cl_kernel), public :: kernel_phase
  type(cl_kernel), public :: dkernel_dot_matrix
  type(cl_kernel), public :: zkernel_dot_matrix
  type(cl_kernel), public :: zkernel_dot_matrix_spinors
  type(cl_kernel), public :: dkernel_dot_vector
  type(cl_kernel), public :: zkernel_dot_vector
  type(cl_kernel), public :: kernel_nrm2_vector
  type(cl_kernel), public :: dzmul
  type(cl_kernel), public :: zzmul

  ! kernels used locally
  type(cl_kernel)         :: set_zero

  interface opencl_create_buffer
    module procedure opencl_create_buffer_4
  end interface

  interface opencl_write_buffer
    module procedure iopencl_write_buffer_1, dopencl_write_buffer_1, zopencl_write_buffer_1
    module procedure iopencl_write_buffer_2, dopencl_write_buffer_2, zopencl_write_buffer_2
    module procedure iopencl_write_buffer_3, dopencl_write_buffer_3, zopencl_write_buffer_3
  end interface

  interface opencl_read_buffer
    module procedure iopencl_read_buffer_1, dopencl_read_buffer_1, zopencl_read_buffer_1
    module procedure iopencl_read_buffer_2, dopencl_read_buffer_2, zopencl_read_buffer_2
    module procedure iopencl_read_buffer_3, dopencl_read_buffer_3, zopencl_read_buffer_3
  end interface

  interface opencl_set_kernel_arg
    module procedure                 &
      opencl_set_kernel_arg_buffer,  &
      iopencl_set_kernel_arg_data,   &
      dopencl_set_kernel_arg_data,   &
      zopencl_set_kernel_arg_data,   &
      opencl_set_kernel_arg_local
  end interface

#endif

  type(profile_t), save :: prof_read, prof_write, prof_kernel_run

  integer, parameter  ::      &
    OPENCL_GPU         = -1,  &
    OPENCL_CPU         = -2,  &
    OPENCL_ACCELERATOR = -3,  &
    OPENCL_DEFAULT     = -4


  integer, parameter  ::      &
    CL_PLAT_INVALID   = -1,   &
    CL_PLAT_AMD       = -2,   &
    CL_PLAT_NVIDIA    = -3,   &
    CL_PLAT_ATI       = -4,   &
    CL_PLAT_INTEL     = -5

  ! a "convenience" public variable
  integer, public :: cl_status

  integer, parameter :: OPENCL_MAX_FILE_LENGTH = 10000

  integer :: buffer_alloc_count
  integer(8) :: allocated_mem

  contains

    pure logical function opencl_is_enabled() result(enabled)
#ifdef HAVE_OPENCL
      enabled = opencl%enabled
#else
      enabled = .false.
#endif
    end function opencl_is_enabled

    ! ------------------------------------------
    
    subroutine opencl_init(base_grp)
      type(mpi_grp_t),  intent(inout) :: base_grp

      logical  :: disable, default
      integer  :: device_type
      integer  :: idevice, iplatform, ndevices, idev, cl_status, ret_devices, nplatforms, iplat
      character(len=256) :: device_name
#ifdef HAVE_OPENCL
      type(cl_program) :: prog
      type(cl_platform_id), allocatable :: allplatforms(:)
      type(cl_device_id), allocatable :: alldevices(:)
      type(profile_t), save :: prof_init
#endif

      PUSH_SUB(opencl_init)

      buffer_alloc_count = 0

      !%Variable DisableOpenCL
      !%Type logical
      !%Default yes
      !%Section Execution::OpenCL
      !%Description
      !% If Octopus was compiled with OpenCL support, it will try to
      !% initialize and use an OpenCL device. By setting this variable
      !% to <tt>yes</tt> you tell Octopus not to use OpenCL.
      !%End

#ifndef HAVE_OPENCL
      default = .true.
#else
      default = .false.
#endif
      call parse_logical(datasets_check('DisableOpenCL'), default, disable)
      opencl%enabled = .not. disable

#ifndef HAVE_OPENCL
      if(opencl%enabled) then
        message(1) = 'Octopus was compiled without OpenCL support.'
        call messages_fatal(1)
      end if
#endif

      if(.not. opencl_is_enabled()) then
        POP_SUB(opencl_init)
        return
      end if

      !%Variable OpenCLPlatform
      !%Type integer
      !%Default 0
      !%Section Execution::OpenCL
      !%Description
      !% This variable selects the OpenCL platform that Octopus will
      !% use. You can give an explicit platform number or use one of
      !% the options that select a particular vendor
      !% implementation. Platform 0 is used by default.
      !%Option amd -2
      !% Use the AMD OpenCL platform.
      !%Option nvidia -3
      !% Use the Nvidia OpenCL platform.
      !%Option ati -4
      !% Use the ATI (old AMD) OpenCL platform.
      !%Option intel -5
      !% Use the Intel OpenCL platform.
      !%End
      call parse_integer(datasets_check('OpenCLPlatform'), 0, iplatform)

      !%Variable OpenCLDevice
      !%Type integer
      !%Default 0
      !%Section Execution::OpenCL
      !%Description
      !% This variable selects the OpenCL device that Octopus will
      !% use. You can specify one of the options below or a numerical
      !% id to select a specific device.
      !%Option gpu -1
      !% If available, Octopus will use a GPU for OpenCL. This is the default.
      !%Option cpu -2
      !% If available, Octopus will use a GPU for OpenCL.
      !%Option accelerator -3
      !% If available, Octopus will use an accelerator for OpenCL.
      !%Option cl_default -4
      !% Octopus will use the default device specified by the OpenCL
      !% implementation.
      !%End
      call parse_integer(datasets_check('OpenCLDevice'), OPENCL_GPU, idevice)

      if(idevice < OPENCL_DEFAULT) then
        message(1) = 'Invalid OpenCLDevice.'
        call messages_fatal(1)
      end if

      call messages_print_stress(stdout, "OpenCL")

#ifdef HAVE_OPENCL
      call profiling_in(prof_init, 'CL_INIT')
      
      call clGetPlatformIDs(nplatforms, cl_status)
      if(cl_status /= CL_SUCCESS) call opencl_print_error(cl_status, "GetPlatformIDs")
 
      SAFE_ALLOCATE(allplatforms(1:nplatforms))

      call clGetPlatformIDs(allplatforms, iplat, cl_status)
      if(cl_status /= CL_SUCCESS) call opencl_print_error(cl_status, "GetPlatformIDs")

      call messages_write('Info: Available CL platforms: ')
      call messages_write(nplatforms)
      call messages_info()

      do iplat = 1, nplatforms

        call clGetPlatformInfo(allplatforms(iplat), CL_PLATFORM_NAME, device_name, cl_status)

        if(iplatform < 0) then
          if(iplatform == get_platform_id(device_name)) iplatform = iplat - 1
        end if

        if(iplatform == iplat - 1) then
          call messages_write('    * Platform ')
        else
          call messages_write('      Platform ')
        end if

        call messages_write(iplat - 1)
        call messages_write(' : '//device_name)
        call clGetPlatformInfo(allplatforms(iplat), CL_PLATFORM_VERSION, device_name, cl_status)
        call messages_write(' ('//trim(device_name)//')')
        call messages_info()
      end do

      call messages_info()

      if(iplatform >= nplatforms .or. iplatform < 0) then
        call messages_write('Requested CL platform does not exist')
        if(iplatform > 0) then 
          call messages_write('(platform = ')
          call messages_write(iplatform)
          call messages_write(').')
        end if
        call messages_fatal()
      end if

      opencl%platform_id = allplatforms(iplatform + 1)

      SAFE_DEALLOCATE_A(allplatforms)

      call clGetDeviceIDs(opencl%platform_id, CL_DEVICE_TYPE_ALL, ndevices, cl_status)

      call messages_write('Info: Available CL devices: ')
      call messages_write(ndevices)
      call messages_info()

      SAFE_ALLOCATE(alldevices(1:ndevices))

      ! list all devices

      call clGetDeviceIDs(opencl%platform_id, CL_DEVICE_TYPE_ALL, alldevices, ret_devices, cl_status)

      do idev = 1, ndevices
        call messages_write('      Device ')
        call messages_write(idev - 1)
        call clGetDeviceInfo(alldevices(idev), CL_DEVICE_NAME, device_name, cl_status)
        call messages_write(' : '//device_name)
        call messages_info()
      end do

      select case(idevice)
        case(OPENCL_GPU)
          device_type = CL_DEVICE_TYPE_GPU
        case(OPENCL_CPU)
          device_type = CL_DEVICE_TYPE_CPU
        case(OPENCL_ACCELERATOR)
          device_type = CL_DEVICE_TYPE_ACCELERATOR
        case(OPENCL_DEFAULT)
          device_type = CL_DEVICE_TYPE_DEFAULT
        case default
          device_type = CL_DEVICE_TYPE_ALL
      end select

      ! now get a list of the selected type
      call clGetDeviceIDs(opencl%platform_id, device_type, alldevices, ret_devices, cl_status)

      ! the number of devices can be smaller
      ndevices = ret_devices

      if(idevice < 0) then
        if(base_grp%size > 1) then
          ! with MPI we have to select the device so multiple GPUs in one
          ! node are correctly distributed
          call select_device(idevice)
        else
          idevice = 0
        end if
      end if

      if(idevice >= ndevices) then
        call messages_write('Requested CL device does not exist (device = ')
        call messages_write(idevice)
        call messages_write(', platform = ')
        call messages_write(iplatform)
        call messages_write(').')
        call messages_fatal()
      end if

      opencl%device = alldevices(idevice + 1)

      if(mpi_grp_is_root(base_grp)) call device_info()

      ! create the context
      opencl%context = clCreateContext(opencl%platform_id, opencl%device, cl_status)
      if(cl_status /= CL_SUCCESS) call opencl_print_error(cl_status, "CreateContext")

      SAFE_DEALLOCATE_A(alldevices)

      opencl%command_queue = clCreateCommandQueue(opencl%context, opencl%device, CL_QUEUE_PROFILING_ENABLE, cl_status)
      if(cl_status /= CL_SUCCESS) call opencl_print_error(cl_status, "CreateCommandQueue")
      
      call clGetDeviceInfo(opencl%device, CL_DEVICE_MAX_WORK_GROUP_SIZE, opencl%max_workgroup_size, cl_status)
      call clGetDeviceInfo(opencl%device, CL_DEVICE_LOCAL_MEM_SIZE, opencl%local_memory_size, cl_status)

      ! now initialize the kernels
      call opencl_build_program(prog, trim(conf%share)//'/opencl/set_zero.cl')
      call opencl_create_kernel(set_zero, prog, "set_zero")
      call opencl_create_kernel(set_zero_part, prog, "set_zero_part")
      call opencl_release_program(prog)
      
      call opencl_build_program(prog, trim(conf%share)//'/opencl/vpsi.cl')
      call opencl_create_kernel(kernel_vpsi, prog, "vpsi")
      call opencl_create_kernel(kernel_vpsi_spinors, prog, "vpsi_spinors")
      call opencl_release_program(prog)
      
      call opencl_build_program(prog, trim(conf%share)//'/opencl/axpy.cl', flags = '-DRTYPE_DOUBLE')
      call opencl_create_kernel(kernel_daxpy, prog, "daxpy")
      call opencl_release_program(prog)

      call opencl_build_program(prog, trim(conf%share)//'/opencl/axpy.cl', flags = '-DRTYPE_COMPLEX')
      call opencl_create_kernel(kernel_zaxpy, prog, "zaxpy")
      call opencl_release_program(prog)

      call opencl_build_program(prog, trim(conf%share)//'/opencl/projector.cl')
      call opencl_create_kernel(kernel_projector_bra, prog, "projector_bra")
      call opencl_create_kernel(kernel_projector_ket, prog, "projector_ket")
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

      call opencl_build_program(prog, trim(conf%share)//'/opencl/mesh_batch.cl')
      call opencl_create_kernel(dkernel_dot_vector, prog, "ddot_vector")
      call opencl_create_kernel(zkernel_dot_vector, prog, "zdot_vector")
      call opencl_create_kernel(dkernel_dot_matrix, prog, "ddot_matrix")
      call opencl_create_kernel(zkernel_dot_matrix, prog, "zdot_matrix")
      call opencl_create_kernel(zkernel_dot_matrix_spinors, prog, "zdot_matrix_spinors")
      call opencl_create_kernel(kernel_nrm2_vector, prog, "nrm2_vector")
      call opencl_release_program(prog)

      call opencl_build_program(prog, trim(conf%share)//'/opencl/mul.cl', flags = '-DRTYPE_DOUBLE')
      call opencl_create_kernel(dzmul, prog, "dzmul")
      call opencl_release_program(prog)

      call opencl_build_program(prog, trim(conf%share)//'/opencl/mul.cl', flags = '-DRTYPE_COMPLEX')
      call opencl_create_kernel(zzmul, prog, "zzmul")
      call opencl_release_program(prog)

#ifdef HAVE_CLAMDBLAS
      call clAmdBlasSetup(cl_status)
      if(cl_status /= clAmdBlasSuccess) call clblas_print_error(cl_status, 'clAmdBlasSetup')
#endif

#ifdef HAVE_CLAMDFFT
      call clAmdFftSetup(cl_status)
      if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clAmdFftSetup')
#endif

      call profiling_out(prof_init)
#endif

      call messages_print_stress(stdout)

      POP_SUB(opencl_init)

    contains
      
      subroutine select_device(idevice)
        integer, intent(inout) :: idevice
#if defined(HAVE_MPI) && defined(HAVE_OPENCL)
        integer :: irank
        character(len=256) :: device_name

        PUSH_SUB(opencl_init.select_device)

        idevice = mod(base_grp%rank, ndevices)

        call MPI_Barrier(base_grp%comm, mpi_err)
        call messages_write('Info: CL device distribution:')
        call messages_info()
        do irank = 0, base_grp%size - 1
          if(irank == base_grp%rank) then
            call clGetDeviceInfo(alldevices(idevice + 1), CL_DEVICE_NAME, device_name, cl_status)
            call messages_write('      MPI node ')
            call messages_write(base_grp%rank)
            call messages_write(' -> CL device ')
            call messages_write(idevice)
            call messages_write(' : '//device_name)
            call messages_info(all_nodes = .true.)
          end if
          call MPI_Barrier(base_grp%comm, mpi_err)
        end do
#endif

        POP_SUB(opencl_init.select_device)
      end subroutine select_device

      subroutine device_info()

#ifdef HAVE_OPENCL
        integer(8) :: val 
        character(len=256) :: val_str

        PUSH_SUB(opencl_init.device_info)

        call messages_new_line()
        call messages_write('Selected CL device:')
        call messages_new_line()

        call clGetDeviceInfo(opencl%device, CL_DEVICE_VENDOR, val_str, cl_status)
        call messages_write('      Device vendor          : '//trim(val_str))
        call messages_new_line()

        call clGetDeviceInfo(opencl%device, CL_DEVICE_NAME, val_str, cl_status)
        call messages_write('      Device name            : '//trim(val_str))
        call messages_new_line()

        call clGetDeviceInfo(opencl%device, CL_DRIVER_VERSION, val_str, cl_status)
        call messages_write('      Driver version         : '//trim(val_str))
        call messages_new_line()

        call clGetDeviceInfo(opencl%device, CL_DEVICE_MAX_COMPUTE_UNITS, val, cl_status)
        call messages_write('      Compute units          :')
        call messages_write(val)
        call messages_new_line()

        call clGetDeviceInfo(opencl%device, CL_DEVICE_MAX_CLOCK_FREQUENCY, val, cl_status)
        call messages_write('      Clock frequency        :')
        call messages_write(val)
        call messages_write(' GHz')
        call messages_new_line()

        call clGetDeviceInfo(opencl%device, CL_DEVICE_GLOBAL_MEM_SIZE, val, cl_status)
        call messages_write('      Device memory          :')
        call messages_write(val/(1024**2))
        call messages_write(' Mb')
        call messages_new_line()

        call clGetDeviceInfo(opencl%device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, val, cl_status)
        call messages_write('      Max alloc size         :')
        call messages_write(val/(1024**2))
        call messages_write(' Mb')
        call messages_new_line()

        call clGetDeviceInfo(opencl%device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, val, cl_status)
        call messages_write('      Device cache           :')
        call messages_write(val/1024)
        call messages_write(' Kb')
        call messages_new_line()

        call clGetDeviceInfo(opencl%device, CL_DEVICE_LOCAL_MEM_SIZE, val, cl_status)
        call messages_write('      Local memory           :')
        call messages_write(val/1024)
        call messages_write(' Kb')
        call messages_new_line()

        call clGetDeviceInfo(opencl%device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, val, cl_status)
        call messages_write('      Constant memory        :')
        call messages_write(val/1024)
        call messages_write(' Kb')
        call messages_new_line()

        call clGetDeviceInfo(opencl%device, CL_DEVICE_MAX_WORK_GROUP_SIZE, val, cl_status)
        call messages_write('      Max. workgroup size    :')
        call messages_write(val)
        call messages_new_line()

        call messages_write('      Extension cl_khr_fp64  :')
        call messages_write(f90_cl_device_has_extension(opencl%device, "cl_khr_fp64"))
        call messages_new_line()

        call messages_write('      Extension cl_amd_fp64  :')
        call messages_write(f90_cl_device_has_extension(opencl%device, "cl_amd_fp64"))
        call messages_new_line()

        call messages_info()

        POP_SUB(opencl_init.device_info)
#endif
      end subroutine device_info

    end subroutine opencl_init

    ! ------------------------------------------

    integer function get_platform_id(platform_name) result(platform_id)
      character(len=*), intent(in) :: platform_name

      platform_id = CL_PLAT_INVALID
      if(index(platform_name, 'AMD') > 0)    platform_id = CL_PLAT_AMD
      if(index(platform_name, 'ATI') > 0)    platform_id = CL_PLAT_ATI
      if(index(platform_name, 'NVIDIA') > 0) platform_id = CL_PLAT_NVIDIA
      if(index(platform_name, 'Intel') > 0)  platform_id = CL_PLAT_INTEL
    end function get_platform_id

    ! ------------------------------------------

    subroutine opencl_end()
      integer :: ierr

      PUSH_SUB(opencl_end)

#ifdef HAVE_CLAMDBLAS
      call clAmdBlasTearDown()
#endif

#ifdef HAVE_CLAMDFFT
      call clAmdFftTearDown()
#endif

      if(opencl_is_enabled()) then
#ifdef HAVE_OPENCL
        call opencl_release_kernel(kernel_vpsi)
        call opencl_release_kernel(kernel_vpsi_spinors)
        call opencl_release_kernel(set_zero)
        call opencl_release_kernel(set_zero_part)
        call opencl_release_kernel(kernel_daxpy)
        call opencl_release_kernel(kernel_zaxpy)
        call opencl_release_kernel(kernel_copy)
        call opencl_release_kernel(kernel_projector_bra)
        call opencl_release_kernel(kernel_projector_ket)
        call opencl_release_kernel(dpack)
        call opencl_release_kernel(zpack)
        call opencl_release_kernel(dunpack)
        call opencl_release_kernel(zunpack)
        call opencl_release_kernel(kernel_subarray_gather)
        call opencl_release_kernel(kernel_density_real)
        call opencl_release_kernel(kernel_density_complex)
        call opencl_release_kernel(kernel_phase)
        call opencl_release_kernel(dkernel_dot_matrix)
        call opencl_release_kernel(zkernel_dot_matrix)
        call opencl_release_kernel(dkernel_dot_vector)
        call opencl_release_kernel(zkernel_dot_vector)
        call opencl_release_kernel(zkernel_dot_matrix_spinors)

        
        call clReleaseCommandQueue(opencl%command_queue, ierr)

        if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "ReleaseCommandQueue")
        call clReleaseContext(opencl%context, cl_status)

        if(buffer_alloc_count /= 0) then
          call messages_write('OpenCL:')
          call messages_write(real(allocated_mem, REAL_PRECISION), fmt = 'f12.1', units = unit_megabytes, align_left = .true.)
          call messages_write(' in ')
          call messages_write(buffer_alloc_count)
          call messages_write(' buffers were not deallocated.')
          call messages_warning()
        end if
#endif
      end if

      POP_SUB(opencl_end)
    end subroutine opencl_end

    ! ------------------------------------------

    integer function opencl_padded_size(nn) result(psize)
      integer,        intent(in) :: nn

#ifdef HAVE_OPENCL
      integer :: modnn, bsize

      bsize = opencl_max_workgroup_size()

      psize = nn
      modnn = mod(nn, bsize)
      if(modnn /= 0) psize = psize + bsize - modnn
#else
      psize = nn
#endif
    end function opencl_padded_size

#ifdef HAVE_OPENCL
    ! ------------------------------------------

    subroutine opencl_create_buffer_4(this, flags, type, size)
      type(opencl_mem_t), intent(inout) :: this
      integer,            intent(in)    :: flags
      type(type_t),       intent(in)    :: type
      integer,            intent(in)    :: size
      
      integer(8) :: fsize
      integer :: ierr

      PUSH_SUB(opencl_create_buffer_4)

      this%type = type
      this%size = size
      fsize = int(size, 8)*types_get_size(type)

      ASSERT(fsize >= 0)

      this%mem = clCreateBuffer(opencl%context, flags, fsize, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clCreateBuffer")

      INCR(buffer_alloc_count, 1)
      INCR(allocated_mem, fsize)

      POP_SUB(opencl_create_buffer_4)
    end subroutine opencl_create_buffer_4

    ! ------------------------------------------

    subroutine opencl_release_buffer(this)
      type(opencl_mem_t), intent(inout) :: this

      integer :: ierr

      PUSH_SUB(opencl_release_buffer)

      call clReleaseMemObject(this%mem, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clReleaseMemObject")

      INCR(buffer_alloc_count, -1)
      INCR(allocated_mem, -int(this%size, 8)*types_get_size(this%type))

      this%size = 0

      POP_SUB(opencl_release_buffer)
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

    subroutine opencl_finish()
      integer :: ierr

      PUSH_SUB(opencl_finish)

      call profiling_in(prof_kernel_run, "CL_KERNEL_RUN")

      call clFinish(opencl%command_queue, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, 'clFinish') 

      call profiling_out(prof_kernel_run)
      POP_SUB(opencl_finish)
    end subroutine opencl_finish

    ! ------------------------------------------

    subroutine opencl_set_kernel_arg_buffer(kernel, narg, buffer)
      type(cl_kernel),    intent(inout) :: kernel
      integer,            intent(in)    :: narg
      type(opencl_mem_t), intent(in)    :: buffer
      
      integer :: ierr

      PUSH_SUB(opencl_set_kernel_arg_buffer)

      call clSetKernelArg(kernel, narg, buffer%mem, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clSetKernelArg_buf")

      POP_SUB(opencl_set_kernel_arg_buffer)

    end subroutine opencl_set_kernel_arg_buffer

    ! ------------------------------------------

    subroutine opencl_set_kernel_arg_local(kernel, narg, type, size)
      type(cl_kernel),    intent(inout) :: kernel
      integer,            intent(in)    :: narg
      type(type_t),       intent(in)    :: type
      integer,            intent(in)    :: size

      integer :: ierr
      integer(8) :: size_in_bytes

      PUSH_SUB(opencl_set_kernel_arg_local)

      size_in_bytes = int(size, 8)*types_get_size(type)
      
      if(size_in_bytes > opencl%local_memory_size) then
        write(message(1), '(a,f12.6,a)') "CL Error: requested local memory: ", dble(size_in_bytes)/1024.0, " Kb"
        write(message(2), '(a,f12.6,a)') "          available local memory: ", dble(opencl%local_memory_size)/1024.0, " Kb"
        call messages_fatal(2)
      else if(size_in_bytes <= 0) then
        write(message(1), '(a,i10)') "CL Error: invalid local memory size: ", size_in_bytes
        call messages_fatal(1)
      end if

      call clSetKernelArgLocal(kernel, narg, size_in_bytes, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "set_kernel_arg_local")

      POP_SUB(opencl_set_kernel_arg_local)

    end subroutine opencl_set_kernel_arg_local

    ! ------------------------------------------

    subroutine opencl_kernel_run(kernel, globalsizes, localsizes)
      type(cl_kernel),    intent(inout) :: kernel
      integer,            intent(in)    :: globalsizes(:)
      integer,            intent(in)    :: localsizes(:)
      
      integer :: dim, ierr
      integer(8) :: gsizes(1:3)
      integer(8) :: lsizes(1:3)

      PUSH_SUB(opencl_kernel_run)
      call profiling_in(prof_kernel_run, "CL_KERNEL_RUN")

      dim = ubound(globalsizes, dim = 1)

      ASSERT(dim == ubound(localsizes, dim = 1))
      ASSERT(all(localsizes <= opencl_max_workgroup_size()))
      ASSERT(all(mod(globalsizes, localsizes) == 0))
     
      gsizes(1:dim) = int(globalsizes(1:dim), 8)
      lsizes(1:dim) = int(localsizes(1:dim), 8)

      call clEnqueueNDRangeKernel(opencl%command_queue, kernel, gsizes(1:dim), lsizes(1:dim), ierr)
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
      type(cl_kernel), intent(inout) :: kernel

      integer(8) :: workgroup_size8
      integer    :: ierr

      call clGetKernelWorkGroupInfo(kernel, opencl%device, CL_KERNEL_WORK_GROUP_SIZE, workgroup_size8, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "EnqueueNDRangeKernel")
      workgroup_size = workgroup_size8

    end function opencl_kernel_workgroup_size

    ! -----------------------------------------------

    subroutine opencl_build_program(prog, filename, flags)
      type(cl_program),           intent(inout) :: prog
      character(len=*),           intent(in)    :: filename
      character(len=*), optional, intent(in)    :: flags
      
      character(len = OPENCL_MAX_FILE_LENGTH) :: string
      integer :: ierr, ierrlog, iunit, irec
      type(profile_t), save :: prof

      PUSH_SUB(opencl_build_program)
      call profiling_in(prof, "CL_COMPILE", exclude = .true.)

      string = ''

      call io_assign(iunit)
      open(unit = iunit, file = trim(filename), access='direct', status = 'old', action = 'read', iostat = ierr, recl = 1)
      irec = 1
      do
        read(unit = iunit, rec = irec, iostat = ierr) string(irec:irec) 
        if (ierr /= 0) exit
        if(irec == OPENCL_MAX_FILE_LENGTH) then
          call messages_write('CL source file is too big: '//trim(filename)//'.')
          call messages_new_line()
          call messages_write("       Increase 'OPENCL_MAX_FILE_LENGTH'.")
          call messages_fatal()
        end if
        irec = irec + 1
      end do

      close(unit = iunit)
      call io_free(iunit)

      call messages_write("Building CL program '"//trim(filename)//"'.")
      call messages_info()

      prog = clCreateProgramWithSource(opencl%context, string, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clCreateProgramWithSource")

      ! build the compilation flags
      string='-w'
      ! full optimization
      string=trim(string)//' -cl-denorms-are-zero'
      string=trim(string)//' -cl-strict-aliasing'
      string=trim(string)//' -cl-mad-enable'
      string=trim(string)//' -cl-unsafe-math-optimizations'
      string=trim(string)//' -cl-finite-math-only'
      string=trim(string)//' -cl-fast-relaxed-math'

      string=trim(string)//' -I'//trim(conf%share)//'/opencl/'
      
      if (f90_cl_device_has_extension(opencl%device, "cl_amd_fp64")) then
        string = trim(string)//' -DEXT_AMD_FP64'
      else if(f90_cl_device_has_extension(opencl%device, "cl_khr_fp64")) then
        string = trim(string)//' -DEXT_KHR_FP64'
      else
        call messages_write('Octopus requires an OpenCL device with double-precision support.')
        call messages_fatal()
      end if

      if(present(flags)) then
        string = trim(string)//' '//trim(flags)
      end if

      if(in_debug_mode) then
        message(1) = "Debug info: compilation flags '"//trim(string)//"'. "
        call messages_info(1)
      end if

      call clBuildProgram(prog, trim(string), ierr)

      call clGetProgramBuildInfo(prog, opencl%device, CL_PROGRAM_BUILD_LOG, string, ierrlog)
      if(ierrlog /= CL_SUCCESS) call opencl_print_error(ierrlog, "clGetProgramBuildInfo")
      
      if(len(trim(string)) > 0) write(stderr, '(a)') trim(string)

      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clBuildProgram")

      call profiling_out(prof)
      POP_SUB(opencl_build_program)
    end subroutine opencl_build_program

    ! -----------------------------------------------

    subroutine opencl_release_program(prog)
      type(cl_program),    intent(inout) :: prog

      integer :: ierr

      PUSH_SUB(opencl_release_program)

      call clReleaseProgram(prog, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clReleaseProgram")

      POP_SUB(opencl_release_program)
    end subroutine opencl_release_program

    ! -----------------------------------------------

    subroutine opencl_release_kernel(prog)
      type(cl_kernel),      intent(inout) :: prog

      integer :: ierr

      PUSH_SUB(opencl_release_kernel)

      call clReleaseKernel(prog, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clReleaseKernel")

      POP_SUB(opencl_release_kernel)
    end subroutine opencl_release_kernel

    ! -----------------------------------------------
    subroutine opencl_create_kernel(kernel, prog, name)
      type(cl_kernel),  intent(inout) :: kernel
      type(cl_program), intent(inout) :: prog
      character(len=*), intent(in)    :: name

      integer :: ierr
      type(profile_t), save :: prof

      PUSH_SUB(opencl_create_kernel)
      call profiling_in(prof, "CL_BUILD_KERNEL", exclude = .true.)

      kernel = clCreateKernel(prog, name, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clCreateKernel")

      call profiling_out(prof)
      POP_SUB(opencl_create_kernel)
    end subroutine opencl_create_kernel

    ! ------------------------------------------------
    
    subroutine opencl_print_error(ierr, name)
      integer,          intent(in) :: ierr
      character(len=*), intent(in) :: name

      character(len=40) :: errcode
    
      PUSH_SUB(opencl_print_error)

      select case(ierr)
      case(CL_SUCCESS); errcode = 'CL_SUCCESS '
      case(CL_DEVICE_NOT_FOUND); errcode = 'CL_DEVICE_NOT_FOUND '
      case(CL_DEVICE_NOT_AVAILABLE); errcode = 'CL_DEVICE_NOT_AVAILABLE '
      case(CL_COMPILER_NOT_AVAILABLE); errcode = 'CL_COMPILER_NOT_AVAILABLE '
      case(CL_MEM_OBJECT_ALLOCATION_FAILURE); errcode = 'CL_MEM_OBJECT_ALLOCATION_FAILURE '
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
      case(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR); errcode = 'CL_INVALID_IMAGE_FORMAT_DESCRIPTOR '
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
      case(CL_PLATFORM_NOT_FOUND_KHR); errcode = 'CL_PLATFORM_NOT_FOUND_KHR'
      case default
        write(errcode, '(i10)') ierr
        errcode = 'UNKNOWN ERROR CODE ('//trim(adjustl(errcode))//')'
      end select

      message(1) = 'OpenCL '//trim(name)//' '//trim(errcode)
      call messages_fatal(1)
  
      POP_SUB(opencl_print_error)
    end subroutine opencl_print_error

    ! ----------------------------------------------------

    subroutine clblas_print_error(ierr, name)
      integer,          intent(in) :: ierr
      character(len=*), intent(in) :: name

      character(len=40) :: errcode
    
      PUSH_SUB(clblas_print_error)
#ifdef HAVE_CLAMDBLAS
      select case(ierr)
      case(clAmdBlasSuccess);                    errcode = 'clAmdBlasSuccess'
      case(clAmdBlasInvalidValue);               errcode = 'clAmdBlasInvalidValue'
      case(clAmdBlasInvalidCommandQueue);        errcode = 'clAmdBlasInvalidCommandQueue'
      case(clAmdBlasInvalidContext);             errcode = 'clAmdBlasInvalidContext'
      case(clAmdBlasInvalidMemObject);           errcode = 'clAmdBlasInvalidMemObject'
      case(clAmdBlasInvalidDevice);              errcode = 'clAmdBlasInvalidDevice'
      case(clAmdBlasInvalidEventWaitList);       errcode = 'clAmdBlasInvalidEventWaitList'
      case(clAmdBlasOutOfResources);             errcode = 'clAmdBlasOutOfResources'
      case(clAmdBlasOutOfHostMemory);            errcode = 'clAmdBlasOutOfHostMemory'
      case(clAmdBlasInvalidOperation);           errcode = 'clAmdBlasInvalidOperation'
      case(clAmdBlasCompilerNotAvailable);       errcode = 'clAmdBlasCompilerNotAvailable'
      case(clAmdBlasBuildProgramFailure );       errcode = 'clAmdBlasBuildProgramFailure'
      case(clAmdBlasNotImplemented);             errcode = 'clAmdBlasNotImplemented'
      case(clAmdBlasNotInitialized);             errcode = 'clAmdBlasNotInitialized'
      case(clAmdBlasInvalidMatA);                errcode = 'clAmdBlasInvalidMatA'
      case(clAmdBlasInvalidMatB);                errcode = 'clAmdBlasInvalidMatB'
      case(clAmdBlasInvalidMatC);                errcode = 'clAmdBlasInvalidMatC'
      case(clAmdBlasInvalidVecX);                errcode = 'clAmdBlasInvalidVecX'
      case(clAmdBlasInvalidVecY);                errcode = 'clAmdBlasInvalidVecY'
      case(clAmdBlasInvalidDim);                 errcode = 'clAmdBlasInvalidDim'
      case(clAmdBlasInvalidLeadDimA);            errcode = 'clAmdBlasInvalidLeadDimA'
      case(clAmdBlasInvalidLeadDimB);            errcode = 'clAmdBlasInvalidLeadDimB'
      case(clAmdBlasInvalidLeadDimC);            errcode = 'clAmdBlasInvalidLeadDimC'
      case(clAmdBlasInvalidIncX);                errcode = 'clAmdBlasInvalidIncX'
      case(clAmdBlasInvalidIncY);                errcode = 'clAmdBlasInvalidIncY'
      case(clAmdBlasInsufficientMemMatA);        errcode = 'clAmdBlasInsufficientMemMatA'
      case(clAmdBlasInsufficientMemMatB);        errcode = 'clAmdBlasInsufficientMemMatB'
      case(clAmdBlasInsufficientMemMatC);        errcode = 'clAmdBlasInsufficientMemMatC'
      case(clAmdBlasInsufficientMemVecX);        errcode = 'clAmdBlasInsufficientMemVecX'
      case(clAmdBlasInsufficientMemVecY);        errcode = 'clAmdBlasInsufficientMemVecY'
      case default
        write(errcode, '(i10)') ierr
        errcode = 'UNKNOWN ERROR CODE ('//trim(adjustl(errcode))//')'
      end select
#endif

      message(1) = 'clAmdBlas '//trim(name)//' '//trim(errcode)
      call messages_fatal(1)
  
      POP_SUB(clblas_print_error)
    end subroutine clblas_print_error

    ! ----------------------------------------------------
    subroutine clfft_print_error(ierr, name)
      integer,          intent(in) :: ierr
      character(len=*), intent(in) :: name

      character(len=40) :: errcode

      PUSH_SUB(clfft_print_error)
#ifdef HAVE_CLAMDFFT
      select case(ierr)
      case(CLFFT_INVALID_GLOBAL_WORK_SIZE);          errcode = 'CLFFT_INVALID_GLOBAL_WORK_SIZE' 
      case(CLFFT_INVALID_MIP_LEVEL);                 errcode = 'CLFFT_INVALID_MIP_LEVEL' 
      case(CLFFT_INVALID_BUFFER_SIZE);               errcode = 'CLFFT_INVALID_BUFFER_SIZE' 
      case(CLFFT_INVALID_GL_OBJECT);                 errcode = 'CLFFT_INVALID_GL_OBJECT' 
      case(CLFFT_INVALID_OPERATION);                 errcode = 'CLFFT_INVALID_OPERATION' 
      case(CLFFT_INVALID_EVENT);                     errcode = 'CLFFT_INVALID_EVENT' 
      case(CLFFT_INVALID_EVENT_WAIT_LIST);           errcode = 'CLFFT_INVALID_EVENT_WAIT_LIST' 
      case(CLFFT_INVALID_GLOBAL_OFFSET);             errcode = 'CLFFT_INVALID_GLOBAL_OFFSET' 
      case(CLFFT_INVALID_WORK_ITEM_SIZE);            errcode = 'CLFFT_INVALID_WORK_ITEM_SIZE' 
      case(CLFFT_INVALID_WORK_GROUP_SIZE);           errcode = 'CLFFT_INVALID_WORK_GROUP_SIZE' 
      case(CLFFT_INVALID_WORK_DIMENSION);            errcode = 'CLFFT_INVALID_WORK_DIMENSION' 
      case(CLFFT_INVALID_KERNEL_ARGS);               errcode = 'CLFFT_INVALID_KERNEL_ARGS' 
      case(CLFFT_INVALID_ARG_SIZE);                  errcode = 'CLFFT_INVALID_ARG_SIZE' 
      case(CLFFT_INVALID_ARG_VALUE);                 errcode = 'CLFFT_INVALID_ARG_VALUE' 
      case(CLFFT_INVALID_ARG_INDEX);                 errcode = 'CLFFT_INVALID_ARG_INDEX' 
      case(CLFFT_INVALID_KERNEL);                    errcode = 'CLFFT_INVALID_KERNEL' 
      case(CLFFT_INVALID_KERNEL_DEFINITION);         errcode = 'CLFFT_INVALID_KERNEL_DEFINITION' 
      case(CLFFT_INVALID_KERNEL_NAME);               errcode = 'CLFFT_INVALID_KERNEL_NAME' 
      case(CLFFT_INVALID_PROGRAM_EXECUTABLE);        errcode = 'CLFFT_INVALID_PROGRAM_EXECUTABLE' 
      case(CLFFT_INVALID_PROGRAM);                   errcode = 'CLFFT_INVALID_PROGRAM' 
      case(CLFFT_INVALID_BUILD_OPTIONS);             errcode = 'CLFFT_INVALID_BUILD_OPTIONS' 
      case(CLFFT_INVALID_BINARY);                    errcode = 'CLFFT_INVALID_BINARY' 
      case(CLFFT_INVALID_SAMPLER);                   errcode = 'CLFFT_INVALID_SAMPLER' 
      case(CLFFT_INVALID_IMAGE_SIZE);                errcode = 'CLFFT_INVALID_IMAGE_SIZE' 
      case(CLFFT_INVALID_IMAGE_FORMAT_DESCRIPTOR);   errcode = 'CLFFT_INVALID_IMAGE_FORMAT_DESCRIPTOR' 
      case(CLFFT_INVALID_MEM_OBJECT);                errcode = 'CLFFT_INVALID_MEM_OBJECT' 
      case(CLFFT_INVALID_HOST_PTR);                  errcode = 'CLFFT_INVALID_HOST_PTR' 
      case(CLFFT_INVALID_COMMAND_QUEUE);             errcode = 'CLFFT_INVALID_COMMAND_QUEUE' 
      case(CLFFT_INVALID_QUEUE_PROPERTIES);          errcode = 'CLFFT_INVALID_QUEUE_PROPERTIES' 
      case(CLFFT_INVALID_CONTEXT);                   errcode = 'CLFFT_INVALID_CONTEXT' 
      case(CLFFT_INVALID_DEVICE);                    errcode = 'CLFFT_INVALID_DEVICE' 
      case(CLFFT_INVALID_PLATFORM);                  errcode = 'CLFFT_INVALID_PLATFORM' 
      case(CLFFT_INVALID_DEVICE_TYPE);               errcode = 'CLFFT_INVALID_DEVICE_TYPE' 
      case(CLFFT_INVALID_VALUE);                     errcode = 'CLFFT_INVALID_VALUE' 
      case(CLFFT_MAP_FAILURE);                       errcode = 'CLFFT_MAP_FAILURE' 
      case(CLFFT_BUILD_PROGRAM_FAILURE);             errcode = 'CLFFT_BUILD_PROGRAM_FAILURE' 
      case(CLFFT_IMAGE_FORMAT_NOT_SUPPORTED);        errcode = 'CLFFT_IMAGE_FORMAT_NOT_SUPPORTED' 
      case(CLFFT_IMAGE_FORMAT_MISMATCH);             errcode = 'CLFFT_IMAGE_FORMAT_MISMATCH' 
      case(CLFFT_MEM_COPY_OVERLAP);                  errcode = 'CLFFT_MEM_COPY_OVERLAP' 
      case(CLFFT_PROFILING_INFO_NOT_AVAILABLE);      errcode = 'CLFFT_PROFILING_INFO_NOT_AVAILABLE' 
      case(CLFFT_OUT_OF_HOST_MEMORY);                errcode = 'CLFFT_OUT_OF_HOST_MEMORY' 
      case(CLFFT_OUT_OF_RESOURCES);                  errcode = 'CLFFT_OUT_OF_RESOURCES' 
      case(CLFFT_MEM_OBJECT_ALLOCATION_FAILURE);     errcode = 'CLFFT_MEM_OBJECT_ALLOCATION_FAILURE' 
      case(CLFFT_COMPILER_NOT_AVAILABLE);            errcode = 'CLFFT_COMPILER_NOT_AVAILABLE' 
      case(CLFFT_DEVICE_NOT_AVAILABLE);              errcode = 'CLFFT_DEVICE_NOT_AVAILABLE' 
      case(CLFFT_DEVICE_NOT_FOUND);                  errcode = 'CLFFT_DEVICE_NOT_FOUND' 
      case(CLFFT_SUCCESS);                           errcode = 'CLFFT_SUCCESS' 
      case(CLFFT_BUGCHECK);                          errcode = 'CLFFT_BUGCHECK' 
      case(CLFFT_NOTIMPLEMENTED);                    errcode = 'CLFFT_NOTIMPLEMENTED' 
      case(CLFFT_FILE_NOT_FOUND);                    errcode = 'CLFFT_FILE_NOT_FOUND' 
      case(CLFFT_FILE_CREATE_FAILURE);               errcode = 'CLFFT_FILE_CREATE_FAILURE' 
      case(CLFFT_VERSION_MISMATCH);                  errcode = 'CLFFT_VERSION_MISMATCH' 
      case(CLFFT_INVALID_PLAN);                      errcode = 'CLFFT_INVALID_PLAN'
      case(CLFFT_DEVICE_NO_DOUBLE);                  errcode = 'CLFFT_DEVICE_NO_DOUBLE' 
      case(CLFFT_ENDSTATUS);                         errcode = 'CLFFT_ENDSTATUS' 
      case default
        write(errcode, '(i10)') ierr
        errcode = 'UNKNOWN ERROR CODE ('//trim(adjustl(errcode))//')'
      end select
#endif

      message(1) = 'clAmdFft '//trim(name)//' '//trim(errcode)
      call messages_fatal(1)

      POP_SUB(clfft_print_error)
    end subroutine clfft_print_error

    ! ----------------------------------------------------

    logical function f90_cl_device_has_extension(device, extension) result(has)
      type(cl_device_id), intent(inout) :: device
      character(len=*),   intent(in)    :: extension

      integer :: cl_status
      character(len=2048) :: all_extensions

      call clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, all_extensions, cl_status)

      has = index(all_extensions, extension) /= 0

    end function f90_cl_device_has_extension

    ! ---------------------------------------------------------
    
    integer pure function opencl_pad(size, blk) result(pad)
      integer, intent(in) :: size
      integer, intent(in) :: blk
      
      integer :: mm
      
      mm = mod(size, blk)
      if(mm == 0) then
        pad = size
      else
        pad = size + blk - mm
      end if
    end function opencl_pad
    
    ! ----------------------------------------------------
    
    subroutine opencl_set_buffer_to_zero(buffer, type, nval)
      type(opencl_mem_t), intent(inout) :: buffer
      type(type_t),       intent(in)    :: type
      integer,            intent(in)    :: nval

      integer :: nval_real, bsize
      
      PUSH_SUB(opencl_set_buffer_to_zero)

      ASSERT(type == TYPE_CMPLX .or. type == TYPE_FLOAT)
      
      nval_real = nval*types_get_size(type)/8

      call opencl_set_kernel_arg(set_zero, 0, nval_real)
      call opencl_set_kernel_arg(set_zero, 1, buffer)

      bsize = opencl_kernel_workgroup_size(set_zero)

      call opencl_kernel_run(set_zero, (/ opencl_pad(nval_real, bsize) /), (/ bsize /))
      call opencl_finish()
      
      POP_SUB(opencl_set_buffer_to_zero)
    end subroutine opencl_set_buffer_to_zero

    ! ----------------------------------------------------
    
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

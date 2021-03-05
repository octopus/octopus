!! Copyright (C) 2010-2016 X. Andrade
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

#if defined(HAVE_OPENCL) && defined(HAVE_CUDA)
#error "Cannot compile with OpenCL and Cuda support at the same time"
#endif

#if defined(HAVE_OPENCL) || defined(HAVE_CUDA)
#define HAVE_ACCEL 1
#endif

module accel_oct_m
  use alloc_cache_oct_m
#ifdef HAVE_OPENCL
  use cl
#endif
#ifdef HAVE_CLBLAS
  use clblas
#endif
  use cuda_oct_m
#ifdef HAVE_CLFFT
  use clfft
#endif
  use global_oct_m
  use iso_c_binding
  use loct_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use types_oct_m
  use parser_oct_m
  use profiling_oct_m
  use unit_system_oct_m

  implicit none 

  private
  
  public ::                       &
    accel_context_t,              &
    accel_device_t,               &
    accel_mem_t,                  &
    accel_kernel_t,               &
    accel_t,                      &
    accel_is_enabled,             &
    accel_allow_CPU_only,         &
    accel_init,                   &
    accel_end,                    &
    accel_padded_size,            &
    accel_mem_nullify,            &
    accel_kernel_start_call,      &
    accel_kernel_build,           &
    accel_create_buffer,          &
    accel_write_buffer,           &
    accel_read_buffer,            &
    accel_release_buffer,         &
    accel_buffer_is_allocated,    &
    accel_finish,                 &
    accel_set_kernel_arg,         &
    accel_max_workgroup_size,     &
    accel_kernel_workgroup_size,  &
    accel_kernel_run,             &
    accel_set_buffer_to_zero,     &
    accel_use_shared_mem,         &
    clblas_print_error,           &
    clfft_print_error,            &
    accel_local_memory_size,      &
    accel_global_memory_size,     &
    accel_max_size_per_dim,       &
    accel_get_device_pointer,     &
    accel_set_stream,             &
    accel_synchronize_all_streams
  
#ifdef HAVE_OPENCL
  integer, public, parameter ::                 &
    ACCEL_MEM_READ_ONLY  = CL_MEM_READ_ONLY,    &
    ACCEL_MEM_READ_WRITE = CL_MEM_READ_WRITE,   &
    ACCEL_MEM_WRITE_ONLY = CL_MEM_WRITE_ONLY
#else
  integer, public, parameter ::                 &
    ACCEL_MEM_READ_ONLY  = 0,                   &
    ACCEL_MEM_READ_WRITE = 1,                   &
    ACCEL_MEM_WRITE_ONLY = 2
#endif

  type accel_context_t
    ! Components are public by default
#ifdef HAVE_OPENCL
    type(cl_context) :: cl_context
#elif defined(HAVE_CUDA)
    type(c_ptr)      :: cuda_context
#else
    integer          :: dummy
#endif
  end type accel_context_t

  type accel_device_t
    ! Components are public by default
#ifdef HAVE_OPENCL
    type(cl_device_id) :: cl_device
#elif defined(HAVE_CUDA)
    type(c_ptr)      :: cuda_device
#else
    integer         :: dummy
#endif
  end type accel_device_t

  type accel_t
    ! Components are public by default
    type(accel_context_t)  :: context
    type(accel_device_t)   :: device
#ifdef HAVE_OPENCL
    type(cl_command_queue) :: command_queue
#endif
    type(c_ptr)            :: cublas_handle
    type(c_ptr)            :: cuda_stream
    type(c_ptr)            :: module_map
    integer                :: max_workgroup_size
    integer(8)             :: local_memory_size
    integer(8)             :: global_memory_size
    logical                :: enabled
    logical                :: allow_CPU_only
    logical                :: shared_mem
    logical                :: cuda_mpi
    integer                :: warp_size
  end type accel_t

  type accel_mem_t
    ! Components are public by default
#ifdef HAVE_OPENCL
    type(cl_mem)           :: mem
#else
    type(c_ptr)            :: mem
#endif
    integer(SIZEOF_SIZE_T) :: size
    type(type_t)           :: type
    integer                :: flags
    logical                :: allocated
  end type accel_mem_t

  type accel_kernel_t
    ! Components are public by default
#ifdef HAVE_OPENCL
    type(cl_kernel)               :: kernel
#endif
#ifdef HAVE_CUDA
    type(c_ptr)                   :: cuda_kernel
    type(c_ptr)                   :: cuda_module
    type(c_ptr)                   :: arguments
#endif
    integer(8)                    :: cuda_shared_mem
    logical                       :: initialized = .false.
    type(accel_kernel_t), pointer :: next
    integer                       :: arg_count
  end type accel_kernel_t

  type(accel_t), public :: accel

  ! the kernels
  type(accel_kernel_t), public, target, save :: kernel_vpsi
  type(accel_kernel_t), public, target, save :: kernel_vpsi_spinors
  type(accel_kernel_t), public, target, save :: kernel_daxpy
  type(accel_kernel_t), public, target, save :: kernel_zaxpy
  type(accel_kernel_t), public, target, save :: kernel_copy
  type(accel_kernel_t), public, target, save :: dpack
  type(accel_kernel_t), public, target, save :: zpack
  type(accel_kernel_t), public, target, save :: dunpack
  type(accel_kernel_t), public, target, save :: zunpack
  type(accel_kernel_t), public, target, save :: kernel_ghost_reorder
  type(accel_kernel_t), public, target, save :: kernel_density_real
  type(accel_kernel_t), public, target, save :: kernel_density_complex
  type(accel_kernel_t), public, target, save :: kernel_density_spinors
  type(accel_kernel_t), public, target, save :: kernel_phase
  type(accel_kernel_t), public, target, save :: kernel_phase_spiral
  type(accel_kernel_t), public, target, save :: dkernel_dot_matrix
  type(accel_kernel_t), public, target, save :: zkernel_dot_matrix
  type(accel_kernel_t), public, target, save :: zkernel_dot_matrix_spinors
  type(accel_kernel_t), public, target, save :: dkernel_batch_axpy
  type(accel_kernel_t), public, target, save :: zkernel_batch_axpy
  type(accel_kernel_t), public, target, save :: dkernel_batch_dotp
  type(accel_kernel_t), public, target, save :: zkernel_batch_dotp
  type(accel_kernel_t), public, target, save :: dzmul
  type(accel_kernel_t), public, target, save :: zzmul
  type(accel_kernel_t), public, target, save :: set_one

  ! kernels used locally
  type(accel_kernel_t), save :: set_zero

  interface accel_create_buffer
    module procedure accel_create_buffer_4, accel_create_buffer_8
  end interface accel_create_buffer

  interface accel_write_buffer
    module procedure iaccel_write_buffer_0, daccel_write_buffer_0, zaccel_write_buffer_0
    module procedure iaccel_write_buffer_1, daccel_write_buffer_1, zaccel_write_buffer_1
    module procedure iaccel_write_buffer_2, daccel_write_buffer_2, zaccel_write_buffer_2
    module procedure iaccel_write_buffer_3, daccel_write_buffer_3, zaccel_write_buffer_3
  end interface accel_write_buffer

  interface accel_read_buffer
    module procedure iaccel_read_buffer_1, daccel_read_buffer_1, zaccel_read_buffer_1
    module procedure iaccel_read_buffer_2, daccel_read_buffer_2, zaccel_read_buffer_2
    module procedure iaccel_read_buffer_3, daccel_read_buffer_3, zaccel_read_buffer_3
  end interface accel_read_buffer

  interface accel_set_kernel_arg
    module procedure                       &
      accel_set_kernel_arg_buffer,  &
      iaccel_set_kernel_arg_data,   &
      daccel_set_kernel_arg_data,   &
      zaccel_set_kernel_arg_data,   &
      accel_set_kernel_arg_local
  end interface accel_set_kernel_arg

  interface accel_get_device_pointer
    module procedure iaccel_get_device_pointer_1
    module procedure iaccel_get_device_pointer_2
    module procedure daccel_get_device_pointer_1, zaccel_get_device_pointer_1
    module procedure daccel_get_device_pointer_2, zaccel_get_device_pointer_2
  end interface accel_get_device_pointer

  type(profile_t), save :: prof_read, prof_write

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

  integer :: buffer_alloc_count
  integer(8) :: allocated_mem
  type(accel_kernel_t), pointer :: head
  type(alloc_cache_t) :: memcache
  
contains

  pure logical function accel_is_enabled() result(enabled)
#ifdef HAVE_ACCEL
    enabled = accel%enabled
#else
    enabled = .false.
#endif
  end function accel_is_enabled

  ! ------------------------------------------

  pure logical function accel_allow_CPU_only() result(allow)
#ifdef HAVE_ACCEL
    allow = accel%allow_CPU_only
#else
    allow = .true.
#endif
  end function accel_allow_CPU_only

  ! ------------------------------------------

  subroutine accel_init(base_grp, namespace)
    type(mpi_grp_t),     intent(inout) :: base_grp
    type(namespace_t),   intent(in)    :: namespace
    
    logical  :: disable, default, run_benchmark
    integer  :: idevice, iplatform
#ifdef HAVE_OPENCL
    integer  :: device_type
    integer :: cl_status, idev
    integer  :: ndevices, ret_devices, nplatforms, iplat
    character(len=256) :: device_name
    type(cl_platform_id) :: platform_id
    type(cl_program) :: prog
    type(cl_platform_id), allocatable :: allplatforms(:)
    type(cl_device_id), allocatable :: alldevices(:)
    type(profile_t), save :: prof_init
#endif

    PUSH_SUB(accel_init)

    buffer_alloc_count = 0

    !%Variable DisableAccel    
    !%Type logical
    !%Default yes
    !%Section Execution::Accel
    !%Description
    !% If Octopus was compiled with OpenCL or CUDA support, it will
    !% try to initialize and use an accelerator device. By setting this
    !% variable to <tt>yes</tt> you force Octopus not to use an accelerator even it is available.
    !%End
    call messages_obsolete_variable(namespace, 'DisableOpenCL', 'DisableAccel')
#ifdef HAVE_ACCEL
    default = .false.
#else
    default = .true.
#endif
    call parse_variable(namespace, 'DisableAccel', default, disable)
    accel%enabled = .not. disable
    
#ifndef HAVE_ACCEL
    if(accel%enabled) then
      message(1) = 'Octopus was compiled without OpenCL or Cuda support.'
      call messages_fatal(1)
    end if
#endif

    if(.not. accel_is_enabled()) then
      POP_SUB(accel_init)
      return
    end if

    !%Variable AccelPlatform
    !%Type integer
    !%Default 0
    !%Section Execution::Accel
    !%Description
    !% This variable selects the OpenCL platform that Octopus will
    !% use. You can give an explicit platform number or use one of
    !% the options that select a particular vendor
    !% implementation. Platform 0 is used by default.
    !%
    !% This variable has no effect for CUDA.
    !%Option amd -2
    !% Use the AMD OpenCL platform.
    !%Option nvidia -3
    !% Use the Nvidia OpenCL platform.
    !%Option ati -4
    !% Use the ATI (old AMD) OpenCL platform.
    !%Option intel -5
    !% Use the Intel OpenCL platform.
    !%End
    call parse_variable(namespace, 'AccelPlatform', 0, iplatform)

    call messages_obsolete_variable(namespace, 'OpenCLPlatform', 'AccelPlatform')
    
    !%Variable AccelDevice
    !%Type integer
    !%Default gpu
    !%Section Execution::Accel
    !%Description
    !% This variable selects the OpenCL or CUDA accelerator device
    !% that Octopus will use. You can specify one of the options below
    !% or a numerical id to select a specific device.
    !%
    !% Values >= 0 select the device to be used. In case of MPI enabled runs
    !% devices are distributed in a round robin fashion, starting at this value.
    !%Option gpu -1
    !% If available, Octopus will use a GPU.
    !%Option cpu -2
    !% If available, Octopus will use a CPU (only for OpenCL).
    !%Option accelerator -3
    !% If available, Octopus will use an accelerator (only for OpenCL).
    !%Option accel_default -4
    !% Octopus will use the default device specified by the implementation.
    !% implementation.
    !%End
    call parse_variable(namespace, 'AccelDevice', OPENCL_GPU, idevice)

    call messages_obsolete_variable(namespace, 'OpenCLDevice', 'AccelDevice')
    
    if(idevice < OPENCL_DEFAULT) then
      call messages_write('Invalid AccelDevice')
      call messages_fatal()
    end if

    call messages_print_stress(stdout, "GPU acceleration")

#ifdef HAVE_CUDA
    if(idevice<0) idevice = 0
    call cuda_init(accel%context%cuda_context, accel%device%cuda_device, accel%cuda_stream, &
      idevice, base_grp%rank)
#ifdef HAVE_MPI
    write(message(1), '(A, I5.5, A, I5.5)') "Rank ", base_grp%rank, " uses device number ", idevice
    call messages_info(1, all_nodes = .true.)
#endif

    ! no shared mem support in our cuda interface (for the moment)
    accel%shared_mem = .true.

    call cublas_init(accel%cublas_handle, accel%cuda_stream)
#endif
    
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

    platform_id = allplatforms(iplatform + 1)

    SAFE_DEALLOCATE_A(allplatforms)

    call clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ALL, ndevices, cl_status)

    call messages_write('Info: Available CL devices: ')
    call messages_write(ndevices)
    call messages_info()

    SAFE_ALLOCATE(alldevices(1:ndevices))

    ! list all devices

    call clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ALL, alldevices, ret_devices, cl_status)

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
    call clGetDeviceIDs(platform_id, device_type, alldevices, ret_devices, cl_status)

    if(ret_devices < 1) then
      ! we didnt find a device of the selected type, we ask for the default device
      call clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, alldevices, ret_devices, cl_status)

      if(ret_devices < 1) then
        ! if this does not work, we ask for all devices
        call clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ALL, alldevices, ret_devices, cl_status)
      end if

      if(ret_devices < 1) then
        call messages_write('Cannot find an OpenCL device')
        call messages_fatal()
      end if
    end if

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

    accel%device%cl_device = alldevices(idevice + 1)

    ! create the context
    accel%context%cl_context = clCreateContext(platform_id, accel%device%cl_device, cl_status)
    if(cl_status /= CL_SUCCESS) call opencl_print_error(cl_status, "CreateContext")

    SAFE_DEALLOCATE_A(alldevices)

    accel%command_queue = clCreateCommandQueue(accel%context%cl_context, accel%device%cl_device, &
      CL_QUEUE_PROFILING_ENABLE, cl_status)
    if(cl_status /= CL_SUCCESS) call opencl_print_error(cl_status, "CreateCommandQueue")

    call clGetDeviceInfo(accel%device%cl_device, CL_DEVICE_TYPE, device_type, cl_status)

    select case(device_type)
    case(CL_DEVICE_TYPE_GPU)
      accel%shared_mem = .true.
    case(CL_DEVICE_TYPE_CPU, CL_DEVICE_TYPE_ACCELERATOR)
      accel%shared_mem = .false.
    case default
      accel%shared_mem = .false.
    end select

#ifdef HAVE_CLBLAS
    call clblasSetup(cl_status)
    if(cl_status /= clblasSuccess) call clblas_print_error(cl_status, 'clblasSetup')
#endif

#ifdef HAVE_CLFFT
    call clfftSetup(cl_status)
    if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clfftSetup')
#endif

    call profiling_out(prof_init)
#endif

    ! Get some device information that we will need later
    
    ! total memory
#ifdef HAVE_OPENCL
    call clGetDeviceInfo(accel%device%cl_device, CL_DEVICE_GLOBAL_MEM_SIZE, accel%global_memory_size, cl_status)
    call clGetDeviceInfo(accel%device%cl_device, CL_DEVICE_LOCAL_MEM_SIZE, accel%local_memory_size, cl_status)
    call clGetDeviceInfo(accel%device%cl_device, CL_DEVICE_MAX_WORK_GROUP_SIZE, accel%max_workgroup_size, cl_status)
    accel%warp_size = 1
#endif
#ifdef HAVE_CUDA
    call cuda_device_total_memory(accel%device%cuda_device, accel%global_memory_size)
    call cuda_device_shared_memory(accel%device%cuda_device, accel%local_memory_size)
    call cuda_device_max_threads_per_block(accel%device%cuda_device, accel%max_workgroup_size)
    call cuda_device_get_warpsize(accel%device%cuda_device, accel%warp_size)
#endif
      
    if(mpi_grp_is_root(base_grp)) call device_info()

    ! initialize the cache used to speed up allocations
    call alloc_cache_init(memcache, nint(CNST(0.25)*accel%global_memory_size, 8))
    
    ! now initialize the kernels
    call accel_kernel_global_init()

    call accel_kernel_start_call(set_zero, 'set_zero.cl', "set_zero")
    call accel_kernel_start_call(set_one, 'set_one.cl', "set_one")
    call accel_kernel_start_call(kernel_vpsi, 'vpsi.cl', "vpsi")
    call accel_kernel_start_call(kernel_vpsi_spinors, 'vpsi.cl', "vpsi_spinors")
    call accel_kernel_start_call(kernel_daxpy, 'axpy.cl', "daxpy", flags = '-DRTYPE_DOUBLE')
    call accel_kernel_start_call(kernel_zaxpy, 'axpy.cl', "zaxpy", flags = '-DRTYPE_COMPLEX')
    call accel_kernel_start_call(dkernel_batch_axpy, 'axpy.cl', "dbatch_axpy_function", flags = '-lineinfo -DRTYPE_DOUBLE')
    call accel_kernel_start_call(zkernel_batch_axpy, 'axpy.cl', "zbatch_axpy_function", flags = '-lineinfo -DRTYPE_COMPLEX')
    call accel_kernel_start_call(dkernel_batch_dotp, 'mesh_batch_single.cl', "dbatch_mf_dotp", flags = '-lineinfo')
    call accel_kernel_start_call(zkernel_batch_dotp, 'mesh_batch_single.cl', "zbatch_mf_dotp", flags = '-lineinfo')
    call accel_kernel_start_call(dpack, 'pack.cl', "dpack")
    call accel_kernel_start_call(zpack, 'pack.cl', "zpack")
    call accel_kernel_start_call(dunpack, 'pack.cl', "dunpack")
    call accel_kernel_start_call(zunpack, 'pack.cl', "zunpack")
    call accel_kernel_start_call(kernel_copy, 'copy.cl', "copy")
    call accel_kernel_start_call(kernel_ghost_reorder, 'ghost.cl', "ghost_reorder")
    call accel_kernel_start_call(kernel_density_real, 'density.cl', "density_real")
    call accel_kernel_start_call(kernel_density_complex, 'density.cl', "density_complex")
    call accel_kernel_start_call(kernel_density_spinors, 'density.cl', "density_spinors")
    call accel_kernel_start_call(kernel_phase, 'phase.cl', "phase")
    call accel_kernel_start_call(dkernel_dot_matrix, 'mesh_batch.cl', "ddot_matrix")
    call accel_kernel_start_call(zkernel_dot_matrix, 'mesh_batch.cl', "zdot_matrix")
    call accel_kernel_start_call(zkernel_dot_matrix_spinors, 'mesh_batch.cl', "zdot_matrix_spinors")

    
    call accel_kernel_start_call(dzmul, 'mul.cl', "dzmul", flags = '-DRTYPE_DOUBLE')
    call accel_kernel_start_call(zzmul, 'mul.cl', "zzmul", flags = '-DRTYPE_COMPLEX')

    !%Variable AccelBenchmark
    !%Type logical
    !%Default no
    !%Section Execution::Accel
    !%Description
    !% If this variable is set to yes, Octopus will run some
    !% routines to benchmark the performance of the accelerator device.
    !%End
    call parse_variable(namespace, 'AccelBenchmark', .false., run_benchmark)

    call messages_obsolete_variable(namespace, 'OpenCLBenchmark', 'AccelBenchmark')
    
    if(run_benchmark) then
      call opencl_check_bandwidth()
    end if

    !%Variable CudaAwareMPI
    !%Type logical
    !%Section Execution::Accel
    !%Description
    !% If Octopus was compiled with CUDA support and MPI support and if the MPI
    !% implementation is CUDA-aware (i.e., it supports communication using device pointers),
    !% this switch can be set to true to use the CUDA-aware MPI features. The advantage
    !% of this approach is that it can do, e.g., peer-to-peer copies between devices without
    !% going through the host memmory.
    !% The default is false, except when the configure switch --enable-cudampi is set, in which
    !% case this variable is set to true.
    !%End
#ifdef HAVE_CUDA_MPI
    default = .true.
#else
    default = .false.
#endif
    call parse_variable(namespace, 'CudaAwareMPI', default, accel%cuda_mpi)
    if(accel%cuda_mpi) then
      call messages_write("Using CUDA-aware MPI.")
      call messages_info()
    end if


    !%Variable AllowCPUonly
    !%Type logical
    !%Section Execution::Accel
    !%Description
    !% In order to prevent waste of resources, the code will normally stop when the GPU is disabled due to 
    !% incomplete implementations or incompatibilities. AllowCPUonly = yes overrides this and allows the 
    !% code execution also in these cases.
    !%End
#if defined (HAVE_ACCEL)
    default = .false.
#else
    default = .true.
#endif
    call parse_variable(namespace, 'AllowCPUonly', default, accel%allow_CPU_only)



    call messages_print_stress(stdout)

    POP_SUB(accel_init)

  contains

#if defined(HAVE_MPI) && defined(HAVE_OPENCL)
    subroutine select_device(idevice)
      integer, intent(inout) :: idevice
      integer :: irank
      character(len=256) :: device_name

      PUSH_SUB(accel_init.select_device)

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

      POP_SUB(accel_init.select_device)
    end subroutine select_device
#endif

    subroutine device_info()
#ifdef HAVE_OPENCL
      integer(8) :: val
#endif
#ifdef HAVE_CUDA
      integer :: version
#endif
      integer :: major, minor
      character(len=256) :: val_str
      
      PUSH_SUB(accel_init.device_info)

      call messages_new_line()
      call messages_write('Selected device:')
      call messages_new_line()

#ifdef HAVE_OPENCL
      call messages_write('      Framework              : OpenCL')
#endif
#ifdef HAVE_CUDA
      call messages_write('      Framework              : CUDA')
#endif
      call messages_info()

#ifdef HAVE_CUDA
      call messages_write('      Device type            : GPU', new_line = .true.)
      call messages_write('      Device vendor          : NVIDIA Corporation', new_line = .true.)
#endif

#ifdef HAVE_OPENCL
      call clGetDeviceInfo(accel%device%cl_device, CL_DEVICE_TYPE, val, cl_status)
      call messages_write('      Device type            :')
      select case(int(val, 4))
      case(CL_DEVICE_TYPE_GPU)
        call messages_write(' GPU')
      case(CL_DEVICE_TYPE_CPU)
        call messages_write(' CPU')
      case(CL_DEVICE_TYPE_ACCELERATOR)
        call messages_write(' accelerator')
      end select
      call messages_new_line()

      call clGetDeviceInfo(accel%device%cl_device, CL_DEVICE_VENDOR, val_str, cl_status)
      call messages_write('      Device vendor          : '//trim(val_str))
      call messages_new_line()
#endif
      
#ifdef HAVE_OPENCL
      call clGetDeviceInfo(accel%device%cl_device, CL_DEVICE_NAME, val_str, cl_status)
#endif
#ifdef HAVE_CUDA
      call cuda_device_name(accel%device%cuda_device, val_str)
#endif
      call messages_write('      Device name            : '//trim(val_str))
      call messages_new_line()
      
#ifdef HAVE_CUDA
      call cuda_device_capability(accel%device%cuda_device, major, minor)
#endif
      call messages_write('      Cuda capabilities      :')
      call messages_write(major, fmt = '(i2)')
      call messages_write('.')
      call messages_write(minor, fmt = '(i1)')
      call messages_new_line()

      ! VERSION
#ifdef HAVE_OPENCL
      call clGetDeviceInfo(accel%device%cl_device, CL_DRIVER_VERSION, val_str, cl_status)
      call messages_write('      Driver version         : '//trim(val_str))
#endif
#ifdef HAVE_CUDA
      call cuda_driver_version(version)
      call messages_write('      Driver version         : ')
      call messages_write(version)
#endif
      call messages_new_line()

      
#ifdef HAVE_OPENCL
      call clGetDeviceInfo(accel%device%cl_device, CL_DEVICE_MAX_COMPUTE_UNITS, val, cl_status)
      call messages_write('      Compute units          :')
      call messages_write(val)
      call messages_new_line()

      call clGetDeviceInfo(accel%device%cl_device, CL_DEVICE_MAX_CLOCK_FREQUENCY, val, cl_status)
      call messages_write('      Clock frequency        :')
      call messages_write(val)
      call messages_write(' GHz')
      call messages_new_line()
#endif

      call messages_write('      Device memory          :')
      call messages_write(accel%global_memory_size, units=unit_megabytes)
      call messages_new_line()

      call messages_write('      Local/shared memory    :')
      call messages_write(accel%local_memory_size, units=unit_kilobytes)
      call messages_new_line()
      
    
#ifdef HAVE_OPENCL
      call clGetDeviceInfo(accel%device%cl_device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, val, cl_status)
      call messages_write('      Max alloc size         :')
      call messages_write(val, units = unit_megabytes)
      call messages_new_line()

      call clGetDeviceInfo(accel%device%cl_device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, val, cl_status)
      call messages_write('      Device cache           :')
      call messages_write(val, units = unit_kilobytes)
      call messages_new_line()

      call clGetDeviceInfo(accel%device%cl_device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, val, cl_status)
      call messages_write('      Constant memory        :')
      call messages_write(val, units = unit_kilobytes)
      call messages_new_line()
#endif

      call messages_write('      Max. group/block size  :')
      call messages_write(accel%max_workgroup_size)
      call messages_new_line()
      

#ifdef HAVE_OPENCL
      call messages_write('      Extension cl_khr_fp64  :')
      call messages_write(f90_cl_device_has_extension(accel%device%cl_device, "cl_khr_fp64"))
      call messages_new_line()

      call messages_write('      Extension cl_amd_fp64  :')
      call messages_write(f90_cl_device_has_extension(accel%device%cl_device, "cl_amd_fp64"))
      call messages_new_line()
#endif
      
      call messages_info()


      POP_SUB(accel_init.device_info)
    end subroutine device_info

  end subroutine accel_init

  ! ------------------------------------------
#ifdef HAVE_OPENCL
  integer function get_platform_id(platform_name) result(platform_id)
    character(len=*), intent(in) :: platform_name

    platform_id = CL_PLAT_INVALID
    if(index(platform_name, 'AMD') > 0)    platform_id = CL_PLAT_AMD
    if(index(platform_name, 'ATI') > 0)    platform_id = CL_PLAT_ATI
    if(index(platform_name, 'NVIDIA') > 0) platform_id = CL_PLAT_NVIDIA
    if(index(platform_name, 'Intel') > 0)  platform_id = CL_PLAT_INTEL
  end function get_platform_id
#endif
  ! ------------------------------------------

  subroutine accel_end()
#ifdef HAVE_OPENCL
    integer :: ierr
#endif
    integer(8) :: hits, misses
    real(8) :: volume_hits, volume_misses
    logical :: found
    type(accel_mem_t) :: tmp

    PUSH_SUB(accel_end)

    if(accel_is_enabled()) then

      do 
        call alloc_cache_get(memcache, ALLOC_CACHE_ANY_SIZE, found, tmp%mem)
        if(.not. found) exit

#ifdef HAVE_OPENCL
        call clReleaseMemObject(tmp%mem, ierr)
        if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clReleaseMemObject")
#endif
#ifdef HAVE_CUDA
        call cuda_mem_free(tmp%mem)
#endif
      end do

      call alloc_cache_end(memcache, hits, misses, volume_hits, volume_misses)

      call messages_print_stress(stdout, "Acceleration-device allocation cache")

      call messages_new_line()
      call messages_write('    Number of allocations    =')
      call messages_write(hits + misses, new_line = .true.)
      call messages_write('    Volume of allocations    =')
      call messages_write(volume_hits + volume_misses, fmt = 'f18.1', units = unit_gigabytes, align_left = .true., &
        new_line = .true.)
      call messages_write('    Hit ratio                =')
      call messages_write(hits/TOFLOAT(hits + misses)*100, fmt='(f6.1)', align_left = .true.)
      call messages_write('%', new_line = .true.)
      call messages_write('    Volume hit ratio         =')
      call messages_write(volume_hits/(volume_hits + volume_misses)*100, fmt='(f6.1)', align_left = .true.)
      call messages_write('%')
      call messages_new_line()
      call messages_info()

      call messages_print_stress(stdout)
    end if
    
    call accel_kernel_global_end()

#ifdef HAVE_CLBLAS
    call clblasTearDown()
#endif

#ifdef HAVE_CLFFT
    call clfftTearDown()
#endif

    if(accel_is_enabled()) then
#ifdef HAVE_CUDA
      call cublas_end(accel%cublas_handle)
      call cuda_end(accel%context%cuda_context, accel%device%cuda_device)
#endif

#ifdef HAVE_OPENCL
      call clReleaseCommandQueue(accel%command_queue, ierr)

      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "ReleaseCommandQueue")
      call clReleaseContext(accel%context%cl_context, cl_status)
#endif
      
      if(buffer_alloc_count /= 0) then
        call messages_write('Accel:')
        call messages_write(TOFLOAT(allocated_mem), fmt = 'f12.1', units = unit_megabytes, align_left = .true.)
        call messages_write(' in ')
        call messages_write(buffer_alloc_count)
        call messages_write(' buffers were not deallocated.')
        call messages_fatal()
      end if

    end if

    POP_SUB(accel_end)
  end subroutine accel_end

  ! ------------------------------------------

  elemental subroutine accel_mem_nullify(this)
    type(accel_mem_t), intent(out) :: this

    !> To be implemented.
    this%size = 0
    this%flags = 0
    this%allocated = .false.
    
  end subroutine accel_mem_nullify

  ! ------------------------------------------

  integer function accel_padded_size(nn) result(psize)
    integer,        intent(in) :: nn

    integer :: modnn, bsize
    
    psize = nn

    if(accel_is_enabled()) then

      bsize = accel_max_workgroup_size()
      
      psize = nn
      modnn = mod(nn, bsize)
      if(modnn /= 0) psize = psize + bsize - modnn

    end if
    
  end function accel_padded_size

  ! ------------------------------------------

  subroutine accel_create_buffer_4(this, flags, type, size)
    type(accel_mem_t),  intent(inout) :: this
    integer,            intent(in)    :: flags
    type(type_t),       intent(in)    :: type
    integer,            intent(in)    :: size

    call accel_create_buffer_8(this, flags, type, int(size, 8))
  end subroutine accel_create_buffer_4

  ! ------------------------------------------

  subroutine accel_create_buffer_8(this, flags, type, size)
    type(accel_mem_t),  intent(inout) :: this
    integer,            intent(in)    :: flags
    type(type_t),       intent(in)    :: type
    integer(8),         intent(in)    :: size

    integer(8) :: fsize
    logical    :: found
#ifdef HAVE_OPENCL
    integer :: ierr
#endif

    PUSH_SUB(accel_create_buffer_8)

    this%type = type
    this%size = size
    this%flags = flags
    fsize = int(size, 8)*types_get_size(type)
    this%allocated = .true.
    
    if(fsize > 0) then

      call alloc_cache_get(memcache, fsize, found, this%mem)

      if(.not. found) then
#ifdef HAVE_OPENCL
        this%mem = clCreateBuffer(accel%context%cl_context, flags, fsize, ierr)
        if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clCreateBuffer")
#endif
#ifdef HAVE_CUDA
        call cuda_mem_alloc(this%mem, fsize)
#endif
      end if
      
      buffer_alloc_count = buffer_alloc_count + 1
      allocated_mem = allocated_mem + fsize

    end if
      
    POP_SUB(accel_create_buffer_8)
  end subroutine accel_create_buffer_8
  
  ! ------------------------------------------

  subroutine accel_release_buffer(this)
    type(accel_mem_t), intent(inout) :: this

#ifdef HAVE_OPENCL
    integer :: ierr
#endif
    logical :: put
    integer(8) :: fsize

    PUSH_SUB(accel_release_buffer)

    if(this%size > 0) then

      fsize = int(this%size, 8)*types_get_size(this%type)
      
      call alloc_cache_put(memcache, fsize, this%mem, put) 

      if(.not. put) then
#ifdef HAVE_OPENCL
        call clReleaseMemObject(this%mem, ierr)
        if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clReleaseMemObject")
#endif
#ifdef HAVE_CUDA
        call cuda_mem_free(this%mem)
#endif
      end if
      
      buffer_alloc_count = buffer_alloc_count - 1
      allocated_mem = allocated_mem + fsize

    end if
    
    this%size = 0
    this%flags = 0

    this%allocated = .false.
    
    POP_SUB(accel_release_buffer)
  end subroutine accel_release_buffer
    
  ! ------------------------------------------
  
  logical pure function accel_buffer_is_allocated(this) result(allocated)
    type(accel_mem_t), intent(in) :: this

    allocated = this%allocated
  end function accel_buffer_is_allocated
    
  ! -----------------------------------------

  subroutine accel_finish()
#ifdef HAVE_OPENCL
    integer :: ierr
#endif

    ! no push_sub, called too frequently
    
    if(accel_is_enabled()) then
#ifdef HAVE_OPENCL
      call clFinish(accel%command_queue, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, 'clFinish')
#endif
#ifdef HAVE_CUDA
      call cuda_context_synchronize()
#endif
    end if
  end subroutine accel_finish

  ! ------------------------------------------

  subroutine accel_set_kernel_arg_buffer(kernel, narg, buffer)
    type(accel_kernel_t), intent(inout) :: kernel
    integer,              intent(in)    :: narg
    type(accel_mem_t),    intent(in)    :: buffer

#ifdef HAVE_OPENCL
    integer :: ierr
#endif

    ASSERT(accel_buffer_is_allocated(buffer))
    
    ! no push_sub, called too frequently
#ifdef HAVE_OPENCL
    call clSetKernelArg(kernel%kernel, narg, buffer%mem, ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clSetKernelArg_buf")
#endif

#ifdef HAVE_CUDA
    call cuda_kernel_set_arg_buffer(kernel%arguments, buffer%mem, narg)
#endif
   
  end subroutine accel_set_kernel_arg_buffer
  
  ! ------------------------------------------

  subroutine accel_set_kernel_arg_local(kernel, narg, type, size)
    type(accel_kernel_t), intent(inout) :: kernel
    integer,              intent(in)    :: narg
    type(type_t),         intent(in)    :: type
    integer,              intent(in)    :: size

#ifdef HAVE_OPENCL
    integer :: ierr
#endif
    integer(8) :: size_in_bytes

    PUSH_SUB(accel_set_kernel_arg_local)

    
    size_in_bytes = int(size, 8)*types_get_size(type)

    if(size_in_bytes > accel%local_memory_size) then
      write(message(1), '(a,f12.6,a)') "CL Error: requested local memory: ", TOFLOAT(size_in_bytes)/1024.0, " Kb"
      write(message(2), '(a,f12.6,a)') "          available local memory: ", TOFLOAT(accel%local_memory_size)/1024.0, " Kb"
      call messages_fatal(2)
    else if(size_in_bytes <= 0) then
      write(message(1), '(a,i10)') "CL Error: invalid local memory size: ", size_in_bytes
      call messages_fatal(1)
    end if

#ifdef HAVE_CUDA
    kernel%cuda_shared_mem = size_in_bytes
#endif

#ifdef HAVE_OPENCL
    call clSetKernelArgLocal(kernel%kernel, narg, size_in_bytes, ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "set_kernel_arg_local")
#endif

    POP_SUB(accel_set_kernel_arg_local)
  end subroutine accel_set_kernel_arg_local

  ! ------------------------------------------

  subroutine accel_kernel_run(kernel, globalsizes, localsizes)
    type(accel_kernel_t), intent(inout) :: kernel
    integer,              intent(in)    :: globalsizes(:)
    integer,              intent(in)    :: localsizes(:)

    integer :: dim
#ifdef HAVE_OPENCL
    integer :: ierr
#endif
    integer(8) :: gsizes(1:3)
    integer(8) :: lsizes(1:3)

    ! no push_sub, called too frequently

    ! cuda needs all dimensions
    gsizes = 1
    lsizes = 1
    
    dim = ubound(globalsizes, dim = 1)

    ASSERT(dim == ubound(localsizes, dim = 1))

    ! if one size is zero, there is nothing to do
    if(any(globalsizes == 0)) return

    ASSERT(all(localsizes > 0))
    ASSERT(all(localsizes <= accel_max_workgroup_size()))
    ASSERT(all(mod(globalsizes, localsizes) == 0))

    gsizes(1:dim) = int(globalsizes(1:dim), 8)
    lsizes(1:dim) = int(localsizes(1:dim), 8)

#ifdef HAVE_OPENCL
    call clEnqueueNDRangeKernel(accel%command_queue, kernel%kernel, gsizes(1:dim), lsizes(1:dim), ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "EnqueueNDRangeKernel")
#endif

#ifdef HAVE_CUDA
    gsizes(1:3) = gsizes(1:3)/lsizes(1:3)

    ASSERT(gsizes(1) < 2_8**31 - 1_8)
    ASSERT(all(gsizes(2:3) <= 65535_8))
    
    call cuda_launch_kernel(kernel%cuda_kernel, gsizes(1), lsizes(1), kernel%cuda_shared_mem, kernel%arguments)

    kernel%cuda_shared_mem = 0    
#endif
    
  end subroutine accel_kernel_run

  ! -----------------------------------------------

  integer pure function accel_max_workgroup_size() result(max_workgroup_size)
    max_workgroup_size = accel%max_workgroup_size
  end function accel_max_workgroup_size

  ! -----------------------------------------------

  integer function accel_kernel_workgroup_size(kernel) result(workgroup_size)
    type(accel_kernel_t), intent(inout) :: kernel

#ifdef HAVE_OPENCL
    integer(8) :: workgroup_size8
    integer :: ierr
#endif

    workgroup_size = 0

#ifdef HAVE_OPENCL
    call clGetKernelWorkGroupInfo(kernel%kernel, accel%device%cl_device, CL_KERNEL_WORK_GROUP_SIZE, workgroup_size8, ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "EnqueueNDRangeKernel")
    workgroup_size = workgroup_size8
#endif

#ifdef HAVE_CUDA
    workgroup_size = accel%max_workgroup_size
#endif

  end function accel_kernel_workgroup_size

  ! -----------------------------------------------

#ifdef HAVE_OPENCL
  subroutine opencl_build_program(prog, filename, flags)
    type(cl_program),           intent(inout) :: prog
    character(len=*),           intent(in)    :: filename
    character(len=*), optional, intent(in)    :: flags

    character(len = 1000) :: string
    character(len = 256) :: share_string
    integer :: ierr, ierrlog, iunit, irec, newlen
    
    PUSH_SUB(opencl_build_program)

    string = '#include "'//trim(filename)//'"'

    if(debug%info) then
      call messages_write("Building CL program '"//trim(filename)//"'.")
      call messages_info()
    end if

    prog = clCreateProgramWithSource(accel%context%cl_context, trim(string), ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clCreateProgramWithSource")

    ! build the compilation flags
    string='-w'
    ! full optimization
    string=trim(string)//' -cl-denorms-are-zero'
    ! The following flag gives an error with the Xeon Phi
    !    string=trim(string)//' -cl-strict-aliasing'
    string=trim(string)//' -cl-mad-enable'
    string=trim(string)//' -cl-unsafe-math-optimizations'
    string=trim(string)//' -cl-finite-math-only'
    string=trim(string)//' -cl-fast-relaxed-math'

    share_string='-I'//trim(conf%share)//'/opencl/'

    if (f90_cl_device_has_extension(accel%device%cl_device, "cl_khr_fp64")) then
      string = trim(string)//' -DEXT_KHR_FP64'
    else if(f90_cl_device_has_extension(accel%device%cl_device, "cl_amd_fp64")) then
      string = trim(string)//' -DEXT_AMD_FP64'
    else
      call messages_write('Octopus requires an OpenCL device with double-precision support.')
      call messages_fatal()
    end if

    if(accel_use_shared_mem()) then
      string = trim(string)//' -DSHARED_MEM'
    end if

    if(present(flags)) then
      string = trim(string)//' '//trim(flags)
    end if

    if(debug%info) then
      call messages_write("Debug info: compilation flags '"//trim(string), new_line = .true.)
      call messages_write('  '//trim(share_string)//"'.")
      call messages_info()
    end if

    string = trim(string)//' '//trim(share_string)

    call clBuildProgram(prog, trim(string), ierr)

    call clGetProgramBuildInfo(prog, accel%device%cl_device, CL_PROGRAM_BUILD_LOG, string, ierrlog)
    if(ierrlog /= CL_SUCCESS) call opencl_print_error(ierrlog, "clGetProgramBuildInfo")

    ! CL_PROGRAM_BUILD_LOG seems to have a useless '\n' in it
    newlen = scan(string, achar(010), back = .true.) - 1
    if(newlen >= 0) string = string(1:newlen)
    
    if(len(trim(string)) > 0) write(stderr, '(a)') trim(string)

    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clBuildProgram")
    
    POP_SUB(opencl_build_program)
  end subroutine opencl_build_program
#endif

  ! -----------------------------------------------
#ifdef HAVE_OPENCL
  subroutine opencl_release_program(prog)
    type(cl_program),    intent(inout) :: prog

    integer :: ierr

    PUSH_SUB(opencl_release_program)

    call clReleaseProgram(prog, ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clReleaseProgram")

    POP_SUB(opencl_release_program)
  end subroutine opencl_release_program
#endif

  ! -----------------------------------------------

#ifdef HAVE_OPENCL
  subroutine opencl_release_kernel(prog)
    type(cl_kernel),      intent(inout) :: prog

    integer :: ierr

    PUSH_SUB(opencl_release_kernel)

#ifdef HAVE_OPENCL
    call clReleaseKernel(prog, ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clReleaseKernel")
#endif

    POP_SUB(opencl_release_kernel)
  end subroutine opencl_release_kernel
#endif

#ifdef HAVE_OPENCL
  ! -----------------------------------------------
  subroutine opencl_create_kernel(kernel, prog, name)
    type(cl_kernel),  intent(inout) :: kernel
    type(cl_program), intent(inout) :: prog
    character(len=*), intent(in)    :: name

    integer :: ierr
    type(profile_t), save :: prof

    PUSH_SUB(opencl_create_kernel)
    call profiling_in(prof, "CL_BUILD_KERNEL", exclude = .true.)
    
#ifdef HAVE_OPENCL
    kernel = clCreateKernel(prog, name, ierr)
    if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "clCreateKernel")
#endif

    call profiling_out(prof)
    POP_SUB(opencl_create_kernel)
  end subroutine opencl_create_kernel
#endif
  
  ! ------------------------------------------------
#ifdef HAVE_OPENCL
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
#endif

  ! ----------------------------------------------------

  subroutine clblas_print_error(ierr, name)
    integer,          intent(in) :: ierr
    character(len=*), intent(in) :: name

    character(len=40) :: errcode

    PUSH_SUB(clblas_print_error)
#ifdef HAVE_CLBLAS
    select case(ierr)
    case(clblasSuccess);                    errcode = 'clblasSuccess'
    case(clblasInvalidValue);               errcode = 'clblasInvalidValue'
    case(clblasInvalidCommandQueue);        errcode = 'clblasInvalidCommandQueue'
    case(clblasInvalidContext);             errcode = 'clblasInvalidContext'
    case(clblasInvalidMemObject);           errcode = 'clblasInvalidMemObject'
    case(clblasInvalidDevice);              errcode = 'clblasInvalidDevice'
    case(clblasInvalidEventWaitList);       errcode = 'clblasInvalidEventWaitList'
    case(clblasOutOfResources);             errcode = 'clblasOutOfResources'
    case(clblasOutOfHostMemory);            errcode = 'clblasOutOfHostMemory'
    case(clblasInvalidOperation);           errcode = 'clblasInvalidOperation'
    case(clblasCompilerNotAvailable);       errcode = 'clblasCompilerNotAvailable'
    case(clblasBuildProgramFailure );       errcode = 'clblasBuildProgramFailure'
    case(clblasNotImplemented);             errcode = 'clblasNotImplemented'
    case(clblasNotInitialized);             errcode = 'clblasNotInitialized'
    case(clblasInvalidMatA);                errcode = 'clblasInvalidMatA'
    case(clblasInvalidMatB);                errcode = 'clblasInvalidMatB'
    case(clblasInvalidMatC);                errcode = 'clblasInvalidMatC'
    case(clblasInvalidVecX);                errcode = 'clblasInvalidVecX'
    case(clblasInvalidVecY);                errcode = 'clblasInvalidVecY'
    case(clblasInvalidDim);                 errcode = 'clblasInvalidDim'
    case(clblasInvalidLeadDimA);            errcode = 'clblasInvalidLeadDimA'
    case(clblasInvalidLeadDimB);            errcode = 'clblasInvalidLeadDimB'
    case(clblasInvalidLeadDimC);            errcode = 'clblasInvalidLeadDimC'
    case(clblasInvalidIncX);                errcode = 'clblasInvalidIncX'
    case(clblasInvalidIncY);                errcode = 'clblasInvalidIncY'
    case(clblasInsufficientMemMatA);        errcode = 'clblasInsufficientMemMatA'
    case(clblasInsufficientMemMatB);        errcode = 'clblasInsufficientMemMatB'
    case(clblasInsufficientMemMatC);        errcode = 'clblasInsufficientMemMatC'
    case(clblasInsufficientMemVecX);        errcode = 'clblasInsufficientMemVecX'
    case(clblasInsufficientMemVecY);        errcode = 'clblasInsufficientMemVecY'
    case default
      write(errcode, '(i10)') ierr
      errcode = 'UNKNOWN ERROR CODE ('//trim(adjustl(errcode))//')'
    end select
#endif

    message(1) = 'clblas '//trim(name)//' '//trim(errcode)
    call messages_fatal(1)

    POP_SUB(clblas_print_error)
  end subroutine clblas_print_error

  ! ----------------------------------------------------
  subroutine clfft_print_error(ierr, name)
    integer,          intent(in) :: ierr
    character(len=*), intent(in) :: name

    character(len=40) :: errcode

    PUSH_SUB(clfft_print_error)
#ifdef HAVE_CLFFT
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

    message(1) = 'clfft '//trim(name)//' '//trim(errcode)
    call messages_fatal(1)

    POP_SUB(clfft_print_error)
  end subroutine clfft_print_error

  ! ----------------------------------------------------

#ifdef HAVE_OPENCL
  logical function f90_cl_device_has_extension(device, extension) result(has)
    type(cl_device_id), intent(inout) :: device
    character(len=*),   intent(in)    :: extension

    integer :: cl_status
    character(len=2048) :: all_extensions

#ifdef HAVE_OPENCL
    call clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, all_extensions, cl_status)
#endif
    
    has = index(all_extensions, extension) /= 0

  end function f90_cl_device_has_extension
#endif
  
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

  subroutine accel_set_buffer_to_zero(buffer, type, nval, offset)
    type(accel_mem_t),  intent(inout) :: buffer
    type(type_t),       intent(in)    :: type
    integer,            intent(in)    :: nval
    integer, optional,  intent(in)    :: offset

    integer :: nval_real, bsize, offset_real

    PUSH_SUB(accel_set_buffer_to_zero)

    ASSERT(type == TYPE_CMPLX .or. type == TYPE_FLOAT)

    if(nval > 0) then
      
      nval_real = nval*(types_get_size(type)/8)
      offset_real = optional_default(offset, 0)*(types_get_size(type)/8)
      
      ASSERT(nval_real > 0)
      
      call accel_set_kernel_arg(set_zero, 0, nval_real)
      call accel_set_kernel_arg(set_zero, 1, offset_real)
      call accel_set_kernel_arg(set_zero, 2, buffer)
      
      bsize = accel_kernel_workgroup_size(set_zero)

           
      call accel_kernel_run(set_zero, (/ opencl_pad(nval_real, bsize) /), (/ bsize /))
      call accel_finish()

    end if
      
    POP_SUB(accel_set_buffer_to_zero)
  end subroutine accel_set_buffer_to_zero

  ! ----------------------------------------------------

  subroutine opencl_check_bandwidth()
    integer :: itime
    integer, parameter :: times = 10
    integer :: size
    FLOAT   :: time, stime
    FLOAT   :: read_bw, write_bw
    type(accel_mem_t) :: buff
    FLOAT, allocatable :: data(:)

    call messages_new_line()
    call messages_write('Info: Benchmarking the bandwidth between main memory and device memory')
    call messages_new_line()
    call messages_info()

    call messages_write(' Buffer size   Read bw  Write bw')
    call messages_new_line()
    call messages_write('       [MiB]   [MiB/s]   [MiB/s]')
    call messages_info()

    size = 15000
    do 
      SAFE_ALLOCATE(data(1:size))
      call accel_create_buffer(buff, ACCEL_MEM_READ_WRITE, TYPE_FLOAT, size)

      stime = loct_clock()
      do itime = 1, times
        call accel_write_buffer(buff, size, data)
        call accel_finish()
      end do
      time = (loct_clock() - stime)/TOFLOAT(times)

      write_bw = TOFLOAT(size)*CNST(8.0)/time

      stime = loct_clock()
      do itime = 1, times
        call accel_read_buffer(buff, size, data)
      end do
      call accel_finish()

      time = (loct_clock() - stime)/TOFLOAT(times)
      read_bw = TOFLOAT(size)*CNST(8.0)/time

      call messages_write(size*CNST(8.0)/CNST(1024.0)**2)
      call messages_write(write_bw/CNST(1024.0)**2, fmt = '(f10.1)')
      call messages_write(read_bw/CNST(1024.0)**2, fmt = '(f10.1)')
      call messages_info()

      call accel_release_buffer(buff)

      SAFE_DEALLOCATE_A(data)

      size = int(size*2.0)

      if(size > 50000000) exit
    end do
  end subroutine opencl_check_bandwidth

  ! ----------------------------------------------------

  logical pure function accel_use_shared_mem() result(use_shared_mem)
    
    use_shared_mem = accel%shared_mem

  end function accel_use_shared_mem

  !------------------------------------------------------------

  subroutine accel_kernel_global_init()
    
    PUSH_SUB(accel_kernel_global_init)

    nullify(head)

    call cuda_module_map_init(accel%module_map)
    
    POP_SUB(accel_kernel_global_init)
  end subroutine accel_kernel_global_init

  !------------------------------------------------------------
  
  subroutine accel_kernel_global_end()
    type(accel_kernel_t), pointer :: next_head

    PUSH_SUB(accel_kernel_global_end)

    do
      if(.not. associated(head)) exit
      next_head => head%next
      call accel_kernel_end(head)
      head => next_head
    end do

    if(accel_is_enabled()) then
      call cuda_module_map_end(accel%module_map)
    end if
    
    POP_SUB(accel_kernel_global_end)
  end subroutine accel_kernel_global_end

  !------------------------------------------------------------

  subroutine accel_kernel_build(this, file_name, kernel_name, flags)
    type(accel_kernel_t),        intent(inout) :: this
    character(len=*),            intent(in)    :: file_name
    character(len=*),            intent(in)    :: kernel_name
    character(len=*), optional,  intent(in)    :: flags

    type(profile_t), save :: prof
#ifdef HAVE_OPENCL
    type(cl_program) :: prog
#endif
#ifdef HAVE_CUDA
    character(len=1000) :: all_flags
    type(c_ptr) :: cuda_module
#endif
   
    PUSH_SUB(accel_kernel_build)

    call profiling_in(prof, "ACCEL_COMPILE", exclude = .true.)

#ifdef HAVE_CUDA
    all_flags = '-I'//trim(conf%share)//'/opencl/'

    if(accel_use_shared_mem()) then
      all_flags = trim(all_flags)//' -DSHARED_MEM'
    end if
    
    if(present(flags)) then
      all_flags = trim(all_flags)//' '//trim(flags)
    end if
    
    call cuda_build_program(accel%module_map, this%cuda_module, accel%device%cuda_device, trim(file_name), trim(all_flags))
    
    call cuda_create_kernel(this%cuda_kernel, this%cuda_module, trim(kernel_name))
    call cuda_alloc_arg_array(this%arguments)

    this%cuda_shared_mem = 0
#endif

#ifdef HAVE_OPENCL
    call opencl_build_program(prog, trim(conf%share)//'/opencl/'//trim(file_name), flags = flags)
    call opencl_create_kernel(this%kernel, prog, trim(kernel_name))
    call opencl_release_program(prog)
#endif

    this%initialized = .true.

    call profiling_out(prof)
    
    POP_SUB(accel_kernel_build)
  end subroutine accel_kernel_build

  !------------------------------------------------------------

  subroutine accel_kernel_end(this)
    type(accel_kernel_t), intent(inout) :: this
#ifdef HAVE_OPENCL
    integer :: ierr
#endif

      PUSH_SUB(accel_kernel_end)

#ifdef HAVE_CUDA
      call cuda_free_arg_array(this%arguments)
      call cuda_release_kernel(this%cuda_kernel)
      ! modules are not released here, since they are not associated to a kernel
#endif
      
#ifdef HAVE_OPENCL
      call clReleaseKernel(this%kernel, ierr)
      if(ierr /= CL_SUCCESS) call opencl_print_error(ierr, "release_kernel")
#endif
      this%initialized = .false.

      POP_SUB(accel_kernel_end)
  end subroutine accel_kernel_end

  !------------------------------------------------------------

  subroutine accel_kernel_start_call(this, file_name, kernel_name, flags)
    type(accel_kernel_t), target, intent(inout) :: this
    character(len=*),             intent(in)    :: file_name
    character(len=*),             intent(in)    :: kernel_name
    character(len=*), optional,   intent(in)    :: flags

    PUSH_SUB(accel_kernel_start_call)

    if(.not. this%initialized) then
      call accel_kernel_build(this, file_name, kernel_name, flags)
      this%next => head
      head => this
    end if

    POP_SUB(accel_kernel_start_call)
  end subroutine accel_kernel_start_call

  !--------------------------------------------------------------

  integer(8) pure function accel_global_memory_size() result(size)

    size = accel%global_memory_size
    
  end function accel_global_memory_size

  !--------------------------------------------------------------
  
  integer(8) pure function accel_local_memory_size() result(size)

    size = accel%local_memory_size
    
  end function accel_local_memory_size

  !--------------------------------------------------------------

  integer pure function accel_max_size_per_dim(dim) result(size)
    integer, intent(in) :: dim

    size = 0
#ifdef HAVE_OPENCL
    size = 2**30
#endif
#ifdef HAVE_CUDA
    if(dim == 1) size = 2**30
    size = 32768
#endif
  end function accel_max_size_per_dim

  ! ------------------------------------------------------

  subroutine accel_set_stream(stream_number)
    integer, intent(in) :: stream_number

    PUSH_SUB(accel_set_stream)

    if(accel_is_enabled()) then
#ifdef HAVE_CUDA
      call cuda_set_stream(accel%cuda_stream, stream_number)
      call cublas_set_stream(accel%cublas_handle, accel%cuda_stream)
#endif
    end if

    POP_SUB(accel_set_stream)
  end subroutine accel_set_stream

  ! ------------------------------------------------------

  subroutine accel_synchronize_all_streams()
    PUSH_SUB(accel_synchronize_all_streams)

#ifdef HAVE_CUDA
    call cuda_synchronize_all_streams()
#endif

    POP_SUB(accel_synchronize_all_streams)
  end subroutine accel_synchronize_all_streams

#include "undef.F90"
#include "real.F90"
#include "accel_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "accel_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "accel_inc.F90"

end module accel_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

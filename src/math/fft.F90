!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2011 J. Alberdi-Rodriguez, P. Garcia Risueño, M. Oliveira
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

!> Fast Fourier Transform module.
!! This module provides a single interface that works with different
!! FFT implementations.
module fft_oct_m
  use accel_oct_m
#ifdef HAVE_OPENCL
  use cl
#ifdef HAVE_CLFFT
  use clfft
#endif
#endif
  use fftw_oct_m
  use fftw_params_oct_m
  use global_oct_m
  use,intrinsic :: iso_c_binding
  use lalg_basic_oct_m
  use loct_math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use nfft_oct_m
#if defined(HAVE_OPENMP) && defined(HAVE_FFTW3_THREADS)
  use omp_lib
#endif
  use parser_oct_m
  use pfft_oct_m
  use pfft_params_oct_m
  use pnfft_oct_m
  use profiling_oct_m
  use types_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::            &
    fft_t,             &
    fft_all_init,      &
    fft_all_end,       &
    fft_init,          &
    fft_init_stage1,   &
    fft_end,           &
    fft_copy,          &
    fft_get_dims,      &
    pad_feq,           &
    dfft_forward,      &
    zfft_forward,      &
    dfft_backward,     &
    zfft_backward,     &
    fft_scaling_factor


  !> global constants
  integer, public, parameter :: &
       FFT_NONE    = 0,         &
       FFT_REAL    = 1,         &
       FFT_COMPLEX = 2

  integer, public, parameter :: &
       FFTLIB_NONE  = 0, &
       FFTLIB_FFTW  = 1, &
       FFTLIB_PFFT  = 2, &
       FFTLIB_ACCEL = 3, &
       FFTLIB_NFFT  = 4, &
       FFTLIB_PNFFT = 5
       
  integer, parameter :: &
    FFT_MAX  = 10, &
    FFT_NULL = -1
  
  type fft_t
    private
    integer         :: slot = 0 !< in which slot do we have this fft

    integer, public :: type    !< is the fft real or complex
    integer, public :: library !< what library are we using

    integer         :: comm           !< MPI communicator
    integer         :: rs_n_global(3) !< total size of the fft in each direction in real space
    integer         :: fs_n_global(3) !< total size of the fft in each direction in fourier space
    integer         :: rs_n(3)        !< local size of the fft in in each direction real space
    integer         :: fs_n(3)        !< local size of the fft in in each direction fourier space
    integer         :: rs_istart(1:3) !< where does the local portion of the function start in real space
    integer         :: fs_istart(1:3) !< where does the local portion of the function start in fourier space

    integer, public :: stride_rs(1:3)
    integer, public :: stride_fs(1:3)

    type(c_ptr) :: planf                  !< plan for forward transform
    type(c_ptr) :: planb                  !< plan for backward transform
    !integer(ptrdiff_t_kind) :: pfft_planf !< PFFT plan for forward transform
    !integer(ptrdiff_t_kind) :: pfft_planb !< PFFT plan for backward transform

    !> The following arrays have to be stored here and allocated in the initialization routine because of PFFT 
    !> These arrays are also used for FFTW, as we want to have aligned memory
    FLOAT, pointer, public :: drs_data(:,:,:) !< array used to store the function in real space.
    CMPLX, pointer, public :: zrs_data(:,:,:) !< array used to store the function in real space.
    CMPLX, pointer, public ::  fs_data(:,:,:) !< array used to store the function in Fourier space.
#ifdef HAVE_CLFFT
    !> data for clfft
    type(clfftPlanHandle) :: cl_plan_fw 
    type(clfftPlanHandle) :: cl_plan_bw !< for real transforms we need a different plan, so we always use 2
#endif
    type(c_ptr)           :: cuda_plan_fw
    type(c_ptr)           :: cuda_plan_bw
    type(nfft_t),  public :: nfft
    type(pnfft_t), public :: pnfft

    logical, public :: aligned_memory
  end type fft_t

  interface dfft_forward
    module procedure dfft_forward_1d, dfft_forward_cl, dfft_forward_3d
  end interface dfft_forward

  interface zfft_forward
    module procedure zfft_forward_1d, zfft_forward_cl, zfft_forward_3d
  end interface zfft_forward

  interface dfft_backward
    module procedure dfft_backward_1d, dfft_backward_cl, dfft_backward_3d
  end interface dfft_backward

  interface zfft_backward
    module procedure zfft_backward_1d, zfft_backward_cl, zfft_backward_3d
  end interface zfft_backward

  logical, save, public :: fft_initialized = .false.
  integer, save         :: fft_refs(FFT_MAX)
  type(fft_t), save     :: fft_array(FFT_MAX)
  logical               :: fft_optimize
  integer, save         :: fft_prepare_plan
  integer, public       :: fft_default_lib = -1
  type(nfft_t), save    :: nfft_options
  type(pnfft_t), save   :: pnfft_options

  integer, parameter ::  &
    CUFFT_R2C = int(z'2a'),   &
    CUFFT_C2R = int(z'2c'),   &
    CUFFT_C2C = int(z'29'),   &
    CUFFT_D2Z = int(z'6a'),   &
    CUFFT_Z2D = int(z'6c'),   &
    CUFFT_Z2Z = int(z'69')

contains

  ! ---------------------------------------------------------
  !> initialize the table
  subroutine fft_all_init(namespace)
    type(namespace_t),      intent(in)   :: namespace
    
    integer :: ii
#if defined(HAVE_OPENMP) && defined(HAVE_FFTW3_THREADS)
    integer :: iret
#endif

    PUSH_SUB(fft_all_init)

    fft_initialized = .true.
    
    !%Variable FFTOptimize
    !%Type logical
    !%Default yes
    !%Section Mesh::FFTs
    !%Description
    !% Should <tt>octopus</tt> optimize the FFT dimensions? 
    !% This means that the mesh to which FFTs are applied is not taken to be as small
    !% as possible: some points may be added to each direction in order to get a "good number"
    !% for the performance of the FFT algorithm.
    !% The best FFT grid dimensions are given by <math>2^a 3^b 5^c 7^d 11^e 13^f</math>
    !% where <math>a,b,c,d</math> are arbitrary and <math>e,f</math> are 0 or 1.
    !% (<a href=http://www.fftw.org/doc/Complex-DFTs.html>ref</a>).
    !% In some cases, namely when using
    !% the split-operator, or Suzuki-Trotter propagators, this option should be turned off.
    !% For spatial FFTs in periodic directions, the grid is never optimized, but a warning will
    !% be written if the number is not good, with a suggestion of a better one to use, so you
    !% can try a different spacing if you want to get a good number.
    !%End
    call parse_variable(namespace, 'FFTOptimize', .true., fft_optimize)
    do ii = 1, FFT_MAX
      fft_refs(ii) = FFT_NULL
    end do

    !%Variable FFTPreparePlan
    !%Type integer
    !%Default fftw_measure
    !%Section Mesh::FFTs
    !%Description
    !% The FFTs are performed in octopus with the help of <a href=http://www.fftw.org>FFTW</a> and similar packages.
    !% Before doing the actual computations, this package prepares a "plan", which means that 
    !% the precise numerical strategy to be followed to compute the FFT is machine/compiler-dependent,
    !% and therefore the software attempts to figure out which is this precise strategy (see the
    !% FFTW documentation for details). This plan preparation, which has to be done for each particular
    !% FFT shape, can be done exhaustively and carefully (slow), or merely estimated. Since this is
    !% a rather critical numerical step, by default it is done carefully, which implies a longer initial
    !% initialization, but faster subsequent computations. You can change this behaviour by changing
    !% this <tt>FFTPreparePlan</tt> variable, and in this way you can force FFTW to do a fast guess or
    !% estimation of which is the best way to perform the FFT.
    !%Option fftw_measure 0
    !% This is the default, and implies a longer initialization, but involves a more careful analysis
    !% of the strategy to follow, and therefore more efficient FFTs.
    !%Option fftw_estimate 64
    !% This is the "fast initialization" scheme, in which the plan is merely guessed from "reasonable"
    !% assumptions.
    !%Option fftw_patient 32
    !% It is like fftw_measure, but considers a wider range of algorithms and often produces a
    !% "more optimal" plan (especially for large transforms), but at the expense of several times
    !% longer planning time (especially for large transforms).
    !%Option fftw_exhaustive 8
    !% It is like fftw_patient, but considers an even wider range of algorithms, 
    !% including many that we think are unlikely to be fast, to produce the most optimal
    !%  plan but with a substantially increased planning time.
    !%End
    call parse_variable(namespace, 'FFTPreparePlan', FFTW_MEASURE, fft_prepare_plan)
    if(.not. varinfo_valid_option('FFTPreparePlan', fft_prepare_plan)) then
      call messages_input_error(namespace, 'FFTPreparePlan')
    end if

    !%Variable FFTLibrary
    !%Type integer
    !%Section Mesh::FFTs
    !%Default fftw
    !%Description
    !% (experimental) You can select the FFT library to use.
    !%Option fftw 1
    !% Uses FFTW3 library.
    !%Option pfft 2
    !% (experimental) Uses PFFT library, which has to be linked.
    !%Option accel 3
    !% (experimental) Uses a GPU accelerated library. This only
    !% works if Octopus was compiled with Cuda or OpenCL support.
    !%End
    call parse_variable(namespace, 'FFTLibrary', FFTLIB_FFTW, fft_default_lib)

    if(fft_default_lib == FFTLIB_ACCEL) then
#if ! (defined(HAVE_CLFFT) || defined(HAVE_CUDA))
      call messages_write('You have selected the Accelerated FFT, but Octopus was compiled', new_line = .true.)
      call messages_write('without clfft (OpenCL) or Cuda support.')
      call messages_fatal()
#endif
      if(.not. accel_is_enabled()) then
        call messages_write('You have selected the accelerated FFT, but acceleration is disabled.')
        call messages_fatal()
      end if
    end if
    
#if defined(HAVE_OPENMP) && defined(HAVE_FFTW3_THREADS)
    if(omp_get_max_threads() > 1) then

      call messages_write('Info: Initializing Multi-threaded FFTW')
      call messages_info()
      
      iret = fftw_init_threads()
      if (iret == 0) then 
        call messages_write('Initialization of FFTW3 threads failed.')
        call messages_fatal()
      end if
      call fftw_plan_with_nthreads(omp_get_max_threads())

    end if
#endif

    call nfft_guru_options(nfft_options, namespace)
    call pnfft_guru_options(pnfft_options, namespace)
    
    POP_SUB(fft_all_init)
  end subroutine fft_all_init


  ! ---------------------------------------------------------
  !> delete all plans
  subroutine fft_all_end()
    integer :: ii

    PUSH_SUB(fft_all_end)

    do ii = 1, FFT_MAX
      if(fft_refs(ii) /= FFT_NULL) then
        call fft_end(fft_array(ii))
      end if
    end do

#ifdef HAVE_PFFT
    call pfft_cleanup()
#endif

#if defined(HAVE_OPENMP) && defined(HAVE_FFTW3_THREADS)
    call fftw_cleanup_threads()
#else
    call fftw_cleanup()
#endif

    fft_initialized = .false.
    
    POP_SUB(fft_all_end)
  end subroutine fft_all_end

  ! ---------------------------------------------------------
  subroutine fft_init(this, nn, dim, type, library, optimize, optimize_parity, mpi_comm, mpi_grp, use_aligned)
    type(fft_t),       intent(inout) :: this     !< FFT data type
    integer,           intent(inout) :: nn(3)    !< Size of the box
    integer,           intent(in)    :: dim      !< Dimensions of the box
    integer,           intent(in)    :: type     !< The type of the FFT; real or complex
    integer,           intent(in)    :: library  !< FFT library to be used. Can be overridden.
    logical,           intent(in)    :: optimize(3) !< whether we should optimize grid in each direction
    integer,           intent(in)    :: optimize_parity(3) !< choose optimized grid in each direction as
                                                 !! even (0), odd (1), or whatever (negative).
    integer, optional, intent(out)   :: mpi_comm !< MPI communicator
    type(mpi_grp_t), optional, intent(in) :: mpi_grp !< the mpi_group we want to use for the parallelization
    logical, optional                :: use_aligned !< For FFTW we can use aligned memory

    integer :: ii, jj, fft_dim, idir, column_size, row_size, alloc_size, n3
    integer :: n_1, n_2, n_3, nn_temp(3)
    integer :: library_
    type(mpi_grp_t) :: mpi_grp_
    integer(8) :: number_points

#ifdef HAVE_CLFFT
    real(8) :: scale
    integer :: status
#endif
#ifdef HAVE_PFFT
    integer :: ierror
#endif 

    PUSH_SUB(fft_init)

    ASSERT(fft_initialized)
    
    ASSERT(type == FFT_REAL .or. type == FFT_COMPLEX)

    mpi_grp_ = mpi_world
    if(present(mpi_grp)) mpi_grp_ = mpi_grp

    this%aligned_memory = optional_default(use_aligned, .false.)

    ! First, figure out the dimensionality of the FFT.
    fft_dim = 0
    do ii = 1, dim
      if(nn(ii) <= 1) exit
      fft_dim = fft_dim + 1
    end do

    if(fft_dim  ==  0) then
      message(1) = "Internal error in fft_init: apparently, a 1x1x1 FFT is required."
      call messages_fatal(1)
    end if

    if(fft_dim > 3) call messages_not_implemented('FFT for dimension > 3')

    library_ = library
    nn_temp(1:fft_dim) = nn(1:fft_dim)

    select case (library_)
    case (FFTLIB_ACCEL)
    
      do ii = 1, fft_dim
        ! the AMD OpenCL FFT only supports sizes 2, 3 and 5, but size
        ! 5 gives an fpe error on the Radeon 7970 (APPML 1.8), so we
        ! only use factors 2 and 3
#ifdef HAVE_CLFFT
        nn_temp(ii) = fft_size(nn(ii), (/2, 3/))
#else
        nn_temp(ii) = fft_size(nn(ii), (/2, 3, 5, 7/))
#endif
        if(fft_optimize .and. optimize(ii)) nn(ii) = nn_temp(ii)
      end do 
      
      ! if we can't optimize, in some cases we can't use the library
      if(any(nn(1:fft_dim) /= nn_temp(1:fft_dim))) then
        call messages_write('Invalid grid size for clfft. FFTW will be used instead.')
        call messages_warning()
        library_ = FFTLIB_FFTW
      end if
      
    case (FFTLIB_NFFT)
            
      do ii = 1, fft_dim
        !NFFT likes even grids
        !The underlying FFT grids are optimized inside the nfft_init routine
        if(int(nn(ii)/2)*2 /= nn(ii) .and. (fft_optimize .and. optimize(ii)) )&
          nn(ii)=nn(ii)+1 
      end do 

    case (FFTLIB_PNFFT)
          
      do ii = 1, fft_dim
        !also PNFFT likes even grids
        if(int(nn(ii)/2)*2 /= nn(ii)) nn(ii)=nn(ii)+1 
      end do 
          
      if(fft_dim < 3) &
          call messages_not_implemented('PNFFT support for dimension < 3')
                    
    case default

      if(fft_dim < 3 .and. library_ == FFTLIB_PFFT) &
           call messages_not_implemented('PFFT support for dimension < 3')


      ! FFT optimization
      if(any(optimize_parity(1:fft_dim) > 1)) then
        message(1) = "Internal error in fft_init: optimize_parity must be negative, 0, or 1."
        call messages_fatal(1)
      end if
      
      do ii = 1, fft_dim
        call loct_fft_optimize(nn_temp(ii), optimize_parity(ii))
        if(fft_optimize .and. optimize(ii)) nn(ii) = nn_temp(ii)
      end do
      
    end select

    ! find out if fft has already been allocated
    jj = 0
    do ii = FFT_MAX, 1, -1
      if(fft_refs(ii) /= FFT_NULL) then
        if(all(nn(1:dim) == fft_array(ii)%rs_n_global(1:dim)) .and. type == fft_array(ii)%type &
             .and. library_ == fft_array(ii)%library .and. library_ /= FFTLIB_NFFT &
             .and. library_ /= FFTLIB_PNFFT &
             .and. this%aligned_memory .eqv. fft_array(ii)%aligned_memory) then
             ! NFFT and PNFFT plans are always allocated from scratch since they 
             ! are very likely to be different
          this = fft_array(ii)              ! return a copy
          fft_refs(ii) = fft_refs(ii) + 1  ! increment the ref count
          if (present(mpi_comm)) mpi_comm = fft_array(ii)%comm ! also return the MPI communicator
          POP_SUB(fft_init)
          return
        end if
      else
        jj = ii
      end if
    end do

    if(jj == 0) then
      message(1) = "Not enough slots for FFTs."
      message(2) = "Please increase FFT_MAX in fft.F90 and recompile."
      call messages_fatal(2)
    end if

    ! jj now contains an empty slot
    fft_refs(jj) = 1
    fft_array(jj)%slot     = jj
    fft_array(jj)%type     = type
    fft_array(jj)%library  = library_
    fft_array(jj)%rs_n_global(1:dim) = nn(1:dim)
    fft_array(jj)%rs_n_global(dim+1:) = 1
    nullify(fft_array(jj)%drs_data)
    nullify(fft_array(jj)%zrs_data)
    nullify(fft_array(jj)%fs_data)

    fft_array(jj)%aligned_memory = this%aligned_memory 

    ! Initialize parallel communicator
    select case (library_)
    case (FFTLIB_PFFT)
#ifdef HAVE_PFFT
      call pfft_init()
 
      call pfft_decompose(mpi_grp_%size, column_size, row_size)
 
      ierror = pfft_create_procmesh_2d(mpi_grp_%comm, column_size, row_size, fft_array(jj)%comm)        
   
      if (ierror /= 0) then
        message(1) = "The number of rows and columns in PFFT processor grid is not equal to "
        message(2) = "the number of processor in the MPI communicator."
        message(3) = "Please check it."
        call messages_fatal(3)
      end if
#endif

    case (FFTLIB_PNFFT)
      call pnfft_init_procmesh(fft_array(jj)%pnfft, mpi_grp_, fft_array(jj)%comm) 

    case default
      fft_array(jj)%comm = -1

    end select
    
    if (present(mpi_comm)) mpi_comm = fft_array(jj)%comm

    ! Get dimentions of arrays
    select case (library_)
    case (FFTLIB_FFTW)
      call fftw_get_dims(fft_array(jj)%rs_n_global, type == FFT_REAL, fft_array(jj)%fs_n_global)
      fft_array(jj)%rs_n = fft_array(jj)%rs_n_global
      fft_array(jj)%fs_n = fft_array(jj)%fs_n_global
      fft_array(jj)%rs_istart = 1
      fft_array(jj)%fs_istart = 1

      if(this%aligned_memory) then
        call fftw_alloc_memory(fft_array(jj)%rs_n_global, type == FFT_REAL, fft_array(jj)%fs_n_global, &
                               fft_array(jj)%drs_data, fft_array(jj)%zrs_data, fft_array(jj)%fs_data)
      end if

    case (FFTLIB_PFFT)
#ifdef HAVE_PFFT     
      call pfft_get_dims(fft_array(jj)%rs_n_global, mpi_comm, type == FFT_REAL, &
           alloc_size, fft_array(jj)%fs_n_global, fft_array(jj)%rs_n, &
           fft_array(jj)%fs_n, fft_array(jj)%rs_istart, fft_array(jj)%fs_istart)
      !write(*,"(6(A,3I4,/),A,I10,/)") "PFFT: rs_n_global = ",fft_array(jj)%rs_n_global,&
      !  "fs_n_global = ",fft_array(jj)%fs_n_global,&
      !  "rs_n        = ",fft_array(jj)%rs_n,&
      !  "fs_n        = ",fft_array(jj)%fs_n,&
      !  "rs_istart   = ",fft_array(jj)%rs_istart,&
      !  "fs_istart   = ",fft_array(jj)%fs_istart,&
      !  "alloc_size  = ",alloc_size
#endif

      ! Allocate memory. Note that PFFT may need extra memory space 
      ! and that in fourier space the function will be transposed
      if (type == FFT_REAL) then
        n_1 = max(1, fft_array(jj)%rs_n(1))
        n_2 = max(1, fft_array(jj)%rs_n(2))
        n_3 = max(1, fft_array(jj)%rs_n(3))

        n3 = ceiling(real(2*alloc_size)/real(n_1*n_2))
        SAFE_ALLOCATE(fft_array(jj)%drs_data(1:n_1, 1:n_2, 1:n3))
      else
        n3 = ceiling(real(alloc_size)/real(fft_array(jj)%rs_n(1)*fft_array(jj)%rs_n(2)))
        SAFE_ALLOCATE(fft_array(jj)%zrs_data(1:fft_array(jj)%rs_n(1), 1:fft_array(jj)%rs_n(2), 1:n3))
      end if

      n_1 = max(1, fft_array(jj)%fs_n(1))
      n_2 = max(1, fft_array(jj)%fs_n(2))
      n_3 = max(1, fft_array(jj)%fs_n(3))

      n3 = ceiling(real(alloc_size)/real(n_3*n_1))
      SAFE_ALLOCATE(fft_array(jj)%fs_data(1:n_3, 1:n_1, 1:n3))

    case(FFTLIB_ACCEL)
      call fftw_get_dims(fft_array(jj)%rs_n_global, (type == FFT_REAL), fft_array(jj)%fs_n_global)
      fft_array(jj)%rs_n = fft_array(jj)%rs_n_global
      fft_array(jj)%fs_n = fft_array(jj)%fs_n_global
      fft_array(jj)%rs_istart = 1
      fft_array(jj)%fs_istart = 1

    case(FFTLIB_NFFT)
      fft_array(jj)%fs_n_global = fft_array(jj)%rs_n_global
      fft_array(jj)%rs_n = fft_array(jj)%rs_n_global
      fft_array(jj)%fs_n = fft_array(jj)%fs_n_global
      fft_array(jj)%rs_istart = 1
      fft_array(jj)%fs_istart = 1
    
    case(FFTLIB_PNFFT)       
      fft_array(jj)%fs_n_global = fft_array(jj)%rs_n_global
      fft_array(jj)%rs_n = fft_array(jj)%rs_n_global
      fft_array(jj)%fs_n = fft_array(jj)%fs_n_global
      fft_array(jj)%rs_istart = 1
      fft_array(jj)%fs_istart = 1
      ! indices partition is performed together with the plan preparation


    end select

    ! Prepare plans
    select case (library_)
    case (FFTLIB_FFTW)
      if(.not. this%aligned_memory) then
        call fftw_prepare_plan(fft_array(jj)%planf, fft_dim, fft_array(jj)%rs_n_global, &
           type == FFT_REAL, FFTW_FORWARD, fft_prepare_plan+FFTW_UNALIGNED)
        call fftw_prepare_plan(fft_array(jj)%planb, fft_dim, fft_array(jj)%rs_n_global, &
           type == FFT_REAL, FFTW_BACKWARD, fft_prepare_plan+FFTW_UNALIGNED)
      else
        if(type == FFT_REAL) then
          call fftw_prepare_plan(fft_array(jj)%planf, fft_dim, fft_array(jj)%rs_n_global, &
             type == FFT_REAL, FFTW_FORWARD, fft_prepare_plan, &
             din_=fft_array(jj)%drs_data, cout_=fft_array(jj)%fs_data)
          call fftw_prepare_plan(fft_array(jj)%planb, fft_dim, fft_array(jj)%rs_n_global, &
             type == FFT_REAL, FFTW_BACKWARD, fft_prepare_plan, &
             din_=fft_array(jj)%drs_data, cout_=fft_array(jj)%fs_data)
        else
          call fftw_prepare_plan(fft_array(jj)%planf, fft_dim, fft_array(jj)%rs_n_global, &
             type == FFT_REAL, FFTW_FORWARD, fft_prepare_plan, &
             cin_=fft_array(jj)%zrs_data, cout_=fft_array(jj)%fs_data)
          call fftw_prepare_plan(fft_array(jj)%planb, fft_dim, fft_array(jj)%rs_n_global, &
             type == FFT_REAL, FFTW_BACKWARD, fft_prepare_plan, &
             cin_=fft_array(jj)%zrs_data, cout_=fft_array(jj)%fs_data)
        end if
      end if

    case(FFTLIB_NFFT)
     call nfft_copy_info(this%nfft,fft_array(jj)%nfft) !copy default parameters set in the calling routine 
     call nfft_init(fft_array(jj)%nfft, nfft_options, fft_array(jj)%rs_n_global, &
                    fft_dim, fft_array(jj)%rs_n_global, optimize = .true.)

    case (FFTLIB_PFFT)
#ifdef HAVE_PFFT     
      if(type == FFT_REAL) then
        call pfft_prepare_plan_r2c(fft_array(jj)%planf, fft_array(jj)%rs_n_global, fft_array(jj)%drs_data, &
             fft_array(jj)%fs_data, FFTW_FORWARD, fft_prepare_plan, mpi_comm)
        call pfft_prepare_plan_c2r(fft_array(jj)%planb, fft_array(jj)%rs_n_global, fft_array(jj)%fs_data, &
             fft_array(jj)%drs_data, FFTW_BACKWARD, fft_prepare_plan, mpi_comm)
      else
        call pfft_prepare_plan_c2c(fft_array(jj)%planf, fft_array(jj)%rs_n_global, fft_array(jj)%zrs_data, &
             fft_array(jj)%fs_data, FFTW_FORWARD, fft_prepare_plan, mpi_comm)
        call pfft_prepare_plan_c2c(fft_array(jj)%planb, fft_array(jj)%rs_n_global, fft_array(jj)%fs_data, &
             fft_array(jj)%zrs_data, FFTW_BACKWARD, fft_prepare_plan, mpi_comm)
      end if
#endif
    case (FFTLIB_PNFFT)
      call pnfft_copy_params(this%pnfft, fft_array(jj)%pnfft) ! pass default parameters like in NFFT

      ! NOTE:
      ! PNFFT (likewise NFFT) breaks the symmetry between real space and Fourier space
      ! by allowing the possibility to have an unstructured grid in rs and by 
      ! using different parallelizations (the rs is transposed w.r.t. fs).
      ! Octopus, in fourier_space_m, uses the convention for which the mapping 
      ! between rs and fs is done with a forward transform (and fs->rs with backward).
      ! This is exactly the opposite of the definitions used by all the libraries 
      ! performing FFTs (PNFFT and NFFT included) [see e.g. M. Frigo, and S. G. Johnson, Proc. 
      ! IEEE 93, 216-231 (2005)].
      ! While this leads to no problem on ordinary ffts where fs and rs can be exchanged 
      ! it does makes a fundamental difference for PNFFT (for some reason I don`t know NFFT 
      ! is still symmetric).
      ! Therefore, in order to perform rs->fs tranforms with PNFFT one should use the 
      ! backward transform.     

      call pnfft_init_plan(fft_array(jj)%pnfft, pnfft_options, mpi_comm, fft_array(jj)%fs_n_global, &
           fft_array(jj)%fs_n, fft_array(jj)%fs_istart, fft_array(jj)%rs_n, fft_array(jj)%rs_istart)

    case(FFTLIB_ACCEL)

      fft_array(jj)%stride_rs(1) = 1
      fft_array(jj)%stride_fs(1) = 1
      do ii = 2, fft_dim
        fft_array(jj)%stride_rs(ii) = fft_array(jj)%stride_rs(ii - 1)*fft_array(jj)%rs_n(ii - 1)
        fft_array(jj)%stride_fs(ii) = fft_array(jj)%stride_fs(ii - 1)*fft_array(jj)%fs_n(ii - 1)
      end do

#ifdef HAVE_CUDA
      call cuda_fft_plan3d(fft_array(jj)%cuda_plan_fw, &
        fft_array(jj)%rs_n_global(3), fft_array(jj)%rs_n_global(2), fft_array(jj)%rs_n_global(1), CUFFT_D2Z)
      call cuda_fft_plan3d(fft_array(jj)%cuda_plan_bw, &
        fft_array(jj)%rs_n_global(3), fft_array(jj)%rs_n_global(2), fft_array(jj)%rs_n_global(1), CUFFT_Z2D)
#endif
      
#ifdef HAVE_CLFFT

      ! create the plans
      call clfftCreateDefaultPlan(fft_array(jj)%cl_plan_fw, accel%context%cl_context, &
        fft_dim, int(fft_array(jj)%rs_n_global, 8), status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftCreateDefaultPlan')

      call clfftCreateDefaultPlan(fft_array(jj)%cl_plan_bw, accel%context%cl_context, &
        fft_dim, int(fft_array(jj)%rs_n_global, 8), status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftCreateDefaultPlan')

      ! set precision

      call clfftSetPlanPrecision(fft_array(jj)%cl_plan_fw, CLFFT_DOUBLE, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanPrecision')

      call clfftSetPlanPrecision(fft_array(jj)%cl_plan_bw, CLFFT_DOUBLE, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanPrecision')

      ! set number of transforms to 1

      call clfftSetPlanBatchSize(fft_array(jj)%cl_plan_fw, 1_8, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanBatchSize')

      call clfftSetPlanBatchSize(fft_array(jj)%cl_plan_bw, 1_8, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanBatchSize')

      ! set the type precision to double

      call clfftSetPlanPrecision(fft_array(jj)%cl_plan_fw, CLFFT_DOUBLE, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanPrecision')

      call clfftSetPlanPrecision(fft_array(jj)%cl_plan_bw, CLFFT_DOUBLE, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanPrecision')


      ! set the layout

      if(type == FFT_REAL) then

        call clfftSetLayout(fft_array(jj)%cl_plan_fw, CLFFT_REAL, CLFFT_HERMITIAN_INTERLEAVED, status)
        if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetLayout')

        call clfftSetLayout(fft_array(jj)%cl_plan_bw, CLFFT_HERMITIAN_INTERLEAVED, CLFFT_REAL, status)
        if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetLayout')

      else

        call clfftSetLayout(fft_array(jj)%cl_plan_fw, CLFFT_COMPLEX_INTERLEAVED, CLFFT_COMPLEX_INTERLEAVED, status)
        if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetLayout')

        call clfftSetLayout(fft_array(jj)%cl_plan_bw, CLFFT_COMPLEX_INTERLEAVED, CLFFT_COMPLEX_INTERLEAVED, status)
        if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetLayout')

      end if

      ! set the plans as at out of place

      call clfftSetResultLocation(fft_array(jj)%cl_plan_fw, CLFFT_OUTOFPLACE, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetResultLocation')

      call clfftSetResultLocation(fft_array(jj)%cl_plan_bw, CLFFT_OUTOFPLACE, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetResultLocation')

      ! the strides
      
      call clfftSetPlanInStride(fft_array(jj)%cl_plan_fw, fft_dim, int(fft_array(jj)%stride_rs, 8), status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanInStride')

      call clfftSetPlanOutStride(fft_array(jj)%cl_plan_fw, fft_dim, int(fft_array(jj)%stride_fs, 8), status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanOutStride')

      call clfftSetPlanInStride(fft_array(jj)%cl_plan_bw, fft_dim, int(fft_array(jj)%stride_fs, 8), status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanInStride')

      call clfftSetPlanOutStride(fft_array(jj)%cl_plan_bw, fft_dim, int(fft_array(jj)%stride_rs, 8), status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanOutStride')
       
      ! set the scaling factors

      scale = 1.0_8/(product(real(fft_array(jj)%rs_n_global(1:fft_dim), 8)))

      call clfftSetPlanScale(fft_array(jj)%cl_plan_fw, CLFFT_FORWARD, 1.0_8, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanScale')

      call clfftSetPlanScale(fft_array(jj)%cl_plan_fw, CLFFT_BACKWARD, scale, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanScale')

      if(type == FFT_REAL) then
        
        call clfftSetPlanScale(fft_array(jj)%cl_plan_bw, CLFFT_FORWARD, 1.0_8, status)
        if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanScale')
        
        call clfftSetPlanScale(fft_array(jj)%cl_plan_bw, CLFFT_BACKWARD, scale, status)
        if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanScale')

      else
        
        call clfftSetPlanScale(fft_array(jj)%cl_plan_bw, CLFFT_FORWARD, scale, status)
        if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanScale')
        
        call clfftSetPlanScale(fft_array(jj)%cl_plan_bw, CLFFT_BACKWARD, 1.0_8, status)
        if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftSetPlanScale')

      end if

      ! now 'bake' the plans, this signals that the plans are ready to use

      call clfftBakePlan(fft_array(jj)%cl_plan_fw, accel%command_queue, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftBakePlan')

      call clfftBakePlan(fft_array(jj)%cl_plan_bw, accel%command_queue, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clfftBakePlan')

#endif

    case default
      call messages_write('Invalid FFT library.')
      call messages_fatal()
    end select
    
    this = fft_array(jj)

    ! Write information
    if (.not. (library_ == FFTLIB_NFFT .or. library_ == FFTLIB_PNFFT)) then
      call messages_write('Info: FFT grid dimensions       =')
      number_points = 1
      do idir = 1, dim
        call messages_write(fft_array(jj)%rs_n_global(idir))
        if(idir < dim) call messages_write(" x ")
        ! do the multiplication in a integer(8) to avoid overflow for large grids
        number_points = number_points * fft_array(jj)%rs_n_global(idir)
      end do
      call messages_new_line()

      call messages_write('      Total grid size           =')
      call messages_write(number_points)
      call messages_write(' (')
      call messages_write(number_points*CNST(8.0), units = unit_megabytes, fmt = '(f6.1)')
      call messages_write(' )')
      if(any(nn(1:fft_dim) /= nn_temp(1:fft_dim))) then
        call messages_new_line()
        call messages_write('      Inefficient FFT grid. A better grid would be: ')
        do idir = 1, fft_dim
          call messages_write(nn_temp(idir))
        end do
      end if
      call messages_info()
    end if
    
    select case (library_)
    case (FFTLIB_PFFT)
      write(message(1),'(a)') "Info: FFT library = PFFT"
      write(message(2),'(a)') "Info: PFFT processor grid"
      write(message(3),'(a, i9)') " No. of processors                = ", mpi_grp_%size
      write(message(4),'(a, i9)') " No. of columns in the proc. grid = ", column_size
      write(message(5),'(a, i9)') " No. of rows    in the proc. grid = ", row_size
      write(message(6),'(a, i9)') " The size of integer is = ", C_INTPTR_T
      call messages_info(6)

    case (FFTLIB_PNFFT)
      call messages_write("Info: FFT library = PNFFT")
      call messages_info()
      call pnfft_write_info(fft_array(jj)%pnfft)
      
    case (FFTLIB_NFFT)
      call messages_write("Info: FFT library = NFFT")
      call messages_info()
      call nfft_write_info(fft_array(jj)%nfft)

    end select

    POP_SUB(fft_init)
  end subroutine fft_init
  
  ! ---------------------------------------------------------
  !> Some fft-libraries (only NFFT for the moment) need an additional 
  !! precomputation stage that depends on the spatial grid whose size 
  !! may change after fft_init
  subroutine fft_init_stage1(this, namespace, XX, nn)
    type(fft_t),       intent(inout) :: this     !< FFT data type
    !> NFFT spatial nodes on x-axis XX(:,1), y-axis XX(:,2),
    !! and z-axis XX(:,3) 
    type(namespace_t), intent(in)    :: namespace
    FLOAT,             intent(in)    :: XX(:,:)  
    integer, optional, intent(in)    :: nn(:)  
 
    integer :: slot

    PUSH_SUB(fft_init_stage1)

    ASSERT(size(XX,2) == 3)

    slot = this%slot
    select case (fft_array(slot)%library)
    case (FFTLIB_FFTW)
    !Do nothing 
    case (FFTLIB_NFFT)
      ASSERT(present(nn))
      call nfft_precompute(fft_array(slot)%nfft, &
          XX(1:nn(1),1), XX(1:nn(2),2), XX(1:nn(3),3)) 

    case (FFTLIB_PFFT)
    !Do nothing 
    case(FFTLIB_ACCEL)
    !Do nothing 
    case(FFTLIB_PNFFT)
      call pnfft_set_sp_nodes(fft_array(slot)%pnfft, namespace, XX)

    case default
      call messages_write('Invalid FFT library.')
      call messages_fatal()
    end select



    POP_SUB(fft_init_stage1)
  end subroutine fft_init_stage1
  ! ---------------------------------------------------------
  subroutine fft_end(this)
    type(fft_t), intent(inout) :: this

    integer :: ii
#ifdef HAVE_CLFFT
    integer :: status
#endif

    PUSH_SUB(fft_end)

    ii = this%slot
    if(fft_refs(ii) == FFT_NULL) then
      message(1) = "Trying to deallocate FFT that has not been allocated."
      call messages_warning(1)
    else
      if(fft_refs(ii) > 1) then
        fft_refs(ii) = fft_refs(ii) - 1
      else
        select case (fft_array(ii)%library)
        case (FFTLIB_FFTW)
          call fftw_destroy_plan(fft_array(ii)%planf)
          call fftw_destroy_plan(fft_array(ii)%planb)

          if(this%aligned_memory) then
            call fftw_free_memory(this%type == FFT_REAL, &
              fft_array(ii)%drs_data, fft_array(ii)%zrs_data, fft_array(ii)%fs_data)
          end if

        case (FFTLIB_PFFT)
#ifdef HAVE_PFFT
          call pfft_destroy_plan(fft_array(ii)%planf)
          call pfft_destroy_plan(fft_array(ii)%planb)
#endif
          SAFE_DEALLOCATE_P(fft_array(ii)%drs_data)
          SAFE_DEALLOCATE_P(fft_array(ii)%zrs_data)
          SAFE_DEALLOCATE_P(fft_array(ii)%fs_data)

        case(FFTLIB_ACCEL)
#ifdef HAVE_CUDA
          call cuda_fft_destroy(fft_array(ii)%cuda_plan_fw)
          call cuda_fft_destroy(fft_array(ii)%cuda_plan_bw)
#endif
#ifdef HAVE_CLFFT
          call clfftDestroyPlan(fft_array(ii)%cl_plan_fw, status)
          call clfftDestroyPlan(fft_array(ii)%cl_plan_bw, status)
#endif

        case(FFTLIB_NFFT)
          call nfft_end(fft_array(ii)%nfft)
          
        case(FFTLIB_PNFFT)
          call pnfft_end(fft_array(ii)%pnfft)
          
        end select
        fft_refs(ii) = FFT_NULL
      end if
    end if
    this%slot = 0

    POP_SUB(fft_end)
  end subroutine fft_end

  ! ---------------------------------------------------------
  subroutine fft_copy(fft_i, fft_o)
    type(fft_t), intent(in)    :: fft_i
    type(fft_t), intent(inout) :: fft_o

    PUSH_SUB(fft_copy)

    if (fft_o%slot > 0) then
      call fft_end(fft_o)
    end if
    ASSERT(fft_i%slot>=1.and.fft_i%slot<=FFT_MAX)
    ASSERT(fft_refs(fft_i%slot) > 0)

    fft_o = fft_i
    fft_refs(fft_i%slot) = fft_refs(fft_i%slot) + 1

    POP_SUB(fft_copy)
  end subroutine fft_copy

  ! ---------------------------------------------------------
  subroutine fft_get_dims(fft, rs_n_global, fs_n_global, rs_n, fs_n, rs_istart, fs_istart)
    type(fft_t), intent(in)  :: fft
    integer,     intent(out) :: rs_n_global(1:3)
    integer,     intent(out) :: fs_n_global(1:3)
    integer,     intent(out) :: rs_n(1:3)
    integer,     intent(out) :: fs_n(1:3)
    integer,     intent(out) :: rs_istart(1:3)
    integer,     intent(out) :: fs_istart(1:3)

    integer :: slot

    PUSH_SUB(fft_get_dims)

    slot = fft%slot
    rs_n_global(1:3) = fft_array(slot)%rs_n_global(1:3)
    fs_n_global(1:3) = fft_array(slot)%fs_n_global(1:3)
    rs_n(1:3) = fft_array(slot)%rs_n(1:3)
    fs_n(1:3) = fft_array(slot)%fs_n(1:3)
    rs_istart(1:3) = fft_array(slot)%rs_istart(1:3)
    fs_istart(1:3) = fft_array(slot)%fs_istart(1:3)

    POP_SUB(fft_get_dims)
  end subroutine fft_get_dims

  ! ---------------------------------------------------------
  !> convert between array index and G-vector
  function pad_feq(ii, nn, mode)
    integer, intent(in) :: ii,nn
    logical, intent(in) :: mode
    integer :: pad_feq

    ! no push_sub: called too frequently

    if(mode) then      ! index to frequency number
      if( ii <= nn/2 + 1 ) then
        pad_feq = ii - 1
      else
        pad_feq = ii - nn - 1
      end if
    else
      if( ii >= 0 ) then
        pad_feq = ii + 1
      else
        pad_feq = ii + nn + 1
      end if
    end if

    return
  end function pad_feq

  ! -------------------------------------------------------

  integer function fft_size(size, factors)
    integer, intent(in) :: size
    integer, intent(in) :: factors(:)

    integer :: nfactors
    integer :: nondiv
    integer, allocatable :: exponents(:)

    PUSH_SUB(fft_size)

    nfactors = ubound(factors, dim = 1)

    SAFE_ALLOCATE(exponents(1:nfactors))

    fft_size = size
    do 
      call get_exponents(fft_size, nfactors, factors, exponents, nondiv)
      if(nondiv == 1) exit
      fft_size = fft_size + 1
    end do

    SAFE_DEALLOCATE_A(exponents)

    POP_SUB(fft_size)
  end function fft_size

  ! -------------------------------------------------------

  subroutine get_exponents(num, nfactors, factors, exponents, nondiv)
    integer, intent(in)  :: num
    integer, intent(in)  :: nfactors
    integer, intent(in)  :: factors(:)
    integer, intent(out) :: exponents(:)
    integer, intent(out) :: nondiv

    integer :: ifactor

    PUSH_SUB(get_exponents)

    nondiv = num
    do ifactor = 1, nfactors
      exponents(ifactor) = 0
      do
        if(mod(nondiv, factors(ifactor)) /= 0) exit
        nondiv = nondiv/factors(ifactor)
        exponents(ifactor) = exponents(ifactor) + 1
      end do
    end do
    
    POP_SUB(get_exponents)
  end subroutine get_exponents


  ! ----------------------------------------------------------

  subroutine fft_operation_count(fft)
    type(fft_t), intent(in)  :: fft

    real(8) :: fullsize

    PUSH_SUB(fft_operation_count)

    fullsize = product(TOFLOAT(fft%fs_n(1:3)))
    call profiling_count_operations(CNST(5.0)*fullsize*log(fullsize)/log(M_TWO))

    POP_SUB(fft_operation_count)
  end subroutine fft_operation_count


  ! ----------------------------------------------------------
  
  !> This function returns the factor required to normalize a function
  !> after a forward and backward transform.  
  FLOAT pure function fft_scaling_factor(fft) result(scaling_factor)
    type(fft_t), intent(in)  :: fft

    ! for the moment this factor is handled by the backwards transform for most libraries
    scaling_factor = M_ONE
    
    select case (fft_array(fft%slot)%library)
    case(FFTLIB_ACCEL)
#ifdef HAVE_CUDA
      scaling_factor = M_ONE/TOFLOAT(fft_array(fft%slot)%rs_n_global(1))
      scaling_factor = scaling_factor/TOFLOAT(fft_array(fft%slot)%rs_n_global(2))
      scaling_factor = scaling_factor/TOFLOAT(fft_array(fft%slot)%rs_n_global(3))
#endif
    end select
  
  end function fft_scaling_factor

#include "undef.F90"
#include "real.F90"
#include "fft_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "fft_inc.F90"

end module fft_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

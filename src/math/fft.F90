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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

#define FFT_MAX 10
#define FFT_NULL -1

!> Fast Fourier Transform module.
!! This module provides a single interface that works with different
!! FFT implementations.
module fft_m
  use c_pointer_m
#ifdef HAVE_OPENCL  
  use cl
#ifdef HAVE_CLAMDFFT
  use clAmdFft
#endif
#endif
  use datasets_m
  use fftw_m
  use global_m
  use lalg_basic_m
  use loct_math_m
  use math_m
  use messages_m
  use mpi_m
#if defined(HAVE_NFFT)
  use nfft_m
#endif
#if defined(HAVE_OPENMP) && defined(HAVE_FFTW3_THREADS)
  use omp_lib
#endif
  use octcl_kernel_m
  use opencl_m
  use parser_m
  use pfft_m
  use pfft_params_m
  use profiling_m
  use types_m
  use varinfo_m

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
    dfft_forward_cl,   &
    zfft_forward_cl,   &
    dfft_backward_cl,  &
    zfft_backward_cl,  &
    dfft_forward1,     &
    zfft_forward1,     &
    dfft_backward1,    &
    zfft_backward1


  !> global constants
  integer, public, parameter :: &
       FFT_NONE    = 0,         &
       FFT_REAL    = 1,         &
       FFT_COMPLEX = 2

  integer, public, parameter :: &
       FFTLIB_NONE  = 0, &
       FFTLIB_FFTW  = 1, &
       FFTLIB_PFFT  = 2, &
       FFTLIB_CLAMD = 3, &
       FFTLIB_NFFT  = 4
       

  type fft_t
    integer     :: slot    !< in which slot do we have this fft

    integer     :: type    !< is the fft real or complex
    integer     :: library !< what library are we using (FFTLIB_FFTW or FFTLIB_PFFT)

    integer     :: comm           !< MPI communicator
    integer     :: rs_n_global(3) !< total size of the fft in each direction in real space
    integer     :: fs_n_global(3) !< total size of the fft in each direction in fourier space
    integer     :: rs_n(3)        !< local size of the fft in in each direction real space
    integer     :: fs_n(3)        !< local size of the fft in in each direction fourier space
    integer     :: rs_istart(1:3) !< where does the local portion of the function start in real space
    integer     :: fs_istart(1:3) !< where does the local portion of the function start in fourier space

    type(c_ptr) :: planf                  !< plan for forward transform
    type(c_ptr) :: planb                  !< plan for backward transform
    integer(ptrdiff_t_kind) :: pfft_planf !< PFFT plan for forward transform
    integer(ptrdiff_t_kind) :: pfft_planb !< PFFT plan for backward transform

    !> The following arrays have to be stored here and allocated in the initialization routine because of PFFT 
    FLOAT, pointer :: drs_data(:,:,:) !< array used to store the function in real space that is passed to PFFT.
    CMPLX, pointer :: zrs_data(:,:,:) !< array used to store the function in real space that is passed to PFFT.
    CMPLX, pointer ::  fs_data(:,:,:) !< array used to store the function in fourier space that is passed to PFFT
    logical        :: cl_use_real
#ifdef HAVE_CLAMDFFT
    !> data for clAmdFft
    type(clAmdFftPlanHandle) :: cl_plan_fw 
    type(clAmdFftPlanHandle) :: cl_plan_bw !< for real transforms we need a different plan, so we always use 2
#endif
#ifdef HAVE_NFFT
    type(nfft_t) :: nfft 
#endif
  end type fft_t

  integer, save     :: fft_refs(FFT_MAX)
  type(fft_t), save :: fft_array(FFT_MAX)
  logical           :: fft_optimize
  integer           :: fft_prepare_plan

contains

  ! ---------------------------------------------------------
  !> initialize the table
  subroutine fft_all_init()
    integer :: ii, iret
    FLOAT   :: time_limit

    PUSH_SUB(fft_all_init)

    !%Variable FFTOptimize
    !%Type logical
    !%Default yes
    !%Section Mesh::FFTs
    !%Description
    !% Should <tt>octopus</tt> optimize the FFT dimensions? 
    !% This means that the mesh to which FFTs are applied is not taken to be as small
    !% as possible: some points may be added to each direction in order to get a "good number"
    !% for the performance of the FFT algorithm.
    !% The best FFT grid dimensions are given by <math>2^a*3^b*5^c*7^d*11^e*13^f</math>
    !% where <math>a,b,c,d</math> are arbitrary and <math>e,f</math> are 0 or 1.
    !% (http://www.fftw.org/doc/Complex-DFTs.html)
    !% In some cases, namely when using
    !% the split-operator, or Suzuki-Trotter propagators, this option should be turned off.
    !% For spatial FFTs in periodic directions, the grid is never optimized, but a warning will
    !% be written if the number is not good, with a suggestion of a better one to use, so you
    !% can try a different spacing if you want to get a good number.
    !%End
    call parse_logical(datasets_check('FFTOptimize'), .true., fft_optimize)
    do ii = 1, FFT_MAX
      fft_refs(ii) = FFT_NULL
    end do

    !%Variable FFTPreparePlan
    !%Type integer
    !%Default fftw_measure
    !%Section Mesh::FFTs
    !%Description
    !% The FFTs are performed in octopus with the help of the FFTW package (http://www.fftw.org).
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
    !%End
    call parse_integer(datasets_check('FFTPreparePlan'), FFTW_MEASURE, fft_prepare_plan)
    if(.not. varinfo_valid_option('FFTPreparePlan', fft_prepare_plan)) call input_error('FFTPreparePlan')
     
!    !%Variable FFTPlanTimeLimit
!    !%Type float
!    !%Default -1
!    !%Section Mesh::FFTs
!    !%Description
!    !% Sometimes the FFTW_MEASURE takes a lot of time to compute
!    !% many different plans. With this variable is possible to limit
!    !% the time (in seconds) that has roughly going to use for the
!    !% creation of the plan. If a negative value (default one) is
!    !% assigned, there is no restriction.
!    !%End   
!    call parse_float(datasets_check('FFTPlanTimeLimit'), -M_ONE, time_limit)    
!    call fftw_set_timelimit(time_limit)

#if defined(HAVE_OPENMP) && defined(HAVE_FFTW3_THREADS)
    if(omp_get_max_threads() > 1) then

      call messages_write('Info: Initializing Multi-threaded FFTW')
      call messages_info()
      
      call fftw_init_threads(iret)
      if (iret == 0) then 
        call messages_write('Initialization of FFTW3 threads failed.')
        call messages_fatal()
      end if
      call fftw_plan_with_nthreads(omp_get_max_threads())

    end if
#endif

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
    call pfft_cleanup
#endif

#if defined(HAVE_OPENMP) && defined(HAVE_FFTW3_THREADS)
    call fftw_cleanup_threads()
#else
    call fftw_cleanup()
#endif

    POP_SUB(fft_all_end)
  end subroutine fft_all_end

  ! ---------------------------------------------------------
  subroutine fft_init(this, nn, dim, type, library, optimize, optimize_parity, mpi_comm, mpi_grp)
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

    integer :: ii, jj, fft_dim, idir, column_size, row_size, alloc_size, ierror, n3
    integer :: n_1, n_2, n_3, nn_temp(3), status
    integer :: library_
    character(len=100) :: str_tmp
    type(mpi_grp_t) :: mpi_grp_

#ifdef HAVE_CLAMDFFT
    integer(8), allocatable :: stride_rs(:), stride_fs(:)
    real(8) :: scale
#endif

    PUSH_SUB(fft_init)

    ASSERT(type == FFT_REAL .or. type == FFT_COMPLEX)

    mpi_grp_ = mpi_world
    if(present(mpi_grp)) mpi_grp_ = mpi_grp

    ! First, figure out the dimensionality of the FFT.
    fft_dim = 0
    do ii = 1, dim
      if(nn(ii) <= 1) exit
      fft_dim = fft_dim + 1
    end do

    if(fft_dim .eq. 0) then
      message(1) = "Internal error in fft_init: apparently, a 1x1x1 FFT is required."
      call messages_fatal(1)
    end if

    if(fft_dim > 3) call messages_not_implemented('FFT for dimension > 3')

    library_ = library

    select case (library_)
    case (FFTLIB_CLAMD)
    
      nn_temp(1:fft_dim) = nn(1:fft_dim)
      do ii = 1, fft_dim
        ! the AMD OpenCL FFT only supports sizes 2, 3 and 5, but size
        ! 5 gives an fpe error on the Radeon 7970 (APPML 1.8), so we
        ! only use factors 2 and 3
        nn_temp(ii) = fft_size(nn(ii), (/2, 3/))
        if(fft_optimize .and. optimize(ii)) nn(ii) = nn_temp(ii)
      end do 
      
      ! if we can't optimize, in some cases we can't use the library
      if(any(nn(1:fft_dim) /= nn_temp(1:fft_dim))) then
        call messages_write('Invalid grid size for clAmdFft. FFTW will be used instead.')
        call messages_warning()
        library_ = FFTLIB_FFTW
      endif
      
    case (FFTLIB_NFFT)
            
      do ii = 1, fft_dim
        !NFFT likes even grids
        !The underlying FFT grids are optimized inside the nfft_init routine
        if(int(nn(ii)/2)*2 .ne. nn(ii)) nn(ii)=nn(ii)+1 
      end do 
          
    case default

      if(fft_dim < 3 .and. library_ == FFTLIB_PFFT) &
           call messages_not_implemented('PFFT support for dimension < 3')

      ! FFT optimization
      if(any(optimize_parity(1:fft_dim) > 1)) then
        message(1) = "Internal error in fft_init: optimize_parity must be negative, 0, or 1."
        call messages_fatal(1)
      endif
      
      nn_temp(1:fft_dim) = nn(1:fft_dim)
      do ii = 1, fft_dim
        call loct_fft_optimize(nn_temp(ii), optimize_parity(ii))
        if(fft_optimize .and. optimize(ii)) nn(ii) = nn_temp(ii)
      end do
      
      if(any(nn(1:fft_dim) /= nn_temp(1:fft_dim))) then
        write(message(1),'(a,9i6)') "Inefficient FFT grid. A better grid would be", nn_temp(1:fft_dim)
        call messages_warning(1)
      endif
    
    end select

    ! find out if fft has already been allocated
    jj = 0
    do ii = FFT_MAX, 1, -1
      if(fft_refs(ii) /= FFT_NULL) then
        if(all(nn(1:dim) == fft_array(ii)%rs_n_global(1:dim)) .and. type == fft_array(ii)%type &
             .and. library_ == fft_array(ii)%library .and. library_ /= FFTLIB_NFFT) then
             !FIXME: For the moment NFFT plans are always allocated from scratch since they 
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
      message(2) = "Please increase FFT_MAX in fftw3.F90 and recompile."
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

    ! Initialize parallel communicator
    if (library_ == FFTLIB_PFFT) then
#ifdef HAVE_PFFT
      call pfft_init()
 
      call pfft_decompose(mpi_grp_%size, column_size, row_size)
 
      call pfft_create_procmesh_2d(ierror, mpi_grp_%comm, column_size, row_size, fft_array(jj)%comm)        
   
      if (ierror .ne. 0) then
        message(1) = "The number of rows and columns in PFFT processor grid is not equal to "
        message(2) = "the number of processor in the MPI communicator."
        message(3) = "Please check it."
        call messages_fatal(3)
      end if
#endif
    else
      fft_array(jj)%comm = -1
    end if
    if (present(mpi_comm)) mpi_comm = fft_array(jj)%comm

    ! Get dimentions of arrays
    select case (library_)
    case (FFTLIB_FFTW)
      call fftw_get_dims(fft_array(jj)%rs_n_global, type == FFT_REAL, fft_array(jj)%fs_n_global)
      fft_array(jj)%rs_n = fft_array(jj)%rs_n_global
      fft_array(jj)%fs_n = fft_array(jj)%fs_n_global
      fft_array(jj)%rs_istart = 1
      fft_array(jj)%fs_istart = 1

    case (FFTLIB_PFFT)
#ifdef HAVE_PFFT     
      call pfft_get_dims(fft_array(jj)%rs_n_global, mpi_comm, type == FFT_REAL, &
           alloc_size, fft_array(jj)%fs_n_global, fft_array(jj)%rs_n, &
           fft_array(jj)%fs_n, fft_array(jj)%rs_istart, fft_array(jj)%fs_istart)
#endif

      ! Allocate memory. Note that PFFT may need extra memory space 
      ! and that in fourier space the function will be transposed
      if (type == FFT_REAL) then
        n_1 = max(1, fft_array(jj)%rs_n(1))
        n_2 = max(1, fft_array(jj)%rs_n(2))
        n_3 = max(1, fft_array(jj)%rs_n(3))

        n3 = ceiling(real(2*alloc_size)/real(n_1*n_2))
        SAFE_ALLOCATE(fft_array(jj)%drs_data(1:n_1, 1:n_2, 1:n3))

        ! For real functions, PFFT increases the size of rs_n(1) by 1, such that rs_n(1) = nn(1) + 1.
        ! The rest of the code does not need to know about this.
        fft_array(jj)%rs_n(1) = fft_array(jj)%rs_n(1) - 1
      else
        n3 = ceiling(real(alloc_size)/real(fft_array(jj)%rs_n(1)*fft_array(jj)%rs_n(2)))
        SAFE_ALLOCATE(fft_array(jj)%zrs_data(1:fft_array(jj)%rs_n(1), 1:fft_array(jj)%rs_n(2), 1:n3))
      end if

      n_1 = max(1, fft_array(jj)%fs_n(1))
      n_2 = max(1, fft_array(jj)%fs_n(2))
      n_3 = max(1, fft_array(jj)%fs_n(3))

      n3 = ceiling(real(alloc_size)/real(n_3*n_1))
      SAFE_ALLOCATE(fft_array(jj)%fs_data(1:n_3, 1:n_1, 1:n3))

    case(FFTLIB_CLAMD)
      fft_array(jj)%cl_use_real = .false.
      fft_array(jj)%cl_use_real = .true.
      call fftw_get_dims(fft_array(jj)%rs_n_global, (fft_array(jj)%cl_use_real .and. type == FFT_REAL), fft_array(jj)%fs_n_global)
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
      
    end select

    ! Prepare plans
    select case (library_)
    case (FFTLIB_FFTW)
      call fftw_prepare_plan(fft_array(jj)%planf, fft_dim, fft_array(jj)%rs_n_global, &
           type == FFT_REAL, FFTW_FORWARD, fft_prepare_plan+FFTW_UNALIGNED)
      call fftw_prepare_plan(fft_array(jj)%planb, fft_dim, fft_array(jj)%rs_n_global, &
           type == FFT_REAL, FFTW_BACKWARD, fft_prepare_plan+FFTW_UNALIGNED)

    case(FFTLIB_NFFT)
#ifdef HAVE_NFFT
     call nfft_copy_info(this%nfft,fft_array(jj)%nfft) !copy default parameters set in the calling routine 
     call nfft_init(fft_array(jj)%nfft, fft_array(jj)%rs_n_global, &
                    fft_dim, fft_array(jj)%rs_n_global(1) , type, optimize = .true.)
#endif

    case (FFTLIB_PFFT)
#ifdef HAVE_PFFT     
      if(type == FFT_REAL) then
        call pfft_prepare_plan_r2c(fft_array(jj)%pfft_planf, fft_array(jj)%rs_n_global, fft_array(jj)%drs_data, &
             fft_array(jj)%fs_data, FFTW_FORWARD, fft_prepare_plan, mpi_comm)
        call pfft_prepare_plan_c2r(fft_array(jj)%pfft_planb, fft_array(jj)%rs_n_global, fft_array(jj)%fs_data, &
             fft_array(jj)%drs_data, FFTW_BACKWARD, fft_prepare_plan, mpi_comm)
      else
        call pfft_prepare_plan_c2c(fft_array(jj)%pfft_planf, fft_array(jj)%rs_n_global, fft_array(jj)%zrs_data, &
             fft_array(jj)%fs_data, FFTW_FORWARD, fft_prepare_plan, mpi_comm)
        call pfft_prepare_plan_c2c(fft_array(jj)%pfft_planb, fft_array(jj)%rs_n_global, fft_array(jj)%fs_data, &
             fft_array(jj)%zrs_data, FFTW_BACKWARD, fft_prepare_plan, mpi_comm)
      end if
#endif
    case(FFTLIB_CLAMD)
#ifdef HAVE_CLAMDFFT

      ! create the plans
      call clAmdFftCreateDefaultPlan(fft_array(jj)%cl_plan_fw, opencl%context, fft_dim, int(fft_array(jj)%rs_n_global, 8), status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftCreateDefaultPlan')

      call clAmdFftCreateDefaultPlan(fft_array(jj)%cl_plan_bw, opencl%context, fft_dim, int(fft_array(jj)%rs_n_global, 8), status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftCreateDefaultPlan')

      ! set precision

      call clAmdFftSetPlanPrecision(fft_array(jj)%cl_plan_fw, CLFFT_DOUBLE, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanPrecision')

      call clAmdFftSetPlanPrecision(fft_array(jj)%cl_plan_bw, CLFFT_DOUBLE, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanPrecision')

      ! set number of transforms to 1

      call clAmdFftSetPlanBatchSize(fft_array(jj)%cl_plan_fw, 1_8, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanBatchSize')

      call clAmdFftSetPlanBatchSize(fft_array(jj)%cl_plan_bw, 1_8, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanBatchSize')

      ! set the type precision to double

      call clAmdFftSetPlanPrecision(fft_array(jj)%cl_plan_fw, CLFFT_DOUBLE, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanPrecision')

      call clAmdFftSetPlanPrecision(fft_array(jj)%cl_plan_bw, CLFFT_DOUBLE, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanPrecision')


      ! set the layout

      if(fft_array(jj)%cl_use_real .and. type == FFT_REAL) then

        call clAmdFftSetLayout(fft_array(jj)%cl_plan_fw, CLFFT_REAL, CLFFT_HERMITIAN_INTERLEAVED, status)
        if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetLayout')

        call clAmdFftSetLayout(fft_array(jj)%cl_plan_bw, CLFFT_HERMITIAN_INTERLEAVED, CLFFT_REAL, status)
        if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetLayout')

      else

        call clAmdFftSetLayout(fft_array(jj)%cl_plan_fw, CLFFT_COMPLEX_INTERLEAVED, CLFFT_COMPLEX_INTERLEAVED, status)
        if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetLayout')

        call clAmdFftSetLayout(fft_array(jj)%cl_plan_bw, CLFFT_COMPLEX_INTERLEAVED, CLFFT_COMPLEX_INTERLEAVED, status)
        if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetLayout')

      end if

      ! set the plans as at out of place

      call clAmdFftSetResultLocation(fft_array(jj)%cl_plan_fw, CLFFT_OUTOFPLACE, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetResultLocation')

      call clAmdFftSetResultLocation(fft_array(jj)%cl_plan_bw, CLFFT_OUTOFPLACE, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetResultLocation')

      
      ! invert the stride for Fortran arrays
      SAFE_ALLOCATE(stride_rs(1:fft_dim))
      SAFE_ALLOCATE(stride_fs(1:fft_dim))

!      call clAmdFftGetPlanInStride(fft_array(jj)%cl_plan_bw, fft_dim, stride_rs, status)
!      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftGetPlanInStride')


!      call clAmdFftGetPlanInStride(fft_array(jj)%cl_plan_fw, fft_dim, stride_fs, status)
!      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftGetPlanInStride')

!      print*, "STRIDE IN     ", stride_rs
!      print*, "STRIDE OUT    ", stride_fs

      stride_rs(1) = 1
      stride_fs(1) = 1
      do ii = 2, fft_dim
        stride_rs(ii) = stride_rs(ii - 1)*fft_array(jj)%rs_n(ii - 1)
        stride_fs(ii) = stride_fs(ii - 1)*fft_array(jj)%fs_n(ii - 1)
      end do

!      print*, "STRIDE NEW IN ", stride_rs
!      print*, "STRIDE NEW OUT", stride_fs

      ! the strides
      
      call clAmdFftSetPlanInStride(fft_array(jj)%cl_plan_fw, fft_dim, stride_rs, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanInStride')

      call clAmdFftSetPlanOutStride(fft_array(jj)%cl_plan_fw, fft_dim, stride_fs, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanOutStride')

      call clAmdFftSetPlanInStride(fft_array(jj)%cl_plan_bw, fft_dim, stride_fs, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanInStride')

      call clAmdFftSetPlanOutStride(fft_array(jj)%cl_plan_bw, fft_dim, stride_rs, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanOutStride')
       
      SAFE_DEALLOCATE_A(stride_rs)
      SAFE_DEALLOCATE_A(stride_fs)

      ! set the scaling factors

      scale = 1.0_8/(product(real(fft_array(jj)%rs_n_global(1:fft_dim), 8)))

      call clAmdFftSetPlanScale(fft_array(jj)%cl_plan_fw, CLFFT_FORWARD, 1.0_8, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanScale')

      call clAmdFftSetPlanScale(fft_array(jj)%cl_plan_fw, CLFFT_BACKWARD, scale, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanScale')

      call clAmdFftSetPlanScale(fft_array(jj)%cl_plan_bw, CLFFT_FORWARD, 1.0_8, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanScale')

      call clAmdFftSetPlanScale(fft_array(jj)%cl_plan_bw, CLFFT_BACKWARD, scale, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftSetPlanScale')

      ! now 'bake' the plans, this signals that the plans are ready to use

      call clAmdFftBakePlan(fft_array(jj)%cl_plan_fw, opencl%command_queue, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftBakePlan')

      call clAmdFftBakePlan(fft_array(jj)%cl_plan_bw, opencl%command_queue, status)
      if(status /= CLFFT_SUCCESS) call clfft_print_error(status, 'clAmdFftBakePlan')

#endif


    case default
      call messages_write('Invalid FFT library.')
      call messages_fatal()
    end select
    
    this = fft_array(jj)
    
    ! Write information
    write(message(1), '(a)') "Info: FFT allocated with size ("
    do idir = 1, dim
      write(str_tmp, '(i7,a)') fft_array(jj)%rs_n_global(idir)
      if(idir == dim) then
        message(1) = trim(message(1)) // trim(str_tmp) // ") in slot "
      else
        message(1) = trim(message(1)) // trim(str_tmp) // ","
      endif
    enddo
    write(str_tmp, '(i2)') jj
    message(1) = trim(message(1)) // trim(str_tmp)
    select case (library_)
    case (FFTLIB_FFTW)
      message(2) = "Info: FFT library = FFTW3"
      call messages_info(2)
    case (FFTLIB_PFFT)
      write(message(2),'(a)') "Info: FFT library = PFFT"
      write(message(3),'(a)') "Info: PFFT processor grid"
      write(message(4),'(a, i9)') " No. of processors                = ", mpi_grp_%size
      write(message(5),'(a, i9)') " No. of columns in the proc. grid = ", column_size
      write(message(6),'(a, i9)') " No. of rows    in the proc. grid = ", row_size
      write(message(7),'(a, i9)') " The size of integer is = ", ptrdiff_t_kind
      call messages_info(7)
    case (FFTLIB_NFFT)
      message(2) = "Info: FFT library = NFFT"
      call messages_info(2)
#ifdef HAVE_NFFT
      call nfft_write_info(fft_array(jj)%nfft)
#endif
    end select

    POP_SUB(fft_init)
  end subroutine fft_init
  
  ! ---------------------------------------------------------
  !> Some fft-libraries (only NFFT for the moment) need an additional 
  !! precomputation stage that depends on the spatial grid whose size 
  !! may change after fft_init
  subroutine fft_init_stage1(this, XX)
    type(fft_t),       intent(inout) :: this     !< FFT data type
    !> NFFT spatial nodes on x-axis XX(:,1), y-axis XX(:,2),
    !! and z-axis XX(:,3) 
    FLOAT, optional,   intent(in)    :: XX(:,:)  
 
    integer :: slot

    PUSH_SUB(fft_init_stage1)

    ASSERT(size(XX,2) == 3)

    slot = this%slot
    select case (fft_array(slot)%library)
    case (FFTLIB_FFTW)
    !Do nothing 
    case (FFTLIB_NFFT)
#ifdef HAVE_NFFT
      call nfft_precompute(fft_array(slot)%nfft, XX(:,1),XX(:,2),XX(:,3)) 
#endif
    case (FFTLIB_PFFT)
    !Do nothing 
    case(FFTLIB_CLAMD)
    !Do nothing 
    case default
      call messages_write('Invalid FFT library.')
      call messages_fatal()
    end select



    POP_SUB(fft_init_stage1)
  end subroutine fft_init_stage1
  ! ---------------------------------------------------------
  subroutine fft_end(this)
    type(fft_t), intent(inout) :: this

    integer :: ii, status

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
        case (FFTLIB_PFFT)
#ifdef HAVE_PFFT
          call pfft_destroy_plan(fft_array(ii)%pfft_planf)
          call pfft_destroy_plan(fft_array(ii)%pfft_planb)
#endif
          SAFE_DEALLOCATE_P(fft_array(ii)%drs_data)
          SAFE_DEALLOCATE_P(fft_array(ii)%zrs_data)
          SAFE_DEALLOCATE_P(fft_array(ii)%fs_data)
#ifdef HAVE_CLAMDFFT
        case(FFTLIB_CLAMD)
          call clAmdFftDestroyPlan(fft_array(ii)%cl_plan_fw, status)
          call clAmdFftDestroyPlan(fft_array(ii)%cl_plan_bw, status)
#endif

#ifdef HAVE_NFFT
        case(FFTLIB_NFFT)
        call nfft_end(fft_array(ii)%nfft)
#endif
        end select
        fft_refs(ii) = FFT_NULL
        write(message(1), '(a,i4)') "Info: FFT deallocated from slot ", ii
        call messages_info(1)
      end if
    end if

    POP_SUB(fft_end)
  end subroutine fft_end

  ! ---------------------------------------------------------
  subroutine fft_copy(fft_i, fft_o)
    type(fft_t), intent(in)  :: fft_i
    type(fft_t), intent(out) :: fft_o

    PUSH_SUB(fft_copy)

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

    nfactors = ubound(factors, dim = 1)

    SAFE_ALLOCATE(exponents(1:nfactors))

    fft_size = size
    do 
      call get_exponents(fft_size, nfactors, factors, exponents, nondiv)
      if(nondiv == 1) exit
      fft_size = fft_size + 1
    end do

    SAFE_DEALLOCATE_A(exponents)

  end function fft_size

  ! -------------------------------------------------------

  subroutine get_exponents(num, nfactors, factors, exponents, nondiv)
    integer, intent(in)  :: num
    integer, intent(in)  :: nfactors
    integer, intent(in)  :: factors(:)
    integer, intent(out) :: exponents(:)
    integer, intent(out) :: nondiv

    integer :: ifactor

    nondiv = num
    do ifactor = 1, nfactors
      exponents(ifactor) = 0
      do
        if(mod(nondiv, factors(ifactor)) /= 0) exit
        nondiv = nondiv/factors(ifactor)
        exponents(ifactor) = exponents(ifactor) + 1
      end do
    end do
    
  end subroutine get_exponents

#include "undef.F90"
#include "real.F90"
#include "fft_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "fft_inc.F90"

end module fft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

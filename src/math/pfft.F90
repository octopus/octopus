!! Copyright (C) 2011 J. Alberdi, P. Garcia Risueño, M. Oliveira
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

#if defined(SINGLE_PRECISION)
#  define PDFFT(x) dpfftf_ ## x
#else
#  define PDFFT(x) dpfft_ ## x
#endif

!> The includes for the PFFT
module pfft_params_m
  use fftw_m
  implicit none

#ifdef HAVE_PFFT
  include "pfft.f"
#else
  integer, parameter :: ptrdiff_t_kind = 8
#endif
end module pfft_params_m
 
!> The low level module to work with the PFFT library.
!! http://www-user.tu-chemnitz.de/~mpip/software.php?lang=en
module pfft_m
  use global_m
  use math_m
  use pfft_params_m
  implicit none

#ifdef HAVE_PFFT
  !> PFFT initialization routines
  interface pfft_init
    subroutine PDFFT(init)
    end subroutine PDFFT(init)
  end interface pfft_init

  interface pfft_create_procmesh_2d
    subroutine PDFFT(create_procmesh_2d) (ierror, mpi_comm, n1, n2, mpi_comm_2d)
      integer, intent(in)  :: mpi_comm, n1, n2
      integer, intent(out) :: ierror, mpi_comm_2d
    end subroutine PDFFT(create_procmesh_2d)
  end interface pfft_create_procmesh_2d

  interface
    !> PFFT basic interface to get the local size of input and output arrays. Real to complex version
    subroutine PDFFT(local_size_dft_r2c_3d)(alloc_size, n, mpi_comm, pfft_flags, &
         local_ni, local_i_start, local_no, local_o_start)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(inout) :: alloc_size    !< Size that has to be allocated
      integer(ptrdiff_t_kind), intent(inout) :: n             !< The size of the global matrix
      integer,                 intent(in)    :: mpi_comm, pfft_flags
      integer(ptrdiff_t_kind), intent(inout) :: local_ni      !< Local number of elements of the input matrix. Rank = 3
      integer(ptrdiff_t_kind), intent(inout) :: local_i_start !< Local start point of the input matrix. Rank = 3
      integer(ptrdiff_t_kind), intent(inout) :: local_no      !< Local number of elements of the output matrix. Rank = 3
      integer(ptrdiff_t_kind), intent(inout) :: local_o_start !< Local start point of the output matrix. Rank = 3
    end subroutine PDFFT(local_size_dft_r2c_3d)
    
    !> PFFT advanced interface to get the local size of input and output arrays. Real to complex version
    subroutine PDFFT(local_size_many_dft_r2c)(alloc_size, rank_n, n, ni, no, howmany, iblock, oblock, &
         mpi_comm, pfft_flags, local_ni, local_i_start, local_no, local_o_start)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(inout) :: alloc_size !< Size that has to be allocated
      integer(ptrdiff_t_kind), intent(in)    :: rank_n     !< The dimensions of the matrices
      integer(ptrdiff_t_kind), intent(inout) :: n       !< The size of the global matrix
      integer(ptrdiff_t_kind), intent(inout) :: ni      !< An integer array with the same length as 'n'. 
                                                        !! It gives the size of the FFT input (non-oversampled)
      integer(ptrdiff_t_kind), intent(inout) :: no      !< An integer array with the same length as 'n'. 
                                                        !! It gives the size of the FFT output
      integer(ptrdiff_t_kind), intent(in)    :: howmany !< More than one FFT could be calculated at once, 
                                                        !! but the input has to be stored interleaved - (Just set howmany to 1)
      integer(ptrdiff_t_kind), intent(in)    :: iblock  !< The block size of the input data distribution (Just use the constant
                                                        !! PFFT_DEFAULT_BLOCKS to set the blocks to their default value)
      integer(ptrdiff_t_kind), intent(in)    :: oblock  !< The block size of the output data distribution (Just use the constant
                                                        !! PFFT_DEFAULT_BLOCKS to set the blocks to their default value)
      integer,                 intent(in)    :: mpi_comm, pfft_flags
      integer(ptrdiff_t_kind), intent(inout) :: local_ni      !< Local number of elements of the input matrix. Rank = 3
      integer(ptrdiff_t_kind), intent(inout) :: local_i_start !< Local start point of the input matrix. Rank = 3
      integer(ptrdiff_t_kind), intent(inout) :: local_no      !< Local number of elements of the output matrix. Rank = 3
      integer(ptrdiff_t_kind), intent(inout) :: local_o_start !< Local start point of the output matrix. Rank = 3
    end subroutine PDFFT(local_size_many_dft_r2c)

    !> PFFT basic interface to get the local size of input and output arrays. Complex to complex version
    subroutine PDFFT(local_size_dft_3d)(alloc_size, n, mpi_comm, pfft_flags, &
         local_ni, local_i_start, local_no, local_o_start)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(inout) :: alloc_size    !< Size that has to be allocated
      integer(ptrdiff_t_kind), intent(inout) :: n             !< The size of the global matrix
      integer,                 intent(in)    :: mpi_comm, pfft_flags
      integer(ptrdiff_t_kind), intent(inout) :: local_ni      !< Local number of elements of the input matrix. Rank = 3
      integer(ptrdiff_t_kind), intent(inout) :: local_i_start !< Local start point of the input matrix. Rank = 3
      integer(ptrdiff_t_kind), intent(inout) :: local_no      !< Local number of elements of the output matrix. Rank = 3
      integer(ptrdiff_t_kind), intent(inout) :: local_o_start !< Local start point of the output matrix. Rank = 3
    end subroutine PDFFT(local_size_dft_3d)

    !> PFFT advanced interface to get the local size of input and output arrays
    subroutine PDFFT(local_size_many_dft)(alloc_size, rank_n, n, ni, no, howmany, iblock, oblock, &
         mpi_comm, pfft_flags, &
         local_ni, local_i_start, local_no, local_o_start)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(inout) :: alloc_size    !< Size that has to be allocated
      integer(ptrdiff_t_kind), intent(in)    :: rank_n        !< The dimensions of the matrices
      integer(ptrdiff_t_kind), intent(inout) :: n       !< The size of the global matrix
      integer(ptrdiff_t_kind), intent(inout) :: ni      !< An integer array with the same length as 'n'. 
                                                        !! It gives the size of the FFT input (non-oversampled)
      integer(ptrdiff_t_kind), intent(inout) :: no      !< An integer array with the same length as 'n'. 
                                                        !! It gives the size of the FFT output
      integer(ptrdiff_t_kind), intent(in)    :: howmany !< More than one FFT could be calculated at once, 
                                                        !! but the input has to be stored interleaved - (Just set howmany to 1)
      integer(ptrdiff_t_kind), intent(in)    :: iblock  !< The block size of the input data distribution (Just use the constant
                                                        !! PFFT_DEFAULT_BLOCKS to set the blocks to their default value)
      integer(ptrdiff_t_kind), intent(in)    :: oblock  !< The block size of the output data distribution(Just use the constant
                                                        !! PFFT_DEFAULT_BLOCKS to set the blocks to their default value)
      integer,                 intent(in)    :: mpi_comm, pfft_flags
      integer(ptrdiff_t_kind), intent(inout) :: local_ni      !< Local number of elements of the input matrix. Rank = 3
      integer(ptrdiff_t_kind), intent(inout) :: local_i_start !< Local start point of the input matrix. Rank = 3
      integer(ptrdiff_t_kind), intent(inout) :: local_no      !< Local number of elements of the output matrix. Rank = 3
      integer(ptrdiff_t_kind), intent(inout) :: local_o_start !< Local start point of the output matrix. Rank = 3
    end subroutine PDFFT(local_size_many_dft)

    !> PFFT simple interface to create plan from real to complex
    subroutine PDFFT(plan_dft_r2c_3d)(plan, n, in, out, mpi_comm, sign, pfft_flags, fftw_flags)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(out)   :: plan    !< The plan that is created by PFFT
      integer(ptrdiff_t_kind), intent(inout) :: n       !< The size of the global matrix
      FLOAT,                   intent(in)    :: in      !< The input matrix that is going to be used to do the transform
      CMPLX,                   intent(in)    :: out     !< The output matrix that is going to be used to do the transform
      integer,                 intent(in)    :: mpi_comm, sign, pfft_flags, fftw_flags
    end subroutine PDFFT(plan_dft_r2c_3d)

    !> PFFT simple interface to create plan from complex to real
    subroutine PDFFT(plan_dft_c2r_3d)(plan, n, in, out, mpi_comm, sign, pfft_flags, fftw_flags)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(out)   :: plan    !< The plan that is created by PFFT
      integer(ptrdiff_t_kind), intent(inout) :: n       !< The size of the global matrix
      CMPLX,                   intent(in)    :: in      !< The input matrix that is going to be used to do the transform
      FLOAT,                   intent(in)    :: out     !< The output matrix that is going to be used to do the transform
      integer,                 intent(in)    :: mpi_comm, sign, pfft_flags, fftw_flags
    end subroutine PDFFT(plan_dft_c2r_3d)

    !> Advanced interface for creating plan of FFT from real to complex
    subroutine PDFFT(plan_many_dft_r2c)(plan, rank_n, n, ni, no, howmany, &
         iblock, oblock, in, out, mpi_comm, sign, pfft_flags, fftw_flags)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(out)   :: plan    !< The plan that is created by PFFT
      integer(ptrdiff_t_kind), intent(in)    :: rank_n  !< The dimensions of the matrices
      integer(ptrdiff_t_kind), intent(inout) :: n       !< The size of the global matrix
      integer(ptrdiff_t_kind), intent(inout) :: ni      !< An integer array with the same length as 'n'. 
                                                        !! It gives the size of the FFT input (non-oversampled)
      integer(ptrdiff_t_kind), intent(inout) :: no      !< An integer array with the same length as 'n'. 
                                                        !! It gives the size of the FFT output
      integer(ptrdiff_t_kind), intent(in)    :: howmany !< More than one FFT could be calculated at once, 
                                                        !! but the input has to be stored interleaved - (Just set howmany to 1)
      integer(ptrdiff_t_kind), intent(in)    :: iblock  !< The block size of the input data distribution (Just use the constant
                                                        !! PFFT_DEFAULT_BLOCKS to set the blocks to their default value)
      integer(ptrdiff_t_kind), intent(in)    :: oblock  !< The block size of the output data distribution (Just use the constant 
                                                        !! PFFT_DEFAULT_BLOCKS to set the blocks to their default value)
      FLOAT,                   intent(in)    :: in      !< The input matrix that is going to be used to do the transform
      CMPLX,                   intent(in)    :: out     !< The output matrix that is going to be used to do the transform
      integer,                 intent(in)    :: mpi_comm, sign, pfft_flags, fftw_flags
    end subroutine PDFFT(plan_many_dft_r2c) 

    !> Advanced interface for creating plan of FFT from complex to real
    subroutine PDFFT(plan_many_dft_c2r)(plan, rank_n, n, ni, no, howmany, &
         iblock, oblock, in, out, mpi_comm, sign, pfft_flags, fftw_flags)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(out)   :: plan    !< The plan that is created by PFFT
      integer(ptrdiff_t_kind), intent(in)    :: rank_n  !< The dimensions of the matrices
      integer(ptrdiff_t_kind), intent(inout) :: n       !< The size of the global matrix
      integer(ptrdiff_t_kind), intent(inout) :: ni      !< An integer array with the same length as 'n'. 
                                                        !! It gives the size of the FFT input (non-oversampled)
      integer(ptrdiff_t_kind), intent(inout) :: no      !< An integer array with the same length as 'n'. 
                                                        !! It gives the size of the FFT output
      integer(ptrdiff_t_kind), intent(in)    :: howmany !< More than one FFT could be calculated at once, 
                                                        !! but the input has to be stored interleaved - (Just set howmany to 1)
      integer(ptrdiff_t_kind), intent(in)    :: iblock  !< The block size of the input data distribution (Just use the constant
                                                        !! PFFT_DEFAULT_BLOCKS to set the blocks to their default value)
      integer(ptrdiff_t_kind), intent(in)    :: oblock  !< The block size of the output data distribution (Just use the constant 
                                                        !! PFFT_DEFAULT_BLOCKS to set the blocks to their default value)
      CMPLX,                   intent(in)    :: in      !< The input matrix that is going to be used to do the transform
      FLOAT,                   intent(in)    :: out     !< The output matrix that is going to be used to do the transform
      integer,                 intent(in)    :: mpi_comm, sign, pfft_flags, fftw_flags
    end subroutine PDFFT(plan_many_dft_c2r)

    !> Advanced interface for creating plan of FFT
    subroutine PDFFT(plan_dft_3d)(plan, n, in, out, mpi_comm, sign, pfft_flags, fftw_flags)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(out)   :: plan    !< The plan that is created by PFFT
      integer(ptrdiff_t_kind), intent(inout) :: n       !< The size of the global matrix
      CMPLX,                   intent(in)    :: in      !< The input matrix that is going to be used to do the transform
      CMPLX,                   intent(in)    :: out     !< The output matrix that is going to be used to do the transform
      integer,                 intent(in)    :: mpi_comm, sign, pfft_flags, fftw_flags
    end subroutine PDFFT(plan_dft_3d)
  end interface

  interface pfft_execute
    subroutine PDFFT(execute)(plan)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(inout) :: plan
    end subroutine PDFFT(execute)
  end interface pfft_execute

  interface pfft_destroy_plan
    subroutine PDFFT(destroy_plan)(plan)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(inout) :: plan
    end subroutine PDFFT(destroy_plan)
  end interface pfft_destroy_plan

  interface pfft_cleanup
    subroutine PDFFT(cleanup)
    end subroutine PDFFT(cleanup)
  end interface pfft_cleanup

contains

  !> Decompose all available processors in 2D processor grid, most equally possible
  subroutine pfft_decompose(n_proc, dim1, dim2)
    integer, intent(in)  :: n_proc !< Number of processors
    integer, intent(out) :: dim1   !< First out dimension
    integer, intent(out) :: dim2   !< Second out dimension

    integer :: np, i
    
    PUSH_SUB(pfft_decompose)

    ASSERT(n_proc > 0)
    
    dim1 = 1
    dim2 = 1
    np = n_proc
    i = n_proc-1
    
    if (is_prime(n_proc)) then
      dim1 = n_proc
    else
      do
        if (i <= 1) exit
        if(mod(np,i) /= 0.or.(.not. is_prime(i))) then
          i=i-1
          cycle
        end if
        np=np/i
        if (dim1 <= dim2) then
          dim1 = dim1*i
        else
          dim2 = dim2*i
        end if
      end do
    end if

    ASSERT(dim1*dim2 == n_proc)
    
    POP_SUB(pfft_decompose)
  end subroutine pfft_decompose

  !> Octopus subroutine to prepare a PFFT plan real to complex
  subroutine pfft_prepare_plan_r2c(plan, n, in, out, sign, flags, mpi_comm)
    integer(ptrdiff_t_kind), intent(out)   :: plan       !< The plan that is created by PFFT
    integer,                 intent(in)    :: n(:)       !< The size of the global matrix
    FLOAT,   pointer,        intent(inout) :: in(:,:,:)  !< The input matrix that is going to be used to do the transform
    CMPLX,   pointer,        intent(inout) :: out(:,:,:) !< The output matrix that is going to be used to do the transform
    integer,                 intent(in)    :: sign       !< Sign flag to decide FFT direction. Has to be FFTW_FORWARD
    integer,                 intent(in)    :: flags      !< Flags for FFT library. Could be changed with the input variable
                                                         !! FFTPreparePlan. Default value is FFTW_MEASURE
    integer,                 intent(in)    :: mpi_comm   !< MPI communicator
    
    integer(ptrdiff_t_kind) :: tmp_n(3)

    PUSH_SUB(pfft_prepare_plan_r2c)

    ASSERT(sign == FFTW_FORWARD)

    tmp_n(1:3) = n(1:3)

    call PDFFT(plan_dft_r2c_3d) (plan, tmp_n(1), in(1,1,1), out(1,1,1), mpi_comm, sign, PFFT_TRANSPOSED_OUT, flags) 

    POP_SUB(pfft_prepare_plan_r2c)
  end subroutine pfft_prepare_plan_r2c

  !> Octopus subroutine to prepare a PFFT plan real to complex
  subroutine pfft_prepare_plan_c2r(plan, n, in, out, sign, flags, mpi_comm)
    integer(ptrdiff_t_kind), intent(out)   :: plan       !< The plan that is created by PFFT
    integer,                 intent(in)    :: n(:)       !< The size of the global matrix
    CMPLX,   pointer,        intent(inout) :: in(:,:,:)  !< The input matrix that is going to be used to do the transform
    FLOAT,   pointer,        intent(inout) :: out(:,:,:) !< The output matrix that is going to be used to do the transform
    integer,                 intent(in)    :: sign       !< Sign flag to decide FFT direction. Has to be FFTW_BACKWARD
    integer,                 intent(in)    :: flags      !< Flags for FFT library. Could be changed with the input variable
                                                         !! FFTPreparePlan. Default value is FFTW_MEASURE
    integer,                 intent(in)    :: mpi_comm   !< MPI communicator
    
    integer(ptrdiff_t_kind) :: tmp_n(3)
    
    PUSH_SUB(pfft_prepare_plan_c2r)

    ASSERT(sign == FFTW_BACKWARD)

    tmp_n(1:3) = n(1:3)

    call PDFFT(plan_dft_c2r_3d) (plan, tmp_n(1), in(1,1,1), out(1,1,1), mpi_comm, sign, PFFT_TRANSPOSED_IN, flags) 

    POP_SUB(pfft_prepare_plan_c2r)
  end subroutine pfft_prepare_plan_c2r

  
  !> Octopus subroutine to prepare a PFFT plan real to complex
  subroutine pfft_prepare_plan_c2c(plan, n, in, out, sign, flags, mpi_comm)
    integer(ptrdiff_t_kind), intent(out)   :: plan       !< The plan that is created by PFFT
    integer,                 intent(in)    :: n(:)       !< The size of the global matrix
    CMPLX,   pointer,        intent(inout) :: in(:,:,:)  !< The input matrix that is going to be used to do the transform
    CMPLX,   pointer,        intent(inout) :: out(:,:,:) !< The output matrix that is going to be used to do the transform
    integer,                 intent(in)    :: sign       !< Sign flag to decide FFT direction. 
                                                         !! Has to be FFTW_FORWARD or FFTW_BACKWARD
    integer,                 intent(in)    :: flags      !< Flags for FFT library. Could be changed with the input variable
                                                         !! FFTPreparePlan. Default value is FFTW_MEASURE
    integer,                 intent(in)    :: mpi_comm   !< MPI communicator

    integer :: pfft_flags
    integer(ptrdiff_t_kind) :: tmp_n(3)

    PUSH_SUB(pfft_prepare_plan_c2c)

    tmp_n(1:3) = n(1:3)

    if (sign == FFTW_FORWARD) then
      pfft_flags = PFFT_TRANSPOSED_OUT
    else
      pfft_flags = PFFT_TRANSPOSED_IN
    end if

    call PDFFT(plan_dft_3d) (plan, tmp_n(1), in(1,1,1), out(1,1,1), mpi_comm, sign, pfft_flags, flags)

    POP_SUB(pfft_prepare_plan_c2c)
  end subroutine pfft_prepare_plan_c2c

  ! ---------------------------------------------------------
  subroutine pfft_get_dims(rs_n_global, mpi_comm, is_real, alloc_size, fs_n_global, rs_n, fs_n, rs_istart, fs_istart)
    integer, intent(in)  :: rs_n_global(1:3) !< The general number of elements in each dimension in real space
    integer, intent(in)  :: mpi_comm         !< MPI comunicator
    logical, intent(in)  :: is_real          !< The input is real or complex
    integer, intent(out) :: alloc_size       !< Number of elements that has to be allocated locally
    integer, intent(out) :: fs_n_global(1:3) !< The general number of elements in each dimension in Fourier space
    integer, intent(out) :: rs_n(1:3)        !< Local number of elements in each direction in real space
    integer, intent(out) :: fs_n(1:3)        !< Local number of elements in each direction in Fourier space
    integer, intent(out) :: rs_istart(1:3)   !< Where does the local portion of the function start in real space
    integer, intent(out) :: fs_istart(1:3)   !< Where does the local portion of the function start in Fourier space

    integer(ptrdiff_t_kind) :: tmp_alloc_size, tmp_n(3)
    integer(ptrdiff_t_kind) :: tmp_rs_n(3), tmp_rs_istart(3)
    integer(ptrdiff_t_kind) :: tmp_fs_n(3), tmp_fs_istart(3)

    PUSH_SUB(pfft_get_dims)

    fs_n_global(1:3) = rs_n_global(1:3)
    tmp_n(1:3) = rs_n_global(1:3)

    if (is_real) then
      call PDFFT(local_size_dft_r2c_3d)(tmp_alloc_size, tmp_n(1), mpi_comm, PFFT_TRANSPOSED_OUT, &
           tmp_rs_n(1), tmp_rs_istart(1), tmp_fs_n(1), tmp_fs_istart(1)) 
      fs_n_global(1) = rs_n_global(1)/2 + 1
    else
      call PDFFT(local_size_dft_3d)(tmp_alloc_size, tmp_n(1), mpi_comm, PFFT_TRANSPOSED_OUT, &
                 tmp_rs_n(1), tmp_rs_istart(1), tmp_fs_n(1), tmp_fs_istart(1))
    end if

    alloc_size       = tmp_alloc_size
    rs_n(1:3)        = tmp_rs_n(1:3)
    fs_n(1:3)        = tmp_fs_n(1:3)
    rs_istart(1:3)   = tmp_rs_istart(1:3)
    fs_istart(1:3)   = tmp_fs_istart(1:3)

    POP_SUB(pfft_get_dims)
  end subroutine pfft_get_dims

#endif

end module pfft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

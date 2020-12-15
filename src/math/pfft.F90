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

!> The includes for the PFFT
module pfft_params_oct_m
  use, intrinsic :: iso_c_binding
  use fftw_params_oct_m
  implicit none

  private

#ifdef HAVE_PFFT
  
  ! make public some symbols from the pfft library
  public :: pfft_cleanup,               &
            pfft_create_procmesh_2d,    &
            pfft_destroy_plan,          &
            pfft_execute,               &
            pfft_init,                  &
            pfft_local_size_dft_r2c_3d, &
            pfft_local_size_dft_3d,     &
            pfft_plan_dft_c2r_3d,       &
            pfft_plan_dft_r2c_3d,       &
            pfft_plan_dft_3d,           &
            PFFT_INT,                   &
            PFFT_PTRDIFF_T,             &
            PFFT_FLOAT,                 &
            PFFT_DOUBLE,                &
            PFFT_LDOUBLE,               &
            PFFT_UNSIGNED,              &
            PFFT_TRANSPOSED_IN,         &
            PFFT_TRANSPOSED_OUT,        &
            PFFT_ESTIMATE,              &
            PFFT_TUNE,                  &
            PFFT_DESTROY_INPUT,         &
            PFFT_MEASURE,               &
            PFFT_FORWARD,               &
            PFFT_BACKWARD

  include "pfft.f03"
#endif

end module pfft_params_oct_m
 
!> The low level module to work with the PFFT library.
!! http://www-user.tu-chemnitz.de/~mpip/software.php?lang=en
module pfft_oct_m
  use global_oct_m
  use, intrinsic :: iso_c_binding
  use math_oct_m
  use messages_oct_m
  use pfft_params_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::                    &
    pfft_decompose,            &
    pfft_prepare_plan_r2c,     &
    pfft_prepare_plan_c2r,     &
    pfft_prepare_plan_c2c,     &
    pfft_get_dims

contains

  ! ---------------------------------------------------------
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

  ! ---------------------------------------------------------
  !> Octopus subroutine to prepare a PFFT plan real to complex
  subroutine pfft_prepare_plan_r2c(plan, n, in, out, sign, flags, mpi_comm)
    type(C_PTR),      intent(out)   :: plan       !< The plan that is created by PFFT
    integer,          intent(in)    :: n(:)       !< The size of the global matrix
    FLOAT,   pointer, intent(inout) :: in(:,:,:)  !< The input matrix that is going to be used to do the transform
    CMPLX,   pointer, intent(inout) :: out(:,:,:) !< The output matrix that is going to be used to do the transform
    integer,          intent(in)    :: sign       !< Sign flag to decide FFT direction. Has to be FFTW_FORWARD
    integer,          intent(in)    :: flags      !< Flags for FFT library. Could be changed with the input variable
                                                  !! FFTPreparePlan. Default value is FFTW_MEASURE
    integer,          intent(in)    :: mpi_comm   !< MPI communicator
    
    integer(C_INTPTR_T) :: tmp_n(3)
    type(profile_t), save   :: prof

    PUSH_SUB(pfft_prepare_plan_r2c)
    call profiling_in(prof,"PFFT_PLAN_R2C")

#ifdef HAVE_PFFT
    ASSERT(sign == PFFT_FORWARD)
#endif

    tmp_n(1:3) = n(3:1:-1)

#ifdef HAVE_PFFT
    plan = pfft_plan_dft_r2c_3d(tmp_n,in,out,mpi_comm,sign,&
      PFFT_TRANSPOSED_OUT + PFFT_MEASURE + PFFT_DESTROY_INPUT + PFFT_TUNE) 
#else
    plan = C_NULL_PTR
#endif

    call profiling_out(prof)
    POP_SUB(pfft_prepare_plan_r2c)
  end subroutine pfft_prepare_plan_r2c

  ! ---------------------------------------------------------
  !> Octopus subroutine to prepare a PFFT plan real to complex
  subroutine pfft_prepare_plan_c2r(plan, n, in, out, sign, flags, mpi_comm)
    type(C_PTR),      intent(out)   :: plan       !< The plan that is created by PFFT
    integer,          intent(in)    :: n(:)       !< The size of the global matrix
    CMPLX,   pointer, intent(inout) :: in(:,:,:)  !< The input matrix that is going to be used to do the transform
    FLOAT,   pointer, intent(inout) :: out(:,:,:) !< The output matrix that is going to be used to do the transform
    integer,          intent(in)    :: sign       !< Sign flag to decide FFT direction. Has to be FFTW_BACKWARD
    integer,          intent(in)    :: flags      !< Flags for FFT library. Could be changed with the input variable
                                                  !! FFTPreparePlan. Default value is FFTW_MEASURE
    integer,          intent(in)    :: mpi_comm   !< MPI communicator
    
    integer(C_INTPTR_T) :: tmp_n(3)
    type(profile_t), save   :: prof
    PUSH_SUB(pfft_prepare_plan_c2r)
    call profiling_in(prof,"PFFT_PLAN_C2R")

#ifdef HAVE_PFFT
    ASSERT(sign == PFFT_BACKWARD)
#endif

    tmp_n(1:3) = n(3:1:-1)

#ifdef HAVE_PFFT
    plan = pfft_plan_dft_c2r_3d(tmp_n,in,out,mpi_comm, sign, &
      PFFT_TRANSPOSED_IN + PFFT_MEASURE + PFFT_DESTROY_INPUT + PFFT_TUNE) 
    !call PDFFT(plan_dft_c2r_3d) (plan, tmp_n(1), in(1,1,1), out(1,1,1), mpi_comm, sign, &
    !     PFFT_TRANSPOSED_IN + PFFT_MEASURE + PFFT_DESTROY_INPUT + PFFT_TUNE) 
#else
    plan = C_NULL_PTR
#endif

    call profiling_out(prof)
    POP_SUB(pfft_prepare_plan_c2r)
  end subroutine pfft_prepare_plan_c2r

  ! ---------------------------------------------------------  
  !> Octopus subroutine to prepare a PFFT plan real to complex
  subroutine pfft_prepare_plan_c2c(plan, n, in, out, sign, flags, mpi_comm)
    type(C_PTR),      intent(out)   :: plan       !< The plan that is created by PFFT
    integer,          intent(in)    :: n(:)       !< The size of the global matrix
    CMPLX,   pointer, intent(inout) :: in(:,:,:)  !< The input matrix that is going to be used to do the transform
    CMPLX,   pointer, intent(inout) :: out(:,:,:) !< The output matrix that is going to be used to do the transform
    integer,          intent(in)    :: sign       !< Sign flag to decide FFT direction. 
                                                         !! Has to be FFTW_FORWARD or FFTW_BACKWARD
    integer,          intent(in)    :: flags      !< Flags for FFT library. Could be changed with the input variable
                                                  !! FFTPreparePlan. Default value is FFTW_MEASURE
    integer,          intent(in)    :: mpi_comm   !< MPI communicator

#ifdef HAVE_PFFT
    integer(C_INT) :: pfft_flags
#endif
    integer(C_INTPTR_T) :: tmp_n(3)

    PUSH_SUB(pfft_prepare_plan_c2c)

    tmp_n(1:3) = n(3:1:-1)

#ifdef HAVE_PFFT
    if (sign == PFFT_FORWARD) then
      pfft_flags = PFFT_TRANSPOSED_OUT
    else
      pfft_flags = PFFT_TRANSPOSED_IN
    end if

    plan = pfft_plan_dft_3d(tmp_n,in,out,mpi_comm, sign, pfft_flags)
    !call PDFFT(plan_dft_3d) (plan, tmp_n(1), in(1,1,1), out(1,1,1), mpi_comm, sign, pfft_flags)
#else
    plan = C_NULL_PTR
#endif

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

    integer(C_INTPTR_T) :: tmp_alloc_size, tmp_n(3)
    integer(C_INTPTR_T) :: tmp_rs_n(3), tmp_rs_istart(3)
    integer(C_INTPTR_T) :: tmp_fs_n(3), tmp_fs_istart(3)

    PUSH_SUB(pfft_get_dims)

    fs_n_global(1:3) = rs_n_global(1:3)
    ! if calling the ISO_C_BINDING routines from PFFT, one must take into account that
    ! the array order is the C array order
    tmp_n(1:3) = rs_n_global(3:1:-1)

    tmp_alloc_size = 0
    if (is_real) then
#ifdef HAVE_PFFT
       tmp_alloc_size = pfft_local_size_dft_r2c_3d(tmp_n, mpi_comm, PFFT_TRANSPOSED_OUT, &
         tmp_rs_n, tmp_rs_istart, tmp_fs_n, tmp_fs_istart) 
#endif
       !call PDFFT(local_size_dft_r2c_3d)(tmp_alloc_size, tmp_n(1), mpi_comm, PFFT_TRANSPOSED_OUT, &
       !       tmp_rs_n(1), tmp_rs_istart(1), tmp_fs_n(1), tmp_fs_istart(1)) 
       fs_n_global(1) = rs_n_global(1)/2 + 1
     else
#ifdef HAVE_PFFT
       tmp_alloc_size = pfft_local_size_dft_3d(tmp_n,mpi_comm, PFFT_TRANSPOSED_OUT, &
         tmp_rs_n,tmp_rs_istart,tmp_fs_n,tmp_fs_istart)
#endif
       !call PDFFT(local_size_dft_3d)(tmp_alloc_size, tmp_n(1), mpi_comm, PFFT_TRANSPOSED_OUT, &
       !             tmp_rs_n(1), tmp_rs_istart(1), tmp_fs_n(1), tmp_fs_istart(1))
    end if

    alloc_size       = int(tmp_alloc_size)
    rs_n(1:3)        = int(tmp_rs_n(3:1:-1))
    fs_n(1:3)        = int(tmp_fs_n(3:1:-1))
    rs_istart(1:3)   = int(tmp_rs_istart(3:1:-1)+1)
    fs_istart(1:3)   = int(tmp_fs_istart(3:1:-1)+1)

    POP_SUB(pfft_get_dims)
  end subroutine pfft_get_dims

end module pfft_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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

module pfft_params_m
  use fftw_m
  implicit none

#ifdef HAVE_PFFT
  include "pfft.f"
#else
  integer, parameter :: ptrdiff_t_kind = 8
#endif
end module pfft_params_m

module pfft_m
  use global_m
  use math_m
  use pfft_params_m
  implicit none

#ifdef HAVE_PFFT

  ! ----------------- init ------------------
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
    subroutine PDFFT(local_size_dft_r2c_3d)(alloc_size, n, mpi_comm, pfft_flags, &
         local_ni, local_i_start, local_no, local_o_start)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(inout) :: alloc_size
      integer(ptrdiff_t_kind), intent(inout) :: n
      integer,                 intent(in)    :: mpi_comm, pfft_flags
      integer(ptrdiff_t_kind), intent(inout) :: local_ni, local_i_start, local_no, local_o_start
    end subroutine PDFFT(local_size_dft_r2c_3d)

    subroutine PDFFT(local_size_dft_3d)(alloc_size, n, mpi_comm, pfft_flags, &
         local_ni, local_i_start, local_no, local_o_start)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(inout) :: alloc_size
      integer(ptrdiff_t_kind), intent(inout) :: n
      integer,                 intent(in)    :: mpi_comm, pfft_flags
      integer(ptrdiff_t_kind), intent(inout) :: local_ni, local_i_start, local_no, local_o_start
    end subroutine PDFFT(local_size_dft_3d)

    subroutine PDFFT(plan_dft_r2c_3d)(plan, n, in, out, mpi_comm, sign, pfft_flags, fftw_flags)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(out)   :: plan
      integer(ptrdiff_t_kind), intent(inout) :: n
      FLOAT,                   intent(in)    :: in
      CMPLX,                   intent(in)    :: out
      integer,                 intent(in)    :: mpi_comm, sign, pfft_flags, fftw_flags
    end subroutine PDFFT(plan_dft_r2c_3d)

    subroutine PDFFT(plan_dft_c2r_3d)(plan, n, in, out, mpi_comm, sign, pfft_flags, fftw_flags)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(in)    :: plan
      integer(ptrdiff_t_kind), intent(inout) :: n
      CMPLX,                   intent(in)    :: in
      FLOAT,                   intent(in)    :: out
      integer,                 intent(in)    :: mpi_comm, sign, pfft_flags, fftw_flags
    end subroutine PDFFT(plan_dft_c2r_3d)

    subroutine PDFFT(plan_dft_3d)(plan, n, in, out, mpi_comm, sign, pfft_flags, fftw_flags)
      use pfft_params_m
      integer(ptrdiff_t_kind), intent(out)   :: plan
      integer(ptrdiff_t_kind), intent(inout) :: n
      CMPLX,                   intent(in)    :: in
      CMPLX,                   intent(in)    :: out
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

  ! ---------------------------------------------------------
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
  subroutine pfft_prepare_plan_r2c(plan, n, in, out, sign, flags, mpi_comm)
    integer(ptrdiff_t_kind), intent(out)   :: plan
    integer,                 intent(in)    :: n(:)
    FLOAT,   pointer,        intent(inout) :: in(:,:,:)
    CMPLX,   pointer,        intent(inout) :: out(:,:,:)
    integer,                 intent(in)    :: sign
    integer,                 intent(in)    :: flags
    integer,                 intent(in)    :: mpi_comm
    
    integer(ptrdiff_t_kind) :: tmp_n(3)

    PUSH_SUB(pfft_prepare_plan_r2c)

    ASSERT(sign == FFTW_FORWARD)

    tmp_n(1:3) = n(1:3)

    call PDFFT(plan_dft_r2c_3d) (plan, tmp_n(1), in(1,1,1), out(1,1,1), mpi_comm, sign, PFFT_TRANSPOSED_OUT, flags)

    POP_SUB(pfft_prepare_plan_r2c)
  end subroutine pfft_prepare_plan_r2c

  ! ---------------------------------------------------------
  subroutine pfft_prepare_plan_c2r(plan, n, in, out, sign, flags, mpi_comm)
    integer(ptrdiff_t_kind), intent(out)   :: plan
    integer,                 intent(in)    :: n(:)
    CMPLX,   pointer,        intent(inout) :: in(:,:,:)
    FLOAT,   pointer,        intent(inout) :: out(:,:,:)
    integer,                 intent(in)    :: sign
    integer,                 intent(in)    :: flags
    integer,                 intent(in)    :: mpi_comm
    
    integer(ptrdiff_t_kind) :: tmp_n(3)

    PUSH_SUB(pfft_prepare_plan_c2r)

    ASSERT(sign == FFTW_BACKWARD)

    tmp_n(1:3) = n(1:3)

    call PDFFT(plan_dft_c2r_3d) (plan, tmp_n(1), in(1,1,1), out(1,1,1), mpi_comm, sign, PFFT_TRANSPOSED_IN, flags)

    POP_SUB(pfft_prepare_plan_c2r)
  end subroutine pfft_prepare_plan_c2r

  ! ---------------------------------------------------------
  subroutine pfft_prepare_plan_c2c(plan, n, in, out, sign, flags, mpi_comm)
    integer(ptrdiff_t_kind), intent(out)   :: plan
    integer,                 intent(in)    :: n(:)
    CMPLX,   pointer,        intent(inout) :: in(:,:,:)
    CMPLX,   pointer,        intent(inout) :: out(:,:,:)
    integer,                 intent(in)    :: sign
    integer,                 intent(in)    :: flags
    integer,                 intent(in)    :: mpi_comm    

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
    integer, intent(in)  :: rs_n_global(1:3)
    integer, intent(in)  :: mpi_comm
    logical, intent(in)  :: is_real
    integer, intent(out) :: alloc_size
    integer, intent(out) :: fs_n_global(1:3)
    integer, intent(out) :: rs_n(1:3)
    integer, intent(out) :: fs_n(1:3)
    integer, intent(out) :: rs_istart(1:3)
    integer, intent(out) :: fs_istart(1:3)

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

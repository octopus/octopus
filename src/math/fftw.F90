!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2011 M. Oliveira
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
#  define DFFTW(x) sfftw_ ## x
#else
#  define DFFTW(x) dfftw_ ## x
#endif

module fftw_m
  use c_pointer_m
  use global_m
  use messages_m
  use profiling_m
  implicit none

  ! fftw constants. this is just a copy from file fftw3.f,
  ! distributed with fftw package.
  integer, parameter ::                &
    FFTW_R2HC                =      0, &
    FFTW_HC2R                =      1, &
    FFTW_DHT                 =      2, &
    FFTW_REDFT00             =      3, &
    FFTW_REDFT01             =      4, &
    FFTW_REDFT10             =      5, &
    FFTW_REDFT11             =      6, &
    FFTW_RODFT00             =      7, &
    FFTW_RODFT01             =      8, &
    FFTW_RODFT10             =      9, &
    FFTW_RODFT11             =     10, &
    FFTW_FORWARD             =     -1, &
    FFTW_BACKWARD            =      1, &
    FFTW_MEASURE             =      0, &
    FFTW_DESTROY_INPUT       =      1, &
    FFTW_UNALIGNED           =      2, &
    FFTW_CONSERVE_MEMORY     =      4, &
    FFTW_EXHAUSTIVE          =      8, &
    FFTW_PRESERVE_INPUT      =     16, &
    FFTW_PATIENT             =     32, &
    FFTW_ESTIMATE            =     64, &
    FFTW_ESTIMATE_PATIENT    =    128, &
    FFTW_BELIEVE_PCOST       =    256, &
    FFTW_DFT_R2HC_ICKY       =    512, &
    FFTW_NONTHREADED_ICKY    =   1024, &
    FFTW_NO_BUFFERING        =   2048, &
    FFTW_NO_INDIRECT_OP      =   4096, &
    FFTW_ALLWO_LARGE_GENERIC =   8192, &
    FFTW_NO_RANK_SPLITS      =  16384, &
    FFTW_NO_VRANK_SPLITS     =  32768, &
    FFTW_NO_VRECURSE         =  65536, &
    FFTW_NO_SIMD             = 131072

  ! ----------------- plan_dft ------------------
  interface 
    subroutine DFFTW(plan_dft_r2c_1d)(plan, n, in, out, flags)
      use c_pointer_m
      type(c_ptr), intent(inout) :: plan
      integer,     intent(in)    :: n
      FLOAT,       intent(in)    :: in  ! in(n)
      CMPLX,       intent(out)   :: out ! out(n/2+1)
      integer,     intent(in)    :: flags
    end subroutine DFFTW(plan_dft_r2c_1d)

    subroutine DFFTW(plan_dft_r2c_2d)(plan, n1, n2, in, out, flags)
      use c_pointer_m
      type(c_ptr), intent(inout) :: plan
      integer,     intent(in)    :: n1, n2
      FLOAT,       intent(in)    :: in  ! in(n1, n2)
      CMPLX,       intent(out)   :: out ! out(n1/2+1, n2)
      integer,     intent(in)    :: flags
    end subroutine DFFTW(plan_dft_r2c_2d)
    
    subroutine DFFTW(plan_dft_r2c_3d)(plan, n1, n2, n3, in, out, flags)
      use c_pointer_m
      type(c_ptr), intent(inout) :: plan
      integer,     intent(in)    :: n1, n2, n3
      FLOAT,       intent(in)    :: in  ! in(n1, n2, n3)
      CMPLX,       intent(out)   :: out ! out(n1/2+1, n2, n3)
      integer,     intent(in)    :: flags
    end subroutine DFFTW(plan_dft_r2c_3d)

    subroutine DFFTW(plan_dft_c2r_1d)(plan, n, in, out, flags)
      use c_pointer_m
      type(c_ptr), intent(inout) :: plan
      integer,     intent(in)    :: n
      CMPLX,       intent(in)    :: in  ! in(n/2+1)
      FLOAT,       intent(out)   :: out ! out(n)
      integer,     intent(in)    :: flags
    end subroutine DFFTW(plan_dft_c2r_1d)

    subroutine DFFTW(plan_dft_c2r_2d)(plan, n1, n2, in, out, flags)
      use c_pointer_m
      type(c_ptr), intent(inout) :: plan
      integer,     intent(in)    :: n1, n2
      CMPLX,       intent(in)    :: in  ! in(n1/2+1, n2)
      FLOAT,       intent(out)   :: out ! out(n1, n2)
      integer,     intent(in)    :: flags
    end subroutine DFFTW(plan_dft_c2r_2d)

    subroutine DFFTW(plan_dft_c2r_3d)(plan, n1, n2, n3, in, out, flags)
      use c_pointer_m
      type(c_ptr), intent(inout) :: plan
      integer,     intent(in)    :: n1, n2, n3
      CMPLX,       intent(in)    :: in  ! in(n1/2+1, n2, n3)
      FLOAT,       intent(out)   :: out ! out(n1, n2, n3)
      integer,     intent(in)    :: flags
    end subroutine DFFTW(plan_dft_c2r_3d)

    subroutine DFFTW(plan_dft_1d)(plan, n, in, out, sign, flags)
      use c_pointer_m
      type(c_ptr), intent(inout) :: plan
      integer,     intent(in)    :: n
      CMPLX,       intent(in)    :: in  ! in(n/2+1)
      CMPLX,       intent(out)   :: out ! out(n)
      integer,     intent(in)    :: sign, flags
    end subroutine DFFTW(plan_dft_1d)

    subroutine DFFTW(plan_dft_2d)(plan, n1, n2, in, out, sign, flags)
      use c_pointer_m
      type(c_ptr), intent(inout) :: plan
      integer,     intent(in)    :: n1, n2
      CMPLX,       intent(in)    :: in  ! in(n1/2+1, n2)
      CMPLX,       intent(out)   :: out ! out(n1, n2)
      integer,     intent(in)    :: sign, flags
    end subroutine DFFTW(plan_dft_2d)

    subroutine DFFTW(plan_dft_3d)(plan, n1, n2, n3, in, out, sign, flags)
      use c_pointer_m
      type(c_ptr), intent(inout) :: plan
      integer,     intent(in)    :: n1, n2, n3
      CMPLX,       intent(in)    :: in  ! in(n1/2+1, n2, n3)
      CMPLX,       intent(out)   :: out ! out(n1, n2, n3)
      integer,     intent(in)    :: sign, flags
    end subroutine DFFTW(plan_dft_3d)
  end interface


  ! ----------------- execute_dft ------------------
  interface fftw_execute_dft
    subroutine DFFTW(execute_dft_r2c)(plan, in, out)
      use c_pointer_m
      type(c_ptr), intent(in)  :: plan
      FLOAT,       intent(in)  :: in
      CMPLX,       intent(out) :: out
    end subroutine DFFTW(execute_dft_r2c)

    subroutine DFFTW(execute_dft_c2r)(plan, in, out)
      use c_pointer_m
      type(c_ptr), intent(in)  :: plan
      CMPLX,       intent(in)  :: in
      FLOAT,       intent(out) :: out
    end subroutine DFFTW(execute_dft_c2r)

    subroutine DFFTW(execute_dft)(plan, in, out)
      use c_pointer_m
      type(c_ptr), intent(in)  :: plan
      CMPLX,       intent(in)  :: in
      CMPLX,       intent(out) :: out
    end subroutine DFFTW(execute_dft)
  end interface fftw_execute_dft


  ! ----------------- destroy_plan ------------------
  interface fftw_destroy_plan
    subroutine DFFTW(destroy_plan)(plan)
      use c_pointer_m
      type(c_ptr), intent(inout) :: plan
    end subroutine DFFTW(destroy_plan)
  end interface fftw_destroy_plan

 
  ! ----------------- cleanup ------------------
  interface fftw_cleanup
    subroutine DFFTW(cleanup)
    end subroutine DFFTW(cleanup)
  end interface fftw_cleanup


  ! ----------------- thread related functions---------------

  interface fftw_init_threads
    subroutine DFFTW(init_threads)(iret)
      integer, intent(out) :: iret
    end subroutine DFFTW(init_threads)
  end interface fftw_init_threads

  interface fftw_plan_with_nthreads
    subroutine DFFTW(plan_with_nthreads)(nthreads)
      integer, intent(in) :: nthreads
    end subroutine DFFTW(plan_with_nthreads)
  end interface fftw_plan_with_nthreads

  interface fftw_cleanup_threads
    subroutine DFFTW(cleanup_threads)
    end subroutine DFFTW(cleanup_threads)
  end interface fftw_cleanup_threads

contains

  ! ---------------------------------------------------------
  subroutine fftw_prepare_plan_r2c(plan, dim, n, sign, flags)
    type(c_ptr), intent(inout) :: plan
    integer,     intent(in)   :: dim
    integer,     intent(in)   :: n(:)
    integer,     intent(in)   :: sign
    integer,     intent(in)   :: flags

    FLOAT, allocatable :: in(:,:,:)
    CMPLX, allocatable :: out(:,:,:)
    
    PUSH_SUB(fftw_prepare_plan_r2c)

    ASSERT(sign == FFTW_FORWARD)

    SAFE_ALLOCATE(in(1:n(1), 1:n(2), 1:n(3)))
    SAFE_ALLOCATE(out(1:n(1)/2+1, 1:n(2), 1:n(3)))

    select case (dim)
    case (1)
      call DFFTW(plan_dft_r2c_1d)(plan, n(1), in(1,1,1), out(1,1,1), flags)
    case (2)
      call DFFTW(plan_dft_r2c_2d)(plan, n(1), n(2), in(1,1,1), out(1,1,1), flags)
    case (3)
      call DFFTW(plan_dft_r2c_3d)(plan, n(1), n(2), n(3), in(1,1,1), out(1,1,1), flags)
    end select

    SAFE_DEALLOCATE_A(in)
    SAFE_DEALLOCATE_A(out)

    POP_SUB(fftw_prepare_plan_r2c)
  end subroutine fftw_prepare_plan_r2c

  ! ---------------------------------------------------------
  subroutine fftw_prepare_plan(plan, dim, n, is_real, sign, flags)
    type(c_ptr), intent(inout) :: plan
    integer,     intent(in)    :: dim
    integer,     intent(in)    :: n(:)
    logical,     intent(in)    :: is_real
    integer,     intent(in)    :: sign
    integer,     intent(in)    :: flags

    FLOAT, allocatable :: rin(:,:,:)
    FLOAT, allocatable :: rout(:,:,:)
    CMPLX, allocatable :: cin(:,:,:)
    CMPLX, allocatable :: cout(:,:,:)
    
    PUSH_SUB(fftw_prepare_plan)

    ASSERT(sign == FFTW_FORWARD .or. sign == FFTW_BACKWARD)

    if (is_real) then
      if (sign == FFTW_FORWARD) then
        SAFE_ALLOCATE(rin(1:n(1), 1:n(2), 1:n(3)))
        SAFE_ALLOCATE(cout(1:n(1)/2+1, 1:n(2), 1:n(3)))

        select case (dim)
        case (1)
          call DFFTW(plan_dft_r2c_1d)(plan, n(1), rin(1,1,1), cout(1,1,1), flags)
        case (2)
          call DFFTW(plan_dft_r2c_2d)(plan, n(1), n(2), rin(1,1,1), cout(1,1,1), flags)
        case (3)
          call DFFTW(plan_dft_r2c_3d)(plan, n(1), n(2), n(3), rin(1,1,1), cout(1,1,1), flags)
        end select

        SAFE_DEALLOCATE_A(rin)
        SAFE_DEALLOCATE_A(cout)
      else
        SAFE_ALLOCATE(cin(1:n(1)/2+1, 1:n(2), 1:n(3)))
        SAFE_ALLOCATE(rout(1:n(1), 1:n(2), 1:n(3)))

        select case (dim)
        case (1)
          call DFFTW(plan_dft_c2r_1d)(plan, n(1), cin(1,1,1), rout(1,1,1), flags)
        case (2)
          call DFFTW(plan_dft_c2r_2d)(plan, n(1), n(2), cin(1,1,1), rout(1,1,1), flags)
        case (3)
          call DFFTW(plan_dft_c2r_3d)(plan, n(1), n(2), n(3), cin(1,1,1), rout(1,1,1), flags)
        end select

        SAFE_DEALLOCATE_A(cin)
        SAFE_DEALLOCATE_A(rout)
      end if
    else
      SAFE_ALLOCATE(cin(1:n(1), 1:n(2), 1:n(3)))
      SAFE_ALLOCATE(cout(1:n(1), 1:n(2), 1:n(3)))

      select case (dim)
      case (1)
        call DFFTW(plan_dft_1d)(plan, n(1), cin(1,1,1), cout(1,1,1), sign, flags)
      case (2)
        call DFFTW(plan_dft_2d)(plan, n(1), n(2), cin(1,1,1), cout(1,1,1), sign, flags)
      case (3)
        call DFFTW(plan_dft_3d)(plan, n(1), n(2), n(3), cin(1,1,1), cout(1,1,1), sign, flags)
      end select

      SAFE_DEALLOCATE_A(cin)
      SAFE_DEALLOCATE_A(cout)
    end if

    POP_SUB(fftw_prepare_plan)
  end subroutine fftw_prepare_plan

  ! ---------------------------------------------------------
  subroutine fftw_get_dims(rs_n, is_real, fs_n)
    integer, intent(in)  :: rs_n(:)
    logical, intent(in)  :: is_real
    integer, intent(out) :: fs_n(:)

    PUSH_SUB(fftw_get_dims)

    fs_n = rs_n
    if (is_real) fs_n(1) = rs_n(1)/2 + 1

    POP_SUB(fftw_get_dims)
  end subroutine fftw_get_dims

end module fftw_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

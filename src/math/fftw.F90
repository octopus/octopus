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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module fftw_params_oct_m
  use, intrinsic :: iso_c_binding
  implicit none

  private

  ! make public some needed symbols from the FFTW library
  public :: fftw_cleanup,         &
            fftw_execute_dft,     &
            fftw_execute_dft_c2r, &
            fftw_execute_dft_r2c, &
            fftw_destroy_plan,    &
            fftw_plan_dft_1d,     &
            fftw_plan_dft_2d,     &
            fftw_plan_dft_3d,     &
            fftw_plan_dft_r2c_1d, &
            fftw_plan_dft_r2c_2d, &
            fftw_plan_dft_r2c_3d, &
            fftw_plan_dft_c2r_1d, &
            fftw_plan_dft_c2r_2d, &
            fftw_plan_dft_c2r_3d, &
            C_FFTW_R2R_KIND,      &
            FFTW_R2HC,            &
            FFTW_HC2R,            &
            FFTW_DHT,             &
            FFTW_REDFT00,         &
            FFTW_REDFT01,         &
            FFTW_REDFT10,         &
            FFTW_REDFT11,         &
            FFTW_RODFT00,         &
            FFTW_RODFT01,         &
            FFTW_RODFT10,         &
            FFTW_RODFT11,         &
            FFTW_FORWARD,         &
            FFTW_BACKWARD,        &
            FFTW_MEASURE,         &
            FFTW_DESTROY_INPUT,   &
            FFTW_UNALIGNED
#ifdef HAVE_FFTW3_MPI
  public :: FFTW_MPI_DEFAULT_BLOCK
#endif
#if defined(HAVE_OPENMP) && defined(HAVE_FFTW3_THREADS)
  public :: fftw_init_threads,   &
            fftw_plan_with_nthreads, &
            fftw_cleanup_threads
#endif

#ifdef HAVE_FFTW3_MPI
  include "fftw3-mpi.f03"
#else
  include "fftw3.f03"
#endif
end module fftw_params_oct_m

module fftw_oct_m
  use fftw_params_oct_m
  use global_oct_m
  use, intrinsic :: iso_c_binding
  use messages_oct_m
  use profiling_oct_m
  implicit none

  private

  public :: &
    fftw_prepare_plan,       &
    fftw_get_dims

contains

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
           plan = fftw_plan_dft_r2c_1d(n(1), rin, cout, flags)
           !call DFFTW(plan_dft_r2c_1d)(plan, n(1), rin(1,1,1), cout(1,1,1), flags)
        case (2)
           plan = fftw_plan_dft_r2c_2d(n(2), n(1),rin, cout, flags)
           !call DFFTW(plan_dft_r2c_2d)(plan, n(1), n(2), rin(1,1,1), cout(1,1,1), flags)
        case (3)
           plan = fftw_plan_dft_r2c_3d(n(3), n(2), n(1), rin, cout, flags)
           !call DFFTW(plan_dft_r2c_3d)(plan, n(1), n(2), n(3), rin(1,1,1), cout(1,1,1), flags)
        end select

        SAFE_DEALLOCATE_A(rin)
        SAFE_DEALLOCATE_A(cout)
      else
        SAFE_ALLOCATE(cin(1:n(1)/2+1, 1:n(2), 1:n(3)))
        SAFE_ALLOCATE(rout(1:n(1), 1:n(2), 1:n(3)))

        select case (dim)
        case (1)
           plan = fftw_plan_dft_c2r_1d(n(1), cin, rout, flags)
           !call DFFTW(plan_dft_c2r_1d)(plan, n(1), cin(1,1,1), rout(1,1,1), flags)
        case (2)
           plan = fftw_plan_dft_c2r_2d(n(2), n(1), cin, rout, flags)
           !call DFFTW(plan_dft_c2r_2d)(plan, n(1), n(2), cin(1,1,1), rout(1,1,1), flags)
        case (3)
           plan = fftw_plan_dft_c2r_3d(n(3), n(2), n(1), cin, rout, flags)
           !call DFFTW(plan_dft_c2r_3d)(plan, n(1), n(2), n(3), cin(1,1,1), rout(1,1,1), flags)
        end select

        SAFE_DEALLOCATE_A(cin)
        SAFE_DEALLOCATE_A(rout)
      end if
    else
      SAFE_ALLOCATE(cin(1:n(1), 1:n(2), 1:n(3)))
      SAFE_ALLOCATE(cout(1:n(1), 1:n(2), 1:n(3)))

      select case (dim)
      case (1)
         plan = fftw_plan_dft_1d(n(1), cin, cout, sign, flags)
         !call DFFTW(plan_dft_1d)(plan, n(1), cin(1,1,1), cout(1,1,1), sign, flags)
      case (2)
         plan = fftw_plan_dft_2d(n(2), n(1), cin, cout, sign, flags)
         !call DFFTW(plan_dft_2d)(plan, n(1), n(2), cin(1,1,1), cout(1,1,1), sign, flags)
      case (3)
         plan = fftw_plan_dft_3d(n(3), n(2), n(1), cin, cout, sign, flags)
         !call DFFTW(plan_dft_3d)(plan, n(1), n(2), n(3), cin(1,1,1), cout(1,1,1), sign, flags)
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

end module fftw_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

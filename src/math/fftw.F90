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
            fftw_alloc_real,      &
            fftw_alloc_complex,   &
            fftw_free,            &
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
    fftw_get_dims,           &
    fftw_alloc_memory,       &
    fftw_free_memory


contains

  ! ---------------------------------------------------------
  subroutine fftw_prepare_plan(plan, dim, n, is_real, sign, flags, din_, cin_, cout_)
    type(c_ptr), intent(inout) :: plan
    integer,     intent(in)    :: dim
    integer,     intent(in)    :: n(:)
    logical,     intent(in)    :: is_real
    integer,     intent(in)    :: sign
    integer,     intent(in)    :: flags
    FLOAT, optional, target, intent(in) :: din_(:,:,:)
    CMPLX, optional, target, intent(in) :: cin_(:,:,:)
    CMPLX, optional, target, intent(in) :: cout_(:,:,:)

    FLOAT, pointer :: rin(:,:,:)
    FLOAT, pointer :: rout(:,:,:)
    CMPLX, pointer :: cin(:,:,:)
    CMPLX, pointer :: cout(:,:,:)
    logical :: aligned_memory
    
    PUSH_SUB(fftw_prepare_plan)

    ASSERT(sign == FFTW_FORWARD .or. sign == FFTW_BACKWARD)

    aligned_memory = .false.
    if(present(din_) .or. present(cin_)) then
      ASSERT(present(cout_))
      ASSERT(present(din_) .neqv. present(cin_))
      aligned_memory = .true.
    end if

    if (is_real) then
      if (sign == FFTW_FORWARD) then
        if(.not.aligned_memory) then
          SAFE_ALLOCATE(rin(1:n(1), 1:n(2), 1:n(3)))
          SAFE_ALLOCATE(cout(1:n(1)/2+1, 1:n(2), 1:n(3)))
        else
          rin => din_
          cout => cout_
        end if

        select case (dim)
        case (1)
           plan = fftw_plan_dft_r2c_1d(n(1), rin, cout, flags)
        case (2)
           plan = fftw_plan_dft_r2c_2d(n(2), n(1),rin, cout, flags)
        case (3)
           plan = fftw_plan_dft_r2c_3d(n(3), n(2), n(1), rin, cout, flags)
        end select

        if(.not.aligned_memory) then
          SAFE_DEALLOCATE_P(rin)
          SAFE_DEALLOCATE_P(cout)
        else
          nullify(rin, cout)
        end if
      else
        if(.not.aligned_memory) then
          SAFE_ALLOCATE(cin(1:n(1)/2+1, 1:n(2), 1:n(3)))
          SAFE_ALLOCATE(rout(1:n(1), 1:n(2), 1:n(3)))
        else
          cin => cout_
          rout => din_
        end if

        select case (dim)
        case (1)
           plan = fftw_plan_dft_c2r_1d(n(1), cin, rout, flags)
        case (2)
           plan = fftw_plan_dft_c2r_2d(n(2), n(1), cin, rout, flags)
        case (3)
           plan = fftw_plan_dft_c2r_3d(n(3), n(2), n(1), cin, rout, flags)
        end select

        if(.not.aligned_memory) then
          SAFE_DEALLOCATE_P(cin)
          SAFE_DEALLOCATE_P(rout)
        else
          nullify(cin,rout)
        end if
      end if
    else
      if(.not.aligned_memory) then
        SAFE_ALLOCATE(cin(1:n(1), 1:n(2), 1:n(3)))
        SAFE_ALLOCATE(cout(1:n(1), 1:n(2), 1:n(3)))
      else
        if (sign == FFTW_FORWARD) then
          cin => cin_
          cout => cout_
        else
          cout => cin_
          cin => cout_
        end if
      end if

      select case (dim)
      case (1)
         plan = fftw_plan_dft_1d(n(1), cin, cout, sign, flags)
      case (2)
         plan = fftw_plan_dft_2d(n(2), n(1), cin, cout, sign, flags)
      case (3)
         plan = fftw_plan_dft_3d(n(3), n(2), n(1), cin, cout, sign, flags)
      end select

      if(.not.aligned_memory) then
        SAFE_DEALLOCATE_P(cin)
        SAFE_DEALLOCATE_P(cout)
      else
        nullify(cin, cout)
      end if
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

  ! ---------------------------------------------------------
  subroutine fftw_alloc_memory(rs_n, is_real, fs_n, p_din, p_zin, p_cout)
    integer, intent(in)    :: rs_n(:)
    logical, intent(in)    :: is_real
    integer, intent(in)    :: fs_n(:)
    type(C_PTR),   intent(inout) :: p_din
    type(C_PTR),   intent(inout) :: p_zin
    type(C_PTR),   intent(inout) :: p_cout


    PUSH_SUB(fftw_alloc_memory)

    if(is_real) then
      p_din  = fftw_alloc_real(int(product(rs_n(1:3)), C_SIZE_T))
    else
      p_zin = fftw_alloc_complex(int(product(rs_n(1:3)), C_SIZE_T))
    end if 
    p_cout = fftw_alloc_complex(int(product(fs_n(1:3)), C_SIZE_T))

    POP_SUB(fftw_alloc_memory)
  end subroutine fftw_alloc_memory

  ! ---------------------------------------------------------
  subroutine fftw_free_memory(is_real, p_din, p_zin, p_cout)
    logical, intent(in)    :: is_real
    type(C_PTR),   intent(inout) :: p_din
    type(C_PTR),   intent(inout) :: p_zin
    type(C_PTR),   intent(inout) :: p_cout

    PUSH_SUB(fftw_free_memory)

    if(is_real) then
      call fftw_free(p_din)
    else
      call fftw_free(p_zin)
    end if
    call fftw_free(p_cout)

    POP_SUB(fftw_free_memory)
  end subroutine fftw_free_memory

end module fftw_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

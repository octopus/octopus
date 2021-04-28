!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module fourier_space_oct_m
  use accel_oct_m
  use cube_oct_m
  use cube_function_oct_m
  use global_oct_m
  use math_oct_m
  use messages_oct_m
  use fft_oct_m
#ifdef HAVE_OPENMP
  use omp_lib
#endif
#ifdef HAVE_PFFT
  use pfft_oct_m
#endif
  use profiling_oct_m
  use types_oct_m

  implicit none
  private
  public ::                     &
    cube_function_alloc_fs,     & 
    cube_function_free_fs,      &
    dcube_function_rs2fs,       &
    zcube_function_rs2fs,       &
    dcube_function_fs2rs,       &
    zcube_function_fs2rs,       &
    fourier_space_op_t,         &
    dfourier_space_op_init,     &
    dfourier_space_op_apply,    &
    zfourier_space_op_init,     &
    zfourier_space_op_apply,    &
    fourier_space_op_end

  type fourier_space_op_t
    private
    FLOAT, allocatable :: dop(:, :, :)
    CMPLX, allocatable :: zop(:, :, :)
    logical :: in_device_memory = .false.
    type(accel_mem_t) :: op_buffer
    logical :: real_op

    !Parameters used to generate this kernel
    FLOAT, public :: qq(1:MAX_DIM) = CNST(1e-5) !We just set a very large q to guaranty that the kernel is always
    FLOAT, public :: singularity = M_ZERO
    FLOAT, public :: mu = M_ZERO !< Range separation for the exchange operator
  end type fourier_space_op_t

contains

  ! ---------------------------------------------------------
  !> Allocates locally the Fourier space grid, if PFFT library is not used.
  !! Otherwise, it assigns the PFFT Fourier space grid to the cube Fourier space grid,
  !! via pointer.
  subroutine cube_function_alloc_fs(cube, cf, force_alloc)
    type(cube_t), target,  intent(in)    :: cube
    type(cube_function_t), intent(inout) :: cf
    logical, optional,     intent(in)    :: force_alloc  
      
    integer :: n1, n2, n3
    logical :: is_allocated
    
    PUSH_SUB(cube_function_alloc_fs)
    
    ASSERT(.not. associated(cf%fs))
    ASSERT(allocated(cube%fft))
    
    cf%forced_alloc = optional_default(force_alloc, .false.)

    n1 = max(1, cube%fs_n(1))
    n2 = max(1, cube%fs_n(2))
    n3 = max(1, cube%fs_n(3))

    is_allocated = .false.
    
    select case(cube%fft%library)
    case(FFTLIB_PFFT)
      if(.not. cf%forced_alloc) then 
        is_allocated = .true.
        if(any(cube%fs_n(1:3) == 0)) then
          cf%fs => cube%fft%fs_data(1:1,1:1,1:1)
        else
          cf%fs => cube%fft%fs_data(1:n3,1:n1,1:n2)
        end if
      else ! force allocate transposed with PFFT  
        is_allocated = .true.
        SAFE_ALLOCATE(cf%fs(1:n3, 1:n1, 1:n2))
      end if
    case(FFTLIB_ACCEL)
      if(cf%in_device_memory) then
        is_allocated = .true.
        call accel_create_buffer(cf%fourier_space_buffer, ACCEL_MEM_READ_WRITE, TYPE_CMPLX, product(cube%fs_n(1:3)))
      end if

    case(FFTLIB_FFTW)
      if(.not. cf%forced_alloc) then
        is_allocated = .true.
        cf%fs => cube%fft%fs_data(1:cube%fs_n(1), 1:cube%fs_n(2), 1:cube%fs_n(3))
      end if
    end select

    if(.not. is_allocated) then
      SAFE_ALLOCATE(cf%fs(1:cube%fs_n(1), 1:cube%fs_n(2), 1:cube%fs_n(3)))
    end if
    
    POP_SUB(cube_function_alloc_fs)
  end subroutine cube_function_alloc_fs

  
  ! ---------------------------------------------------------
  !> Deallocates the Fourier space grid
  subroutine cube_function_free_fs(cube, cf)
    type(cube_t),          intent(in)    :: cube
    type(cube_function_t), intent(inout) :: cf
    
    logical :: deallocated

    PUSH_SUB(cube_function_free_fs)

    ASSERT(allocated(cube%fft))

    deallocated = .false.

    select case(cube%fft%library)
    case(FFTLIB_PFFT)
      if(.not. cf%forced_alloc) then
        deallocated = .true.
        nullify(cf%fs)
      end if
    case(FFTLIB_ACCEL)
      if(cf%in_device_memory) then
        deallocated = .true.
        call accel_release_buffer(cf%fourier_space_buffer)
      end if
    case(FFTLIB_FFTW)
      if(.not. cf%forced_alloc) then
        deallocated = .true.
        nullify(cf%fs)
      end if
    end select

    if(.not. deallocated) then
      ASSERT(associated(cf%fs))
      SAFE_DEALLOCATE_P(cf%fs)
    end if
    
    POP_SUB(cube_function_free_fs)
  end subroutine cube_function_free_fs
  

  ! ---------------------------------------------------------  
  subroutine fourier_space_op_end(this)
    type(fourier_space_op_t), intent(inout) :: this
    
    PUSH_SUB(fourier_space_op_end)

    if(this%in_device_memory) then
      call accel_release_buffer(this%op_buffer)
      this%in_device_memory = .false.
    end if
    SAFE_DEALLOCATE_A(this%dop)
    SAFE_DEALLOCATE_A(this%zop)
    

    POP_SUB(fourier_space_op_end)
  end subroutine fourier_space_op_end

#include "undef.F90"
#include "real.F90"
#include "fourier_space_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "fourier_space_inc.F90"

end module fourier_space_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

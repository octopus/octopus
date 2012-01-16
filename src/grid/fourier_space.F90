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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module fourier_space_m
  use cube_m
  use cube_function_m
  use global_m
  use mesh_m
  use messages_m
  use fft_m
#ifdef HAVE_OPENMP
  use omp_lib
#endif
#ifdef HAVE_PFFT
  use pfft_m
#endif
  use profiling_m
  use simul_box_m

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
    fourier_space_op_end,       &
    dfourier_to_mesh,           &
    zfourier_to_mesh

  type fourier_space_op_t
    private
    FLOAT, pointer :: dop(:, :, :)
    CMPLX, pointer :: zop(:, :, :)
  end type fourier_space_op_t

contains

  ! ---------------------------------------------------------

  !> Allocates locally the Fourier space grid, if PFFT library is not used.
  !! Otherwise, it assigns the PFFT Fourier space grid to the cube Fourier space grid,
  !! via pointer.
  subroutine cube_function_alloc_fs(cube, cf)
    type(cube_t),          intent(in)    :: cube
    type(cube_function_t), intent(inout) :: cf
    
    integer :: n1, n2, n3

    PUSH_SUB(cube_function_alloc_fs)
    
    ASSERT(.not.associated(cf%fs))

    n1 = max(1, cube%fs_n(1))
    n2 = max(1, cube%fs_n(2))
    n3 = max(1, cube%fs_n(3))

    if (cube%fft_library /= FFTLIB_PFFT) then
      SAFE_ALLOCATE(cf%fs(1:cube%fs_n(1), 1:cube%fs_n(2), 1:cube%fs_n(3)))
    else
      ASSERT(associated(cube%fft))  
      if(any(cube%fs_n(1:3) == 0)) then
        cf%fs => cube%fft%fs_data(1:1,1:1,1:1)
      else
        cf%fs => cube%fft%fs_data(1:n3,1:n1,1:n2)
      end if
    end if
    
    POP_SUB(cube_function_alloc_fs)
  end subroutine cube_function_alloc_fs

  
  ! ---------------------------------------------------------
  !> Deallocates the Fourier space grid
  subroutine cube_function_free_fs(cube, cf)
    type(cube_t),          intent(in)    :: cube
    type(cube_function_t), intent(inout) :: cf
    
    PUSH_SUB(cube_function_free_fs)
    
    if (cube%fft_library /= FFTLIB_PFFT) then
      SAFE_DEALLOCATE_P(cf%fs)
    else
      nullify(cf%fs)
    end if
    
    POP_SUB(cube_function_free_fs)
  end subroutine cube_function_free_fs
  

  ! ---------------------------------------------------------  
  subroutine fourier_space_op_end(this)
    type(fourier_space_op_t), intent(inout) :: this
    
    PUSH_SUB(fourier_space_op_end)

    SAFE_DEALLOCATE_P(this%dop)
    SAFE_DEALLOCATE_P(this%zop)

    POP_SUB(fourier_space_op_end)
  end subroutine fourier_space_op_end

#include "undef.F90"
#include "real.F90"
#include "fourier_space_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "fourier_space_inc.F90"

end module fourier_space_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

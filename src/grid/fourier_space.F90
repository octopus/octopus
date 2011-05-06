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
  use cube_function_m
  use global_m
  use mesh_m
  use messages_m
  use fft_m
#ifdef HAVE_PFFT
  use pfft_m
#endif
  use profiling_m
  use simul_box_m

  implicit none
  private
  public ::                     &
    dcube_function_alloc_rs,       &
    zcube_function_alloc_rs,       &
    dcube_function_alloc_fs,       & 
    zcube_function_alloc_fs,       &
    dcube_function_free_fs,        &
    zcube_function_free_fs,        &
    dcube_function_fft_init,       &
    zcube_function_fft_init,       &
    cube_function_fft_end,         &
    dcube_function_pfft_init,      &
    zcube_function_pfft_init,      &
    dcube_function_RS2FS,          &
    zcube_function_RS2FS,          &
    dcube_function_FS2RS,          &
    zcube_function_FS2RS,          &
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

  subroutine cube_function_fft_end(cf)
    type(cube_function_t),     intent(inout) :: cf

    PUSH_SUB(cube_function_fft_end)

    if (cf%fft_library == PFFT_LIB) then
#ifdef HAVE_PFFT
      if(associated(cf%pfft)) then
        call pfft_end(cf%pfft)
        SAFE_DEALLOCATE_P(cf%pfft)
      end if
#else
    else
      if(associated(cf%fft)) then
        call fft_end(cf%fft)
        SAFE_DEALLOCATE_P(cf%fft)
      end if
#endif
    end if

    POP_SUB(cube_function_fft_end)

  end subroutine cube_function_fft_end

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

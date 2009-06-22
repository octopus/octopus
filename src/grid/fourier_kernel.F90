!! Copyright (C) 2009 X. Andrade
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
!! $Id: cube_function.F90 4907 2009-01-31 11:21:32Z xavier $

#include "global.h"

module fourier_kernel_m
  use cubic_mesh_m
  use global_m
  use mesh_m
  use messages_m
  use profiling_m
  use simul_box_m
  
  implicit none
  private

  public ::                           &
       fourier_kernel_t,              &
       fourier_kernel_init,           &
       fourier_kernel_end,            &
       fourier_kernel_size_real,      &
       fourier_kernel_size_fourier,   &
       dfourier_kernel_set_value
  
  type fourier_kernel_t
    private
    type(cubic_mesh_t) :: kernel
    integer            :: size_real(1:MAX_DIM)
    integer            :: size_fourier(1:MAX_DIM)
  end type fourier_kernel_t

contains

  subroutine fourier_kernel_init(this, mesh, sizes)
    type(fourier_kernel_t), intent(out) :: this
    type(mesh_t),           intent(in)  :: mesh
    integer,                intent(in)  :: sizes(1:MAX_DIM)

    this%size_real = sizes
    this%size_fourier = sizes
    this%size_fourier(1) = this%size_fourier(1)/2

    call cubic_mesh_init(this%kernel, mesh, sizes)

  end subroutine fourier_kernel_init

  !------------------------------------------------------------------------------

  subroutine fourier_kernel_end(this)
    type(fourier_kernel_t), intent(inout) :: this

    call cubic_mesh_end(this%kernel)
    
  end subroutine fourier_kernel_end

  !------------------------------------------------------------------------------
  function fourier_kernel_size_real(this) result(sizes)
    type(fourier_kernel_t), intent(out) :: this
    integer                             :: sizes(1:MAX_DIM)

    sizes = this%size_real
  end function fourier_kernel_size_real

  !------------------------------------------------------------------------------

  function fourier_kernel_size_fourier(this) result(sizes)
    type(fourier_kernel_t), intent(out) :: this
    integer                             :: sizes(1:MAX_DIM)

    sizes = this%size_fourier
  end function fourier_kernel_size_fourier

  !------------------------------------------------------------------------------

  subroutine dfourier_kernel_set_value(this, ix, iy, iz, val)
    type(fourier_kernel_t), intent(inout) :: this
    integer,                intent(in)    :: ix
    integer,                intent(in)    :: iy
    integer,                intent(in)    :: iz
    FLOAT,                  intent(in)    :: val

    call dcubic_mesh_set_point(this%kernel, ix, iy, iz, val)

  end subroutine dfourier_kernel_set_value
  
end module fourier_kernel_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

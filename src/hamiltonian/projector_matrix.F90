!! Copyright (C) 2010 X. Andrade
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
!! $Id: projector_matrix.F90 6489 2010-04-10 14:59:34Z xavier $

#include "global.h"

module projector_matrix_m
  use global_m
  use messages_m
#ifdef HAVE_OPENCL
  use opencl_m
#endif
  use types_m
  use profiling_m

  implicit none

  private

  public ::                     &
    projector_matrix_t,         &
    projector_matrix_nullify,   &
    projector_matrix_allocate,  &
    projector_matrix_deallocate

  type projector_matrix_t
    integer, pointer :: map(:)
    FLOAT,   pointer :: projectors(:, :)
    FLOAT,   pointer :: scal(:)
    integer          :: npoints
    integer          :: nprojs
#ifdef HAVE_OPENCL
    type(opencl_mem_t) :: buff_map
    type(opencl_mem_t) :: buff_projectors
    type(opencl_mem_t) :: buff_scal
    logical            :: buffers_allocated
#endif    
  end type projector_matrix_t

contains

  elemental subroutine projector_matrix_nullify(this)
    type(projector_matrix_t), intent(out) :: this

    nullify(this%map)
    nullify(this%projectors)
    nullify(this%scal)
    this%buffers_allocated = .false.

  end subroutine projector_matrix_nullify
  
  ! -------------------------------------------------

  subroutine projector_matrix_allocate(this, npoints, nprojs)
    type(projector_matrix_t), intent(out) :: this
    integer,                  intent(in)  :: npoints
    integer,                  intent(in)  :: nprojs

    this%npoints = npoints
    this%nprojs = nprojs

    SAFE_ALLOCATE(this%map(1:npoints))
    SAFE_ALLOCATE(this%projectors(1:npoints, 1:nprojs))
    SAFE_ALLOCATE(this%scal(1:nprojs))
#ifdef HAVE_OPENCL
    if(opencl_is_enabled()) then
      this%buffers_allocated = .true.
      call opencl_create_buffer(this%buff_map, CL_MEM_READ_ONLY, TYPE_INTEGER, npoints)
      call opencl_create_buffer(this%buff_scal, CL_MEM_READ_ONLY, TYPE_FLOAT, nprojs)
      call opencl_create_buffer(this%buff_projectors, CL_MEM_READ_ONLY, TYPE_FLOAT, nprojs*npoints)
    end if
#endif
  end subroutine projector_matrix_allocate

  ! -------------------------------------------------

  subroutine projector_matrix_deallocate(this)
    type(projector_matrix_t), intent(out) :: this

    SAFE_DEALLOCATE_P(this%map)
    SAFE_DEALLOCATE_P(this%projectors)
    SAFE_DEALLOCATE_P(this%scal)
#ifdef HAVE_OPENCL
    if(opencl_is_enabled() .and. this%buffers_allocated) then
      this%buffers_allocated = .false.
      call opencl_release_buffer(this%buff_map)
      call opencl_release_buffer(this%buff_scal)
      call opencl_release_buffer(this%buff_projectors)
    end if
#endif
  end subroutine projector_matrix_deallocate

  ! -------------------------------------------------

end module projector_matrix_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:


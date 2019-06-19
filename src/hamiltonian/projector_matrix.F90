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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module projector_matrix_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::                     &
    projector_matrix_t,         &
    projector_matrix_nullify,   &
    projector_matrix_allocate,  &
    projector_matrix_deallocate

  type projector_matrix_t
    ! Components are public by default
    integer, pointer :: map(:)
    FLOAT,   pointer :: projectors(:, :)
    FLOAT,   pointer :: scal(:)
    FLOAT,   allocatable :: position(:, :)
    integer          :: npoints
    integer          :: nprojs
    FLOAT,   allocatable :: mix(:, :)
  end type projector_matrix_t

contains

  elemental subroutine projector_matrix_nullify(this)
    type(projector_matrix_t), intent(out) :: this

    nullify(this%map)
    nullify(this%projectors)
    nullify(this%scal)

  end subroutine projector_matrix_nullify
  
  ! -------------------------------------------------

  subroutine projector_matrix_allocate(this, npoints, nprojs, has_mix_matrix)
    type(projector_matrix_t), intent(out) :: this
    integer,                  intent(in)  :: npoints
    integer,                  intent(in)  :: nprojs
    logical,                  intent(in)  :: has_mix_matrix

    PUSH_SUB(projector_matrix_allocate)

    this%npoints = npoints
    this%nprojs = nprojs

    SAFE_ALLOCATE(this%map(1:npoints))
    SAFE_ALLOCATE(this%projectors(1:npoints, 1:nprojs))
    SAFE_ALLOCATE(this%scal(1:nprojs))
    SAFE_ALLOCATE(this%position(1:3, 1:npoints))

    if(has_mix_matrix) then
      SAFE_ALLOCATE(this%mix(1:nprojs, 1:nprojs))
    end if

    POP_SUB(projector_matrix_allocate)
  end subroutine projector_matrix_allocate

  ! -------------------------------------------------

  subroutine projector_matrix_deallocate(this)
    type(projector_matrix_t), intent(inout) :: this

    PUSH_SUB(projector_matrix_deallocate)

    SAFE_DEALLOCATE_P(this%map)
    SAFE_DEALLOCATE_P(this%projectors)
    SAFE_DEALLOCATE_P(this%scal)
    SAFE_DEALLOCATE_A(this%position)
    SAFE_DEALLOCATE_A(this%mix)

    POP_SUB(projector_matrix_deallocate)
  end subroutine projector_matrix_deallocate

  ! -------------------------------------------------

end module projector_matrix_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

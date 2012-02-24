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
  use profiling_m
  use types_m

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
  end type projector_matrix_t

contains

  elemental subroutine projector_matrix_nullify(this)
    type(projector_matrix_t), intent(out) :: this
    !
    this%map=>null()
    this%projectors=>null()
    this%scal=>null()
    return
  end subroutine projector_matrix_nullify
  
  ! -------------------------------------------------

  subroutine projector_matrix_allocate(this, npoints, nprojs)
    type(projector_matrix_t), intent(out) :: this
    integer,                  intent(in)  :: npoints
    integer,                  intent(in)  :: nprojs

    PUSH_SUB(projector_matrix_allocate)

    this%npoints = npoints
    this%nprojs = nprojs

    SAFE_ALLOCATE(this%map(1:npoints))
    SAFE_ALLOCATE(this%projectors(1:npoints, 1:nprojs))
    SAFE_ALLOCATE(this%scal(1:nprojs))

    POP_SUB(projector_matrix_allocate)
  end subroutine projector_matrix_allocate

  ! -------------------------------------------------

  subroutine projector_matrix_deallocate(this)
    type(projector_matrix_t), intent(out) :: this

    PUSH_SUB(projector_matrix_deallocate)

    SAFE_DEALLOCATE_P(this%map)
    SAFE_DEALLOCATE_P(this%projectors)
    SAFE_DEALLOCATE_P(this%scal)

    POP_SUB(projector_matrix_deallocate)
  end subroutine projector_matrix_deallocate

  ! -------------------------------------------------

end module projector_matrix_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

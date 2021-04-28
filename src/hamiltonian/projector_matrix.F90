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
    projector_matrix_allocate,  &
    projector_matrix_deallocate

  type projector_matrix_t
    ! Components are public by default
    integer, allocatable :: map(:)
    FLOAT,   allocatable :: dprojectors(:, :)
    CMPLX,   allocatable :: zprojectors(:, :)
    FLOAT,   allocatable :: scal(:)
    FLOAT,   allocatable :: position(:, :)
    integer              :: npoints
    integer              :: nprojs
    FLOAT,   allocatable :: dmix(:, :)
    CMPLX,   allocatable :: zmix(:, :, :)
    logical              :: is_cmplx = .false.
  end type projector_matrix_t

contains

  ! -------------------------------------------------
  subroutine projector_matrix_allocate(this, npoints, nprojs, has_mix_matrix, is_cmplx)
    type(projector_matrix_t), intent(out) :: this
    integer,                  intent(in)  :: npoints
    integer,                  intent(in)  :: nprojs
    logical,                  intent(in)  :: has_mix_matrix
    logical, optional,        intent(in)  :: is_cmplx

    PUSH_SUB(projector_matrix_allocate)

    this%npoints = npoints
    this%nprojs = nprojs

    this%is_cmplx = optional_default(is_cmplx, .false.)

    SAFE_ALLOCATE(this%map(1:npoints))
    if(this%is_cmplx) then
      SAFE_ALLOCATE(this%zprojectors(1:npoints, 1:nprojs))
    else
      SAFE_ALLOCATE(this%dprojectors(1:npoints, 1:nprojs))
    end if
    SAFE_ALLOCATE(this%scal(1:nprojs))
    SAFE_ALLOCATE(this%position(1:3, 1:npoints))

    if(has_mix_matrix) then
      if(this%is_cmplx) then
        SAFE_ALLOCATE(this%zmix(1:nprojs, 1:nprojs, 1:4))
      else
        SAFE_ALLOCATE(this%dmix(1:nprojs, 1:nprojs))
      end if
    end if

    POP_SUB(projector_matrix_allocate)
  end subroutine projector_matrix_allocate

  ! -------------------------------------------------

  subroutine projector_matrix_deallocate(this)
    type(projector_matrix_t), intent(inout) :: this

    PUSH_SUB(projector_matrix_deallocate)

    SAFE_DEALLOCATE_A(this%map)
    SAFE_DEALLOCATE_A(this%dprojectors)
    SAFE_DEALLOCATE_A(this%zprojectors)
    SAFE_DEALLOCATE_A(this%scal)
    SAFE_DEALLOCATE_A(this%position)
    SAFE_DEALLOCATE_A(this%dmix)
    SAFE_DEALLOCATE_A(this%zmix)

    POP_SUB(projector_matrix_deallocate)
  end subroutine projector_matrix_deallocate

  ! -------------------------------------------------

end module projector_matrix_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

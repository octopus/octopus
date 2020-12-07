!! Copyright (C) 2008 X. Andrade
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

module stencil_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                        &
    stencil_t,                     &
    stencil_allocate,              &
    stencil_copy,                  &
    stencil_end,                   &
    stencil_init_center

  type stargeneral_arms_t
    ! Components are public by default
    integer          :: arms(1:3,1:3)
    integer          :: narms  
  end type stargeneral_arms_t


  type stencil_t
    ! Components are public by default
    integer          :: center
    integer          :: size
    integer          :: npoly
    integer, pointer :: points(:, :) 
    
    ! The stargeneral arms
    type(stargeneral_arms_t) :: stargeneral 
  end type stencil_t

contains

  !-------------------------------------------------------  
  subroutine stencil_allocate(this, size)
    type(stencil_t), intent(out) :: this
    integer,         intent(in)  :: size

    PUSH_SUB(stencil_allocate)

    this%size = size
    this%npoly = size

    SAFE_ALLOCATE(this%points(1:MAX_DIM, 1:size))

    this%points = 0

    POP_SUB(stencil_allocate)
  end subroutine stencil_allocate

  !-------------------------------------------------------  
  subroutine stencil_copy(input, output)
    type(stencil_t), intent(in)  :: input
    type(stencil_t), intent(out) :: output
    
    PUSH_SUB(stencil_copy)

    call stencil_allocate(output, input%size)
    output%points(1:MAX_DIM, 1:output%size) = input%points(1:MAX_DIM, 1:output%size)
    output%center = input%center
    output%npoly = input%npoly

    POP_SUB(stencil_copy)
  end subroutine stencil_copy


  !-------------------------------------------------------  
  subroutine stencil_end(this)
    type(stencil_t), intent(inout) :: this

    PUSH_SUB(stencil_end)

    SAFE_DEALLOCATE_P(this%points)

    POP_SUB(stencil_end)
  end subroutine stencil_end

  
  !-------------------------------------------------------  
  subroutine stencil_init_center(this)
    type(stencil_t), intent(inout) :: this

    integer :: ii

    PUSH_SUB(stencil_init_center)

    this%center = -1

    do ii = 1, this%size
      if(all(this%points(1:MAX_DIM, ii) == 0)) this%center = ii
    end do
    
    POP_SUB(stencil_init_center)
  end subroutine stencil_init_center

end module stencil_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

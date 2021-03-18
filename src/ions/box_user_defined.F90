!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2021 M. Oliveira
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

module box_user_defined_oct_m
  use box_oct_m
  use box_parallelepiped_oct_m
  use box_shape_oct_m
  use global_oct_m
  use messages_oct_m
  use parser_oct_m
  use profiling_oct_m
  use string_oct_m

  implicit none

  public :: box_user_defined_t

  !> Class implementing a box defined by a mathematical expression.
  !! This box needs to be inside a parallelepiped in order to avoid problems
  !! when the mathematical expression does not define a closed box in one or
  !! more directions. For example, in 2D, the expression x > y is not bounded,
  !! but if we intersect it with a square box, it defines a triangle.
  type, extends(box_shape_t) :: box_user_defined_t
    private
    type(box_parallelepiped_t), pointer :: outer_box  !< the outer parallelepiped containing the user-defined box
    character(len=1024)                 :: expression !< the mathematical expression defining th box boundaries
  contains
    procedure :: contains_points => box_user_defined_contains_points
    procedure :: write_info => box_user_defined_write_info
    procedure :: short_info => box_user_defined_short_info
    final     :: box_user_defined_finalize
  end type box_user_defined_t

  interface box_user_defined_t
    procedure box_user_defined_constructor
  end interface box_user_defined_t

contains

  !--------------------------------------------------------------
  function box_user_defined_constructor(dim, center, expression, length) result(box)
    integer,             intent(in) :: dim
    FLOAT,               intent(in) :: center(dim)
    character(len=1024), intent(in) :: expression
    FLOAT,               intent(in) :: length(dim) !< length of the parallelepiped along each Cartesian direction
    class(box_user_defined_t), pointer :: box

    PUSH_SUB(box_user_defined_constructor)

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    call box_shape_init(box, dim, center)
    box%expression = expression
    call conv_to_C_string(box%expression)
    box%outer_box => box_parallelepiped_t(dim, center, length)

    POP_SUB(box_user_defined_constructor)
  end function box_user_defined_constructor

  !--------------------------------------------------------------
  subroutine box_user_defined_finalize(this)
    type(box_user_defined_t), intent(inout) :: this

    PUSH_SUB(box_user_defined_finalize)

    call box_shape_end(this)
    SAFE_DEALLOCATE_P(this%outer_box)

    POP_SUB(box_user_defined_finalize)
  end subroutine box_user_defined_finalize

  !--------------------------------------------------------------
  function box_user_defined_contains_points(this, nn, xx) result(contained)
    class(box_user_defined_t), intent(in)  :: this
    integer,                   intent(in)  :: nn
    FLOAT,                     intent(in)  :: xx(:,:)
    logical :: contained(1:nn)

    integer :: ip
    FLOAT :: re, im, rr, xx_centered(this%dim)

    contained = this%outer_box%contains_points(nn, xx)

    do ip = 1, nn
      if (.not. contained(ip)) then
        ! Skip points that are outside the outer box
        cycle
      end if

      xx_centered = xx(ip, :) - this%center
      rr = sqrt(sum(xx_centered**2))
      call parse_expression(re, im, this%dim, xx_centered, rr, M_ZERO, this%expression)
      contained(ip) = re /= M_ZERO .neqv. this%is_inside_out()
    end do

  end function box_user_defined_contains_points

  !--------------------------------------------------------------
  subroutine box_user_defined_write_info(this, iunit)
    class(box_user_defined_t), intent(in) :: this
    integer,                   intent(in) :: iunit

    PUSH_SUB(box_user_defined_write_info)

    write(message(1),'(2x,a)') 'Type = user-defined'
    call messages_info(1, iunit)

    POP_SUB(box_user_defined_write_info)
  end subroutine box_user_defined_write_info

  !--------------------------------------------------------------
  character(len=BOX_INFO_LEN) function box_user_defined_short_info(this, unit_length) result(info)
    class(box_user_defined_t), intent(in) :: this
    type(unit_t),              intent(in) :: unit_length

    PUSH_SUB(box_user_defined_short_info)

    write(info,'(3a)') 'BoxShape = user_defined; BoxShapeUsDef = "', trim(this%expression), '"'

    POP_SUB(box_user_defined_short_info)
  end function box_user_defined_short_info

end module box_user_defined_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

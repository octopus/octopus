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

module box_sphere_oct_m
  use box_oct_m
  use box_shape_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  public :: box_sphere_t

  !> Class implementing a spherical box.
  type, extends(box_shape_t) :: box_sphere_t
    FLOAT   :: radius !< the radius of the sphere
  contains
    procedure :: contains_points => box_sphere_contains_points
    procedure :: write_info => box_sphere_write_info
    procedure :: short_info => box_sphere_short_info
    final     :: box_sphere_finalize
  end type box_sphere_t

  interface box_sphere_t
    procedure box_sphere_constructor
  end interface box_sphere_t

contains

  !--------------------------------------------------------------
  function box_sphere_constructor(dim, center, radius) result(box)
    integer, intent(in) :: dim
    FLOAT,   intent(in) :: center(dim)
    FLOAT,   intent(in) :: radius
    class(box_sphere_t), pointer :: box

    PUSH_SUB(box_sphere_constructor)

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    call box_shape_init(box, dim, center)
    box%radius = radius

    POP_SUB(box_sphere_constructor)
  end function box_sphere_constructor

  !--------------------------------------------------------------
  subroutine box_sphere_finalize(this)
    type(box_sphere_t), intent(inout) :: this

    PUSH_SUB(box_sphere_finalize)

    call box_shape_end(this)

    POP_SUB(box_sphere_finalize)
  end subroutine box_sphere_finalize

  !--------------------------------------------------------------
  function box_sphere_contains_points(this, nn, xx) result(contained)
    class(box_sphere_t), intent(in)  :: this
    integer,             intent(in)  :: nn
    FLOAT,               intent(in)  :: xx(:,:)
    logical :: contained(1:nn)

    integer :: ip

    do ip = 1, nn
      contained(ip) = sum((xx(ip, 1:this%dim) - this%center(1:this%dim))**2) <= (this%radius + BOX_BOUNDARY_DELTA)**2 &
        .neqv. this%is_inside_out()
    end do

  end function box_sphere_contains_points

  !--------------------------------------------------------------
  subroutine box_sphere_write_info(this, iunit)
    class(box_sphere_t), intent(in) :: this
    integer,             intent(in) :: iunit

    PUSH_SUB(box_sphere_write_info)

    write(message(1),'(2x,a)') 'Type = sphere'
    write(message(2),'(2x,3a,f7.3)') 'Radius  [', trim(units_abbrev(units_out%length)), '] = ', &
      units_from_atomic(units_out%length, this%radius)
    call messages_info(2, iunit)

    POP_SUB(box_sphere_write_info)
  end subroutine box_sphere_write_info

  !--------------------------------------------------------------
  character(len=BOX_INFO_LEN) function box_sphere_short_info(this, unit_length) result(info)
    class(box_sphere_t), intent(in) :: this
    type(unit_t),        intent(in) :: unit_length

    PUSH_SUB(box_sphere_short_info)

    write(info,'(a,f11.6,a,a)') 'BoxShape = sphere; Radius =', units_from_atomic(unit_length, this%radius), ' ', &
      trim(units_abbrev(unit_length))

    POP_SUB(box_sphere_short_info)
  end function box_sphere_short_info

end module box_sphere_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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

module box_cylinder_oct_m
  use box_oct_m
  use box_shape_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  public :: box_cylinder_t

  !> Class implementing a cylinder box. Currently this is restricted to right
  !! circular cylinders.
  type, extends(box_shape_t) :: box_cylinder_t
    private
    integer :: dir         !< axis along which the cylinder lies (for the moment
                           !! must be one of the Cartesian axes, but this could
                           !! be generalized)
    FLOAT, public :: radius      !< the radius of the cylinder
    FLOAT, public :: half_length !< half the length of the cylinder

    logical :: periodic_boundaries = .false. !< are the bases of the cylinder to be treated as periodic?
  contains
    procedure :: contains_points => box_cylinder_contains_points
    procedure :: write_info => box_cylinder_write_info
    procedure :: short_info => box_cylinder_short_info
    final     :: box_cylinder_finalize
  end type box_cylinder_t

  interface box_cylinder_t
    procedure box_cylinder_constructor
  end interface box_cylinder_t

contains

  !--------------------------------------------------------------
  function box_cylinder_constructor(dim, center, radius, dir, length, namespace, periodic_boundaries) result(box)
    integer,            intent(in) :: dim
    FLOAT,              intent(in) :: center(dim)
    FLOAT,              intent(in) :: radius              !< cylinder radius
    integer,            intent(in) :: dir                 !< cartesian direction along which the cylinder lies
    FLOAT,              intent(in) :: length              !< lenght of the cylinder along the axis
    type(namespace_t),  intent(in) :: namespace
    logical, optional,  intent(in) :: periodic_boundaries !< are the bases of the cylinder to be treated as periodic?
    class(box_cylinder_t), pointer :: box

    PUSH_SUB(box_cylinder_constructor)

    ! Sanity checks
    if (dim <= 2) then
      message(1) = "Cannot create a cylinder in 1D or 2D. Use sphere if you want a circle."
      call messages_fatal(1, namespace=namespace)
    end if
    if (dir > dim) then
      message(1) = "Direction for cylinder axis cannot be larger than the dimension of the space in which the cylinder lives."
      call messages_fatal(1, namespace=namespace)
    end if

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    call box_shape_init(box, dim, center)
    box%radius = radius
    box%dir = dir
    box%half_length = M_HALF*length

    if (present(periodic_boundaries)) then
      box%periodic_boundaries = periodic_boundaries
    end if
    if (box%periodic_boundaries .and. dir /= 1) then
      ! Currently the rest of the code expects the first periodic dimension to
      ! be always x and cylinders can only be periodic along the axis. Therefore
      ! we can only have periodic cylinders along the x axis.
      message(1) = "Cylinder can only have periodic boundaries if the axis is along the x direction."
      call messages_fatal(1, namespace=namespace)
    end if

    POP_SUB(box_cylinder_constructor)
  end function box_cylinder_constructor

  !--------------------------------------------------------------
  subroutine box_cylinder_finalize(this)
    type(box_cylinder_t), intent(inout) :: this

    PUSH_SUB(box_cylinder_finalize)

    call box_shape_end(this)

    POP_SUB(box_cylinder_finalize)
  end subroutine box_cylinder_finalize

  !--------------------------------------------------------------
  recursive function box_cylinder_contains_points(this, nn, xx) result(contained)
    class(box_cylinder_t), intent(in)  :: this
    integer,               intent(in)  :: nn
    FLOAT,                 intent(in)  :: xx(:,:)
    logical :: contained(1:nn)

    integer :: ip, idim
    FLOAT :: rr2, vv(this%dim)

    do ip = 1, nn
      vv = xx(ip, 1:this%dim) - this%center(1:this%dim)

      ! First check if we are "inside" along the axis direction. If not, do not bother checking the other directions.
      if (.not. this%periodic_boundaries) then
        contained(ip) = abs(vv(this%dir)) <= this%half_length + BOX_BOUNDARY_DELTA .neqv. this%is_inside_out()
      else
        ! When periodic, we exclude one of the faces from the box. Also, we
        ! never consider the box to be inside out along the periodic dimension.
        contained(ip) = abs(vv(this%dir) + BOX_BOUNDARY_DELTA)  <= this%half_length
      end if
      if (.not. contained(ip)) cycle

      ! Check if we are inside along the directions perpendicular to the axis
      rr2 = M_ZERO
      do idim = 1, this%dim
        if (idim == this%dir) cycle
        rr2 = rr2 + vv(idim)**2
      end do
      contained(ip) = rr2 <= (this%radius + BOX_BOUNDARY_DELTA)**2 .neqv. this%is_inside_out()
    end do

  end function box_cylinder_contains_points

  !--------------------------------------------------------------
  subroutine box_cylinder_write_info(this, iunit)
    class(box_cylinder_t), intent(in) :: this
    integer,               intent(in) :: iunit

    PUSH_SUB(box_cylinder_write_info)

    write(message(1),'(2x,a)') 'Type = cylinder'
    write(message(2),'(2x,3a,f7.3)') 'Radius  [', trim(units_abbrev(units_out%length)), '] = ', &
      units_from_atomic(units_out%length, this%radius)
    write(message(3),'(2x,3a,f7.3)') 'Xlength [', trim(units_abbrev(units_out%length)), '] = ', &
      units_from_atomic(units_out%length, this%half_length)
    call messages_info(3, iunit)

    POP_SUB(box_cylinder_write_info)
  end subroutine box_cylinder_write_info

  !--------------------------------------------------------------
  character(len=BOX_INFO_LEN) function box_cylinder_short_info(this, unit_length) result(info)
    class(box_cylinder_t), intent(in) :: this
    type(unit_t),          intent(in) :: unit_length

    PUSH_SUB(box_cylinder_short_info)

    write(info, '(a,f11.6,a,a,a,f11.6,a,a)') 'BoxShape = cylinder, Radius =', units_from_atomic(unit_length, this%radius), ' ', &
      trim(units_abbrev(unit_length)), '; Xlength =', units_from_atomic(unit_length, this%half_length), ' ', &
      trim(units_abbrev(unit_length))

    POP_SUB(box_cylinder_short_info)
  end function box_cylinder_short_info

end module box_cylinder_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

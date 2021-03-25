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

module box_parallelepiped_oct_m
  use box_oct_m
  use box_shape_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  public :: box_parallelepiped_t

  !> Class implementing a parallelepiped box. Currently this is restricted to a
  !! rectangular cuboid (all the faces must be rectangular) and the vectors
  !! generating the parallelepiped must be along the Cartesian axes.
  type, extends(box_shape_t) :: box_parallelepiped_t
    private
    FLOAT, allocatable, public :: half_length(:) !< half the length of the parallelepiped in each direction.

    integer, public :: n_periodic_boundaries = 0 !< in how many directions the parallelepiped boundaries are periodic
  contains
    procedure :: contains_points => box_parallelepiped_contains_points
    procedure :: write_info => box_parallelepiped_write_info
    procedure :: short_info => box_parallelepiped_short_info
    final     :: box_parallelepiped_finalize
  end type box_parallelepiped_t

  interface box_parallelepiped_t
    procedure box_parallelepiped_constructor
  end interface box_parallelepiped_t

contains

  !--------------------------------------------------------------
  function box_parallelepiped_constructor(dim, center, length, n_periodic_boundaries) result(box)
    integer,            intent(in) :: dim
    FLOAT,              intent(in) :: center(dim)
    FLOAT,              intent(in) :: length(dim)           !< length of the parallelepiped along each Cartesian direction
    integer, optional,  intent(in) :: n_periodic_boundaries !< in how many directions the parallelepiped boundaries are periodic
    class(box_parallelepiped_t), pointer :: box

    PUSH_SUB(box_parallelepiped_constructor)

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    call box_shape_init(box, dim, center)
    SAFE_ALLOCATE(box%half_length(1:dim))
    box%half_length = M_HALF*length
    if (present(n_periodic_boundaries)) then
      box%n_periodic_boundaries = n_periodic_boundaries
    end if

    POP_SUB(box_parallelepiped_constructor)
  end function box_parallelepiped_constructor

  !--------------------------------------------------------------
  subroutine box_parallelepiped_finalize(this)
    type(box_parallelepiped_t), intent(inout) :: this

    PUSH_SUB(box_parallelepiped_finalize)

    call box_shape_end(this)
    SAFE_DEALLOCATE_A(this%half_length)

    POP_SUB(box_parallelepiped_finalize)
  end subroutine box_parallelepiped_finalize

  !--------------------------------------------------------------
  function box_parallelepiped_contains_points(this, nn, xx) result(contained)
    class(box_parallelepiped_t), intent(in)  :: this
    integer,                     intent(in)  :: nn
    FLOAT,                       intent(in)  :: xx(:,:)
    logical :: contained(1:nn)

    integer :: ip, idir
    FLOAT :: ulimit(this%dim), llimit(this%dim)

    llimit = this%center - this%half_length - BOX_BOUNDARY_DELTA
    do idir = 1, this%dim
      if (idir <= this%n_periodic_boundaries) then
        ! When periodic, we exclude one of the faces from the box.
        ulimit(idir) = this%center(idir) + this%half_length(idir) - BOX_BOUNDARY_DELTA
      else
        ulimit(idir) = this%center(idir) + this%half_length(idir) + BOX_BOUNDARY_DELTA
      end if
    end do

    do ip = 1, nn
      contained(ip) = .true.
      do idir = 1, this%dim
        contained(ip) = contained(ip) .and. xx(ip, idir) >= llimit(idir) .and. xx(ip, idir) <= ulimit(idir)
        if (idir > this%n_periodic_boundaries) then
          ! We only consider the box to be inside out along the non-periodic directions.
          contained(ip) = contained(ip) .neqv. this%is_inside_out()
        end if
      end do
    end do

  end function box_parallelepiped_contains_points

  !--------------------------------------------------------------
  subroutine box_parallelepiped_write_info(this, iunit)
    class(box_parallelepiped_t), intent(in) :: this
    integer,                     intent(in) :: iunit

    integer :: idir

    PUSH_SUB(box_parallelepiped_write_info)

    write(message(1),'(2x,a)') 'Type = parallelepiped'
    write(message(2),'(2x,3a, 99(f8.3,a))') 'Lengths [', trim(units_abbrev(units_out%length)), '] = (', &
      (units_from_atomic(units_out%length, M_TWO*this%half_length(idir)), ',', idir = 1, this%dim - 1), &
      units_from_atomic(units_out%length, M_TWO*this%half_length(this%dim)), ')'
    call messages_info(2, iunit)

    POP_SUB(box_parallelepiped_write_info)
  end subroutine box_parallelepiped_write_info

  !--------------------------------------------------------------
  character(len=BOX_INFO_LEN) function box_parallelepiped_short_info(this, unit_length) result(info)
    class(box_parallelepiped_t), intent(in) :: this
    type(unit_t),                intent(in) :: unit_length

    integer :: idir

    PUSH_SUB(box_parallelepiped_short_info)

    write(info, '(a,a,a,99(f11.6,a))') 'BoxShape = parallelepiped; Lengths [', trim(units_abbrev(unit_length)),'] = [', &
      (units_from_atomic(unit_length, M_TWO*this%half_length(idir)), ',', idir = 1, this%dim - 1), &
      units_from_atomic(unit_length, M_TWO*this%half_length(this%dim)), ']'

    POP_SUB(box_parallelepiped_short_info)
  end function box_parallelepiped_short_info

end module box_parallelepiped_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

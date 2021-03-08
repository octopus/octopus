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

module box_hypercube_oct_m
  use box_parallelepiped_oct_m
  use box_shape_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  public :: box_hypercube_t

  !> Class implementing an hypercube box. Currently this is restricted to a
  !! rectangular cuboid (all the faces must be rectangular) and the vectors
  !! generating the hypercube must be along the Cartesian axes.
  type, extends(box_parallelepiped_t) :: box_hypercube_t
    private
  contains
    procedure :: write_info => box_hypercube_write_info
    procedure :: write_short_info => box_hypercube_write_short_info
    final     :: box_hypercube_finalize
  end type box_hypercube_t

  interface box_hypercube_t
    procedure box_hypercube_constructor
  end interface box_hypercube_t

contains

  !--------------------------------------------------------------
  function box_hypercube_constructor(dim, center, length, n_periodic_boundaries) result(box)
    integer,            intent(in) :: dim
    FLOAT,              intent(in) :: center(dim)
    FLOAT,              intent(in) :: length(dim)           !< length of the hypercube along each Cartesian direction
    integer, optional,  intent(in) :: n_periodic_boundaries !< in how many directions the hypercube boundaries are periodic
    class(box_hypercube_t), pointer :: box

    PUSH_SUB(box_hypercube_constructor)

    ! Allocate memory
    SAFE_ALLOCATE(box)

    ! Initialize box
    call box_shape_init(box, dim, center)
    SAFE_ALLOCATE(box%half_length(1:dim))
    box%half_length = M_HALF*length
    if (present(n_periodic_boundaries)) then
      box%n_periodic_boundaries = n_periodic_boundaries
    end if

    POP_SUB(box_hypercube_constructor)
  end function box_hypercube_constructor

  !--------------------------------------------------------------
  subroutine box_hypercube_finalize(this)
    type(box_hypercube_t), intent(inout) :: this

    PUSH_SUB(box_hypercube_finalize)

    call box_shape_end(this)
    SAFE_DEALLOCATE_A(this%half_length)

    POP_SUB(box_hypercube_finalize)
  end subroutine box_hypercube_finalize

  !--------------------------------------------------------------
  subroutine box_hypercube_write_info(this, iunit)
    class(box_hypercube_t), intent(in) :: this
    integer,                     intent(in) :: iunit

    integer :: idir

    PUSH_SUB(box_hypercube_write_info)

    write(message(1),'(2x,a)') 'Type = hypercube'
    write(message(2),'(2x,3a, 99(f8.3,a))') 'Lengths [', trim(units_abbrev(units_out%length)), '] = (', &
      (units_from_atomic(units_out%length, M_TWO*this%half_length(idir)), ',', idir = 1, this%dim - 1), &
      units_from_atomic(units_out%length, M_TWO*this%half_length(this%dim)), ')'
    call messages_info(2, iunit)

    POP_SUB(box_hypercube_write_info)
  end subroutine box_hypercube_write_info

  !--------------------------------------------------------------
  subroutine box_hypercube_write_short_info(this, iunit)
    class(box_hypercube_t), intent(in) :: this
    integer,                intent(in) :: iunit

    PUSH_SUB(box_hypercube_write_short_info)

    write(iunit, '(a)') 'BoxShape = hypercube'  ! add parameters?

    POP_SUB(box_hypercube_write_short_info)
  end subroutine box_hypercube_write_short_info

end module box_hypercube_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

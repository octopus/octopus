!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module space_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m

  implicit none

  private

  public ::                   &
    space_t,                  &
    space_init

  type space_t
    ! Components are public by default
    integer :: dim
    integer :: periodic_dim = 0

  contains
    procedure :: is_periodic => space_is_periodic
    procedure :: write_info => space_write_info
    procedure :: short_info => space_short_info
  end type space_t

contains

  ! ---------------------------------------------------------
  subroutine space_init(this, namespace)
    type(space_t),     intent(inout) :: this
    type(namespace_t), intent(in)    :: namespace

    integer, parameter :: default_ndim = 3

    PUSH_SUB(space_init)
    
    !%Variable Dimensions
    !%Type integer
    !%Section System
    !%Default 3
    !%Description
    !% <tt>Octopus</tt> can run in 1, 2 or 3 dimensions, depending on the value of this
    !% variable (or more, if configured with <tt>--with-max-dim=4</tt> or higher).
    !% Note that not all input variables may be available in all cases.
    !%End
    call parse_variable(namespace, 'Dimensions', default_ndim, this%dim)
    if((this%dim>MAX_DIM).or.(this%dim<1)) call messages_input_error(namespace, 'Dimensions')

    !%Variable PeriodicDimensions
    !%Type integer
    !%Default 0
    !%Section System
    !%Description
    !% Define how many directions are to be considered periodic. It has to be a number
    !% between zero and <tt>Dimensions</tt>.
    !%Option 0
    !% No direction is periodic (molecule).
    !%Option 1
    !% The <i>x</i> direction is periodic.
    !%Option 2
    !% The <i>x</i> and <i>y</i> directions are periodic.
    !%Option 3
    !% The <i>x</i>, <i>y</i>, and <i>z</i> directions are periodic.
    !%End
    call parse_variable(namespace, 'PeriodicDimensions', 0, this%periodic_dim)

    if ((this%periodic_dim < 0) .or. (this%periodic_dim > MAX_DIM) .or. (this%periodic_dim > this%dim)) then
      call messages_input_error(namespace, 'PeriodicDimensions')
    end if

    POP_SUB(space_init)
  end subroutine space_init

  !--------------------------------------------------------------
  logical pure function space_is_periodic(this)
    class(space_t), intent(in) :: this

    space_is_periodic = this%periodic_dim > 0

  end function space_is_periodic

  !--------------------------------------------------------------
  subroutine space_write_info(this, iunit)
    class(space_t), intent(in) :: this
    integer,        intent(in) :: iunit

    PUSH_SUB(space_write_info)

    call messages_print_stress(iunit, "Space")

    write(message(1), '(a,i1,a)') 'Octopus will run in ', this%dim, ' dimension(s).'
    write(message(2), '(a,i1,a)') 'Octopus will treat the system as periodic in ', this%periodic_dim, ' dimension(s).'
    call messages_info(2, iunit)

    call messages_print_stress(iunit)

    POP_SUB(space_write_info)
  end subroutine space_write_info

  !--------------------------------------------------------------
  character(len=40) function space_short_info(this) result(info)
    class(space_t), intent(in) :: this

    PUSH_SUB(space_short_info)

    write(info, '(a,i1,a,i1)') 'Dimensions = ', this%dim, '; PeriodicDimensions = ', this%periodic_dim

    POP_SUB(space_short_info)
  end function space_short_info


end module space_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:


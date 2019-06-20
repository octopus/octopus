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
  use parser_oct_m

  implicit none

  private
  public ::       &
    operator(==), &
    operator(/=)

  public ::                   &
    space_t,                  &
    space_init,               &
    space_copy,               &
    space_end

  integer, parameter :: default_ndim = 3

  type space_t
    ! Components are public by default
    integer :: dim
  end type space_t

  interface operator(==)
    module procedure space_equal
  end interface operator(==)

  interface operator(/=)
    module procedure space_not_equal
  end interface operator(/=)

contains

  ! ---------------------------------------------------------
  subroutine space_init(this, dim)
    type(space_t),     intent(inout) :: this
    integer, optional, intent(in)    :: dim

    PUSH_SUB(space_init_simple)
    
    if(present(dim))then
      this%dim=dim
    else
      !%Variable Dimensions
      !%Type integer
      !%Section System
      !%Default 3
      !%Description
      !% <tt>Octopus</tt> can run in 1, 2 or 3 dimensions, depending on the value of this
      !% variable (or more, if configured with <tt>--with-max-dim=4</tt> or higher).
      !% Note that not all input variables may be available in all cases.
      !%End
      call parse_variable(parser, 'Dimensions', default_ndim, this%dim)
    end if
    if((this%dim>MAX_DIM).or.(this%dim<1)) call messages_input_error('Dimensions')

    POP_SUB(space_init_simple)
  end subroutine space_init

  ! ---------------------------------------------------------
  elemental subroutine space_copy(this_out, this_in)
    type(space_t), intent(inout) :: this_out
    type(space_t), intent(in)    :: this_in

    this_out%dim=this_in%dim
  end subroutine space_copy

  ! -----------------------------------------------------
  elemental function space_equal(sa, sb) result(eqv)
    type(space_t), intent(in) :: sa
    type(space_t), intent(in) :: sb

    logical :: eqv

    eqv=(sa%dim==sb%dim)
  end function space_equal
  
  ! -----------------------------------------------------
  elemental function space_not_equal(sa, sb) result(neqv)
    type(space_t), intent(in) :: sa
    type(space_t), intent(in) :: sb
    
    logical :: neqv

    neqv=(sa%dim/=sb%dim)
  end function space_not_equal

  ! ---------------------------------------------------------
  elemental subroutine space_end(this)
    type(space_t), intent(inout) :: this

    this%dim=0
  end subroutine space_end

end module space_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:


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
  subroutine space_init(this)
    type(space_t),     intent(inout) :: this

    PUSH_SUB(space_init)
    
    this%dim = 3
    call messages_obsolete_variable('Dimensions')

    POP_SUB(space_init)
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


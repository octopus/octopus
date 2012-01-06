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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: space.F90 6963 2010-08-25 09:10:42Z xavier $

#include "global.h"

module space_m
  use datasets_m
  use global_m
  use messages_m
  use parser_m

  implicit none

  private
  public ::                &
    space_t,               &
    space_init,            &
    space_init_from_dump,  &
    space_dump,            &
    space_copy,            &
    space_end,             &
    operator(==),          &
    operator(/=)

  interface operator(==)
    module procedure space_equal
  end interface

  interface operator(/=)
    module procedure space_not_equal
  end interface

  type space_t
    integer :: dim
  end type space_t

contains

  ! ---------------------------------------------------------
  subroutine space_init(this)
    type(space_t), intent(out) :: this

    !%Variable Dimensions
    !%Type integer
    !%Section System
    !%Default 3
    !%Description
    !% <tt>Octopus</tt> can run in 1, 2 or 3 dimensions, depending on the value of this
    !% variable. Note that not all input variables may be available in all cases.
    !%End
    call parse_integer(datasets_check('Dimensions'), 3, this%dim)
    if( this%dim > MAX_DIM .or. this%dim < 1) call input_error('Dimensions')
    return
  end subroutine space_init

  ! ---------------------------------------------------------
  subroutine space_init_from_dump(this, iunit)
    type(space_t), intent(inout) :: this
    integer,       intent(in)    :: iunit
    !
    integer :: gb
    !
    read(iunit) gb
    ASSERT(gb==GUARD_BITS)
    read(iunit) this
    read(iunit) gb
    ASSERT(gb==GUARD_BITS)
    return
  end subroutine space_init_from_dump

  ! ---------------------------------------------------------
  subroutine space_dump(this, iunit)
    type(space_t), intent(in) :: this
    integer,       intent(in) :: iunit
    !
    write(iunit) GUARD_BITS
    write(iunit) this
    write(iunit) GUARD_BITS
    return
  end subroutine space_dump

  ! ---------------------------------------------------------
  subroutine space_copy(this_out, this_in)
    type(space_t), intent(inout) :: this_out
    type(space_t), intent(in)    :: this_in
    !
    this_out%dim=this_in%dim
    return
  end subroutine space_copy

  ! -----------------------------------------------------
  elemental function space_equal(sa, sb) result(eqv)
    type(space_t), intent(in) :: sa
    type(space_t), intent(in) :: sb
    !
    logical :: eqv
    !
    eqv=(sa%dim==sb%dim)
    return
  end function space_equal
  
  ! -----------------------------------------------------
  elemental function space_not_equal(sa, sb) result(neqv)
    type(space_t), intent(in) :: sa
    type(space_t), intent(in) :: sb
    !
    logical :: neqv
    !
    neqv=(sa%dim/=sb%dim)
    return
  end function space_not_equal

  ! ---------------------------------------------------------
  subroutine space_end(this)
    type(space_t), intent(inout) :: this
    !
    this%dim=0
    return
  end subroutine space_end

end module space_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

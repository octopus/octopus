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
!! $Id: space.F90 6963 2010-08-25 09:10:42Z xavier $

#include "global.h"

module space_m

  use datasets_m
  use global_m
  use messages_m
  use parser_m

  use json_m, only: JSON_OK, json_object_t, json_init, json_set, json_get

  implicit none

  private
  public ::       &
    operator(==), &
    operator(/=)

  public ::                   &
    space_t,                  &
    space_init,               &
    space_create_data_object, &
    space_copy,               &
    space_end

  integer, parameter :: default_ndim = 3

  type space_t
    integer :: dim = 0
  end type space_t

  interface operator(==)
    module procedure space_equal
  end interface operator(==)

  interface operator(/=)
    module procedure space_not_equal
  end interface operator(/=)

  interface space_init
    module procedure space_init_simple
    module procedure space_init_data_object
  end interface space_init

contains

  ! ---------------------------------------------------------
  subroutine space_init_simple(this, dim)
    type(space_t),     intent(inout) :: this
    integer, optional, intent(in)    :: dim
    !
    if(present(dim))then
      this%dim=dim
    else
      !%Variable Dimensions
      !%Type integer
      !%Section System
      !%Default 3
      !%Description
      !% <tt>Octopus</tt> can run in 1, 2 or 3 dimensions, depending on the value of this
      !% variable. Note that not all input variables may be available in all cases.
      !%End
      call parse_integer(datasets_check('Dimensions'), default_ndim, this%dim)
    end if
    if((this%dim>MAX_DIM).or.(this%dim<1)) call input_error('Dimensions')
    return
  end subroutine space_init_simple

  ! ---------------------------------------------------------
  subroutine space_init_data_object(this, config)
    type(space_t),       intent(out) :: this
    type(json_object_t), intent(in)  :: config
    !
    integer :: ndim, ierr
    !
    call json_get(config, "dimensions", ndim, ierr=ierr)
    if(ierr/=JSON_OK)ndim=default_ndim
    call space_init_simple(this, ndim)
    return
  end subroutine space_init_data_object

  ! ---------------------------------------------------------
  subroutine space_create_data_object(this, config)
    type(space_t),       intent(in)  :: this
    type(json_object_t), intent(out) :: config
    !
    call json_init(config)
    call json_set(config, "dimensions", this%dim)
    return
  end subroutine space_create_data_object

  ! ---------------------------------------------------------
  elemental subroutine space_copy(this_out, this_in)
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
  elemental subroutine space_end(this)
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


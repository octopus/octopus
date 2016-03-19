!! Copyright (C) 2015 X. Andrade
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
!! $Id$

#include "global.h"

module element_oct_m
  use global_oct_m
  use io_oct_m
  use loct_oct_m
  use messages_oct_m
  use profiling_oct_m
  
  implicit none

  private
  public ::                    &
    element_t,                 &
    element_init,              &
    element_end,               &
    element_mass,              &
    element_vdw_radius,        &
    element_valid,             &
    element_atomic_number

  type element_t
    logical          :: valid
    integer          :: atomic_number
    character(len=3) :: symbol
    FLOAT            :: mass
    FLOAT            :: vdw_radius
  end type element_t
  
contains

  ! -----------------------------------

  subroutine element_init(this, label)
    type(element_t),   intent(out)   :: this
    character(len=*),  intent(in)    :: label

    integer :: iunit, nelement, ii, ilend
    character(len=MAX_PATH_LEN) :: fname

    PUSH_SUB(element_init)
    
    do ilend = 1, len(label)
      if( iachar(label(ilend:ilend)) >= iachar('a') .and. iachar(label(ilend:ilend)) <= iachar('z') ) cycle
      if( iachar(label(ilend:ilend)) >= iachar('A') .and. iachar(label(ilend:ilend)) <= iachar('Z') ) cycle
      exit
    end do
    ilend = ilend - 1
          
    this%valid = .false.

    fname = trim(conf%share)//'/pseudopotentials/elements'
    
    nelement = max(0, loct_number_of_lines(fname) - 3) 

    iunit = io_open(trim(conf%share)//'/pseudopotentials/elements', action='read', status='old', die=.false.)

    ! skip comment lines
    read(iunit, *)
    read(iunit, *)
    read(iunit, *)

    do ii = 1, nelement
      read(iunit, *) this%symbol, this%atomic_number, this%mass, this%vdw_radius
      if(trim(this%symbol) == label(1:ilend)) then
        this%valid = .true.
        exit
      end if
    end do

    call io_close(iunit)

    POP_SUB(element_init)
  end subroutine element_init
  
  ! -----------------------------------

  subroutine element_end(this)
    type(element_t),   intent(inout) :: this

    PUSH_SUB(element_end)
    
    this%valid = .false.

    POP_SUB(element_end)
  end subroutine element_end

  ! ------------------------------------

  pure FLOAT function element_mass(this) result(mass)
    type(element_t),   intent(in)    :: this

    mass = this%mass
  end function element_mass

  ! ------------------------------------

  pure FLOAT function element_vdw_radius(this) result(vdw_radius)
    type(element_t),   intent(in)    :: this

    vdw_radius = this%vdw_radius
  end function element_vdw_radius

  ! ------------------------------------

  pure logical function element_valid(this) result(valid)
    type(element_t),   intent(in)    :: this

    valid = this%valid
  end function element_valid

  ! ------------------------------------

  pure integer function element_atomic_number(this) result(atomic_number)
    type(element_t),   intent(in)    :: this

    atomic_number = this%atomic_number
  end function element_atomic_number
  
end module element_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

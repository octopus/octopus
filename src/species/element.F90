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

module element_m
  use global_m
  use io_m
  use loct_m
  use messages_m
  use profiling_m
  
  implicit none

  private
  public ::             &
    element_t

  type element_t
    private
    logical          :: valid_
    integer          :: atomic_number
    character(len=3) :: symbol
    FLOAT            :: mass_
  contains
    procedure        :: init
    procedure        :: end
    procedure        :: valid
    procedure        :: mass
  end type element_t
  
contains

  ! -----------------------------------

  subroutine init(this, label)
    class(element_t),   intent(out)   :: this
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
          
    this%valid_ = .false.

    fname = trim(conf%share)//'/pseudopotentials/elements'
    
    nelement = max(0, loct_number_of_lines(fname) - 2) 

    iunit = io_open(trim(conf%share)//'/pseudopotentials/elements', action='read', status='old', die=.false.)

    ! skip comment lines
    read(iunit, *)
    read(iunit, *)

    do ii = 1, nelement
      read(iunit, *) this%symbol, this%atomic_number, this%mass_
      if(trim(this%symbol) == label(1:ilend)) then
        this%valid_ = .true.
        exit
      end if
    end do
    
    call io_close(iunit)

    POP_SUB(element_init)
  end subroutine init
  
  ! -----------------------------------

  subroutine end(this)
    class(element_t),   intent(inout) :: this

    PUSH_SUB(element_end)
    
    this%valid_ = .false.

    POP_SUB(element_end)
  end subroutine end

  ! ------------------------------------

  pure FLOAT function mass(this)
    class(element_t),   intent(in)    :: this

    mass = this%mass_
  end function mass

  ! ------------------------------------

  pure logical function valid(this)
    class(element_t),   intent(in)    :: this

    valid = this%valid_
  end function valid
  
end module element_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

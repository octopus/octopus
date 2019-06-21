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

#include "global.h"

module element_oct_m
  
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
    integer(8) :: dummy
  end type element_t

  interface
    
    ! -------------------------------------------------
    
    subroutine element_init(self, symbol)
      import :: element_t
      implicit none
      
      type(element_t),  intent(out)   :: self
      character(len=*), intent(in)    :: symbol
    end subroutine element_init

    ! -------------------------------------------------
    
    subroutine element_end(self)
      import :: element_t
      implicit none

      type(element_t),   intent(inout) :: self
    end subroutine element_end

    ! ------------------------------------

    real(8) function element_mass(self)
      import :: element_t
      implicit none
      
      type(element_t),   intent(in)    :: self
    end function element_mass

    ! ------------------------------------

    real(8) function element_vdw_radius(self)
      import :: element_t
      implicit none
      
      type(element_t),   intent(in)    :: self
    end function element_vdw_radius
    
    ! ------------------------------------

    integer function element_atomic_number(self)
      import :: element_t
      implicit none
      
      type(element_t),   intent(in)    :: self
    end function element_atomic_number

  end interface
  

contains

  logical function element_valid(self) result(valid)
    type(element_t),   intent(in)    :: self

    interface
      integer function element_valid_low(self)
        import :: element_t
        implicit none
        
        type(element_t),   intent(in)    :: self
      end function element_valid_low
    end interface
    
    valid = element_valid_low(self) /= 0
  end function element_valid

end module element_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

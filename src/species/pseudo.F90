!! Copyright (C) 2018 X. Andrade
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

module pseudo_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::                     &
    pseudo_t,                   &
    pseudo_init,                &
    pseudo_end,                 &
    pseudo_type,                &
    pseudo_valence_charge,      &
    pseudo_mesh_spacing,        &
    pseudo_mass,                &
    pseudo_lmax,                &
    pseudo_llocal

  !these values have to match with those on base.hpp
  integer, parameter, public ::               &
    PSEUDO_TYPE_ULTRASOFT         = 30,       &
    PSEUDO_TYPE_NORM_CONSERVING   = 31,       &
    PSEUDO_TYPE_KLEINMAN_BYLANDER = 32
  
  type pseudo_t
    private
    integer(8) :: dummy
  end type pseudo_t

  interface
    subroutine pseudo_init(pseudo, filename, ierr)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(out)   :: pseudo
      character(len=*), intent(in)    :: filename
      integer,          intent(out)   :: ierr
    end subroutine pseudo_init

    subroutine pseudo_end(pseudo)
      import :: pseudo_t
      implicit none

      type(pseudo_t),   intent(out)   :: pseudo
    end subroutine pseudo_end

    integer function pseudo_type(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(out)   :: pseudo
    end function pseudo_type
    
    real(8) function pseudo_valence_charge(pseudo)
      import :: pseudo_t
      implicit none

      type(pseudo_t),   intent(out)   :: pseudo
    end function pseudo_valence_charge

    real(8) function pseudo_mesh_spacing(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(out)   :: pseudo
    end function pseudo_mesh_spacing
    
    real(8) function pseudo_mass(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(out)   :: pseudo
    end function pseudo_mass

    integer function pseudo_lmax(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(out)   :: pseudo
    end function pseudo_lmax

    integer function pseudo_llocal(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(out)   :: pseudo
    end function pseudo_llocal
    
  end interface
  
contains


  
end module pseudo_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

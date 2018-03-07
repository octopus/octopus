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
    pseudo_mesh_size,           &
    pseudo_mass,                &
    pseudo_lmax,                &
    pseudo_llocal,              &
    pseudo_nchannels

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

      type(pseudo_t),   intent(inout) :: pseudo
    end subroutine pseudo_end

    integer function pseudo_type(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_type
    
    real(8) function pseudo_valence_charge(pseudo)
      import :: pseudo_t
      implicit none

      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_valence_charge

    real(8) function pseudo_mesh_spacing(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_mesh_spacing

    integer function pseudo_mesh_size(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_mesh_size
    
    real(8) function pseudo_mass(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_mass

    integer function pseudo_lmax(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_lmax

    integer function pseudo_llocal(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_llocal

    integer function pseudo_nchannels(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_nchannels

    subroutine pseudo_local_potential(pseudo, local_potential)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      real(8),          intent(in)    :: local_potential(:)
    end subroutine pseudo_local_potential

    subroutine pseudo_projector(pseudo, l, ic, projector)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      integer,          intent(in)    :: l
      integer,          intent(in)    :: ic
      real(8),          intent(in)    :: projector(:)
    end subroutine pseudo_projector

    real(8) function pseudo_dij(pseudo, l, ic)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      integer,          intent(in)    :: l
      integer,          intent(in)    :: ic
    end function pseudo_dij

    subroutine pseudo_radial_potential(pseudo, l, ic, radial_potential)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      integer,          intent(in)    :: l
      integer,          intent(in)    :: ic
      real(8),          intent(in)    :: radial_potential(:)
    end subroutine pseudo_radial_potential

    subroutine pseudo_radial_function(pseudo, l, ic, radial_function)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      integer,          intent(in)    :: l
      integer,          intent(in)    :: ic
      real(8),          intent(in)    :: radial_function(:)
    end subroutine pseudo_radial_function

  end interface
  
contains

  logical function pseudo_has_projectors(pseudo, l)
    type(pseudo_t),   intent(in)      :: pseudo
    integer,          intent(in)      :: l
    
    interface
      integer function pseudo_has_projectors_low(pseudo, l)
        import :: pseudo_t
        implicit none
        
        type(pseudo_t),   intent(in)      :: pseudo
        integer,          intent(in)      :: l
      end function pseudo_has_projectors_low
    end interface

    pseudo_has_projectors = (pseudo_has_projectors_low(pseudo, l) /= 0)
    
  end function pseudo_has_projectors

end module pseudo_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

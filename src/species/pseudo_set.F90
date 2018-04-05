!! Copyright (C) 2018 X. Andrade
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the Lesser GNU General Public License as published by
!! the Free Software Foundation; either version 3, or (at your option)
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

module pseudo_set_oct_m
  use element_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::                                  &
    pseudo_set_t,                            &
    pseudo_set_init,                         &
    pseudo_set_end,                          &
    pseudo_set_has,                          &
    pseudo_set_file_path,                    &
    pseudo_set_lmax,                         &
    pseudo_set_llocal,                       &
    pseudo_set_spacing,                      &
    pseudo_set_radius
  
  type pseudo_set_t
    private
    integer(8) :: dummy
  end type pseudo_set_t

  interface
    
    ! -------------------------------------------------
    
    subroutine pseudo_set_init(pseudo_set, filename, ierr)
      import :: pseudo_set_t
      implicit none
      
      type(pseudo_set_t), intent(out)   :: pseudo_set
      character(len=*),   intent(in)    :: filename
      integer,            intent(out)   :: ierr
    end subroutine pseudo_set_init

    ! -------------------------------------------------
    
    subroutine pseudo_set_end(pseudo_set)
      import :: pseudo_set_t
      implicit none

      type(pseudo_set_t),   intent(inout) :: pseudo_set
    end subroutine pseudo_set_end

    ! -------------------------------------------------

    integer function pseudo_set_lmax(pseudo_set, element)
      use element_oct_m
      import :: pseudo_set_t
      implicit none
      
      type(pseudo_set_t), intent(in)    :: pseudo_set
      type(element_t),    intent(in)    :: element
    end function pseudo_set_lmax
    
    ! -------------------------------------------------

    integer function pseudo_set_llocal(pseudo_set, element)
      use element_oct_m
      import :: pseudo_set_t
      implicit none
      
      type(pseudo_set_t), intent(in)    :: pseudo_set
      type(element_t),    intent(in)    :: element
    end function pseudo_set_llocal

    ! -------------------------------------------------

    real(8) function pseudo_set_spacing(pseudo_set, element)
      use element_oct_m
      import :: pseudo_set_t
      implicit none
      
      type(pseudo_set_t), intent(in)    :: pseudo_set
      type(element_t),    intent(in)    :: element
    end function pseudo_set_spacing
    
    ! -------------------------------------------------

    real(8) function pseudo_set_radius(pseudo_set, element)
      use element_oct_m
      import :: pseudo_set_t
      implicit none
      
      type(pseudo_set_t), intent(in)    :: pseudo_set
      type(element_t),    intent(in)    :: element
    end function pseudo_set_radius

  end interface


contains

  logical function pseudo_set_has(pseudo_set, element)
    type(pseudo_set_t), intent(in)    :: pseudo_set
    type(element_t),    intent(in)    :: element
    
    interface
      integer function pseudo_set_has_low(pseudo_set, element)
        use element_oct_m
        import :: pseudo_set_t
        implicit none
      
        type(pseudo_set_t), intent(in)    :: pseudo_set
        type(element_t),    intent(in)    :: element
      end function pseudo_set_has_low
    end interface

    pseudo_set_has = pseudo_set_has_low(pseudo_set, element) /= 0
    
  end function pseudo_set_has

  ! -----------------------------------------------------------------
  
  character(len=MAX_PATH_LEN) function pseudo_set_file_path(pseudo_set, element)
    type(pseudo_set_t), intent(in)    :: pseudo_set
    type(element_t),    intent(in)    :: element
    
    interface
      subroutine pseudo_set_file_path_low(pseudo_set, element, path)
        use element_oct_m
        import :: pseudo_set_t
        implicit none
      
        type(pseudo_set_t), intent(in)    :: pseudo_set
        type(element_t),    intent(in)    :: element
        character(len=*),   intent(out)   :: path
      end subroutine pseudo_set_file_path_low
    end interface

    call pseudo_set_file_path_low(pseudo_set, element, pseudo_set_file_path)
    
  end function pseudo_set_file_path
  
end module pseudo_set_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

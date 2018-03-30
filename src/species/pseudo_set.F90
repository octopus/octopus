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
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::                                  &
    pseudo_set_t,                            &
    pseudo_set_init,                         &
    pseudo_set_end
  
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
    
  end interface
  
end module pseudo_set_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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
!! $Id: global.F90 3587 2007-11-22 16:43:00Z xavier $

#include "global.h"

module types_m

  implicit none
  
  private

  type, public :: c_pointer_t
    private
    integer, pointer :: p
  end type c_pointer_t

  type, public :: block_t
    private
    integer, pointer :: p
  end type block_t

end module types_m


module c_pointer_m

  use types_m

  implicit none 

  private

  public :: &
       c_pointer_t,      &
       set_null,         &
       is_null

  interface
    subroutine set_null(ptr)
      use types_m
      type(c_pointer_t), intent(out) :: ptr
    end subroutine set_null
  end interface

  interface
    function is_null_int(ptr)
      use types_m
      type(c_pointer_t), intent(in) :: ptr
      integer :: is_null_int
    end function is_null_int
  end interface
  
contains

  logical function is_null(ptr)
    type(c_pointer_t), intent(in) :: ptr
    
    is_null = (is_null_int(ptr) /= 0)
  end function is_null
  
end module c_pointer_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

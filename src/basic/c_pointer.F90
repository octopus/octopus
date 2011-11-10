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

#ifndef F2003_C_PTR 
module c_pointer_types_m 

  implicit none 
  
  type, public :: c_ptr 
    private 
    integer, pointer :: p=>null()
  end type c_ptr
  
end module c_pointer_types_m
#endif

module c_pointer_m
#ifdef F2003_C_PTR
  use iso_c_binding
#else
  use c_pointer_types_m
#endif

  ! This module must be public, because the Sun compiler cannot
  ! declare c_ptr as public.

  implicit none 

#ifndef F2003_C_PTR
  interface
    subroutine set_null(ptr)
      use c_pointer_types_m
      type(c_ptr), intent(out) :: ptr
    end subroutine set_null
  end interface

  interface
    function is_null_int(ptr)
      use c_pointer_types_m
      type(c_ptr), intent(in) :: ptr
      integer :: is_null_int
    end function is_null_int
  end interface
#endif
  
contains

#ifndef F2003_C_PTR
  logical function c_associated(ptr)
    type(c_ptr), intent(in) :: ptr
    c_associated = (is_null_int(ptr) == 0)
  end function c_associated
#else
  subroutine set_null(ptr)
      type(c_ptr), intent(out) :: ptr
      ptr = c_null_ptr
  end subroutine set_null
#endif
end module c_pointer_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

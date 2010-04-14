!! Copyright (C) 2010 X. Andrade
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
!! $Id: opencl_f.F90 3587 2007-11-22 16:43:00Z xavier $

#include "global.h"

module opencl_m
  use c_pointer_m

  implicit none 

  private

  public ::             &
    opencl_init,        &
    opencl_end

  type(c_ptr), public :: opencl

  interface
    subroutine opencl_init(this)
      use c_pointer_m

      type(c_ptr), intent(out) :: this
    end subroutine opencl_init

    subroutine opencl_end(this)
      use c_pointer_m

      type(c_ptr), intent(inout) :: this
    end subroutine opencl_end
  end interface

end module opencl_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

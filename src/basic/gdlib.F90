!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

module gdlib_oct_m
  use iso_c_binding
  use string_oct_m
  implicit none

  private
#if defined(HAVE_GDLIB)
  public ::                    &
    gdlib_image_create_from,   &
    gdlib_image_sx,            &
    gdlib_image_sy,            &
    gdlib_image_get_pixel_rgb, &
    gdlib_imagedestroy

  interface
    integer(c_int) function gdlib_image_sx(im) bind(c)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in) :: im
    end function gdlib_image_sx

    integer(c_int) function gdlib_image_sy(im) bind(c)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(in) :: im
    end function gdlib_image_sy

    subroutine gdlib_image_get_pixel_rgb(im, x, y, r, g, b) bind(c)
      use iso_c_binding
      implicit none
      type(c_ptr),    intent(in)  :: im
      integer(c_int), intent(in)  :: x, y
      integer(c_int), intent(out) :: r, g, b
    end subroutine gdlib_image_get_pixel_rgb

    subroutine gdlib_imagedestroy(im) bind(c, name="gdImageDestroy")
      use iso_c_binding
      implicit none
      type(c_ptr), value :: im
    end subroutine gdlib_imagedestroy
  end interface
#endif

contains

#if defined(HAVE_GDLIB)
  type(c_ptr) function gdlib_image_create_from(filename)
    character(len=*), intent(in) :: filename
    interface
      type(c_ptr) function cgdlib_image_create_from(filename) bind(c, name="gdlib_image_create_from")
        use iso_c_binding
        implicit none
        character(kind=c_char), intent(in) :: filename(*)
      end function cgdlib_image_create_from
    end interface

    gdlib_image_create_from = cgdlib_image_create_from(string_f_to_c(filename))

  end function gdlib_image_create_from
#endif

end module gdlib_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

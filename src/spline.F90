!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

module lib_oct_gsl_spline

  implicit none

  ! Define the which routines can be seen from the outside
  private
  public :: loct_spline_type, loct_spline_init, loct_spline_end, loct_spline_fit, loct_splint

  ! the basic spline datatype
  type loct_spline_type
    integer(POINTER_SIZE) :: spl, acc
  end type loct_spline_type

  integer, parameter :: ntbmax = 200

  ! some oct interfaces
  interface
    subroutine oct_spline_end(spl, acc)
      integer(POINTER_SIZE), intent(inout) :: spl, acc
    end subroutine oct_spline_end

    subroutine oct_spline_fit(nrc, x, y, spl, acc)
      integer, intent(in) :: nrc
      real(8), intent(in) :: x, y
      integer(POINTER_SIZE), intent(inout) :: spl, acc
    end subroutine oct_spline_fit
    
    real(8) function oct_spline_eval(x, spl, acc)
      real(8), intent(in) :: x
      integer(POINTER_SIZE), intent(in) :: spl, acc
    end function oct_spline_eval
  end interface

  interface loct_spline_fit
    module procedure spline_fit4
    module procedure spline_fit8
  end interface

  interface loct_splint
    module procedure splint4
    module procedure splint8
  end interface

contains

  subroutine loct_spline_init(spl)
    type(loct_spline_type), intent(out) :: spl
    
    spl%spl = 0; spl%acc = 0
  end subroutine loct_spline_init

  subroutine loct_spline_end(spl)
    type(loct_spline_type), intent(inout) :: spl
    
    if(spl%spl.ne.0 .and. spl%acc.ne.0) then
      call oct_spline_end(spl%spl, spl%acc)
    end if
    spl%spl = 0; spl%acc = 0
  end subroutine loct_spline_end

  subroutine spline_fit8(nrc, rofi, ffit, spl)
    integer, intent(in) :: nrc
    real(8), intent(in) :: ffit(nrc), rofi(nrc)
    type(loct_spline_type), intent(out) :: spl
    
    call oct_spline_fit(nrc, rofi(1), ffit(1), spl%spl, spl%acc)
  end subroutine spline_fit8

  subroutine spline_fit4(nrc, rofi, ffit, spl)
    integer, intent(in) :: nrc
    real(4), intent(IN) :: rofi(nrc), ffit(nrc)
    type(loct_spline_type), intent(out) :: spl
    
    real(8), allocatable :: rofi8(:), ffit8(:)

    allocate(rofi8(nrc), ffit8(nrc))
    rofi8 = real(rofi, kind=8)
    ffit8 = real(ffit, kind=8)
    call oct_spline_fit(nrc, rofi8(1), ffit8(1), spl%spl, spl%acc)
    deallocate(rofi8, ffit8)
  end subroutine spline_fit4

  real(8) function splint8(spl, x)
    type(loct_spline_type), intent(in) :: spl
    real(8), intent(in) :: x
  
    splint8 = oct_spline_eval(x, spl%spl, spl%acc)
  end function splint8

  real(4) function splint4(spl, x)
    type(loct_spline_type), intent(in) :: spl
    real(4), intent(in) :: x
  
    splint4 = real(oct_spline_eval(real(x, kind=8), spl%spl, spl%acc), kind=4)
  end function splint4

end module lib_oct_gsl_spline

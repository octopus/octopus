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

#include "config_F90.h"

module spline
use global

implicit none

type spline_type
  integer(POINTER_SIZE) :: spl, acc
end type spline_type

integer, parameter :: ntbmax = 200

interface
  subroutine oct_spline_init(nrc, spl, acc)
    integer, intent(in) :: nrc
    integer(POINTER_SIZE), intent(out) :: spl, acc
  end subroutine oct_spline_init
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

contains

subroutine spline_init(spl)
  type(spline_type), intent(out) :: spl

  spl%spl = 0; spl%acc = 0
end subroutine spline_init

subroutine spline_end(spl)
  type(spline_type), intent(inout) :: spl

  if(spl%spl.ne.0 .and. spl%acc.ne.0) then
    call oct_spline_end(spl%spl, spl%acc)
  end if
  spl%spl = 0; spl%acc = 0
end subroutine spline_end

subroutine spline_fit(nrc, rofi, ffit, spl)
  integer, intent(in) :: nrc
  real(r8), intent(IN) :: ffit(nrc), rofi(nrc)
  type(spline_type), intent(out) :: spl

  call oct_spline_fit(nrc, rofi(1), ffit(1), spl%spl, spl%acc)

  return
end subroutine spline_fit

function splint(spl, x)
  real(r8) :: splint

  type(spline_type), intent(IN) :: spl
  real(r8), intent(IN) :: x
  
  splint = oct_spline_eval(x, spl%spl, spl%acc)

  return
end function splint

end module spline

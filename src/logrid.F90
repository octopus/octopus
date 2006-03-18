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
!!
!! $Id$

#include "global.h"

module logrid_m

  use global_m
  use messages_m

  implicit none

  type logrid_t
    FLOAT :: a,b
    integer  :: nrval
    FLOAT, pointer :: rofi(:), drdi(:), s(:)
  end type logrid_t

contains

  ! ---------------------------------------------------------
  subroutine logrid_init(g, a, b, nrval)
    type(logrid_t), intent(out) :: g
    FLOAT, intent(in)           :: a,b
    integer, intent(in)         :: nrval

    FLOAT :: rpb, ea
    integer  :: ir

    g%a = a; g%b = b; g%nrval = nrval

    ALLOCATE(g%rofi(nrval), nrval)
    ALLOCATE(g%drdi(nrval), nrval)
    ALLOCATE(g%s(nrval),    nrval)

    rpb = b; ea = exp(a)
    do ir = 1, nrval
      g%drdi(ir) = a*rpb
      g%s(ir) = sqrt(a*rpb)
      rpb = rpb*ea
      g%rofi(ir) = b * ( exp( a * (ir - 1) ) - M_ONE )
    end do

  end subroutine logrid_init


  ! ---------------------------------------------------------
  subroutine logrid_end(g)
    type(logrid_t), intent(inout) :: g
    deallocate(g%rofi, g%drdi, g%s)
  end subroutine logrid_end


  ! ---------------------------------------------------------
  subroutine derivate_in_log_grid(g, f, dfdr)
    type(logrid_t)     :: g
    FLOAT, intent(in)  :: f(g%nrval)
    FLOAT, intent(out) :: dfdr(g%nrval)

    FLOAT :: x, y, a, b
    integer :: i, nrval

    a = g%a; b = g%b; nrval = g%nrval

    x = M_ONE - exp(-2*a)
    y = M_ONE - exp(-a)

    dfdr(1) = (1/(y*b))*exp(-a)*(f(2)-f(1))
    do i = 2, nrval-1
      dfdr(i) = (1/(x*b))*exp(-i*a)*(f(i+1)-f(i-1))
    end do
    dfdr(nrval) = (1/(y*b))*exp(-(nrval-1)*a)*(f(nrval)-f(nrval-1))

    return
  end subroutine derivate_in_log_grid

end module logrid_m

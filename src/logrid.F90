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

module logrid

use global

implicit none

type logrid_type
  real(r8) :: a,b
  integer  :: nrval
  real(r8), pointer :: rofi(:), drdi(:), s(:)
end type

contains

subroutine init_logrid(g, a, b, nrval)
  type(logrid_type), intent(out) :: g
  real(r8), intent(in)           :: a,b
  integer, intent(in)            :: nrval

  real(r8) :: rpb, ea
  integer  :: ir

  g%a = a; g%b = b; g%nrval = nrval

  allocate(g%rofi(nrval), g%drdi(nrval), g%s(nrval))
  
  rpb = b; ea = exp(a)
  do ir = 1, nrval
    g%drdi(ir) = a*rpb
    g%s(ir) = sqrt(a*rpb)
    rpb = rpb*ea
    g%rofi(ir) = b * ( exp( a * (ir - 1) ) - 1.0_r8 )
  end do

end subroutine init_logrid

subroutine kill_logrid(g)
  type(logrid_type), intent(inout) :: g
  deallocate(g%rofi, g%drdi, g%s)
end subroutine kill_logrid

subroutine derivate_in_log_grid(g, f, dfdr)
  type(logrid_type)     :: g
  real(r8), intent(IN)  :: f(g%nrval)
  real(r8), intent(out) :: dfdr(g%nrval)

  real(r8) :: x, y, a, b
  integer :: i, nrval

  a = g%a; b = g%b; nrval = g%nrval

  x = 1.0_r8 - exp(-2*a)
  y = 1.0_r8 - exp(-a)

  dfdr(1) = (1/(y*b))*exp(-a)*(f(2)-f(1))  
  do i = 2, nrval-1
    dfdr(i) = (1/(x*b))*exp(-i*a)*(f(i+1)-f(i-1))
  enddo
  dfdr(nrval) = (1/(y*b))*exp(-(nrval-1)*a)*(f(nrval)-f(nrval-1))

  return
end subroutine derivate_in_log_grid

end module logrid

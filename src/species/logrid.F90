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
!! $Id$

#include "global.h"

module logrid_m
  use global_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::         &
    logrid_t,       &
    logrid_init,    &
    logrid_end,     &
    logrid_copy,    &
    logrid_index,   &
    derivate_in_log_grid

  integer, parameter, public :: &
    LOGRID_PSF  = 1, & ! log grid used in Troullier-Martins code
    LOGRID_CPI  = 2    ! log grid used in FHI code

  type logrid_t
    integer  :: flavor

    FLOAT    :: a, b
    integer  :: nrval

    FLOAT, pointer :: rofi(:) ! r value of the point i
    FLOAT, pointer :: r2ofi(:) ! r value of the point i
    FLOAT, pointer :: drdi(:) ! jacobian, i.e., the derivative of r in terms of i
    FLOAT, pointer :: s(:)    ! sqrt of drdi
  end type logrid_t

contains

  ! ---------------------------------------------------------
  subroutine logrid_init(g, flavor, a, b, nrval)
    type(logrid_t), intent(out) :: g
    integer,        intent(in)  :: flavor
    FLOAT,          intent(in)  :: a, b
    integer,        intent(in)  :: nrval

    FLOAT :: rpb, ea
    integer  :: ir

    ASSERT(flavor==LOGRID_PSF.or.flavor==LOGRID_CPI)

    g%flavor = flavor
    g%a = a; g%b = b; g%nrval = nrval

    ALLOCATE(g%rofi(nrval), nrval)
    ALLOCATE(g%r2ofi(nrval), nrval)
    ALLOCATE(g%drdi(nrval), nrval)
    ALLOCATE(g%s(nrval),    nrval)

    select case(g%flavor)
    case(LOGRID_PSF)
      rpb = b
      ea  = exp(a)
      do ir = 1, nrval
        g%drdi(ir) = a*rpb
        rpb        = rpb*ea
        g%rofi(ir) = b*(exp(a*(ir-1)) - M_ONE)
      end do

    case(LOGRID_CPI)
      g%rofi(1) = M_ZERO
      g%drdi(1) = M_ZERO

      rpb = log(a)
      g%rofi(2) = b
      g%drdi(2) = b*rpb
      do ir = 3, g%nrval
        g%rofi(ir) = g%rofi(ir-1)*a
        g%drdi(ir) = g%rofi(ir)*rpb
      end do
    end select

    ! calculate sqrt(drdi)
    do ir = 1, g%nrval
      g%s(ir)    = sqrt(g%drdi(ir))
      g%r2ofi(ir) = g%rofi(ir)**2
    end do


  end subroutine logrid_init


  ! ---------------------------------------------------------
  subroutine logrid_end(g)
    type(logrid_t), intent(inout) :: g

    SAFE_DEALLOCATE_P(g%rofi)
    SAFE_DEALLOCATE_P(g%r2ofi)
    SAFE_DEALLOCATE_P(g%drdi)
    SAFE_DEALLOCATE_P(g%s)

  end subroutine logrid_end


  ! ---------------------------------------------------------
  subroutine logrid_copy(gi, go)
    type(logrid_t), intent(in)  :: gi
    type(logrid_t), intent(out) :: go

    go%flavor = gi%flavor
    go%a      = gi%a
    go%b      = gi%b
    go%nrval  = gi%nrval

    ALLOCATE(go%rofi(go%nrval), go%nrval)
    ALLOCATE(go%r2ofi(go%nrval), go%nrval)
    ALLOCATE(go%drdi(go%nrval), go%nrval)
    ALLOCATE(go%s(go%nrval), go%nrval)

    go%rofi(:) = gi%rofi(:)
    go%r2ofi(:) = gi%r2ofi(:)
    go%drdi(:) = gi%drdi(:)
    go%s(:)    = gi%s(:)

  end subroutine logrid_copy


  ! ---------------------------------------------------------
  integer function logrid_index(g, rofi) result(ii)
    type(logrid_t), intent(in) :: g
    FLOAT,          intent(in) :: rofi

    integer :: ir

    ii = 0
    do ir = 1, g%nrval-1

      if(rofi >= g%rofi(ir).and.rofi < g%rofi(ir+1)) then
        if(abs(rofi-g%rofi(ir)) < abs(rofi-g%rofi(ir+1))) then
          ii = ir
        else
          ii = ir + 1
        end if
        exit
      end if

    end do
  end function logrid_index


  ! ---------------------------------------------------------
  subroutine derivate_in_log_grid(g, f, dfdr)
    type(logrid_t), intent(in)   :: g
    FLOAT,          intent(in)   :: f(:)
    FLOAT,          intent(out)  :: dfdr(:)

    integer :: i

    dfdr(1) = (f(2) - f(1))/(g%rofi(2) - g%rofi(1))
    do i = 2, g%nrval-1
      dfdr(i) = (f(i+1) - f(i-1))/(g%rofi(i+1) - g%rofi(i-1))
    end do
    dfdr(g%nrval) = (f(g%nrval) - f(g%nrval-1))/(g%rofi(g%nrval) - g%rofi(g%nrval-1))

  end subroutine derivate_in_log_grid

end module logrid_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

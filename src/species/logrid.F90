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
    logrid_radius,  &
    derivative_in_log_grid

  integer, parameter, public :: &
    LOGRID_PSF  = 1, & !< log grid used in Troullier-Martins code
    LOGRID_CPI  = 2    !< log grid used in FHI code

  type logrid_t
    integer  :: flavor

    FLOAT    :: a, b
    integer  :: nrval

    FLOAT, pointer :: rofi(:)  !< r value of the point i
    FLOAT, pointer :: r2ofi(:) !< r value of the point i
    FLOAT, pointer :: drdi(:)  !< Jacobian, i.e., the derivative of r in terms of i
    FLOAT, pointer :: s(:)     !< sqrt of drdi
  end type logrid_t

contains

  ! ---------------------------------------------------------
  subroutine logrid_init(grid, flavor, aa, bb, nrval)
    type(logrid_t), intent(out) :: grid
    integer,        intent(in)  :: flavor
    FLOAT,          intent(in)  :: aa, bb
    integer,        intent(in)  :: nrval

    FLOAT :: rpb, ea
    integer  :: ir

    PUSH_SUB(logrid_init)

    ASSERT(flavor==LOGRID_PSF.or.flavor==LOGRID_CPI)

    grid%flavor = flavor
    grid%a = aa
    grid%b = bb
    grid%nrval = nrval

    SAFE_ALLOCATE(grid%rofi(1:nrval))
    SAFE_ALLOCATE(grid%r2ofi(1:nrval))
    SAFE_ALLOCATE(grid%drdi(1:nrval))
    SAFE_ALLOCATE(grid%s(1:nrval))

    select case(grid%flavor)
    case(LOGRID_PSF)
      rpb = bb
      ea  = exp(aa)
      do ir = 1, nrval
        grid%drdi(ir) = aa*rpb
        rpb           = rpb*ea
        grid%rofi(ir) = bb*(exp(aa*(ir-1)) - M_ONE)
      end do

    case(LOGRID_CPI)
      grid%rofi(1) = M_ZERO
      grid%drdi(1) = M_ZERO

      rpb = log(aa)
      grid%rofi(2) = bb
      grid%drdi(2) = bb*rpb
      do ir = 3, grid%nrval
        grid%rofi(ir) = grid%rofi(ir-1)*aa
        grid%drdi(ir) = grid%rofi(ir)*rpb
      end do
    end select

    ! calculate sqrt(drdi)
    do ir = 1, grid%nrval
      grid%s(ir)     = sqrt(grid%drdi(ir))
      grid%r2ofi(ir) = grid%rofi(ir)**2
    end do

    POP_SUB(logrid_init)
  end subroutine logrid_init


  ! ---------------------------------------------------------
  subroutine logrid_end(grid)
    type(logrid_t), intent(inout) :: grid

    PUSH_SUB(logrid_end)

    SAFE_DEALLOCATE_P(grid%rofi)
    SAFE_DEALLOCATE_P(grid%r2ofi)
    SAFE_DEALLOCATE_P(grid%drdi)
    SAFE_DEALLOCATE_P(grid%s)

    POP_SUB(logrid_end)
  end subroutine logrid_end


  ! ---------------------------------------------------------
  subroutine logrid_copy(grid_in, grid_out)
    type(logrid_t), intent(in)  :: grid_in
    type(logrid_t), intent(out) :: grid_out

    PUSH_SUB(logrid_copy)

    grid_out%flavor = grid_in%flavor
    grid_out%a      = grid_in%a
    grid_out%b      = grid_in%b
    grid_out%nrval  = grid_in%nrval

    SAFE_ALLOCATE(grid_out%rofi (1:grid_out%nrval))
    SAFE_ALLOCATE(grid_out%r2ofi(1:grid_out%nrval))
    SAFE_ALLOCATE(grid_out%drdi (1:grid_out%nrval))
    SAFE_ALLOCATE(grid_out%s    (1:grid_out%nrval))

    grid_out%rofi(:)  = grid_in%rofi(:)
    grid_out%r2ofi(:) = grid_in%r2ofi(:)
    grid_out%drdi(:)  = grid_in%drdi(:)
    grid_out%s(:)     = grid_in%s(:)

    POP_SUB(logrid_copy)
  end subroutine logrid_copy


  ! ---------------------------------------------------------
  integer function logrid_index(grid, rofi) result(ii)
    type(logrid_t), intent(in) :: grid
    FLOAT,          intent(in) :: rofi

    integer :: ir

    PUSH_SUB(logrid_index)

    ii = 0
    do ir = 1, grid%nrval-1

      if(rofi >= grid%rofi(ir).and.rofi < grid%rofi(ir+1)) then
        if(abs(rofi-grid%rofi(ir)) < abs(rofi-grid%rofi(ir+1))) then
          ii = ir
        else
          ii = ir + 1
        end if
        exit
      end if

    end do

    POP_SUB(logrid_index)
  end function logrid_index


  ! ---------------------------------------------------------
  subroutine derivative_in_log_grid(grid, ff, dfdr)
    type(logrid_t), intent(in)   :: grid    
    FLOAT,          intent(in)   :: ff(:)
    FLOAT,          intent(out)  :: dfdr(:)

    integer :: ii

    PUSH_SUB(derivative_in_log_grid)

    dfdr(1) = (ff(2) - ff(1))/(grid%rofi(2) - grid%rofi(1))
    do ii = 2, grid%nrval-1
      dfdr(ii) = (ff(ii+1) - ff(ii-1))/(grid%rofi(ii+1) - grid%rofi(ii-1))
    end do
    dfdr(grid%nrval) = (ff(grid%nrval) - ff(grid%nrval-1))/(grid%rofi(grid%nrval) - grid%rofi(grid%nrval-1))

    POP_SUB(derivative_in_log_grid)
  end subroutine derivative_in_log_grid

  ! ----------------------------------------------------------
  FLOAT pure function logrid_radius(grid) result(radius)
    type(logrid_t), intent(in)   :: grid

    radius = grid%rofi(grid%nrval)
  end function logrid_radius

end module logrid_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

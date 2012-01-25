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

module stencil_starplus_m
  use global_m
  use messages_m
  use profiling_m
  use stencil_m

  private
  public ::                     &
    stencil_starplus_size_lapl, &
    stencil_starplus_extent,    &
    stencil_starplus_get_lapl,  &
    stencil_starplus_pol_lapl,  &
    stencil_starplus_size_grad, &
    stencil_starplus_get_grad,  &
    stencil_starplus_pol_grad

contains

  ! ---------------------------------------------------------
  integer function stencil_starplus_size_lapl(dim, order) result(n)
    integer, intent(in) :: dim
    integer, intent(in) :: order

    PUSH_SUB(stencil_starplus_size_lapl)

    n = 2*dim*order + 1
    if(dim == 2) n = n + 12
    if(dim == 3) n = n + 44

    POP_SUB(stencil_starplus_size_lapl)
  end function stencil_starplus_size_lapl


  ! ---------------------------------------------------------
  !> Returns maximum extension of the stencil in spatial direction
  !! dir = 1, 2, 3 for a given discretization order.
  integer function stencil_starplus_extent(dir, order)
    integer, intent(in) :: dir
    integer, intent(in) :: order

    integer :: extent

    PUSH_SUB(stencil_starplus_extent)

    extent = 0
    if(dir.ge.1.or.dir.le.3) then
      if(order.le.2) then
        extent = 2
      else
        extent = order
      end if
    end if
    stencil_starplus_extent = extent

    POP_SUB(stencil_starplus_extent)
  end function stencil_starplus_extent


  ! ---------------------------------------------------------
  integer function stencil_starplus_size_grad(dim, order) result(n)
    integer, intent(in) :: dim
    integer, intent(in) :: order

    PUSH_SUB(stencil_starplus_size_grad)

    n = 2*order + 1
    if(dim == 2) n = n + 2
    if(dim == 3) n = n + 4

    POP_SUB(stencil_starplus_size_grad)
  end function stencil_starplus_size_grad


  ! ---------------------------------------------------------
  subroutine stencil_starplus_get_lapl(this, dim, order)
    type(stencil_t), intent(out) :: this
    integer,         intent(in)  :: dim
    integer,         intent(in)  :: order

    integer :: i, j, n

    PUSH_SUB(stencil_starplus_get_lapl)

    call stencil_allocate(this, stencil_starplus_size_lapl(dim, order))

    n = 1
    select case(dim)
    case(1)
      n = 1
      do i = 1, dim
        do j = -order, order
          if(j == 0) cycle
          n = n + 1
          this%points(i, n) = j
        end do
      end do
    case(2)
      n = 1
      do i = 1, dim
        do j = -order, order
          if(j == 0) cycle
          n = n + 1
          this%points(i, n) = j
        end do
      end do
      n = n + 1; this%points(1:2, n) = (/ -2,  1 /)
      n = n + 1; this%points(1:2, n) = (/ -2, -1 /)
      n = n + 1; this%points(1:2, n) = (/ -1,  2 /)
      n = n + 1; this%points(1:2, n) = (/ -1,  1 /)
      n = n + 1; this%points(1:2, n) = (/ -1, -1 /)
      n = n + 1; this%points(1:2, n) = (/ -1, -2 /)
      n = n + 1; this%points(1:2, n) = (/  1,  2 /)
      n = n + 1; this%points(1:2, n) = (/  1,  1 /)
      n = n + 1; this%points(1:2, n) = (/  1, -1 /)
      n = n + 1; this%points(1:2, n) = (/  1, -2 /)
      n = n + 1; this%points(1:2, n) = (/  2,  1 /)
      n = n + 1; this%points(1:2, n) = (/  2, -1 /)
    case(3)
      n = 1
      do i = 1, dim
        do j = -order, order
          if(j == 0) cycle
          n = n + 1
          this%points(i, n) = j
        end do
      end do
      n = n + 1; this%points(1:3, n) = (/ -2,  1, 0 /)
      n = n + 1; this%points(1:3, n) = (/ -2, -1, 0 /)
      n = n + 1; this%points(1:3, n) = (/ -1,  2, 0 /)
      n = n + 1; this%points(1:3, n) = (/ -1,  1, 0 /)
      n = n + 1; this%points(1:3, n) = (/ -1, -1, 0 /)
      n = n + 1; this%points(1:3, n) = (/ -1, -2, 0 /)
      n = n + 1; this%points(1:3, n) = (/  1,  2, 0 /)
      n = n + 1; this%points(1:3, n) = (/  1,  1, 0 /)
      n = n + 1; this%points(1:3, n) = (/  1, -1, 0 /)
      n = n + 1; this%points(1:3, n) = (/  1, -2, 0 /)
      n = n + 1; this%points(1:3, n) = (/  2,  1, 0 /)
      n = n + 1; this%points(1:3, n) = (/  2, -1, 0 /)

      n = n + 1; this%points(1:3, n) = (/ -2, 0,  1 /)
      n = n + 1; this%points(1:3, n) = (/ -2, 0, -1 /)
      n = n + 1; this%points(1:3, n) = (/ -1, 0,  2 /)
      n = n + 1; this%points(1:3, n) = (/ -1, 0,  1 /)
      n = n + 1; this%points(1:3, n) = (/ -1, 0, -1 /)
      n = n + 1; this%points(1:3, n) = (/ -1, 0, -2 /)
      n = n + 1; this%points(1:3, n) = (/  1, 0,  2 /)
      n = n + 1; this%points(1:3, n) = (/  1, 0,  1 /)
      n = n + 1; this%points(1:3, n) = (/  1, 0, -1 /)
      n = n + 1; this%points(1:3, n) = (/  1, 0, -2 /)
      n = n + 1; this%points(1:3, n) = (/  2, 0,  1 /)
      n = n + 1; this%points(1:3, n) = (/  2, 0, -1 /)

      n = n + 1; this%points(1:3, n) = (/ 0, -2,  1 /)
      n = n + 1; this%points(1:3, n) = (/ 0, -2, -1 /)
      n = n + 1; this%points(1:3, n) = (/ 0, -1,  2 /)
      n = n + 1; this%points(1:3, n) = (/ 0, -1,  1 /)
      n = n + 1; this%points(1:3, n) = (/ 0, -1, -1 /)
      n = n + 1; this%points(1:3, n) = (/ 0, -1, -2 /)
      n = n + 1; this%points(1:3, n) = (/ 0,  1,  2 /)
      n = n + 1; this%points(1:3, n) = (/ 0,  1,  1 /)
      n = n + 1; this%points(1:3, n) = (/ 0,  1, -1 /)
      n = n + 1; this%points(1:3, n) = (/ 0,  1, -2 /)
      n = n + 1; this%points(1:3, n) = (/ 0,  2,  1 /)
      n = n + 1; this%points(1:3, n) = (/ 0,  2, -1 /)

      n = n + 1; this%points(1:3, n) = (/ -1, -1, -1 /)
      n = n + 1; this%points(1:3, n) = (/ -1, -1,  1 /)
      n = n + 1; this%points(1:3, n) = (/ -1,  1, -1 /)
      n = n + 1; this%points(1:3, n) = (/ -1,  1,  1 /)
      n = n + 1; this%points(1:3, n) = (/  1, -1, -1 /)
      n = n + 1; this%points(1:3, n) = (/  1, -1,  1 /)
      n = n + 1; this%points(1:3, n) = (/  1,  1, -1 /)
      n = n + 1; this%points(1:3, n) = (/  1,  1,  1 /)

    end select

    call stencil_init_center(this)

    POP_SUB(stencil_starplus_get_lapl)
  end subroutine stencil_starplus_get_lapl


  ! ---------------------------------------------------------
  subroutine stencil_starplus_get_grad(this, dim, dir, order)
    type(stencil_t), intent(out) :: this
    integer,         intent(in)  :: dim
    integer,         intent(in)  :: dir
    integer,         intent(in)  :: order

    integer :: i, n, j

    PUSH_SUB(stencil_star_get_grad)

    call stencil_allocate(this, stencil_starplus_size_grad(dim, order))

    n = 1
    do i = -order, order
      this%points(dir, n) = i
      n = n + 1
    end do
    do j = 1, dim
      if(j==dir) cycle
      this%points(j, n) = -1
      n = n + 1
      this%points(j, n) =  1
      n = n + 1
    end do

    call stencil_init_center(this)

    POP_SUB(stencil_star_get_grad)
  end subroutine stencil_starplus_get_grad


  ! ---------------------------------------------------------
  subroutine stencil_starplus_pol_lapl(dim, order, pol)
    integer, intent(in)  :: dim
    integer, intent(in)  :: order
    integer, intent(out) :: pol(:,:) !< pol(dim, order)

    integer :: i, j, n

    PUSH_SUB(stencil_starplus_pol_lapl)

    n = 1
    select case(dim)
    case(1)
      n = 1
      pol(:,:) = 0
      do i = 1, dim
        do j = 1, 2*order
          n = n + 1
          pol(i, n) = j
        end do
      end do
    case(2)
      n = 1
      pol(:,:) = 0
      do i = 1, dim
        do j = 1, 2*order
          n = n + 1
          pol(i, n) = j
        end do
      end do
      n = n + 1; pol(1:2, n) = (/ 1, 1 /)
      n = n + 1; pol(1:2, n) = (/ 1, 2 /)
      n = n + 1; pol(1:2, n) = (/ 1, 3 /)
      n = n + 1; pol(1:2, n) = (/ 1, 4 /)
      n = n + 1; pol(1:2, n) = (/ 2, 1 /)
      n = n + 1; pol(1:2, n) = (/ 2, 2 /)
      n = n + 1; pol(1:2, n) = (/ 2, 3 /)
      n = n + 1; pol(1:2, n) = (/ 2, 4 /)
      n = n + 1; pol(1:2, n) = (/ 3, 1 /)
      n = n + 1; pol(1:2, n) = (/ 3, 2 /)
      n = n + 1; pol(1:2, n) = (/ 4, 1 /)
      n = n + 1; pol(1:2, n) = (/ 4, 2 /)
    case(3)
      n = 1
      pol(:,:) = 0
      do i = 1, dim
        do j = 1, 2*order
          n = n + 1
          pol(i, n) = j
        end do
      end do

      n = n + 1; pol(1:3, n) = (/ 1, 1, 0 /)
      n = n + 1; pol(1:3, n) = (/ 1, 2, 0 /)
      n = n + 1; pol(1:3, n) = (/ 1, 3, 0 /)
      n = n + 1; pol(1:3, n) = (/ 1, 4, 0 /)
      n = n + 1; pol(1:3, n) = (/ 2, 1, 0 /)
      n = n + 1; pol(1:3, n) = (/ 2, 2, 0 /)
      n = n + 1; pol(1:3, n) = (/ 2, 3, 0 /)
      n = n + 1; pol(1:3, n) = (/ 2, 4, 0 /)
      n = n + 1; pol(1:3, n) = (/ 3, 1, 0 /)
      n = n + 1; pol(1:3, n) = (/ 3, 2, 0 /)
      n = n + 1; pol(1:3, n) = (/ 4, 1, 0 /)
      n = n + 1; pol(1:3, n) = (/ 4, 2, 0 /)

      n = n + 1; pol(1:3, n) = (/ 1, 0, 1 /)
      n = n + 1; pol(1:3, n) = (/ 1, 0, 2 /)
      n = n + 1; pol(1:3, n) = (/ 1, 0, 3 /)
      n = n + 1; pol(1:3, n) = (/ 1, 0, 4 /)
      n = n + 1; pol(1:3, n) = (/ 2, 0, 1 /)
      n = n + 1; pol(1:3, n) = (/ 2, 0, 2 /)
      n = n + 1; pol(1:3, n) = (/ 2, 0, 3 /)
      n = n + 1; pol(1:3, n) = (/ 2, 0, 4 /)
      n = n + 1; pol(1:3, n) = (/ 3, 0, 1 /)
      n = n + 1; pol(1:3, n) = (/ 3, 0, 2 /)
      n = n + 1; pol(1:3, n) = (/ 4, 0, 1 /)
      n = n + 1; pol(1:3, n) = (/ 4, 0, 2 /)

      n = n + 1; pol(1:3, n) = (/ 0, 1, 1 /)
      n = n + 1; pol(1:3, n) = (/ 0, 1, 2 /)
      n = n + 1; pol(1:3, n) = (/ 0, 1, 3 /)
      n = n + 1; pol(1:3, n) = (/ 0, 1, 4 /)
      n = n + 1; pol(1:3, n) = (/ 0, 2, 1 /)
      n = n + 1; pol(1:3, n) = (/ 0, 2, 2 /)
      n = n + 1; pol(1:3, n) = (/ 0, 2, 3 /)
      n = n + 1; pol(1:3, n) = (/ 0, 2, 4 /)
      n = n + 1; pol(1:3, n) = (/ 0, 3, 1 /)
      n = n + 1; pol(1:3, n) = (/ 0, 3, 2 /)
      n = n + 1; pol(1:3, n) = (/ 0, 4, 1 /)
      n = n + 1; pol(1:3, n) = (/ 0, 4, 2 /)

      n = n + 1; pol(1:3, n) = (/ 1, 1, 1 /)
      n = n + 1; pol(1:3, n) = (/ 1, 1, 2 /)
      n = n + 1; pol(1:3, n) = (/ 1, 2, 1 /)
      n = n + 1; pol(1:3, n) = (/ 1, 2, 2 /)
      n = n + 1; pol(1:3, n) = (/ 2, 1, 1 /)
      n = n + 1; pol(1:3, n) = (/ 2, 1, 2 /)
      n = n + 1; pol(1:3, n) = (/ 2, 2, 1 /)
      n = n + 1; pol(1:3, n) = (/ 2, 2, 2 /)

    end select

    POP_SUB(stencil_starplus_pol_lapl)
  end subroutine stencil_starplus_pol_lapl


  ! ---------------------------------------------------------
  subroutine stencil_starplus_pol_grad(dim, dir, order, pol)
    integer, intent(in)  :: dim
    integer, intent(in)  :: dir
    integer, intent(in)  :: order
    integer, intent(out) :: pol(:,:) !< pol(dim, order)

    integer :: j, n

    PUSH_SUB(stencil_starplus_pol_grad)

    pol(:,:) = 0
    do j = 0, 2*order
      pol(dir, j+1) = j
    end do
    n = 2*order + 1

    select case(dim)
    case(2)
      select case(dir)
      case(1)
        n = n + 1; pol(1:2, n) = (/ 0, 1 /)
        n = n + 1; pol(1:2, n) = (/ 0, 2 /)
      case(2)
        n = n + 1; pol(1:2, n) = (/ 1, 0 /)
        n = n + 1; pol(1:2, n) = (/ 2, 0 /)
      end select
    case(3)
      select case(dir)
      case(1)
        n = n + 1; pol(1:3, n) = (/ 0, 1, 0 /)
        n = n + 1; pol(1:3, n) = (/ 0, 2, 0 /)
        n = n + 1; pol(1:3, n) = (/ 0, 0, 1 /)
        n = n + 1; pol(1:3, n) = (/ 0, 0, 2 /)
      case(2)
        n = n + 1; pol(1:3, n) = (/ 1, 0, 0 /)
        n = n + 1; pol(1:3, n) = (/ 2, 0, 0 /)
        n = n + 1; pol(1:3, n) = (/ 0, 0, 1 /)
        n = n + 1; pol(1:3, n) = (/ 0, 0, 2 /)
      case(3)
        n = n + 1; pol(1:3, n) = (/ 1, 0, 0 /)
        n = n + 1; pol(1:3, n) = (/ 2, 0, 0 /)
        n = n + 1; pol(1:3, n) = (/ 0, 1, 0 /)
        n = n + 1; pol(1:3, n) = (/ 0, 2, 0 /)
      end select
    end select

    POP_SUB(stencil_starplus_pol_grad)
  end subroutine stencil_starplus_pol_grad

end module stencil_starplus_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

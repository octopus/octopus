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

module stencil_cube_m
  use global_m
  use messages_m
  use profiling_m
  use stencil_m

  implicit none

  private
  public ::                        &
    stencil_cube_size_lapl,        &
    stencil_cube_extent,           &
    stencil_cube_get_lapl,         &
    stencil_cube_polynomials_lapl, &
    stencil_cube_size_grad,        &
    stencil_cube_get_grad,         &
    stencil_cube_polynomials_grad


contains

  ! ---------------------------------------------------------
  integer function stencil_cube_size_lapl(dim, order)
    integer, intent(in) :: dim
    integer, intent(in) :: order

    PUSH_SUB(stencil_cube_size_lapl)

    stencil_cube_size_lapl = (2*order+1)**dim

    POP_SUB(stencil_cube_size_lapl)
  end function stencil_cube_size_lapl


  ! ---------------------------------------------------------
  !> Returns maximum extension of the stencil in spatial direction
  !! dir = 1, 2, 3 for a given discretization order.
  integer function stencil_cube_extent(dir, order)
    integer, intent(in) :: dir
    integer, intent(in) :: order

    integer :: extent

    PUSH_SUB(stencil_cube_extent)

    extent = 0
    if(dir.ge.1.or.dir.le.3) then
      extent = order
    end if
    stencil_cube_extent = extent

    POP_SUB(stencil_cube_extent)
  end function stencil_cube_extent


  ! ---------------------------------------------------------
  subroutine stencil_cube_get_lapl(this, dim, order)
    type(stencil_t), intent(out) :: this
    integer,         intent(in)  :: dim
    integer,         intent(in)  :: order

    integer :: i, j, k, n

    PUSH_SUB(stencil_cube_get_lapl)

    call stencil_allocate(this, stencil_cube_size_lapl(dim, order))

    n = 1
    select case(dim)
    case(1)
      do i = -order, order
        this%points(1, n) = i
        n = n + 1
      end do
    case(2)
      do i = -order, order
        do j = -order, order
          this%points(1, n) = i
          this%points(2, n) = j
          n = n + 1
        end do
      end do
    case(3)
      do i = -order, order
        do j = -order, order
          do k = -order, order
            this%points(1, n) = i
            this%points(2, n) = j
            this%points(3, n) = k
            n = n + 1
          end do
        end do
      end do
    end select

    call stencil_init_center(this)

    POP_SUB(stencil_cube_get_lapl)
  end subroutine stencil_cube_get_lapl


  ! ---------------------------------------------------------
  subroutine stencil_cube_polynomials_lapl(dim, order, pol)
    integer, intent(in)  :: dim
    integer, intent(in)  :: order
    integer, intent(out) :: pol(:,:) !< pol(dim, order)

    integer :: i, j, k, n

    PUSH_SUB(stencil_cube_polynomials_lapl)

    n = 1
    select case(dim)
    case(1)
      do i = 0, 2*order
        pol(1, n) = i
        n = n + 1
      end do
    case(2)
      do i = 0, 2*order
        do j = 0, 2*order
          pol(1, n) = i
          pol(2, n) = j
          n = n + 1
        end do
      end do
    case(3)
      do i = 0, 2*order
        do j = 0, 2*order
          do k = 0, 2*order
            pol(1, n) = i
            pol(2, n) = j
            pol(3, n) = k
            n = n + 1
          end do
        end do
      end do
    end select

    POP_SUB(stencil_cube_polynomials_lapl)
  end subroutine stencil_cube_polynomials_lapl


  !> Now come the gradient routines. As this stencil is the same for
  !! the laplacian and the gradient, these routines just call the
  !! corresponding ones for the laplacian

  ! ---------------------------------------------------------
  integer function stencil_cube_size_grad(dim, order)
    integer, intent(in) :: dim
    integer, intent(in) :: order

    PUSH_SUB(stencil_cube_size_grad)

    stencil_cube_size_grad = stencil_cube_size_lapl(dim, order)

    POP_SUB(stencil_cube_size_grad)
  end function stencil_cube_size_grad


  ! ---------------------------------------------------------
  subroutine stencil_cube_get_grad(this, dim, order)
    type(stencil_t), intent(out) :: this
    integer,         intent(in)  :: dim
    integer,         intent(in)  :: order

    PUSH_SUB(stencil_cube_get_grad)

    call stencil_cube_get_lapl(this, dim, order)

    POP_SUB(stencil_cube_get_grad)
  end subroutine stencil_cube_get_grad


  ! ---------------------------------------------------------
  subroutine stencil_cube_polynomials_grad(dim, order, pol)
    integer, intent(in)  :: dim
    integer, intent(in)  :: order
    integer, intent(out) :: pol(:,:) !< pol(sb%dim, order)

    PUSH_SUB(stencil_cube_polynomials_grad)

    call stencil_cube_polynomials_lapl(dim, order, pol)

    POP_SUB(stencil_cube_polynomials_grad)
  end subroutine stencil_cube_polynomials_grad

end module stencil_cube_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

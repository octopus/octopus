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

module stencil_star_m
  use global_m
  use math_m
  use messages_m
  use nl_operator_m
  use profiling_m
  use stencil_m

  implicit none

  private
  public ::                        &
    stencil_star_size_lapl,        &
    stencil_star_extent,           &
    stencil_star_get_lapl,         &
    stencil_star_polynomials_lapl, &
    stencil_star_coeff_lapl,       &
    stencil_star_size_grad,        &
    stencil_star_get_grad,         &
    stencil_star_polynomials_grad, &
    stencil_star_coeff_grad

contains

  ! ---------------------------------------------------------
  integer function stencil_star_size_lapl(dim, order)
    integer, intent(in) :: dim
    integer, intent(in) :: order

    PUSH_SUB(stencil_star_size_lapl)

    stencil_star_size_lapl = 2*dim*order + 1

    POP_SUB(stencil_star_size_lapl)
  end function stencil_star_size_lapl


  ! ---------------------------------------------------------
  !> Returns maximum extension of the stencil in spatial direction
  !! dir = 1, 2, 3 for a given discretization order.
  integer function stencil_star_extent(dir, order)
    integer, intent(in) :: dir
    integer, intent(in) :: order

    PUSH_SUB(stencil_star_extent)

    stencil_star_extent = order

    POP_SUB(stencil_star_extent)
  end function stencil_star_extent
  

  ! ---------------------------------------------------------
  subroutine stencil_star_get_lapl(this, dim, order)
    type(stencil_t), intent(out) :: this
    integer,         intent(in)  :: dim
    integer,         intent(in)  :: order

    integer :: ii, jj, nn
    logical :: got_center

    PUSH_SUB(stencil_star_get_lapl)

    call stencil_allocate(this, stencil_star_size_lapl(dim, order))

    got_center = .false.

    nn = 0
    do ii = 1, dim
      do jj = -order, order

        if(jj == 0) then
          if(got_center) then
            cycle
          else
            got_center = .true.
          end if
        end if

        nn = nn + 1
        this%points(ii, nn) = jj
      end do
    end do

    call stencil_init_center(this)

    POP_SUB(stencil_star_get_lapl)
  end subroutine stencil_star_get_lapl


  ! ---------------------------------------------------------
  subroutine stencil_star_polynomials_lapl(dim, order, pol)
    integer, intent(in)  :: dim
    integer, intent(in)  :: order
    integer, intent(out) :: pol(:,:) !< pol(dim, order)

    integer :: i, j, n

    PUSH_SUB(stencil_star_polynomials_lapl)

    n = 1
    pol(:,:) = 0
    do i = 1, dim
      do j = 1, 2*order
        n = n + 1
        pol(i, n) = j
      end do
    end do

    POP_SUB(stencil_star_polynomials_lapl)
  end subroutine stencil_star_polynomials_lapl


  ! ---------------------------------------------------------
  subroutine stencil_star_coeff_lapl(dim, order, h, lapl)
    integer,                intent(in)    :: dim
    integer,                intent(in)    :: order
    FLOAT,                  intent(in)    :: h(:)   !< h(dim)
    type(nl_operator_t),    intent(inout) :: lapl

    integer :: k, i, j, morder
    FLOAT, allocatable :: cc(:,:,:)

    PUSH_SUB(stencil_star_coeff_lapl)

    ASSERT(order >= 1)

    morder = 2*order
    SAFE_ALLOCATE(cc(0:morder, 0:morder, 0:2))
    call weights(2, morder, cc)
    lapl%w_re(1,:) = cc(0, morder, 2)*sum(1/h(1:dim)**2)

    k = 1
    do i = 1, dim
      do j = -order, -1
        k = k + 1
        lapl%w_re(k,:) = cc(-2*j-1, morder, 2) / h(i)**2
      end do

      do j = 1, order
        k = k + 1
        lapl%w_re(k,:) = cc( 2*j,   morder, 2) / h(i)**2
      end do
    end do

    SAFE_DEALLOCATE_A(cc)

    POP_SUB(stencil_star_coeff_lapl)
  end subroutine stencil_star_coeff_lapl


  !> now come the gradient routines

  ! ---------------------------------------------------------
  integer function stencil_star_size_grad(order)
    integer, intent(in) :: order

    PUSH_SUB(stencil_star_size_grad)

    stencil_star_size_grad = 2*order + 1

    POP_SUB(stencil_star_size_grad)
  end function stencil_star_size_grad


  ! ---------------------------------------------------------
  subroutine stencil_star_get_grad(this, dir, order)
    type(stencil_t), intent(out) :: this
    integer, intent(in)  :: dir
    integer, intent(in)  :: order

    integer :: i, n

    PUSH_SUB(stencil_star_get_grad)

    call stencil_allocate(this, stencil_star_size_grad(order))

    n = 1
    do i = -order, order
      this%points(dir, n) = i
      n = n + 1
    end do

    call stencil_init_center(this)

    POP_SUB(stencil_star_get_grad)
  end subroutine stencil_star_get_grad


  ! ---------------------------------------------------------
  subroutine stencil_star_polynomials_grad(dir, order, pol)
    integer, intent(in)  :: dir
    integer, intent(in)  :: order
    integer, intent(out) :: pol(:,:) !< pol(dim, order)

    integer :: j

    PUSH_SUB(stencil_star_polynomials_grad)

    pol(:,:) = 0
    do j = 0, 2*order
      pol(dir, j+1) = j
    end do

    POP_SUB(stencil_star_polynomials_grad)
  end subroutine stencil_star_polynomials_grad


  ! ---------------------------------------------------------
  subroutine stencil_star_coeff_grad(order, h, grad)
    integer,             intent(in)    :: order
    FLOAT,               intent(in)    :: h
    type(nl_operator_t), intent(inout) :: grad

    integer :: j, k, morder
    FLOAT, allocatable :: cc(:,:,:)

    PUSH_SUB(stencil_star_coeff_grad)

    ASSERT(order >= 1)

    morder = 2*order
    SAFE_ALLOCATE(cc(0:morder, 0:morder, 0:1))
    call weights(1, morder, cc)

    k = 1
    do j = -order, -1
      grad%w_re(k,:) = cc(-2*j-1, morder, 1) / h
      k = k + 1
    end do

    grad%w_re(k,:) = cc(0, morder, 1) / h

    do j = 1, order
      k = k + 1
      grad%w_re(k,:) = cc(2*j, morder, 1) / h
    end do

    SAFE_DEALLOCATE_A(cc)

    POP_SUB(stencil_star_coeff_grad)
  end subroutine stencil_star_coeff_grad

end module stencil_star_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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

module stencil_star
  use global
  use math, only: weights
  use nl_operator

  implicit none

contains

  ! ---------------------------------------------------------
  integer function stencil_star_size_lapl(order)
    integer, intent(in) :: order
    
    call push_sub('stencil_star_size_lapl')

    stencil_star_size_lapl = 2*conf%dim*order + 1

    call pop_sub()
  end function stencil_star_size_lapl


  ! ---------------------------------------------------------
  subroutine stencil_star_get_lapl(order, stencil)
    integer, intent(in)  :: order
    integer, intent(out) :: stencil(:,:)

    integer :: i, j, n
    
    call push_sub('stencil_star_get_lapl')

    stencil(:,:) = 0
    n = 1
    do i = 1, conf%dim
      do j = -order, order
        if(j == 0) cycle
        n = n + 1
        stencil(i, n) = j
      end do
    end do
    
    call pop_sub()
  end subroutine stencil_star_get_lapl


  ! ---------------------------------------------------------
  subroutine stencil_star_polynomials_lapl(pol, order)
    integer, intent(out) :: pol(:,:) ! pol(conf%dim, order)
    integer, intent(in)  :: order

    integer :: i, j, n
    
    call push_sub('stencil_star_polynomials_lapl')

    n = 1
    pol(:,:) = 0
    do i = 1, conf%dim
      do j = 1, 2*order
        n = n + 1
        pol(i, n) = j
      end do
    end do

    call pop_sub()
  end subroutine stencil_star_polynomials_lapl


  ! ---------------------------------------------------------
  subroutine stencil_star_coeff_lapl(h, order, lapl)
    FLOAT,                  intent(in)    :: h(:)   ! h(conf%dim)
    integer,                intent(in)    :: order
    type(nl_operator_type), intent(inout) :: lapl

    integer :: k, i, j, morder
    FLOAT, allocatable :: cc(:,:,:)

    call push_sub('stencil_star_coeff_lapl')

    ASSERT(order >= 1)

    morder = 2*order
    allocate(cc(0:morder, 0:morder, 0:2))
    call weights(2, morder, cc)
    lapl%w(1,1) = cc(0, morder, 2)*sum(1/h(1:conf%dim)**2)
      
    k = 1
    do i = 1, conf%dim
      do j = -order, -1
        k = k + 1
        lapl%w(k,1) = cc(-2*j-1, morder, 2) / h(i)**2
      end do
      
      do j = 1, order
        k = k + 1
        lapl%w(k,1) = cc( 2*j,   morder, 2) / h(i)**2
      end do
    end do
      
    deallocate(cc)

    call pop_sub()
  end subroutine stencil_star_coeff_lapl

  
  ! now come the gradient routines


  ! ---------------------------------------------------------
  integer function stencil_star_size_grad(order)
    integer, intent(in) :: order
    
    call push_sub('stencil_star_size_grad')

    stencil_star_size_grad = 2*order + 1

    call pop_sub()
  end function stencil_star_size_grad


  ! ---------------------------------------------------------
  subroutine stencil_star_get_grad(dir, order, stencil)
    integer, intent(in)  :: dir
    integer, intent(in)  :: order
    integer, intent(out) :: stencil(:,:)

    integer :: i, n
      
    call push_sub('stencil_star_get_grad')

    stencil(:,:) = 0
    n = 1
    do i = -order, order
      stencil(dir, n) = i
      n = n + 1
    end do
    
    call pop_sub()
  end subroutine stencil_star_get_grad


  ! ---------------------------------------------------------
  subroutine stencil_star_polynomials_grad(dir, pol, order)
    integer, intent(in)  :: dir
    integer, intent(out) :: pol(:,:) ! pol(conf%dim, order)
    integer, intent(in)  :: order

    integer :: j
      
    call push_sub('stencil_star_polynomials_grad')

    pol(:,:) = 0
    do j = 0, 2*order
      pol(dir, j+1) = j
    end do
    
    call pop_sub()
  end subroutine stencil_star_polynomials_grad

  ! ---------------------------------------------------------
  subroutine stencil_star_coeff_grad(h, order, grad)
    FLOAT,                  intent(in)    :: h
    integer,                intent(in)    :: order
    type(nl_operator_type), intent(inout) :: grad

    integer :: j, k, morder
    FLOAT, allocatable :: cc(:,:,:)

    call push_sub('stencil_star_coeff_grad')

    ASSERT(order >= 1)

    morder = 2*order
    allocate(cc(0:morder, 0:morder, 0:1))
    call weights(1, morder, cc)
  
    k = 1
    do j = -order, -1
      grad%w(k,1) = cc(-2*j-1, morder, 1) / h
      k = k + 1
    end do

    grad%w(k,1) = cc(0, morder, 1) / h

    do j = 1, order
      k = k + 1
      grad%w(k,1) = cc(2*j, morder, 1) / h
    end do

    deallocate(cc)

    call pop_sub()
  end subroutine stencil_star_coeff_grad

end module stencil_star

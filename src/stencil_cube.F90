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

module stencil_cube
  use global
  use nl_operator

  implicit none

contains

  ! ---------------------------------------------------------
  integer function stencil_cube_size_lapl(order)
    integer, intent(in) :: order
    
    call push_sub('stencil_cube_size_lapl')

    stencil_cube_size_lapl = (2*order+1)**conf%dim

    call pop_sub()
  end function stencil_cube_size_lapl


  ! ---------------------------------------------------------
  subroutine stencil_cube_get_lapl(order, stencil)
    integer, intent(in)  :: order
    integer, intent(out) :: stencil(:,:)

    integer :: i, j, k, n

    call push_sub('stencil_cube_get_lapl')

    n = 1
    select case(conf%dim)
    case(1)
      do i = -order, order
        stencil(1, n) = i
        n = n + 1
      end do
    case(2)
      do i = -order, order
        do j = -order, order
          stencil(1, n) = i
          stencil(2, n) = j
          n = n + 1
        end do
      end do
    case(3)
      do i = -order, order
        do j = -order, order
          do k = -order, order
            stencil(1, n) = i
            stencil(2, n) = j
            stencil(3, n) = k
            n = n + 1
          end do
        end do
      end do
    end select
    
    call pop_sub()
  end subroutine stencil_cube_get_lapl
  

  ! ---------------------------------------------------------
  subroutine stencil_cube_polynomials_lapl(pol, order)
    integer, intent(out) :: pol(:,:) ! pol(conf%dim, order)
    integer, intent(in)  :: order

    integer :: i, j, k, n

    call push_sub('stencil_cube_polynomials_lapl')

    n = 1
    select case(conf%dim)
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
    
    call pop_sub()
  end subroutine stencil_cube_polynomials_lapl


  ! Now come the gradient routines. As this stencil is the same for
  ! the laplacian and the gradient, these routines just call the
  ! corresponding ones for the laplacian

  ! ---------------------------------------------------------
  integer function stencil_cube_size_grad(order)
    integer, intent(in) :: order
    
    stencil_cube_size_grad = stencil_cube_size_lapl(order)
  end function stencil_cube_size_grad


  ! ---------------------------------------------------------
  subroutine stencil_cube_get_grad(dir, order, stencil)
    integer, intent(in)  :: dir
    integer, intent(in)  :: order
    integer, intent(out) :: stencil(:,:)

    call stencil_cube_get_lapl(order, stencil)
  end subroutine stencil_cube_get_grad


  ! ---------------------------------------------------------
  subroutine stencil_cube_polynomials_grad(dir, pol, order)
    integer, intent(in)  :: dir
    integer, intent(out) :: pol(:,:) ! pol(conf%dim, order)
    integer, intent(in)  :: order

    call stencil_cube_polynomials_lapl(pol, order)
  end subroutine stencil_cube_polynomials_grad

end module stencil_cube

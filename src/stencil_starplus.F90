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

module stencil_starplus
  use global
  use messages
  use nl_operator

  private
  public :: stencil_starplus_size_lapl, &
            stencil_starplus_get_lapl, &
            stencil_starplus_pol_lapl, &
            stencil_starplus_size_grad, &
            stencil_starplus_get_grad, &
            stencil_starplus_pol_grad

  contains

  integer function stencil_starplus_size_lapl(order) result(n)
    integer, intent(in) :: order

    call push_sub('stencil_starplus_size_lapl')

    n = 2*conf%dim*order + 1
    if(conf%dim == 2) n = n + 4

    call pop_sub()
  end function stencil_starplus_size_lapl

  subroutine stencil_starplus_get_lapl(order, stencil)
    integer, intent(in)  :: order
    integer, intent(out) :: stencil(:,:)

    integer :: i, j, k, n

    call push_sub('stencil_starplus_get_lapl')

    stencil = 0
    n = 1
    select case(conf%dim)
    case(1)
       stencil(:,:) = 0
       n = 1
       do i = 1, conf%dim
          do j = -order, order
             if(j == 0) cycle
             n = n + 1
             stencil(i, n) = j
         end do
       end do
    case(2)
       stencil(:,:) = 0
       n = 1
       do i = 1, conf%dim
          do j = -order, order
             if(j == 0) cycle
             n = n + 1
             stencil(i, n) = j
         end do
       end do
       n = n + 1
       stencil(1, n) = -1
       stencil(2, n) = -1
       n = n + 1
       stencil(1, n) = 1
       stencil(2, n) = -1
       n = n + 1
       stencil(1, n) = -1
       stencil(2, n) = 1
       n = n + 1
       stencil(1, n) = 1
       stencil(2, n) = 1
    case(3)
      stop 'Not yet implemented.'
    end select
    
    call pop_sub()
  end subroutine stencil_starplus_get_lapl

  subroutine stencil_starplus_pol_lapl(pol, order)
    integer, intent(out) :: pol(:,:) ! pol(conf%dim, order)
    integer, intent(in)  :: order

    integer :: i, j, k, n

    call push_sub('stencil_starplus_pol_lapl')

    n = 1
    select case(conf%dim)
    case(1)
      n = 1
      pol(:,:) = 0
      do i = 1, conf%dim
         do j = 1, 2*order
            n = n + 1
            pol(i, n) = j
         end do
      end do
    case(2)
      n = 1
      pol(:,:) = 0
      do i = 1, conf%dim
         do j = 1, 2*order
            n = n + 1
            pol(i, n) = j
         end do
      end do
      n = n + 1
      pol(1, n) = 1
      pol(2, n) = 1
      n = n + 1
      pol(1, n) = 1
      pol(2, n) = 2
      n = n + 1
      pol(1, n) = 2
      pol(2, n) = 1
      n = n + 1
      pol(1, n) = 2
      pol(2, n) = 2
    case(3)
      stop 'Not yet implemented.'
    end select
    
    call pop_sub()
  end subroutine stencil_starplus_pol_lapl


  integer function stencil_starplus_size_grad(order)
    integer, intent(in) :: order
    stencil_starplus_size_grad = stencil_starplus_size_lapl(order)
  end function stencil_starplus_size_grad

  subroutine stencil_starplus_get_grad(order, stencil)
    integer, intent(in)  :: order
    integer, intent(out) :: stencil(:,:)
    call stencil_starplus_get_lapl(order, stencil)
  end subroutine stencil_starplus_get_grad

  subroutine stencil_starplus_pol_grad(pol, order)
    integer, intent(out) :: pol(:,:) ! pol(conf%dim, order)
    integer, intent(in)  :: order
    call stencil_starplus_pol_lapl(pol, order)
  end subroutine stencil_starplus_pol_grad


end module stencil_starplus

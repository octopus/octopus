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

  integer function stencil_starplus_size_lapl(dim, order) result(n)
    integer, intent(in) :: dim
    integer, intent(in) :: order

    call push_sub('stencil_starplus_size_lapl')

    n = 2*dim*order + 1
    if(dim == 2) n = n + 12
    if(dim == 3) n = n + 44

    call pop_sub()
  end function stencil_starplus_size_lapl


  subroutine stencil_starplus_get_lapl(dim, order, stencil)
    integer, intent(in)  :: dim
    integer, intent(in)  :: order
    integer, intent(out) :: stencil(:,:)

    integer :: i, j, k, n

    call push_sub('stencil_starplus_get_lapl')

    stencil = 0
    n = 1
    select case(dim)
    case(1)
       stencil(:,:) = 0
       n = 1
       do i = 1, dim
          do j = -order, order
             if(j == 0) cycle
             n = n + 1
             stencil(i, n) = j
         end do
       end do
    case(2)
       stencil(:,:) = 0
       n = 1
       do i = 1, dim
          do j = -order, order
             if(j == 0) cycle
             n = n + 1
             stencil(i, n) = j
         end do
       end do
       n = n + 1; stencil(1:2, n) = (/ -2,  1 /)
       n = n + 1; stencil(1:2, n) = (/ -2, -1 /)
       n = n + 1; stencil(1:2, n) = (/ -1,  2 /)
       n = n + 1; stencil(1:2, n) = (/ -1,  1 /)
       n = n + 1; stencil(1:2, n) = (/ -1, -1 /)
       n = n + 1; stencil(1:2, n) = (/ -1, -2 /)
       n = n + 1; stencil(1:2, n) = (/  1,  2 /)
       n = n + 1; stencil(1:2, n) = (/  1,  1 /)
       n = n + 1; stencil(1:2, n) = (/  1, -1 /)
       n = n + 1; stencil(1:2, n) = (/  1, -2 /)
       n = n + 1; stencil(1:2, n) = (/  2,  1 /)
       n = n + 1; stencil(1:2, n) = (/  2, -1 /)
    case(3)
       stencil(:,:) = 0
       n = 1
       do i = 1, dim
          do j = -order, order
             if(j == 0) cycle
             n = n + 1
             stencil(i, n) = j
         end do
       end do
       n = n + 1; stencil(1:3, n) = (/ -2,  1, 0 /)
       n = n + 1; stencil(1:3, n) = (/ -2, -1, 0 /)
       n = n + 1; stencil(1:3, n) = (/ -1,  2, 0 /)
       n = n + 1; stencil(1:3, n) = (/ -1,  1, 0 /)
       n = n + 1; stencil(1:3, n) = (/ -1, -1, 0 /)
       n = n + 1; stencil(1:3, n) = (/ -1, -2, 0 /)
       n = n + 1; stencil(1:3, n) = (/  1,  2, 0 /)
       n = n + 1; stencil(1:3, n) = (/  1,  1, 0 /)
       n = n + 1; stencil(1:3, n) = (/  1, -1, 0 /)
       n = n + 1; stencil(1:3, n) = (/  1, -2, 0 /)
       n = n + 1; stencil(1:3, n) = (/  2,  1, 0 /)
       n = n + 1; stencil(1:3, n) = (/  2, -1, 0 /)

       n = n + 1; stencil(1:3, n) = (/ -2, 0,  1 /)
       n = n + 1; stencil(1:3, n) = (/ -2, 0, -1 /)
       n = n + 1; stencil(1:3, n) = (/ -1, 0,  2 /)
       n = n + 1; stencil(1:3, n) = (/ -1, 0,  1 /)
       n = n + 1; stencil(1:3, n) = (/ -1, 0, -1 /)
       n = n + 1; stencil(1:3, n) = (/ -1, 0, -2 /)
       n = n + 1; stencil(1:3, n) = (/  1, 0,  2 /)
       n = n + 1; stencil(1:3, n) = (/  1, 0,  1 /)
       n = n + 1; stencil(1:3, n) = (/  1, 0, -1 /)
       n = n + 1; stencil(1:3, n) = (/  1, 0, -2 /)
       n = n + 1; stencil(1:3, n) = (/  2, 0,  1 /)
       n = n + 1; stencil(1:3, n) = (/  2, 0, -1 /)

       n = n + 1; stencil(1:3, n) = (/ 0, -2,  1 /)
       n = n + 1; stencil(1:3, n) = (/ 0, -2, -1 /)
       n = n + 1; stencil(1:3, n) = (/ 0, -1,  2 /)
       n = n + 1; stencil(1:3, n) = (/ 0, -1,  1 /)
       n = n + 1; stencil(1:3, n) = (/ 0, -1, -1 /)
       n = n + 1; stencil(1:3, n) = (/ 0, -1, -2 /)
       n = n + 1; stencil(1:3, n) = (/ 0,  1,  2 /)
       n = n + 1; stencil(1:3, n) = (/ 0,  1,  1 /)
       n = n + 1; stencil(1:3, n) = (/ 0,  1, -1 /)
       n = n + 1; stencil(1:3, n) = (/ 0,  1, -2 /)
       n = n + 1; stencil(1:3, n) = (/ 0,  2,  1 /)
       n = n + 1; stencil(1:3, n) = (/ 0,  2, -1 /)

       n = n + 1; stencil(1:3, n) = (/ -1, -1, -1 /)
       n = n + 1; stencil(1:3, n) = (/ -1, -1,  1 /)
       n = n + 1; stencil(1:3, n) = (/ -1,  1, -1 /)
       n = n + 1; stencil(1:3, n) = (/ -1,  1,  1 /)
       n = n + 1; stencil(1:3, n) = (/  1, -1, -1 /)
       n = n + 1; stencil(1:3, n) = (/  1, -1,  1 /)
       n = n + 1; stencil(1:3, n) = (/  1,  1, -1 /)
       n = n + 1; stencil(1:3, n) = (/  1,  1,  1 /)

    end select
    
    call pop_sub()
  end subroutine stencil_starplus_get_lapl

  subroutine stencil_starplus_pol_lapl(dim, order, pol)
    integer, intent(in)  :: dim
    integer, intent(in)  :: order
    integer, intent(out) :: pol(:,:) ! pol(dim, order)

    integer :: i, j, k, n

    call push_sub('stencil_starplus_pol_lapl')

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
    
    call pop_sub()
  end subroutine stencil_starplus_pol_lapl


  integer function stencil_starplus_size_grad(dim, order)
    integer, intent(in) :: dim
    integer, intent(in) :: order

    stencil_starplus_size_grad = stencil_starplus_size_lapl(dim, order)
  end function stencil_starplus_size_grad


  subroutine stencil_starplus_get_grad(dim, order, stencil)
    integer, intent(in)  :: dim
    integer, intent(in)  :: order
    integer, intent(out) :: stencil(:,:)

    call stencil_starplus_get_lapl(dim, order, stencil)
  end subroutine stencil_starplus_get_grad


  subroutine stencil_starplus_pol_grad(dim, order, pol)
    integer, intent(in)  :: dim
    integer, intent(in)  :: order
    integer, intent(out) :: pol(:,:) ! pol(sb%dim, order)

    call stencil_starplus_pol_lapl(dim, order, pol)
  end subroutine stencil_starplus_pol_grad


end module stencil_starplus

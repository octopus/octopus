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

! ----------------------------------------
! Implements the variational discretization of the laplacian
! as proposed by P. Maragakis, J. Soler, and E. Kaxiras, PRB 64, 193101 (2001)
! ----------------------------------------

module stencil_variational
  use global
  use nl_operator

  implicit none

contains


  ! ---------------------------------------------------------
  subroutine stencil_variational_coeff_lapl(h, order, lapl)
    FLOAT,                  intent(in)    :: h(:)   ! h(conf%dim)
    integer,                intent(in)    :: order
    type(nl_operator_type), intent(inout) :: lapl

    integer :: i, j, k
    FLOAT, parameter   :: pi2 = M_PI**2
    FLOAT, allocatable :: fp(:)

    call push_sub('stencil_variational_coeff_lapl')

    allocate(fp(order+1))
    select case(order)
    case(1) 
      fp(1) = -pi2/M_TWO
      fp(2) =  pi2/M_FOUR
    case(2) 
      fp(1) = -M_HALF-M_THREE*pi2/M_EIGHT
      fp(2) =  pi2/M_FOUR
      fp(3) =  M_ONE/M_FOUR - pi2/CNST(16.0)
    case(3)
      fp(1) = -M_FIVE/M_SIX - M_FIVE*pi2/CNST(16.0)
      fp(2) =  M_ONE/CNST(12.0) + CNST(15.0)*pi2/CNST(64.0)
      fp(3) =  M_FIVE/CNST(12.0) - M_THREE*pi2/CNST(32.0)
      fp(4) = -M_ONE/CNST(12.0) + pi2/CNST(64.0)
    case(4)
      fp(1) = -CNST(77.0)/CNST(72.0) - CNST(35.0)*pi2/CNST(128.0)
      fp(2) =  M_EIGHT/CNST(45.0) + M_SEVEN*pi2/CNST(32.0)
      fp(3) =  CNST(23.0)/CNST(45.0) - M_SEVEN*pi2/CNST(64.0)
      fp(4) = -M_EIGHT/CNST(45.0) + pi2/CNST(32.0)
      fp(5) =  CNST(17.0)/CNST(720.0) - pi2/CNST(256.0)
    case(5)
      fp(1) = -CNST(449.0)/CNST(360.0) - CNST(63.0)*pi2/CNST(256.0)
      fp(2) =  M_FOUR/CNST(15.0) + CNST(105.0)*pi2/CNST(512.0)
      fp(3) =  CNST(59.0)/CNST(105.0) - CNST(15.0)*pi2/CNST(128.0)
      fp(4) = -CNST(82.0)/CNST(315.0) + CNST(45.0)*pi2/CNST(1024.0)
      fp(5) =  CNST(311.0)/CNST(5040.0) - M_FIVE*pi2/CNST(512.0)
      fp(6) = -M_TWO/CNST(315.0) + pi2/CNST(1024.0)
    case(6)
      fp(1) = -CNST(2497.0)/CNST(1800.0) - CNST(231.0)*pi2/CNST(1024.0)
      fp(2) =  CNST(26.0)/CNST(75.0) + CNST(99.0)*pi2/CNST(512.0)
      fp(3) =  CNST(493.0)/CNST(840.0) - CNST(495.0)*pi2/CNST(4096.0)
      fp(4) = -CNST(103.0)/CNST(315.0) + CNST(55.0)*pi2/CNST(1024.0)
      fp(5) =  CNST(2647.0)/CNST(25200.0) - CNST(33.0)*pi2/CNST(2048.0)
      fp(6) = -CNST(31.0)/CNST(1575.0) + M_THREE*pi2/CNST(1024.0)
      fp(7) =  M_ONE/CNST(600.0) - pi2/CNST(4096.0)
    end select
    
    lapl%w(1,1) = fp(1)*sum(1/h(1:conf%dim)**2)

    k = 1
    do i = 1, conf%dim
      do j = -order, -1
        k = k + 1
        lapl%w(k,1) = fp(-j+1) / h(i)**2
      end do

      do j = 1, order
        k = k + 1
        lapl%w(k,1) = fp( j+1) / h(i)**2
      end do
    end do
    
    deallocate(fp)

    call pop_sub()
  end subroutine stencil_variational_coeff_lapl

end module stencil_variational

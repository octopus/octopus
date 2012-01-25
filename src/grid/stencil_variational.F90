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
! ---------------------------------------------------------
!> Implements the variational discretization of the Laplacian
!! as proposed by P. Maragakis, J. Soler, and E. Kaxiras, PRB 64, 193101 (2001)
!!
!! However, we have introduced a possible variation: incorporating into
!! the expression of the Laplacian a high-frequency filter. This in fact
!! goes against the spirit of the work of Maragakis et al., which
!! attempts to increase the weight of the high frequencies over the
!! conventional finite-difference scheme. But the mathematical
!! machinery used in that work to generate the coefficient is ideal
!! to add a frequency filter.
!! The filter is decided by the optional parameter alpha: The
!! highest frequency allowed by the filter will be alpha*k, where
!! k is the Nyquist frequency of the grid. Thus alpha = 1 means
!! no filter at all.
!!
!! \todo The coefficients should be calculated, and not hard-coded, in
!! a subroutine similar to weights in the math module. It should also
!! allos to get the coefficients to the gradient.
!! \todo This module should look like stencil_star, allowing for
!! coefficients on non-uniform grids.
! ---------------------------------------------------------
module stencil_variational_m
  use global_m
  use messages_m
  use profiling_m
  use nl_operator_m

  implicit none

  private
  public ::                         &
    stencil_variational_coeff_lapl, &
    stencil_variational_extent

contains

  ! ---------------------------------------------------------
  subroutine stencil_variational_coeff_lapl(dim, order, h, lapl, alpha)
    integer,                intent(in)    :: dim
    integer,                intent(in)    :: order
    FLOAT,                  intent(in)    :: h(:)   ! h(dim)
    type(nl_operator_t), intent(inout) :: lapl
    FLOAT, optional,        intent(in)    :: alpha

    integer :: i, j, k
    FLOAT :: alpha_, kmax
    FLOAT, parameter   :: pi2 = M_PI**2
    FLOAT, allocatable :: fp(:)

    PUSH_SUB(stencil_variational_coeff_lapl)

    alpha_ = M_ONE
    if(present(alpha)) alpha_ = alpha
    kmax = M_PI**2 * alpha_

    SAFE_ALLOCATE(fp(1:order+1))
    select case(order)
    case(1)
      fp(1) = -kmax/M_TWO
      fp(2) =  kmax/M_FOUR
    case(2)
      fp(1) = -M_HALF-M_THREE*kmax/M_EIGHT
      fp(2) =  kmax/M_FOUR
      fp(3) =  M_ONE/M_FOUR - kmax/CNST(16.0)
    case(3)
      fp(1) = -M_FIVE/M_SIX - M_FIVE*kmax/CNST(16.0)
      fp(2) =  M_ONE/CNST(12.0) + CNST(15.0)*kmax/CNST(64.0)
      fp(3) =  M_FIVE/CNST(12.0) - M_THREE*kmax/CNST(32.0)
      fp(4) = -M_ONE/CNST(12.0) + kmax/CNST(64.0)
    case(4)
      fp(1) = -CNST(77.0)/CNST(72.0) - CNST(35.0)*kmax/CNST(128.0)
      fp(2) =  M_EIGHT/CNST(45.0) + M_SEVEN*kmax/CNST(32.0)
      fp(3) =  CNST(23.0)/CNST(45.0) - M_SEVEN*kmax/CNST(64.0)
      fp(4) = -M_EIGHT/CNST(45.0) + kmax/CNST(32.0)
      fp(5) =  CNST(17.0)/CNST(720.0) - kmax/CNST(256.0)
    case(5)
      fp(1) = -CNST(449.0)/CNST(360.0) - CNST(63.0)*kmax/CNST(256.0)
      fp(2) =  M_FOUR/CNST(15.0) + CNST(105.0)*kmax/CNST(512.0)
      fp(3) =  CNST(59.0)/CNST(105.0) - CNST(15.0)*kmax/CNST(128.0)
      fp(4) = -CNST(82.0)/CNST(315.0) + CNST(45.0)*kmax/CNST(1024.0)
      fp(5) =  CNST(311.0)/CNST(5040.0) - M_FIVE*kmax/CNST(512.0)
      fp(6) = -M_TWO/CNST(315.0) + kmax/CNST(1024.0)
    case(6)
      fp(1) = -CNST(2497.0)/CNST(1800.0) - CNST(231.0)*kmax/CNST(1024.0)
      fp(2) =  CNST(26.0)/CNST(75.0) + CNST(99.0)*kmax/CNST(512.0)
      fp(3) =  CNST(493.0)/CNST(840.0) - CNST(495.0)*kmax/CNST(4096.0)
      fp(4) = -CNST(103.0)/CNST(315.0) + CNST(55.0)*kmax/CNST(1024.0)
      fp(5) =  CNST(2647.0)/CNST(25200.0) - CNST(33.0)*kmax/CNST(2048.0)
      fp(6) = -CNST(31.0)/CNST(1575.0) + M_THREE*kmax/CNST(1024.0)
      fp(7) =  M_ONE/CNST(600.0) - kmax/CNST(4096.0)
    end select

    lapl%w_re(1,:) = fp(1)*sum(1/h(1:dim)**2)

    k = 1
    do i = 1, dim
      do j = -order, -1
        k = k + 1
        lapl%w_re(k,:) = fp(-j+1) / h(i)**2
      end do

      do j = 1, order
        k = k + 1
        lapl%w_re(k,:) = fp( j+1) / h(i)**2
      end do
    end do

    SAFE_DEALLOCATE_A(fp)

    POP_SUB(stencil_variational_coeff_lapl)
  end subroutine stencil_variational_coeff_lapl


  ! ---------------------------------------------------------
  !> Returns maximum extension of the stencil in spatial direction
  !! dir = 1, 2, 3 for a given discretization order.
  integer function stencil_variational_extent(dir, order)
    integer, intent(in) :: dir
    integer, intent(in) :: order

    PUSH_SUB(stencil_variational_extent)

    stencil_variational_extent = order

    POP_SUB(stencil_variational_extent)
  end function stencil_variational_extent
end module stencil_variational_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

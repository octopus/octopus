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

module curvlinear
  use geometry

  implicit none

!  private
  public :: jacobian

  FLOAT, parameter :: alpha(3) = (/1.5, 0.5, 2.0/)

contains
  subroutine x_to_chi(geo, x, chi)
    type(geometry_type), intent(in) :: geo
    FLOAT, intent(in)  :: x(:)    ! x(conf%dim)
    FLOAT, intent(out) :: chi(:)  ! chi(conf%dim)
    
    integer :: i
    FLOAT :: r, f_alpha

    chi = x
    do i = 1, geo%natoms
      r = max(sqrt(sum((x - geo%atom(i)%x)**2)), CNST(1e-6))
      f_alpha = alpha(1)*alpha(2)/r*tanh(r/alpha(2))*exp(-(r/alpha(3))**2)

      chi = chi + (x - geo%atom(i)%x) * f_alpha
    end do
  end subroutine x_to_chi

  subroutine jacobian(geo, x, chi, J)
    type(geometry_type), intent(in) :: geo
    FLOAT, intent(in)  :: x(:)    ! x(conf%dim)
    FLOAT, intent(out) :: chi(:)  ! chi(conf%dim)
    FLOAT, intent(out) :: J(:,:)  ! J(conf%dim,conf%dim), the Jacobian

    integer :: i, ix, iy
    FLOAT :: r, f_alpha, df_alpha
    FLOAT :: th, ex, ar

    J(1:conf%dim,1:conf%dim) = M_ZERO
    do ix = 1, conf%dim
      J(ix, ix) = M_ONE
      chi(ix)   = x(ix)
    end do

    do i = 1, geo%natoms
      r = max(sqrt(sum((x - geo%atom(i)%x)**2)), CNST(1e-6))
      
      ar = alpha(1)*alpha(2)/r
      th = tanh(r/alpha(2))
      ex = exp(-(r/alpha(3))**2)
      
      f_alpha  = ar * th * ex
      df_alpha = ar*(-th*ex/r + ex/(alpha(2)*cosh(r/alpha(2))**2) - th*M_TWO*r*ex/alpha(3)**2)
      
      do ix = 1, conf%dim
        chi(ix) = chi(ix) + f_alpha*(x(ix)-geo%atom(i)%x(ix))
        
        J(ix, ix) = J(ix, ix) + f_alpha
        do iy = 1, conf%dim
          J(ix, iy) = J(ix, iy) + (x(ix)-geo%atom(i)%x(ix))*(x(iy)-geo%atom(i)%x(iy))/r * df_alpha
        end do
      end do
    end do

  end subroutine jacobian

  subroutine chi_to_x(geo, chi, x)
    type(geometry_type), intent(in) :: geo
    FLOAT, intent(in)  :: chi(:)  ! chi(conf%dim)
    FLOAT, intent(out) :: x(:)    ! x(conf%dim)
    
    integer :: iter
    FLOAT, allocatable :: f(:,:), delta(:,:), J(:,:), chi2(:)
    logical :: conv
  
    allocate(f(conf%dim, 1), delta(conf%dim, 1), J(conf%dim, conf%dim), chi2(conf%dim))

    x(1:conf%dim) = chi(1:conf%dim)
    conv          = .false.

    do iter = 1, 500
      call jacobian(geo, x, chi2, J)
      f(:,1) = chi(1:conf%dim) - chi2(:)

      if(sum(f(:,1)**2) < CNST(1e-6)**2) then
        conv = .true.
        exit
      end if

      call lalg_linsyssolve(conf%dim, 1, J, f, delta)
      x(1:conf%dim) = x(1:conf%dim) + delta(1:conf%dim, 1)
    end do

    if(.not.conv) then
      message(1) = "Newton-Raphson method did not converge for point"
      write(message(2), '(3es14.5)') x(1:conf%dim)
      call write_warning(2)
    end if

    deallocate(f, delta, J, chi2)
  end subroutine chi_to_x

end module curvlinear

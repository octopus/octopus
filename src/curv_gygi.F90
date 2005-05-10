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

! This module implements the curvilinear coordinates given in
! F. Gygi and G. Galli, PRB 52 R2229 (1996).

module curv_gygi
  use lib_oct_parser
  use global
  use units
  use geometry
  use lib_adv_alg

  implicit none

  type curv_gygi_type
    FLOAT :: A       ! local reduction in grid spacing is 1/(1+A)
    FLOAT :: alpha   ! range of enhacement of the resolution
    FLOAT :: beta    ! distance over which Euclidian coordinates are recovered
  end type curv_gygi_type

contains
  
  !-------------------------------------
  subroutine curv_gygi_init(cv)
    type(curv_gygi_type), intent(out) :: cv

    call loct_parse_float(check_inp('CurvGygiA'), M_ONE, cv%A)
    call loct_parse_float(check_inp('CurvGygiAlpha'), M_TWO/units_inp%length%factor, cv%alpha)
    call loct_parse_float(check_inp('CurvGygiBeta'),  M_TWO/units_inp%length%factor, cv%beta)
    
    cv%alpha = cv%alpha*units_inp%length%factor
    cv%beta  = cv%beta *units_inp%length%factor

    if(cv%A<=M_ZERO.or.cv%alpha<=M_ZERO.or.cv%beta<=M_ZERO) then
      message(1) = 'The parameters CurvGygiA, CurvGygiAlpha, and CurvGygiBeta'
      message(2) = 'must all be larger than zero'
      call write_fatal(2)
    end if

  end subroutine curv_gygi_init


  !-------------------------------------
  subroutine curv_gygi_chi2x(cv, geo, chi, x)
    type(curv_gygi_type), intent(in)  :: cv
    type(geometry_type),  intent(in)  :: geo
    FLOAT,                intent(in)  :: chi(:)  ! chi(conf%dim)
    FLOAT,                intent(out) :: x(:)    ! x(conf%dim)

    ! parameters
    integer, parameter :: max_iter = 500
    FLOAT,   parameter :: x_conv   = CNST(1e-6)

    ! local variables
    integer :: iter
    FLOAT, allocatable :: f(:,:), delta(:,:), J(:,:), chi2(:)
    logical :: conv
    
    allocate(f(conf%dim, 1), delta(conf%dim, 1), J(conf%dim, conf%dim), chi2(conf%dim))
    
    x(1:conf%dim) = chi(1:conf%dim)
    conv          = .false.
    
    do iter = 1, max_iter
      call curv_gygi_jacobian(cv, geo, x, chi2, J)
      f(:,1) = chi(1:conf%dim) - chi2(:)
      
      if(sum(f(:,1)**2) < x_conv**2) then
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
    
    ! clean up
    deallocate(f, delta, J, chi2)

  end subroutine curv_gygi_chi2x


  !-------------------------------------
  subroutine curv_gygi_jacobian(cv, geo, x, chi, J)
    type(curv_gygi_type), intent(in)  :: cv
    type(geometry_type),  intent(in)  :: geo
    FLOAT,                intent(in)  :: x(:)    ! x(conf%dim)
    FLOAT,                intent(out) :: chi(:)  ! chi(conf%dim)
    FLOAT,                intent(out) :: J(:,:)  ! J(conf%dim,conf%dim), the Jacobian
    
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
      
      ar = cv%A*cv%alpha/r
      th = tanh(r/cv%alpha)
      ex = exp(-(r/cv%beta)**2)
      
      f_alpha  = ar * th * ex
      df_alpha = ar*(-th*ex/r + ex/(cv%alpha*cosh(r/cv%alpha)**2) - th*M_TWO*r*ex/cv%beta**2)
      
      do ix = 1, conf%dim
        chi(ix) = chi(ix) + f_alpha*(x(ix)-geo%atom(i)%x(ix))
        
        J(ix, ix) = J(ix, ix) + f_alpha
        do iy = 1, conf%dim
          J(ix, iy) = J(ix, iy) + (x(ix)-geo%atom(i)%x(ix))*(x(iy)-geo%atom(i)%x(iy))/r * df_alpha
        end do
      end do
    end do

  end subroutine curv_gygi_jacobian

end module curv_gygi

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
! E.L. Briggs, D.J. Sullivan, and J. Bernholc, PRB 54 14362 (1996)
!
! The local refinement was changed for a simple exponential.
! I believe that the recipe given by the authors is too complicated
! for me to sort it out.

module curv_modine
  use lib_oct_parser
  use global
  use units
  use geometry
  use lib_adv_alg

  implicit none

  type curv_modine_type
    FLOAT :: L(3)      ! size of the box
    FLOAT :: xbar      ! size of central flat region (in units of L)
    FLOAT :: Jbar      ! increase in density of points is 1/J
    FLOAT :: Jlocal    ! local (around the atoms) refinement
    FLOAT :: Jrange    ! local refinement range
  end type curv_modine_type

contains

  !-------------------------------------
  subroutine curv_modine_init(l, cv)
    FLOAT,                  intent(in)  :: l(:)  ! l(1:conf%dim)
    type(curv_modine_type), intent(out) :: cv

    call loct_parse_float(check_inp('CurvModineXBar'), M_ONE/M_THREE, cv%xbar)
    call loct_parse_float(check_inp('CurvModineJBar'), M_HALF, cv%Jbar)
    
    cv%L = M_ZERO
    cv%L(1:conf%dim) = l(1:conf%dim)/ cv%Jbar

    if(cv%xbar<M_ZERO.or.cv%xbar>M_ONE) then
      message(1) = 'The parameter "CurvModineXBar" must lie between 0 and 1.'
      call write_fatal(1)
    end if

    call loct_parse_float(check_inp('CurvModineJlocal'), CNST(0.5), cv%Jlocal)
    call loct_parse_float(check_inp('CurvModineJrange'), M_TWO/units_inp%length%factor, cv%Jrange)
    
    if(cv%Jlocal<M_ZERO.or.cv%Jlocal>M_ONE) then
      message(1) = 'The parameter "CurvModineJlocal" must lie between 0 and 1.'
      call write_fatal(1)
    end if

  end subroutine curv_modine_init


  !-------------------------------------
  subroutine curv_modine_chi2x(cv, geo, chi_, x)
    type(curv_modine_type), intent(in)  :: cv
    type(geometry_type),    intent(in)  :: geo
    FLOAT,                  intent(in)  :: chi_(:)  ! chi_(conf%dim)
    FLOAT,                  intent(out) :: x(:)     !   x (conf%dim)
 
    integer, parameter :: q = 3

    FLOAT :: chibar(conf%dim), r, chi
    logical :: neg
    integer :: i
    
    chibar = cv%xbar*cv%L(:)

    do i = 1, conf%dim
      neg = (chi_(i) < 0)
      chi = abs(chi_(i))

      x(i) = cv%Jbar * chi
      if(chi > chibar(i)) then
        r = (chi-chibar(i))/(cv%L(i)-chibar(i))
        x(i) = x(i) + cv%L(i)/M_TWO*(1-cv%Jbar) * r**q *   &
           (q + M_ONE - (q - M_ONE)*r)
      end if

      if(neg) x(i) = -x(i)
    end do

    do i = 1, geo%natoms
      r = max(sqrt(sum((x - geo%atom(i)%x)**2)), CNST(1e-6))
      x(:) = x(:) - cv%Jlocal*(x - geo%atom(i)%x)*exp(-r**2/(M_TWO*cv%Jrange**2))
    end do

  end subroutine curv_modine_chi2x


  !-------------------------------------
  subroutine curv_modine_jacobian_inv(cv, geo, chi_, J)
    type(curv_modine_type), intent(in)  :: cv
    type(geometry_type),    intent(in)  :: geo
    FLOAT,                  intent(in)  :: chi_(:)  ! chi(conf%dim)
    FLOAT,                  intent(out) :: J(:,:)   ! J(conf%dim,conf%dim), the Jacobian
 
    integer, parameter :: q = 3

    FLOAT :: chibar(conf%dim), r, f, chi, J2(conf%dim), x(conf%dim)
    logical :: neg
    integer :: i, ix, iy
    
    chibar = cv%xbar*cv%L(:)

    J2(:) = M_ZERO
    do i = 1, conf%dim
      neg = (chi_(i) < 0)
      chi = abs(chi_(i))

      x(i)  = cv%Jbar * chi
      J2(i) = cv%Jbar
      
      if(chi > chibar(i)) then
        r = (chi-chibar(i))/(cv%L(i)-chibar(i))

        x(i)  = x(i) + cv%L(i)/M_TWO*(1-cv%Jbar) * r**q *   &
           (q + M_ONE - (q - M_ONE)*r)

        J2(i) = J2(i) + cv%L(i)/M_TWO*(1-cv%Jbar) * r**(q-1)/(cv%L(i)-chibar(i)) *   &
           (q*(q+1) - (q**2-1)*r)
      end if

      if(neg) x(i) = -x(i)
    end do

    J(:,:) = M_ZERO
    do i = 1, geo%natoms
      r = max(sqrt(sum((x - geo%atom(i)%x)**2)), CNST(1e-6))
      f = exp(-r**2/(M_TWO*cv%Jrange**2))

      do ix = 1, conf%dim
        J(ix, ix) = M_ONE - cv%Jlocal*f
        do iy = 1, conf%dim
          J(ix, iy) = J(ix, iy) + cv%Jlocal*(x(ix)-geo%atom(i)%x(ix))*(x(iy)-geo%atom(i)%x(iy)) * &
             M_TWO/(M_TWO*cv%Jrange**2) * f
        end do
      end do
    end do

    do ix = 1, conf%dim
      J(ix,:) = J(ix,:)*J2(:)
    end do

  end subroutine curv_modine_jacobian_inv

end module curv_modine

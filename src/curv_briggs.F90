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

! This module implements the curvilinear coordinates given in
! E.L. Briggs, D.J. Sullivan, and J. Bernholc, PRB 54 14362 (1996)
!
! It assumes that the Oxiggen atom is located at x0=0 (see Eq. (12))

module curv_briggs
  use lib_oct_parser
  use global
  use messages
  use syslabels
  use units
  use geometry
  use lib_adv_alg

  implicit none

  type curv_briggs_type
    FLOAT :: L(3)   ! size of the box
    FLOAT :: beta   ! adjustable parameter between 0 and 1 that controls the degree of scaling
  end type curv_briggs_type

contains

  !-------------------------------------
  subroutine curv_briggs_init(l, cv)
    FLOAT,                   intent(in)  :: l(:)  ! l(1:conf%dim)
    type(curv_briggs_type), intent(out) :: cv

    cv%L = M_ZERO
    cv%L(1:conf%dim) = l(1:conf%dim)

    call loct_parse_float(check_inp('CurvBriggsBeta'), M_HALF, cv%beta)
    
    if(cv%beta<M_ZERO.or.cv%beta>M_ONE) then
      message(1) = 'The parameter "CurvBriggsBeta" must lie between 0 and 1.'
      call write_fatal(1)
    end if

  end subroutine curv_briggs_init


  !-------------------------------------
  subroutine curv_briggs_chi2x(cv, chi, x)
    type(curv_briggs_type), intent(in)  :: cv
    FLOAT,                  intent(in)  :: chi(:)  ! chi(conf%dim)
    FLOAT,                  intent(out) :: x(:)    ! x(conf%dim)
 
    integer :: i

    do i = 1, conf%dim
      x(i) = chi(i) - cv%L(i)*cv%beta/(M_TWO*M_PI)*   &
         sin(M_TWO*M_PI*chi(i)/cv%L(i))
    end do

  end subroutine curv_briggs_chi2x


  !-------------------------------------
  subroutine curv_briggs_jacobian_inv(cv, chi, J)
    type(curv_briggs_type), intent(in)  :: cv
    FLOAT,                  intent(in)  :: chi(:)  ! chi(conf%dim)
    FLOAT,                  intent(out) :: J(:,:)  ! J(conf%dim,conf%dim), the Jacobian
 
    integer :: i

    J(:,:) = M_ZERO
    do i = 1, conf%dim
      J(i,i) = M_ONE - cv%beta*cos(M_TWO*M_PI*chi(i)/cv%L(i))
    end do

  end subroutine curv_briggs_jacobian_inv

end module curv_briggs

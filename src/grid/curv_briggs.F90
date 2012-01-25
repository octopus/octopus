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

!> This module implements the curvilinear coordinates given in
!! E.L. Briggs, D.J. Sullivan, and J. Bernholc, PRB 54 14362 (1996)
!!
!! It assumes that the Oxygen atom is located at x0=0 (see Eq. (12))

module curv_briggs_m
  use datasets_m
  use global_m
  use parser_m
  use messages_m
  use simul_box_m

  implicit none

  private
  public ::                     &
    curv_briggs_t,              &
    curv_briggs_init,           &
    curv_briggs_copy,           &
    curv_briggs_chi2x,          &
    curv_briggs_jacobian_inv

  type curv_briggs_t
    FLOAT :: L(MAX_DIM)  !< size of the box
    FLOAT :: beta        !< adjustable parameter between 0 and 1 that controls the degree of scaling
  end type curv_briggs_t

contains

  ! ---------------------------------------------------------
  subroutine curv_briggs_init(cv, sb)
    type(curv_briggs_t), intent(out) :: cv
    type(simul_box_t),   intent(in)  :: sb

    cv%L = M_ZERO
    cv%L(1:sb%dim) = sb%lsize(1:sb%dim)

    call parse_float(datasets_check('CurvBriggsBeta'), M_HALF, cv%beta)

    if(cv%beta<M_ZERO.or.cv%beta>M_ONE) then
      message(1) = 'The parameter "CurvBriggsBeta" must lie between 0 and 1.'
      call messages_fatal(1)
    end if

  end subroutine curv_briggs_init

  ! ---------------------------------------------------------
  subroutine curv_briggs_copy(this_out, this_in)
    type(curv_briggs_t), intent(inout) :: this_out
    type(curv_briggs_t), intent(in)    :: this_in
    !
    PUSH_SUB(curv_briggs_copy)
    this_out%L=this_in%L
    this_out%beta=this_in%beta
    POP_SUB(curv_briggs_copy)
    return
  end subroutine curv_briggs_copy

  ! ---------------------------------------------------------
  subroutine curv_briggs_chi2x(sb, cv, chi, x)
    type(simul_box_t),   intent(in)  :: sb
    type(curv_briggs_t), intent(in)  :: cv
    FLOAT,               intent(in)  :: chi(:)  !< chi(sb%dim)
    FLOAT,               intent(out) :: x(:)    !< x(sb%dim)

    integer :: i

    do i = 1, sb%dim
      x(i) = chi(i) - cv%L(i)*cv%beta/(M_TWO*M_PI)*   &
        sin(M_TWO*M_PI*chi(i)/cv%L(i))
    end do

  end subroutine curv_briggs_chi2x


  ! ---------------------------------------------------------
  subroutine curv_briggs_jacobian_inv(sb, cv, chi, J)
    type(simul_box_t),   intent(in)  :: sb
    type(curv_briggs_t), intent(in)  :: cv
    FLOAT,               intent(in)  :: chi(:)  !< chi(sb%dim)
    FLOAT,               intent(out) :: J(:,:)  !< J(sb%dim,sb%dim), the Jacobian

    integer :: i

    J(:,:) = M_ZERO
    do i = 1, sb%dim
      J(i,i) = M_ONE - cv%beta*cos(M_TWO*M_PI*chi(i)/cv%L(i))
    end do

  end subroutine curv_briggs_jacobian_inv

end module curv_briggs_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

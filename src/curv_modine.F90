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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

! This module implements the curvilinear coordinates given in
! N. A. Modine, G. Zumbach, and E. Kaxiras, Phys. Rev. B 55, 10289-10301 (1997) 
!
! The local refinement was changed for a simple exponential.
! I believe that the recipe given by the authors is too complicated
! for me to sort it out.

module curv_modine_m
  use datasets_m
  use geometry_m
  use geometry_m
  use global_m
  use lib_oct_parser_m
  use messages_m
  use simul_box_m
  use units_m

  implicit none

  private
  public ::                  &
    curv_modine_t,           &
    curv_modine_init,        &
    curv_modine_end,         &
    curv_modine_chi2x,       &
    curv_modine_jacobian_inv

  type curv_modine_t
    FLOAT :: L(MAX_DIM)    ! size of the box
    FLOAT :: xbar          ! size of central flat region (in units of L)
    FLOAT :: Jbar          ! increase in density of points is 1/J

    type(geometry_t), pointer :: geo        ! shortcut to geometry information
    FLOAT,            pointer :: Jlocal(:)  ! local (around the atoms) refinement
    FLOAT,            pointer :: Jrange(:)  ! local refinement range

    integer, pointer :: chi_atoms(:,:)
  end type curv_modine_t

  integer, parameter :: qq = 3

contains

  ! ---------------------------------------------------------
  subroutine curv_modine_init(sb, geo, cv)
    type(simul_box_t),   intent(in)  :: sb
    type(geometry_t),    intent(in)  :: geo
    type(curv_modine_t), intent(out) :: cv

    call loct_parse_float(check_inp('CurvModineXBar'), M_ONE/M_THREE, cv%xbar)
    call loct_parse_float(check_inp('CurvModineJBar'), M_HALF, cv%Jbar)

    cv%L = M_ZERO
    cv%L(1:sb%dim) = sb%lsize(1:sb%dim) / cv%Jbar

    if(cv%xbar<M_ZERO.or.cv%xbar>M_ONE) then
      message(1) = 'The parameter "CurvModineXBar" must lie between 0 and 1.'
      call write_fatal(1)
    end if

    ALLOCATE(cv%Jlocal(geo%natoms), geo%natoms)
    ALLOCATE(cv%Jrange(geo%natoms), geo%natoms)
    ALLOCATE(cv%chi_atoms(sb%dim, geo%natoms), sb%dim*geo%natoms)

    ! WARNING: the reading has to be done for each atom kind
    call loct_parse_float(check_inp('CurvModineJlocal'), CNST(0.5), cv%Jlocal(1))
    call loct_parse_float(check_inp('CurvModineJrange'), M_TWO/units_inp%length%factor, cv%Jrange(1))

    if(cv%Jlocal(1)<M_ZERO.or.cv%Jlocal(1)>M_ONE) then
      message(1) = 'The parameter "CurvModineJlocal" must lie between 0 and 1.'
      call write_fatal(1)
    end if

    cv%Jlocal(:) = cv%Jlocal(1)
    cv%Jrange(:) = cv%Jrange(1)

  end subroutine curv_modine_init


  ! ---------------------------------------------------------
  subroutine curv_modine_end(cv)
    type(curv_modine_t), intent(inout) :: cv

    deallocate(cv%Jlocal, cv%Jrange)
    deallocate(cv%chi_atoms)

  end subroutine curv_modine_end


  ! ---------------------------------------------------------
  subroutine curv_modine_chi2x(sb, geo, cv, chi_, x)
    type(simul_box_t),   intent(in)  :: sb
    type(geometry_t),    intent(in)  :: geo
    type(curv_modine_t), intent(in)  :: cv
    FLOAT,               intent(in)  :: chi_(:)  ! chi_(sb%dim)
    FLOAT,               intent(out) :: x(:)     !   x (sb%dim)

    FLOAT :: chibar(MAX_DIM), r, chi
    logical :: neg
    integer :: i

    chibar(1:sb%dim) = cv%xbar*cv%L(1:sb%dim)

    do i = 1, sb%dim
      neg = (chi_(i) < 0)
      chi = abs(chi_(i))

      x(i) = cv%Jbar * chi
      if(chi > chibar(i)) then
        r = (chi-chibar(i))/(cv%L(i)-chibar(i))
        x(i) = x(i) + cv%L(i)/M_TWO*(1-cv%Jbar) * r**qq *   &
           (qq + M_ONE - (qq - M_ONE)*r)
      end if

      if(neg) x(i) = -x(i)
    end do

    do i = 1, geo%natoms
      r = max(sqrt(sum((x - geo%atom(i)%x)**2)), CNST(1e-6))
      x(:) = x(:) - cv%Jlocal(i)*(x - geo%atom(i)%x)*exp(-r**2/(M_TWO*cv%Jrange(i)**2))
    end do

  end subroutine curv_modine_chi2x


  ! ---------------------------------------------------------
  subroutine curv_modine_jacobian_inv(sb, geo, cv, chi_, J)
    type(simul_box_t),   intent(in)  :: sb
    type(geometry_t),    intent(in)  :: geo
    type(curv_modine_t), intent(in)  :: cv
    FLOAT,               intent(in)  :: chi_(:)  ! chi(sb%dim)
    FLOAT,               intent(out) :: J(:,:)   ! J(sb%dim,sb%dim), the Jacobian

    FLOAT :: chibar(MAX_DIM), r, f, chi, J2(MAX_DIM), x(MAX_DIM)
    logical :: neg
    integer :: i, ix, iy

    chibar(1:sb%dim) = cv%xbar*cv%L(1:sb%dim)

    J2 = M_ZERO
    do i = 1, sb%dim
      neg = (chi_(i) < 0)
      chi = abs(chi_(i))

      x(i)  = cv%Jbar * chi
      J2(i) = cv%Jbar

      if(chi > chibar(i)) then
        r = (chi-chibar(i))/(cv%L(i)-chibar(i))

        x(i)  = x(i) + cv%L(i)/M_TWO*(1-cv%Jbar) * r**qq *   &
           (qq + M_ONE - (qq - M_ONE)*r)

        J2(i) = J2(i) + cv%L(i)/M_TWO*(1-cv%Jbar) * r**(qq-1)/(cv%L(i)-chibar(i)) *   &
           (qq*(qq+1) - (qq**2-1)*r)
      end if

      if(neg) x(i) = -x(i)
    end do

    J(:,:) = M_ZERO
    do i = 1, geo%natoms
      r = max(sqrt(sum((x - geo%atom(i)%x)**2)), CNST(1e-6))
      f = exp(-r**2/(M_TWO*cv%Jrange(i)**2))

      do ix = 1, sb%dim
        J(ix, ix) = M_ONE - cv%Jlocal(i)*f
        do iy = 1, sb%dim
          J(ix, iy) = J(ix, iy) + cv%Jlocal(i)*(x(ix)-geo%atom(i)%x(ix))*(x(iy)-geo%atom(i)%x(iy)) * &
             M_TWO/(M_TWO*cv%Jrange(i)**2) * f
        end do
      end do
    end do

    do ix = 1, sb%dim
      J(ix,:) = J(ix,:)*J2(:)
    end do

  end subroutine curv_modine_jacobian_inv

end module curv_modine_m

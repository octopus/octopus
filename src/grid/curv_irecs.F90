!! Copyright (C) 2015 U. De Giovannini
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

!> This module implements the infinite range complex scaling 
!! coordinate transformation as of 
!! M. Weinmuller, M. Weinmuller, J. Rohland, and A. Scrinzi, arXiv:1509.04947 (2015).

module curv_irecs_oct_m
  use geometry_oct_m
  use global_oct_m
  use loct_pointer_oct_m
  use parser_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                   &
    curv_irecs_t,              &
    curv_irecs_init,           &
    curv_irecs_copy,           &
    curv_irecs_end,            &
    curv_irecs_chi2x,          &
    curv_irecs_x2chi,          &
    curv_irecs_jacobian

  type curv_irecs_t
    FLOAT :: R0            !< the scaling region radius
    FLOAT :: theta         !< the scaling angle in complex plane
    FLOAT :: gamma         !< the exponential samlping parameter 
  end type curv_irecs_t


contains

  ! ---------------------------------------------------------
  subroutine curv_irecs_init(cv, sb)
    type(curv_irecs_t), intent(out) :: cv
    type(simul_box_t),  intent(in)  :: sb

    integer :: ipos, idir

    PUSH_SUB(curv_irecs_init)


    !%Variable CurvIrECSRadius
    !%Type float
    !%Default 10 
    !%Section Mesh::Curvilinear::irECS
    !%Description
    !% The complex scaling radius identifying the value after which 
    !% the scaling transformation is applied 
    !%End
    call parse_variable('CurvIrECSRadius', CNST(10.), cv%R0)

    !%Variable CurvIrECSTheta
    !%Type float
    !%Default 0.3
    !%Section Mesh::Curvilinear::irECS
    !%Description
    !% The complex scaling parameter. 
    !% It will have to be merged with the complex scaling module.
    !%End
    call parse_variable('CurvIrECSTheta', CNST(0.3), cv%theta)

    !%Variable CurvIrECSGamma
    !%Type float
    !%Default 2.0 a.u.
    !%Section Mesh::Curvilinear::irECS
    !%Description
    !% This number determines the region over which the grid is enhanced (range of
    !%End

    call parse_variable('CurvIrECSGamma', CNST(1.0), cv%gamma, unit_one/units_inp%length)

    
    POP_SUB(curv_irecs_init)
  end subroutine curv_irecs_init

  ! ---------------------------------------------------------
  subroutine curv_irecs_copy(this_out, this_in)
    type(curv_irecs_t), intent(inout) :: this_out
    type(curv_irecs_t), intent(in)    :: this_in
    !
    PUSH_SUB(curv_irecs_copy)
    this_out%theta=this_in%theta
    this_out%gamma=this_in%gamma
    this_out%R0=this_in%R0
    POP_SUB(curv_irecs_copy)
    return
  end subroutine curv_irecs_copy

  ! ---------------------------------------------------------
  subroutine curv_irecs_end(cv)
    type(curv_irecs_t), intent(inout) :: cv

    PUSH_SUB(curv_irecs_end)


    POP_SUB(curv_irecs_end)
  end subroutine curv_irecs_end



  ! ---------------------------------------------------------
  subroutine curv_irecs_chi2x(sb, cv, chi, x, phx)
    type(simul_box_t), target, intent(in)  :: sb
    type(curv_irecs_t), target, intent(in)  :: cv
    FLOAT,                     intent(in)  :: chi(:)  !< chi(sb%dim)
    FLOAT,                     intent(out) :: x(:)    !< x(sb%dim)
    CMPLX, optional,           intent(out) :: phx(:)  !< the complex phase of x

    integer :: i
    FLOAT :: r

    ! no push_sub, called too frequently
    
    r = sqrt(sum(chi(1:sb%dim)**2))
    
    if (r < cv%R0) then
      x(sb%dim) =  chi(sb%dim)
      if (present(phx)) phx(sb%dim) = M_ONE 
    else
      do i=1, sb%dim
        x(i)   = chi(i)/(r*cv%gamma) * exp(cv%gamma * (r - cv%R0))
! only real scaling for now        
!         x(i)   = chi(i)/(r*cv%gamma) * exp(cv%gamma * (r - cv%R0) + M_ZI*cv%theta)
        if (present(phx)) phx(i) = exp(M_ZI*cv%theta) 
      end do
    end if

  end subroutine curv_irecs_chi2x


  ! ---------------------------------------------------------
  subroutine curv_irecs_x2chi(sb, cv, x, chi)
    type(simul_box_t), intent(in)  :: sb
    type(curv_irecs_t), intent(in)  :: cv
    FLOAT,             intent(in)  :: x(:)    ! x(sb%dim)
    FLOAT,             intent(out) :: chi(:)  ! chi(sb%dim)

    integer :: i, ia
    FLOAT   :: r, ar, th, ex

    PUSH_SUB(curv_irecs_x2chi)


    POP_SUB(curv_irecs_x2chi)
  end subroutine curv_irecs_x2chi


  ! ---------------------------------------------------------
  subroutine curv_irecs_jacobian(sb, cv, x, chi, J)
    type(simul_box_t), intent(in)  :: sb
    type(curv_irecs_t), intent(in)  :: cv
    FLOAT,             intent(in)  :: x(:)    !< x(sb%dim)
    FLOAT,             intent(out) :: chi(:)  !< chi(sb%dim)
    FLOAT,             intent(out) :: J(:,:)  !< J(sb%dim,sb%dim), the Jacobian

    integer :: i, ix, iy, natoms_
    FLOAT :: r, f_alpha, df_alpha
    FLOAT :: th, ex, ar

    ! no push_sub, called too frequently

    do i =1, sb%dim
      J(i,i) = M_ONE/sb%dim 
    end do

  end subroutine curv_irecs_jacobian

end module curv_irecs_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

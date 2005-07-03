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

module curvlinear
  use lib_oct_parser
  use simul_box
  use geometry
  use curv_gygi
  use curv_briggs
  use curv_modine

  implicit none

  type curvlinear_type
    integer :: method

    type(curv_gygi_type)   :: gygi
    type(curv_briggs_type) :: briggs
    type(curv_modine_type) :: modine
  end type curvlinear_type

  integer, parameter, public :: &
     CURV_METHOD_UNIFORM = 1,  &
     CURV_METHOD_GYGI    = 2,  &
     CURV_METHOD_BRIGGS  = 3,  &
     CURV_METHOD_MODINE  = 4

contains

  !-------------------------------------
  subroutine curvlinear_init(sb, cv)
    type(simul_box_type),  intent(in)  :: sb
    type(curvlinear_type), intent(out) :: cv

    call push_sub('curvlinear_init')

    call loct_parse_int(check_inp('CurvMethod'), CURV_METHOD_UNIFORM, cv%method)
    if(cv%method<CURV_METHOD_UNIFORM.or.cv%method>CURV_METHOD_MODINE) then
      write(message(1), '(a,i2,a)') 'Do not have a "CurvMethod = ', cv%method, '"'
      call write_fatal(1)
    end if

    select case(cv%method)
    case(CURV_METHOD_GYGI)
      call curv_gygi_init(cv%gygi)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_init(sb, cv%briggs)
    case(CURV_METHOD_MODINE)
      call curv_modine_init(sb, cv%modine)
    end select

    call pop_sub()
  end subroutine curvlinear_init


  !-------------------------------------
  subroutine curvlinear_chi2x(sb, geo, cv, chi, x)
    type(simul_box_type),  intent(in)  :: sb
    type(geometry_type),   intent(in)  :: geo
    type(curvlinear_type), intent(in)  :: cv
    FLOAT,                 intent(in)  :: chi(:)  ! chi(conf%dim)
    FLOAT,                 intent(out) :: x(:)    ! x(conf%dim)

    select case(cv%method)
    case(CURV_METHOD_UNIFORM)
      x(1:sb%dim) = chi(1:sb%dim)
    case(CURV_METHOD_GYGI)
      call curv_gygi_chi2x(sb, geo, cv%gygi, chi, x)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_chi2x(sb, cv%briggs, chi, x)
    case(CURV_METHOD_MODINE)
      call curv_modine_chi2x(sb, geo, cv%modine, chi, x)
    end select

  end subroutine curvlinear_chi2x
  

  !-------------------------------------
  FLOAT function curvlinear_det_Jac(sb, geo, cv, x, chi) result(jdet)
    type(simul_box_type),  intent(in)  :: sb
    type(geometry_type),   intent(in)  :: geo
    type(curvlinear_type), intent(in)  :: cv
    FLOAT,                 intent(in)  :: x(:)    !   x(sb%dim)
    FLOAT,                 intent(in)  :: chi(:)  ! chi(sb%dim)

    FLOAT :: dummy(sb%dim), Jac(sb%dim, sb%dim)
    integer :: i

    select case(cv%method)
    case(CURV_METHOD_UNIFORM)
      jdet = M_ONE
    case(CURV_METHOD_GYGI)
      call curv_gygi_jacobian(sb, geo, cv%gygi, x, dummy, Jac)
      jdet = M_ONE/lalg_det(Jac, sb%dim)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_jacobian_inv(sb, cv%briggs, chi, Jac)
      jdet = M_ONE
      do i = 1, sb%dim
        jdet = jdet * Jac(i,i) ! Jacobian is diagonal in this method
      end do
    case(CURV_METHOD_MODINE)
      call curv_modine_jacobian_inv(sb, geo, cv%modine, chi, Jac)
      jdet = M_ONE*lalg_det(Jac, sb%dim)
    end select

  end function curvlinear_det_Jac

end module curvlinear

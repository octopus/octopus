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

module curvilinear_m
  use curv_briggs_m
  use curv_gygi_m
  use curv_modine_m
  use datasets_m
  use geometry_m
  use global_m
  use lalg_adv_m
  use parser_m
  use math_m
  use messages_m
  use profiling_m
  use simul_box_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public ::                     &
    curvilinear_t,               &
    curvilinear_init,            &
    curvilinear_end,             &
    curvilinear_chi2x,           &
    curvilinear_x2chi,           &
    curvilinear_det_Jac,         &
    curvilinear_write_info

  integer, parameter, public :: &
    CURV_METHOD_UNIFORM = 1,    &
    CURV_METHOD_GYGI    = 2,    &
    CURV_METHOD_BRIGGS  = 3,    &
    CURV_METHOD_MODINE  = 4

  type curvilinear_t
    integer :: method
    type(curv_gygi_t)   :: gygi
    type(curv_briggs_t) :: briggs
    type(curv_modine_t) :: modine
  end type curvilinear_t

  character(len=23), parameter :: dump_tag = '*** curvilinear_dump **'

contains

  ! ---------------------------------------------------------
  subroutine curvilinear_init(cv, sb, geo, spacing)
    type(curvilinear_t), intent(out) :: cv
    type(simul_box_t),   intent(in)  :: sb
    type(geometry_t),    intent(in)  :: geo
    FLOAT,               intent(in)  :: spacing(:)

    PUSH_SUB(curvilinear_init)

    !%Variable CurvMethod
    !%Type integer
    !%Default curv_uniform
    !%Section Mesh::Curvilinear
    !%Description
    !% The relevant functions in octopus are represented on a mesh in real space.
    !% This mesh may be an evenly spaced regular rectangular grid (standard mode),
    !% or else an *adaptive* or *curvilinear grid*. We have implemented
    !% three kinds of adaptive meshes, although only one is currently working,
    !% the one invented by F. Gygi (<tt>curv_gygi</tt>). The code will stop if any of
    !% the other two is invoked.
    !%Option curv_uniform 1
    !% Regular, uniform rectangular grid. The default.
    !%Option curv_gygi 2
    !% The deformation of the grid is done according to the scheme described by
    !% F. Gygi [F. Gygi and G. Galli, <i>Phys. Rev. B</i> <b>52</b>, R2229 (1995)].
    !%Option curv_briggs 3
    !% The deformation of the grid is done according to the scheme described by
    !% Briggs [E.L. Briggs, D.J. Sullivan, and J. Bernholc, <i>Phys. Rev. B</i> <b>54</b> 14362 (1996)]
    !% (NOT WORKING).
    !%Option curv_modine 4
    !% The deformation of the grid is done according to the scheme described by
    !% Modine [N.A. Modine, G. Zumbach and E. Kaxiras, <i>Phys. Rev. B</i> <b>55</b>, 10289 (1997)]
    !% (NOT WORKING).
    !%End
    call parse_integer(datasets_check('CurvMethod'), CURV_METHOD_UNIFORM, cv%method)
    if(.not.varinfo_valid_option('CurvMethod', cv%method)) call input_error('CurvMethod')
    call messages_print_var_option(stdout, "CurvMethod", cv%method)

    ! FIXME: The other two methods are apparently not working
    if(cv%method > CURV_METHOD_GYGI) call messages_experimental('Selected curvilinear coordinates method')

    select case(cv%method)
    case(CURV_METHOD_GYGI)
      call curv_gygi_init(cv%gygi, sb, geo)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_init(cv%briggs, sb)
    case(CURV_METHOD_MODINE)
      call curv_modine_init(cv%modine, sb, geo, spacing)
    end select

    POP_SUB(curvilinear_init)
  end subroutine curvilinear_init


  ! ---------------------------------------------------------
  subroutine curvilinear_end(cv)
    type(curvilinear_t), intent(inout) :: cv

    PUSH_SUB(curvilinear_end)

    select case(cv%method)
    case(CURV_METHOD_GYGI)
      call curv_gygi_end(cv%gygi)
    case(CURV_METHOD_BRIGGS)
    case(CURV_METHOD_MODINE)
      call curv_modine_end(cv%modine)
    end select

    POP_SUB(curvilinear_end)
  end subroutine curvilinear_end
  

  ! ---------------------------------------------------------
  subroutine curvilinear_chi2x(sb, cv, chi, x)
    type(simul_box_t),   intent(in)  :: sb
    type(curvilinear_t), intent(in)  :: cv
    FLOAT,               intent(in)  :: chi(MAX_DIM)  ! chi(conf%dim)
    FLOAT,               intent(out) :: x(MAX_DIM)    ! x(conf%dim)

    ! no push_sub because called too frequently
    x = M_ZERO

    select case(cv%method)
    case(CURV_METHOD_UNIFORM)
      x(1:sb%dim) = matmul(sb%rlattice_primitive(1:sb%dim,1:sb%dim), chi(1:sb%dim))
    case(CURV_METHOD_GYGI)
      call curv_gygi_chi2x(sb, cv%gygi, chi, x)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_chi2x(sb, cv%briggs, chi, x)
    case(CURV_METHOD_MODINE)
      call curv_modine_chi2x(sb, cv%modine, chi, x)
    end select

  end subroutine curvilinear_chi2x


  ! ---------------------------------------------------------
  subroutine curvilinear_x2chi(sb, cv, x, chi)
    type(simul_box_t),   intent(in)  :: sb
    type(curvilinear_t), intent(in)  :: cv
    FLOAT,               intent(in)  :: x(MAX_DIM)
    FLOAT,               intent(out) :: chi(MAX_DIM)

    PUSH_SUB(curvilinear_x2chi)
    
    chi = M_ZERO

    select case(cv%method)
    case(CURV_METHOD_UNIFORM)
      chi = matmul(x, sb%klattice_primitive)
    case(CURV_METHOD_GYGI)
      call curv_gygi_x2chi(sb, cv%gygi, x, chi)
    case(CURV_METHOD_BRIGGS, CURV_METHOD_MODINE)
      message(1) = "Internal error in curvilinear_x2chi"
      call messages_fatal(1)
    end select

    POP_SUB(curvilinear_x2chi)
  end subroutine curvilinear_x2chi


  ! ---------------------------------------------------------
  FLOAT function curvilinear_det_Jac(sb, cv, x, chi) result(jdet)
    type(simul_box_t),   intent(in)  :: sb
    type(curvilinear_t), intent(in)  :: cv
    FLOAT,               intent(in)  :: x(:)    !   x(sb%dim)
    FLOAT,               intent(in)  :: chi(:)  ! chi(sb%dim)

    FLOAT :: dummy(MAX_DIM)
    FLOAT, allocatable :: Jac(:,:)
    integer :: i

    ! No PUSH_SUB, called too often

    if(cv%method.ne.CURV_METHOD_UNIFORM) then
      SAFE_ALLOCATE(Jac(1:sb%dim, 1:sb%dim))
    end if

    select case(cv%method)
    case(CURV_METHOD_UNIFORM)
      jdet = sb%volume_element
    case(CURV_METHOD_GYGI)
      call curv_gygi_jacobian(sb, cv%gygi, x, dummy, Jac)
      jdet = M_ONE/lalg_determinant(sb%dim, Jac, invert = .false.)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_jacobian_inv(sb, cv%briggs, chi, Jac)
      jdet = M_ONE
      do i = 1, sb%dim
        jdet = jdet * Jac(i,i) ! Jacobian is diagonal in this method
      end do
    case(CURV_METHOD_MODINE)
      call curv_modine_jacobian_inv(sb, cv%modine, chi, dummy, Jac)
      jdet = M_ONE*lalg_determinant(sb%dim, Jac, invert = .false.)
    end select

    if(cv%method.ne.CURV_METHOD_UNIFORM) then
      SAFE_DEALLOCATE_A(Jac)
    end if

  end function curvilinear_det_Jac

  ! ---------------------------------------------------------
  subroutine curvilinear_write_info(cv, unit)
    type(curvilinear_t), intent(in) :: cv
    integer,            intent(in) :: unit

    PUSH_SUB(curvilinear_write_info)

    select case(cv%method)
    case(CURV_METHOD_GYGI)
      write(message(1), '(a)')  '  Curvilinear Method = gygi'
      write(message(2), '(a)')  '  Gygi Parameters:'
      write(message(3), '(4x,a,f6.3)')  'A = ', cv%gygi%a
      write(message(4), '(4x,3a,f6.3)') 'alpha [', &
                                      trim(units_abbrev(units_out%length)), '] = ', &
                                      units_from_atomic(units_out%length, cv%gygi%alpha)
      write(message(5), '(4x,3a,f6.3)') 'beta  [', &
                                      trim(units_abbrev(units_out%length)), '] = ', &
                                      units_from_atomic(units_out%length, cv%gygi%beta)
      call messages_info(5, unit)

    case(CURV_METHOD_BRIGGS)
      write(message(1), '(a)') '  Curvilinear Method = briggs'
      call messages_info(1, unit)

    case(CURV_METHOD_MODINE)
      write(message(1), '(a)') ' Curvilinear  Method = modine'
      call messages_info(1, unit)

    end select

    POP_SUB(curvilinear_write_info)
  end subroutine curvilinear_write_info

  ! ---------------------------------------------------------

end module curvilinear_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

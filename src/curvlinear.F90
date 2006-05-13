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

module curvlinear_m
  use lib_oct_parser_m
  use simul_box_m
  use geometry_m
  use curv_gygi_m
  use curv_briggs_m
  use curv_modine_m
  use varinfo_m
  use datasets_m
  use math_m

  implicit none

  private
  public ::                     &
    curvlinear_t,               &
    curvlinear_init,            &
    curvlinear_chi2x,           &
    curvlinear_det_Jac,         &
    curvlinear_write_info,      &
    curvlinear_dump,            &
    curvlinear_init_from_file,  &
    operator(.eq.)

  integer, parameter, public :: &
    CURV_METHOD_UNIFORM = 1,    &
    CURV_METHOD_GYGI    = 2,    &
    CURV_METHOD_BRIGGS  = 3,    &
    CURV_METHOD_MODINE  = 4

  type curvlinear_t
    integer :: method
    type(curv_gygi_t)   :: gygi
    type(curv_briggs_t) :: briggs
    type(curv_modine_t) :: modine
  end type curvlinear_t

  interface operator(.eq.)
    module procedure curv_is_equal
  end interface


contains

  ! ---------------------------------------------------------
  subroutine curvlinear_init(sb, cv)
    type(simul_box_t),  intent(in)  :: sb
    type(curvlinear_t), intent(out) :: cv

    call push_sub('curvlinear.curvlinear_init')

    !%Variable CurvMethod
    !%Type integer
    !%Default curv_uniform
    !%Section Mesh::Curvilinear
    !%Description
    !% The relevant functions in octopus are represented on a mesh in real space.
    !% This mesh may be an evenly spaced regular rectangular grid (standard mode),
    !% or else an *adaptive* or *curvilinear grid*. We have implemented
    !% three kinds of adative meshes, although only one is currently working,
    !% the one invented by F. Gygi ("curv_gygi"). The code will stop if any of
    !% the other two is invoked.
    !%Option curv_uniform 1
    !% Regular, uniform rectangular grid. The default.
    !%Option curv_gygi 2
    !% The deformation of the grid is done according to the scheme described by
    !% F. Gygi [F. Gygi and G. Galli, Phys. Rev. B 52, R2229 (1995)]
    !%Option curv_briggs 3
    !% The deformation of the grid is done according to the scheme described by
    !% Briggs [E.L. Briggs, D.J. Sullivan, and J. Bernholc, Phys. Rev. B 54 14362 (1996)]
    !% (NOT WORKING)
    !%Option curv_modine 4
    !% The deformation of the grid is done according to the scheme described by
    !% Modine [N.A. Modine, G. Zumbach and E. Kaxiras, Phys. Rev. B 55, 10289 (1997)]
    !% (NOT WORKING)
    !%End
    call loct_parse_int(check_inp('CurvMethod'), CURV_METHOD_UNIFORM, cv%method)
    if(.not.varinfo_valid_option('CurvMethod', cv%method)) call input_error('CurvMethod')
    call messages_print_var_option(stdout, "CurvMethod", cv%method)

    ! FIXME: The other two methods are apparently not working
    if(cv%method > CURV_METHOD_GYGI.and.(.not.conf%devel_version)) call input_error('CurvMethod')

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


  ! ---------------------------------------------------------
  subroutine curvlinear_chi2x(sb, geo, cv, chi, x)
    type(simul_box_t),  intent(in)  :: sb
    type(geometry_t),   intent(in)  :: geo
    type(curvlinear_t), intent(in)  :: cv
    FLOAT,              intent(in)  :: chi(MAX_DIM)  ! chi(conf%dim)
    FLOAT,              intent(out) :: x(MAX_DIM)    ! x(conf%dim)

    select case(cv%method)
    case(CURV_METHOD_UNIFORM)
      x = chi
    case(CURV_METHOD_GYGI)
      call curv_gygi_chi2x(sb, geo, cv%gygi, chi, x)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_chi2x(sb, cv%briggs, chi, x)
    case(CURV_METHOD_MODINE)
      call curv_modine_chi2x(sb, geo, cv%modine, chi, x)
    end select

  end subroutine curvlinear_chi2x


  ! ---------------------------------------------------------
  FLOAT function curvlinear_det_Jac(sb, geo, cv, x, chi) result(jdet)
    type(simul_box_t),  intent(in)  :: sb
    type(geometry_t),   intent(in)  :: geo
    type(curvlinear_t), intent(in)  :: cv
    FLOAT,              intent(in)  :: x(:)    !   x(sb%dim)
    FLOAT,              intent(in)  :: chi(:)  ! chi(sb%dim)

    FLOAT :: dummy(MAX_DIM), Jac(MAX_DIM, MAX_DIM)
    integer :: i

    select case(cv%method)
    case(CURV_METHOD_UNIFORM)
      jdet = M_ONE
    case(CURV_METHOD_GYGI)
      call curv_gygi_jacobian(sb, geo, cv%gygi, x, dummy(1:sb%dim), Jac(1:sb%dim, 1:sb%dim))
      jdet = M_ONE/lalg_determinant(sb%dim, Jac, invert = .false.)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_jacobian_inv(sb, cv%briggs, chi, Jac(1:sb%dim, 1:sb%dim))
      jdet = M_ONE
      do i = 1, sb%dim
        jdet = jdet * Jac(i,i) ! Jacobian is diagonal in this method
      end do
    case(CURV_METHOD_MODINE)
      call curv_modine_jacobian_inv(sb, geo, cv%modine, chi, Jac(1:sb%dim, 1:sb%dim))
      jdet = M_ONE*lalg_determinant(sb%dim, Jac(1:sb%dim, 1:sb%dim), invert = .false.)
    end select

  end function curvlinear_det_Jac

  ! ---------------------------------------------------------
  subroutine curvlinear_write_info(cv, unit)
    type(curvlinear_t), intent(in) :: cv
    integer,            intent(in) :: unit

    call push_sub('curvlinear.curvlinear_write_info')

    select case(cv%method)
    case(CURV_METHOD_GYGI)
      write(message(1), '(a)')  '  Curvilinear Method = gygi'
      write(message(2), '(a)')  '  Gygi Parameters:'
      write(message(3), '(4x,a,f6.3)')  'A = ', cv%gygi%a
      write(message(4), '(4x,3a,f6.3)') 'alpha [', &
                                      trim(units_out%length%abbrev), '] = ', &
                                      cv%gygi%alpha/units_out%length%factor
      write(message(5), '(4x,3a,f6.3)') 'beta  [', &
                                      trim(units_out%length%abbrev), '] = ', &
                                      cv%gygi%beta/units_out%length%factor
      call write_info(5, unit)

    case(CURV_METHOD_BRIGGS)
      write(message(1), '(a)') '  Curvilinear Method = briggs'
      call write_info(1, unit)

    case(CURV_METHOD_MODINE)
      write(message(1), '(a)') ' Curvilinear  Method = modine'
      call write_info(1, unit)

    end select

    call pop_sub()
  end subroutine curvlinear_write_info

  ! ---------------------------------------------------------
  subroutine curvlinear_dump(cv, iunit)
    type(curvlinear_t), intent(in) :: cv
    integer,            intent(in) :: iunit

    write(iunit, '(a20,i4)')        'method=              ', cv%method
    select case(cv%method)
    case(CURV_METHOD_GYGI)
      write(iunit, '(a20,e22.14)')  'a=                   ', cv%gygi%a
      write(iunit, '(a20,e22.14)')  'alpha=               ', cv%gygi%alpha
      write(iunit, '(a20,e22.14)')  'beta=                ', cv%gygi%beta
    case(CURV_METHOD_BRIGGS)
      write(iunit, '(a20,3e22.14)') 'l=                   ', cv%briggs%l(1:3)
      write(iunit, '(a20,e22.14)')  'beta=                ', cv%briggs%beta
    case(CURV_METHOD_MODINE)
      write(iunit, '(a20,3e22.14)') 'l=                   ', cv%modine%l(1:3)
      write(iunit, '(a20,e22.14)')  'xbar=                ', cv%modine%xbar
      write(iunit, '(a20,e22.14)')  'jbar=                ', cv%modine%jbar
      write(iunit, '(a20,e22.14)')  'jlocal=              ', cv%modine%jlocal
      write(iunit, '(a20,e22.14)')  'jrange=              ', cv%modine%jrange
    end select

  end subroutine curvlinear_dump


  ! ---------------------------------------------------------
  subroutine curvlinear_init_from_file(cv, iunit)
    type(curvlinear_t), intent(out) :: cv
    integer,            intent(in)  :: iunit

    character(len=20) :: str

    read(iunit, *) str, cv%method
    select case(cv%method)
    case(CURV_METHOD_GYGI)
      read(iunit, *) str, cv%gygi%a
      read(iunit, *) str, cv%gygi%alpha
      read(iunit, *) str, cv%gygi%beta
    case(CURV_METHOD_BRIGGS)
      read(iunit, *) str, cv%briggs%l(1:3)
      read(iunit, *) str, cv%briggs%beta
    case(CURV_METHOD_MODINE)
      read(iunit, *) str, cv%modine%l(1:3)
      read(iunit, *) str, cv%modine%xbar
      read(iunit, *) str, cv%modine%jbar
      read(iunit, *) str, cv%modine%jlocal
      read(iunit, *) str, cv%modine%jrange
    end select

  end subroutine curvlinear_init_from_file


  ! ---------------------------------------------------------
  logical function curv_is_equal(cv1, cv2) result(res)
    type(curvlinear_t), intent(in) :: cv1, cv2

    res = .false.
    if(cv1%method .ne. cv2%method)                           return
    select case(cv1%method)
    case(CURV_METHOD_GYGI)
      if(.not. (cv1%gygi%a .app.cv2%gygi%a) )                return
      if(.not. (cv1%gygi%alpha .app.cv2%gygi%alpha) )        return
      if(.not. (cv1%gygi%beta .app.cv2%gygi%beta) )          return
    case(CURV_METHOD_BRIGGS)
      if(.not. (cv1%briggs%l .app.cv2%briggs%l) )            return
      if(.not. (cv1%briggs%beta .app.cv2%briggs%beta) )      return
    case(CURV_METHOD_MODINE)
      if(.not. (cv1%modine%l .app. cv2%modine%l) )           return
      if(.not. (cv1%modine%xbar .app. cv2%modine%xbar) )     return
      if(.not. (cv1%modine%jbar .app. cv2%modine%jbar) )     return
      if(.not. (cv1%modine%jlocal .app. cv2%modine%jlocal) ) return
      if(.not. (cv1%modine%jrange .app. cv2%modine%jrange) ) return
    end select
    res = .true.

  end function curv_is_equal

end module curvlinear_m

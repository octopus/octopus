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

! This module implements the curvilinear coordinates given in
! F. Gygi and G. Galli, PRB 52 R2229 (1996).

module curv_gygi_m
  use datasets_m
  use geometry_m
  use global_m
  use parser_m
  use messages_m
  use profiling_m
  use root_solver_m
  use simul_box_m
  use unit_m
  use unit_system_m

  implicit none

  private
  public ::            &
    curv_gygi_t,       &
    curv_gygi_init,    &
    curv_gygi_end,     &
    curv_gygi_chi2x,   &
    curv_gygi_x2chi,   &
    curv_gygi_jacobian

  type curv_gygi_t
    FLOAT :: A             ! local reduction in grid spacing is 1/(1+A)
    FLOAT :: alpha         ! range of enhancement of the resolution
    FLOAT :: beta          ! distance over which Euclidian coordinates are recovered
    FLOAT, pointer :: pos(:, :)
    integer :: npos
  end type curv_gygi_t

  type(simul_box_t), pointer  :: sb_p
  type(curv_gygi_t), pointer  :: cv_p
  integer :: i_p
  FLOAT :: chi_p(MAX_DIM)

contains

  ! ---------------------------------------------------------
  subroutine curv_gygi_init(cv, sb, geo)
    type(curv_gygi_t), intent(out) :: cv
    type(simul_box_t), intent(in)  :: sb
    type(geometry_t),  intent(in)  :: geo

    integer :: ipos, idir

    PUSH_SUB(curv_gygi_init)

    !%Variable CurvGygiA
    !%Type float
    !%Section Mesh::Curvilinear::Gygi
    !%Description
    !% The grid spacing is reduced locally around each atom, and the reduction is
    !% given by 1/(1+A), where A is specified by this variable, CurvGygiA. So, if
    !% A=1/2 (the default), the grid spacing is reduced to two thirds = 1/(1+1/2).
    !% [This is the <math>A_{\alpha}</math> variable in Eq. 2 of F. Gygi and G. Galli, <i>Phys.
    !% Rev. B</i> <b>52</b>, R2229 (1995)]
    !% It must be larger than zero.
    !%End
    call parse_float(datasets_check('CurvGygiA'), M_HALF, cv%A)
    !%Variable CurvGygiAlpha
    !%Type float
    !%Section Mesh::Curvilinear::Gygi
    !%Description
    !% This number determines the region over which the grid is enhanced (range of
    !% enhancement of the resolution). That is, the grid is enhanced on a sphere
    !% around each atom, whose radius is given by this variable. [This is the <math>a_{\alpha}</math>
    !% variable in Eq. 2 of F. Gygi and G. Galli, <i>Phys. Rev. B</i> <b>52M</b>, R2229 (1995)].
    !% The default is two atomic units.
    !% It must be larger than zero.
    !%End
    call parse_float(datasets_check('CurvGygiAlpha'), units_from_atomic(units_inp%length, M_TWO), cv%alpha)
    !%Variable CurvGygiBeta
    !%Type float
    !%Section Mesh::Curvilinear::Gygi
    !%Description
    !% This number determines the distance over which Euclidean coordinates are
    !% recovered. [This is the <math>b_{\alpha}</math> variable in Eq. 2 of F. Gygi and G. Galli,
    !% <i>Phys. Rev. B</i> <b>52</b>, R2229 (1995)]. The default is four atomic units.
    !% It must be larger than zero.
    !%End
    call parse_float(datasets_check('CurvGygiBeta'),  units_from_atomic(units_inp%length, M_FOUR), cv%beta)

    if(cv%a<=M_ZERO)     call input_error('CurvGygiA')
    if(cv%alpha<=M_ZERO) call input_error('CurvGygiAlpha')
    if(cv%beta<=M_ZERO)  call input_error('CurvGygiBeta')

    cv%alpha = units_to_atomic(units_inp%length, cv%alpha)
    cv%beta  = units_to_atomic(units_inp%length, cv%beta)

    cv%npos = geo%natoms
    SAFE_ALLOCATE(cv%pos(1:cv%npos, 1:sb%dim))
    forall(ipos = 1:cv%npos, idir = 1:sb%dim) cv%pos(ipos, idir) = geo%atom(ipos)%x(idir)
    
    POP_SUB(curv_gygi_init)
  end subroutine curv_gygi_init

  ! ---------------------------------------------------------
  subroutine curv_gygi_end(cv)
    type(curv_gygi_t), intent(inout) :: cv

    PUSH_SUB(curv_gygi_end)

    SAFE_DEALLOCATE_P(cv%pos)

    POP_SUB(curv_gygi_end)
  end subroutine curv_gygi_end

  ! ---------------------------------------------------------
  subroutine getf(y, f, jf)
    FLOAT, intent(in)    :: y(:)
    FLOAT, intent(out)   :: f(:), jf(:, :)

    PUSH_SUB(getf)

    call curv_gygi_jacobian(sb_p, cv_p, y, f, jf, i_p)
    f(1:sb_p%dim) = f(1:sb_p%dim) - chi_p(1:sb_p%dim)

    POP_SUB(getf)
  end subroutine getf 


  ! ---------------------------------------------------------
  subroutine curv_gygi_chi2x(sb, cv, chi, x)
    type(simul_box_t), target, intent(in)  :: sb
    type(curv_gygi_t), target, intent(in)  :: cv
    FLOAT,                     intent(in)  :: chi(:)  ! chi(sb%dim)
    FLOAT,                     intent(out) :: x(:)    ! x(sb%dim)

    ! local variables
    integer :: i
    logical :: conv
    type(root_solver_t) :: rs

    PUSH_SUB(curv_gygi_chi2x)

    call root_solver_init(rs, sb%dim,  &
      solver_type = ROOT_NEWTON, maxiter = 500, abs_tolerance = CNST(1.0e-10))

    sb_p            => sb
    cv_p            => cv
    i_p             =  cv%npos
    chi_p(1:sb%dim) =  chi(1:sb%dim)

    call droot_solver_run(rs, getf, x, conv, startval = chi)

    if(.not.conv) then
      do i = 1, cv%npos
        conv = .false.
        i_p = i
        call droot_solver_run(rs, getf, x, conv, startval = x(1:sb%dim))
      end do
    end if

    nullify(sb_p); nullify(cv_p)

    if(.not.conv) then
      message(1) = "During the construction of the adaptive grid, the Newton-Raphson"
      message(2) = "method did not converge for point:"
      write(message(3),'(9f14.6)') x(1:sb%dim)
      message(4) = "Try varying the Gygi parameters -- usually reducing CurvGygiA or"
      message(5) = "CurvGygiAlpha (or both) solves the problem."
      call messages_fatal(5)
    end if

    POP_SUB(curv_gygi_chi2x)
  end subroutine curv_gygi_chi2x


  ! ---------------------------------------------------------
  subroutine curv_gygi_x2chi(sb, cv, x, chi)
    type(simul_box_t), intent(in)  :: sb
    type(curv_gygi_t), intent(in)  :: cv
    FLOAT,             intent(in)  :: x(:)    ! x(sb%dim)
    FLOAT,             intent(out) :: chi(:)  ! chi(sb%dim)

    integer :: i, ia
    FLOAT   :: r, ar, th, ex

    PUSH_SUB(curv_gygi_x2chi)

    chi(1:sb%dim) = x(1:sb%dim)
    do ia = 1, cv%npos
      r = max(sqrt(sum((x(1:sb%dim) - cv%pos(ia, 1:sb%dim))**2)), CNST(1e-6))
      ar = cv%A*cv%alpha/r
      th = tanh(r/cv%alpha)
      ex = exp(-(r/cv%beta)**2)
      do i = 1, sb%dim
        chi(i) = chi(i) + (x(i) - cv%pos(ia, i)) * cv%a * ar * th * ex
      end do
    end do

    POP_SUB(curv_gygi_x2chi)
  end subroutine curv_gygi_x2chi


  ! ---------------------------------------------------------
  subroutine curv_gygi_jacobian(sb, cv, x, chi, J, natoms)
    type(simul_box_t), intent(in)  :: sb
    type(curv_gygi_t), intent(in)  :: cv
    FLOAT,             intent(in)  :: x(:)    ! x(sb%dim)
    FLOAT,             intent(out) :: chi(:)  ! chi(sb%dim)
    FLOAT,             intent(out) :: J(:,:)  ! J(sb%dim,sb%dim), the Jacobian
    integer, optional, intent(in)  :: natoms

    integer :: i, ix, iy, natoms_
    FLOAT :: r, f_alpha, df_alpha
    FLOAT :: th, ex, ar

    PUSH_SUB(curv_gygi_jacobian)

    J(1:sb%dim,1:sb%dim) = M_ZERO
    do ix = 1, sb%dim
      J(ix, ix) = M_ONE
      chi(ix)   = x(ix)
    end do

    natoms_ = cv%npos
    if(present(natoms)) natoms_ = natoms

    do i = 1, natoms_
      r = max(sqrt(sum((x(1:sb%dim) - cv%pos(i, 1:sb%dim))**2)), CNST(1e-6))

      ar = cv%A*cv%alpha/r
      th = tanh(r/cv%alpha)
      ex = exp(-(r/cv%beta)**2)

      f_alpha  = ar * th * ex
      df_alpha = ar*(-th*ex/r + ex/(cv%alpha*cosh(r/cv%alpha)**2) - th*M_TWO*r*ex/cv%beta**2)

      do ix = 1, sb%dim
        chi(ix) = chi(ix) + f_alpha*(x(ix)-cv%pos(i, ix))

        J(ix, ix) = J(ix, ix) + f_alpha
        do iy = 1, sb%dim
          J(ix, iy) = J(ix, iy) + (x(ix)-cv%pos(i, ix))*(x(iy)-cv%pos(i, iy))/r*df_alpha
        end do
      end do
    end do

    POP_SUB(curv_gygi_jacobian)
  end subroutine curv_gygi_jacobian

end module curv_gygi_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

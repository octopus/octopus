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
! F. Gygi and G. Galli, PRB 52 R2229 (1996).

module curv_gygi
  use lib_oct_parser
  use global
  use messages
  use syslabels
  use units
  use geometry
  use simul_box
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

    !%Variable CurvGygiA
    !%Type float
    !%Section 4 Mesh
    !%Description
    !% The grid spacing is reduced locally around each atom, and the reduction is
    !% given by 1/(1+A), where A is specified by this variable, CurvGygiA. So, if
    !% A=1 (the default), the grid spacing is reduced to one half = 1/(1+1).
    !% [This is the A_{\alpha} variable in Eq. 2 of F. Gygi and G. Galli, Phys.
    !% Rev. B 52, R2229 (1995)]
    !% It must be larger than zero.
    !%End
    call loct_parse_float(check_inp('CurvGygiA'), M_HALF, cv%A)
    !%Variable CurvGygiAlpha
    !%Type float
    !%Section 4 Mesh
    !%Description
    !% This number determines the region over which the grid is enhanced (range of
    !% enhancement of the resolution). That is, the grid is enhanced on a sphere
    !% around each atom, whose radius is given by this variable. [This is the a_{\alpha} 
    !% variable in Eq. 2 of F. Gygi and G. Galli, Phys. Rev. B 52, R2229 (1995)].
    !% The default is two atomic units. 
    !% It must be larger than zero.
    !%End
    call loct_parse_float(check_inp('CurvGygiAlpha'), M_TWO/units_inp%length%factor, cv%alpha)
    !%Variable CurvGygiBeta
    !%Type float
    !%Section 4 Mesh
    !%Description
    !% This number determines the distance over which Euclidean coordinates are
    !% recovered. [This is the b_{\alpha} variable in Eq. 2 of F. Gygi and G. Galli, 
    !% Phys. Rev. B 52, R2229 (1995)]. The default is four atomic units.
    !% It must be larger than zero.
    !%End
    call loct_parse_float(check_inp('CurvGygiBeta'),  M_FOUR/units_inp%length%factor, cv%beta)

    if(cv%a<=M_ZERO)     call input_error('CurvGygiA')
    if(cv%alpha<=M_ZERO) call input_error('CurvGygiAlpha')
    if(cv%beta<=M_ZERO)  call input_error('CurvGygiBeta')

    cv%alpha = cv%alpha*units_inp%length%factor
    cv%beta  = cv%beta *units_inp%length%factor

  end subroutine curv_gygi_init


  !-------------------------------------
  subroutine curv_gygi_chi2x(sb, geo, cv, chi, x)
    type(simul_box_type), intent(in)  :: sb
    type(geometry_type),  intent(in)  :: geo
    type(curv_gygi_type), intent(in)  :: cv
    FLOAT,                intent(in)  :: chi(:)  ! chi(sb%dim)
    FLOAT,                intent(out) :: x(:)    ! x(sb%dim)

    ! parameters
    integer, parameter :: max_iter = 500
    FLOAT,   parameter :: x_conv   = CNST(1e-10)

    ! local variables
    integer :: iter, i
    FLOAT, allocatable :: f(:,:), delta(:,:), J(:,:), chi2(:)
    logical :: conv

    allocate(f(sb%dim, 1), delta(sb%dim, 1), J(sb%dim, sb%dim), chi2(sb%dim))

    x(1:sb%dim) = chi(1:sb%dim)


    conv          = .false.
    do iter = 1, max_iter
       call curv_gygi_jacobian(sb, geo, cv, x, chi2, J)
       f(:,1) = chi(1:sb%dim) - chi2(:)
       if(sum(f(:,1)**2) < x_conv) then
           conv = .true.
           exit
       end if

       call lalg_linsyssolve(sb%dim, 1, J, f, delta)
       x(1:sb%dim) = x(1:sb%dim) + delta(1:sb%dim, 1)
    end do

    if(.not.conv) then
      x(1:sb%dim) = chi(1:sb%dim)
      do i = 1, geo%natoms

       conv          = .false.
       do iter = 1, max_iter
          call curv_gygi_jacobian(sb, geo, cv, x, chi2, J, i)
          f(:,1) = chi(1:sb%dim) - chi2(:)
          if(sum(f(:,1)**2) < x_conv) then
             conv = .true.
             exit
          end if

          call lalg_linsyssolve(sb%dim, 1, J, f, delta)
          x(1:sb%dim) = x(1:sb%dim) + delta(1:sb%dim, 1)
       end do

      enddo
    endif

    if(.not.conv) then
      message(1) = "During the construction of the adaptive grid, the Newton-Raphson"
      message(2) = "method did not converge for point:"
      write(message(3),'(3f14.6)') x(1:sb%dim)
      message(4) = "Try varying the Gygi parameters -- usually reducing CurvGygiA or"
      message(5) = "CurvGygiAlpha (or both) solves the problem."
      call write_fatal(5)
    endif

    ! clean up
    deallocate(f, delta, J, chi2)

  end subroutine curv_gygi_chi2x


  !-------------------------------------
  subroutine curv_gygi_jacobian(sb, geo, cv, x, chi, J, natoms)
    type(simul_box_type), intent(in)  :: sb
    type(geometry_type),  intent(in)  :: geo
    type(curv_gygi_type), intent(in)  :: cv
    FLOAT,                intent(in)  :: x(:)    ! x(sb%dim)
    FLOAT,                intent(out) :: chi(:)  ! chi(sb%dim)
    FLOAT,                intent(out) :: J(:,:)  ! J(sb%dim,sb%dim), the Jacobian
    integer, optional, intent(in) :: natoms

    FLOAT :: natoms_
    integer :: i, ix, iy
    FLOAT :: r, f_alpha, df_alpha
    FLOAT :: th, ex, ar

    J(1:sb%dim,1:sb%dim) = M_ZERO
    do ix = 1, sb%dim
      J(ix, ix) = M_ONE
      chi(ix)   = x(ix)
    end do

    natoms_ = geo%natoms
    if(present(natoms)) natoms_ = natoms
    do i = 1, natoms_
      r = max(sqrt(sum((x(1:sb%dim) - geo%atom(i)%x(1:sb%dim))**2)), CNST(1e-6))

      ar = cv%A*cv%alpha/r
      th = tanh(r/cv%alpha)
      ex = exp(-(r/cv%beta)**2)

      f_alpha  = ar * th * ex
      df_alpha = ar*(-th*ex/r + ex/(cv%alpha*cosh(r/cv%alpha)**2) - th*M_TWO*r*ex/cv%beta**2)

      do ix = 1, sb%dim
        chi(ix) = chi(ix) + f_alpha*(x(ix)-geo%atom(i)%x(ix))

        J(ix, ix) = J(ix, ix) + f_alpha
        do iy = 1, sb%dim
          J(ix, iy) = J(ix, iy) + (x(ix)-geo%atom(i)%x(ix))*(x(iy)-geo%atom(i)%x(iy))/r * df_alpha
        end do
      end do
    end do

  end subroutine curv_gygi_jacobian

end module curv_gygi

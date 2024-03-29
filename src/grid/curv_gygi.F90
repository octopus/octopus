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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

!> This module implements the curvilinear coordinates given in
!! F. Gygi and G. Galli, PRB 52 R2229 (1996).

module curv_gygi_oct_m
  use global_oct_m
  use ions_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use root_solver_oct_m
  use simul_box_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                   &
    curv_gygi_t,              &
    curv_gygi_init,           &
    curv_gygi_copy,           &
    curv_gygi_end,            &
    curv_gygi_chi2x,          &
    curv_gygi_x2chi,          &
    curv_gygi_jacobian

  type curv_gygi_t
    private
    FLOAT, public :: A             !< local reduction in grid spacing is 1/(1+A)
    FLOAT, public :: alpha         !< range of enhancement of the resolution
    FLOAT, public :: beta          !< distance over which Euclidian coordinates are recovered
    FLOAT, allocatable :: pos(:, :)
    integer :: npos
  end type curv_gygi_t

  type(simul_box_t), pointer  :: sb_p
  type(curv_gygi_t), pointer  :: cv_p
  integer :: i_p
  FLOAT :: chi_p(MAX_DIM)

contains

  ! ---------------------------------------------------------
  subroutine curv_gygi_init(cv, namespace, sb, ions, min_scaling_product)
    type(curv_gygi_t), intent(out) :: cv
    type(namespace_t), intent(in)  :: namespace
    type(simul_box_t), intent(in)  :: sb
    type(ions_t),      intent(in)  :: ions
    FLOAT,             intent(out) :: min_scaling_product

    PUSH_SUB(curv_gygi_init)

    !%Variable CurvGygiA
    !%Type float
    !%Default 0.5
    !%Section Mesh::Curvilinear::Gygi
    !%Description
    !% The grid spacing is reduced locally around each atom, and the reduction is
    !% given by 1/(1+<i>A</i>), where <i>A</i> is specified by this variable. So, if
    !% <i>A</i>=1/2 (the default), the grid spacing is reduced to two thirds = 1/(1+1/2).
    !% [This is the <math>A_{\alpha}</math> variable in Eq. 2 of F. Gygi and G. Galli, <i>Phys.
    !% Rev. B</i> <b>52</b>, R2229 (1995)]. It must be larger than zero.
    !%End
    call parse_variable(namespace, 'CurvGygiA', M_HALF, cv%A)

    !%Variable CurvGygiAlpha
    !%Type float
    !%Default 2.0 a.u.
    !%Section Mesh::Curvilinear::Gygi
    !%Description
    !% This number determines the region over which the grid is enhanced (range of
    !% enhancement of the resolution). That is, the grid is enhanced on a sphere
    !% around each atom, whose radius is given by this variable. [This is the <math>a_{\alpha}</math>
    !% variable in Eq. 2 of F. Gygi and G. Galli, <i>Phys. Rev. B</i> <b>52</b>, R2229 (1995)].
    !% It must be larger than zero.
    !%End

    call parse_variable(namespace, 'CurvGygiAlpha', M_TWO, cv%alpha, units_inp%length)
    !%Variable CurvGygiBeta
    !%Type float
    !%Default 4.0 a.u.
    !%Section Mesh::Curvilinear::Gygi
    !%Description
    !% This number determines the distance over which Euclidean coordinates are
    !% recovered. [This is the <math>b_{\alpha}</math> variable in Eq. 2 of F. Gygi and G. Galli,
    !% <i>Phys. Rev. B</i> <b>52</b>, R2229 (1995)]. It must be larger than zero.
    !%End
    call parse_variable(namespace, 'CurvGygiBeta', M_FOUR, cv%beta, units_inp%length)

    if(cv%a<=M_ZERO)     call messages_input_error(namespace, 'CurvGygiA')
    if(cv%alpha<=M_ZERO) call messages_input_error(namespace, 'CurvGygiAlpha')
    if(cv%beta<=M_ZERO)  call messages_input_error(namespace, 'CurvGygiBeta')

    cv%npos = ions%natoms
    SAFE_ALLOCATE(cv%pos(1:sb%dim, 1:cv%npos))
    cv%pos = ions%pos

    call curv_gygi_min_scaling(sb, cv, min_scaling_product)

    POP_SUB(curv_gygi_init)
  end subroutine curv_gygi_init

  ! ---------------------------------------------------------
  subroutine curv_gygi_copy(this_out, this_in)
    type(curv_gygi_t), intent(inout) :: this_out
    type(curv_gygi_t), intent(in)    :: this_in
    !
    PUSH_SUB(curv_gygi_copy)
    this_out%A=this_in%A
    this_out%alpha=this_in%alpha
    this_out%beta=this_in%beta
    SAFE_ALLOCATE_SOURCE_A(this_out%pos, this_in%pos)
    this_out%npos=this_in%npos
    POP_SUB(curv_gygi_copy)
    return
  end subroutine curv_gygi_copy

  ! ---------------------------------------------------------
  subroutine curv_gygi_end(cv)
    type(curv_gygi_t), intent(inout) :: cv

    PUSH_SUB(curv_gygi_end)

    SAFE_DEALLOCATE_A(cv%pos)

    POP_SUB(curv_gygi_end)
  end subroutine curv_gygi_end

  ! ---------------------------------------------------------
  subroutine getf(y, f, jf)
    FLOAT, intent(in)    :: y(:)
    FLOAT, intent(out)   :: f(:), jf(:, :)

    ! no push_sub, called too frequently

    call curv_gygi_jacobian(sb_p, cv_p, y, f, jf, i_p)
    f(1:sb_p%dim) = f(1:sb_p%dim) - chi_p(1:sb_p%dim)

  end subroutine getf 


  ! ---------------------------------------------------------
  subroutine curv_gygi_chi2x(sb, cv, rs, chi, x)
    type(simul_box_t), target, intent(in)   :: sb
    type(curv_gygi_t), target, intent(in)   :: cv
    type(root_solver_t), intent(in)         :: rs
    FLOAT,                     intent(in)  :: chi(:)  !< chi(sb%dim)
    FLOAT,                     intent(out) :: x(:)    !< x(sb%dim)

    integer :: i
    logical :: conv

    ! no push_sub, called too frequently

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
      r = max(norm2(x(1:sb%dim) - cv%pos(1:sb%dim, ia)), CNST(1e-6))
      ar = cv%A*cv%alpha/r
      th = tanh(r/cv%alpha)
      ex = exp(-(r/cv%beta)**2)
      do i = 1, sb%dim
        chi(i) = chi(i) + (x(i) - cv%pos(i, ia)) * cv%a * ar * th * ex
      end do
    end do

    POP_SUB(curv_gygi_x2chi)
  end subroutine curv_gygi_x2chi


  ! ---------------------------------------------------------
  subroutine curv_gygi_jacobian(sb, cv, x, chi, J, natoms)
    type(simul_box_t), intent(in)  :: sb
    type(curv_gygi_t), intent(in)  :: cv
    FLOAT,             intent(in)  :: x(:)    !< x(sb%dim)
    FLOAT,             intent(out) :: chi(:)  !< chi(sb%dim)
    FLOAT,             intent(out) :: J(:,:)  !< J(sb%dim,sb%dim), the Jacobian
    integer, optional, intent(in)  :: natoms

    integer :: i, ix, iy, natoms_
    FLOAT :: r, f_alpha, df_alpha
    FLOAT :: th, ex, ar

    ! no push_sub, called too frequently

    J(1:sb%dim,1:sb%dim) = M_ZERO
    do ix = 1, sb%dim
      J(ix, ix) = M_ONE
      chi(ix)   = x(ix)
    end do

    natoms_ = cv%npos
    if(present(natoms)) natoms_ = natoms

    do i = 1, natoms_
      r = max(norm2(x(1:sb%dim) - cv%pos(1:sb%dim, i)), CNST(1e-6))

      ar = cv%A*cv%alpha/r
      th = tanh(r/cv%alpha)
      ex = exp(-(r/cv%beta)**2)

      f_alpha  = ar * th * ex
      df_alpha = ar*(-th*ex/r + ex/(cv%alpha*cosh(r/cv%alpha)**2) - th*M_TWO*r*ex/cv%beta**2)

      do ix = 1, sb%dim
        chi(ix) = chi(ix) + f_alpha*(x(ix) - cv%pos(ix, i))

        J(ix, ix) = J(ix, ix) + f_alpha
        do iy = 1, sb%dim
          J(ix, iy) = J(ix, iy) + (x(ix) - cv%pos(ix, i))*(x(iy) - cv%pos(iy, i))/r*df_alpha
        end do
      end do
    end do

  end subroutine curv_gygi_jacobian


  ! ---------------------------------------------------------
  subroutine curv_gygi_min_scaling(sb, cv, min_scaling_product)
    type(simul_box_t), intent(in)  :: sb
    type(curv_gygi_t), intent(in)  :: cv
    FLOAT,             intent(out) :: min_scaling_product
    
    min_scaling_product = (1.0 / (1.0 + cv%A))**sb%dim
  end subroutine curv_gygi_min_scaling

end module curv_gygi_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

!! Copyright (C) 2005-2006 Heiko Appel
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

module root_solver_m
  use datasets_m
  use global_m
  use lalg_adv_m
  use messages_m
  use ode_solver_m
  use parser_m
  use profiling_m

  implicit none

  private
  public ::                               &
    root_solver_t,                        &
    root_solver_init,                     &
    root_solver_read,                     &
    droot_solver_run,                     &
    droot_solver_end,                     &
    zroot_solver_run,                     &
    zroot_solver_end,                     &
    zroot_watterstrom
 


  integer, public, parameter ::           &
    ROOT_BISECTION   =  1,                &
    ROOT_BRENT       =  2,                &
    ROOT_NEWTON      =  3,                &
    ROOT_LAGUERRE    =  4,                &
    ROOT_WATTERSTROM =  5,                &
    ROOT_MINVAL      =  ROOT_BISECTION,   &
    ROOT_MAXVAL      =  ROOT_WATTERSTROM

  type root_solver_t
    private
    integer :: solver_type    !< what solver to use (see ROOT_* variables above)_m
    integer :: dim            !< dimensionality of the problem
    integer :: maxiter        !< maximal number of iterations
    integer :: usediter       !< number of iterations actually performed
    FLOAT   :: abs_tolerance
    FLOAT   :: rel_tolerance
    FLOAT   :: ws_radius      !< radius of circle in complex plane; used for initial values
    logical :: have_polynomial
    integer :: poly_order
  end type root_solver_t

  !> a few variables which we have to define globally
  !! for this module
  CMPLX, allocatable :: gbase_coeff(:), gcoeff(:)
  integer            :: gorder

contains

  ! ---------------------------------------------------------
  subroutine root_solver_init(rs, dimensionality, solver_type, maxiter, rel_tolerance, &
    abs_tolerance, have_polynomial, ws_radius)
    type(root_solver_t), intent(out) :: rs
    integer,             intent(in)  :: dimensionality
    integer, optional,   intent(in)  :: solver_type, maxiter
    FLOAT, optional,     intent(in)  :: rel_tolerance, abs_tolerance, ws_radius
    logical, optional,   intent(in)  :: have_polynomial

    ! no push_sub, called too often

    ! Fill in the defaults
    rs%dim             = dimensionality
    rs%solver_type     = ROOT_NEWTON
    rs%maxiter         = 100
    rs%rel_tolerance   = CNST(1.0e-8)
    rs%abs_tolerance   = CNST(1.0e-8)
    rs%have_polynomial = .false.
    rs%ws_radius       = CNST(1.0)

    if(present(solver_type))     rs%solver_type     = solver_type
    if(present(maxiter))         rs%maxiter         = maxiter
    if(present(rel_tolerance))   rs%rel_tolerance   = rel_tolerance
    if(present(abs_tolerance))   rs%abs_tolerance   = abs_tolerance
    if(present(have_polynomial)) rs%have_polynomial = have_polynomial
    if(present(ws_radius))       rs%ws_radius       = ws_radius

  end subroutine root_solver_init


  ! ---------------------------------------------------------
  subroutine root_solver_read(rs)
    type(root_solver_t), intent(out) :: rs

    PUSH_SUB(root_solver_read)

    !%Variable RootSolver
    !%Type integer
    !%Default root_newton
    !%Section Math::General
    !%Description
    !% Specifies what kind of root solver will be used.
    !%Option root_bisection 1
    !% Bisection method.
    !%Option root_brent 2
    !% Brent method.
    !%Option root_newton 3
    !% Newton method.
    !%Option root_laguerre 4
    !% Laguerre method.
    !%Option root_watterstrom 5
    !% Watterstrom method.
    !%End
    call parse_integer(datasets_check('RootSolver'), ROOT_NEWTON, rs%solver_type)
    if( rs%solver_type.lt.ROOT_MINVAL.or.rs%solver_type.gt.ROOT_MAXVAL ) then
      call input_error(datasets_check('RootSolver'))
    end if

    !%Variable RootSolverMaxIter
    !%Type integer
    !%Default 100
    !%Section Math::General
    !%Description
    !% In case of an iterative root solver, this variable determines the maximum number
    !% of iteration steps.
    !%End
    call parse_integer(datasets_check('RootSolverMaxIter'), 100, rs%maxiter)

    !%Variable RootSolverRelTolerance
    !%Type float
    !%Default 1e-8
    !%Section Math::General
    !%Description
    !% Relative tolerance for the root-finding process.
    !%End
    call parse_float(datasets_check('RootSolverRelTolerance'), CNST(1e-8), rs%rel_tolerance)

    !%Variable RootSolverAbsTolerance
    !%Type float
    !%Default 1e-8
    !%Section Math::General
    !%Description
    !% Relative tolerance for the root-finding process.
    !%End
    call parse_float(datasets_check('RootSolverAbsTolerance'), CNST(1e-8), rs%abs_tolerance)

    !%Variable RootSolverHavePolynomial
    !%Type logical
    !%Default no
    !%Section Math::General
    !%Description
    !%  If set to yes, the coefficients of the polynomial have to be passed to
    !%  the root solver.
    !%End
    call parse_logical(datasets_check('RootSolverHavePolynomial'), .false., rs%have_polynomial)

    !%Variable RootSolverWSRadius
    !%Type float
    !%Default 1.0
    !%Section Math::General
    !%Description
    !% Radius of circle in the complex plane. If <tt>RootSolverWSRadius = 1.0</tt>,
    !% the unit roots of an <i>n</i>th-order polynomial are taken as initial values.
    !%End
    call parse_float(datasets_check('RootSolverWSRadius'), CNST( 1.0), rs%ws_radius)

    POP_SUB(root_solver_read)
  end subroutine root_solver_read

  ! ---------------------------------------------------------
  subroutine droot_bisection
    PUSH_SUB(droot_bisection)

    call messages_not_implemented('root bisection')

    POP_SUB(droot_bisection)
  end subroutine droot_bisection


  ! ---------------------------------------------------------
  !> Implementation of J. Comp. Phys., 8, (1971), p. 304-308
  subroutine zroot_watterstrom(rs, roots, coeff)
    type(root_solver_t), intent(in)  :: rs
    CMPLX,                  intent(out) :: roots(:)    !< roots we are searching
    CMPLX,                  intent(in)  :: coeff(:)    !< polynomial coefficients

    type(ode_solver_t) :: os
    CMPLX, allocatable    :: base_roots(:)
    FLOAT   :: theta
    integer :: order, j

    PUSH_SUB(zroot_watterstrom)

    order  = rs%poly_order
    gorder = order

    SAFE_ALLOCATE(gbase_coeff(1:order+1))
    SAFE_ALLOCATE(gcoeff     (1:order+1))
    SAFE_ALLOCATE(base_roots (1:order))

    ! normalize polynomial
    do j = 1, order+1
      gcoeff(j) = coeff(j)/coeff(order+1)
    end do

    gbase_coeff = M_ZERO
    gbase_coeff(1)       = (rs%ws_radius)**order
    gbase_coeff(order+1) = M_ONE

    do j = 1, order
      theta = (M_TWO*j-M_ONE)*M_PI/order
      base_roots(j) = exp(M_zI*theta)*(rs%ws_radius)
    end do

    !%Variable WatterstromODESolver
    !%Type integer
    !%Default ode_pd89
    !%Section Math::General
    !%Description
    !% The Watterstrom method (<i>J. Comp. Phys.</i> <b>8</b>, 304-308 (1971)) transforms
    !% finding roots for <i>n</i>th-order polynomials into the solution of <i>n</i> uncoupled 
    !% ODEs. This variable specifies the solver that should be used for the ODE 
    !% stepping. Valid solver types are those allowed for the <tt>ODESolver</tt> variable.
    !%Option ode_rk4 1
    !% Standard 4th-order Runge-Kutta.
    !%Option ode_fb78 2
    !% Fehlberg solver.
    !%Option ode_vr89 3
    !% Verner solver.
    !%Option ode_pd89 4
    !% Prince-Dormand solver.
    !%End
    call parse_integer(datasets_check('WatterstromODESolver'), ODE_PD89, os%solver_type)

    !%Variable WatterstromODESolverNSteps
    !%Type integer
    !%Default 400
    !%Section Math::General
    !%Description
    !% Number of steps which the chosen ODE solver should perform
    !% in the integration interval [<i>a</i>, <i>b</i>] of the Watterstrom ODE.
    !%End
    call parse_integer(datasets_check('WatterstromODESolverNSteps'), 400, os%nsteps)

    ! set up ODE solver
    os%nsize       = order
    os%tmin        = M_ZERO
    os%tmax        = M_ONE
    call zode_solver_create(os)
    call zode_solver_run(os, func_ws, base_roots, roots)

    SAFE_DEALLOCATE_A(gbase_coeff)
    SAFE_DEALLOCATE_A(gcoeff)
    SAFE_DEALLOCATE_A(base_roots)

    POP_SUB(zroot_watterstrom)

  end subroutine zroot_watterstrom


  ! ---------------------------------------------------------
  subroutine func_ws(size, t, z, res)
    integer, intent(in)  :: size
    FLOAT,   intent(in)  :: t
    CMPLX,   intent(in)  :: z(:)
    CMPLX,   intent(out) :: res(:)

    CMPLX, allocatable   :: numerator(:), denominator(:)
    integer :: j

    PUSH_SUB(func_ws)

    SAFE_ALLOCATE(  numerator(1:size))
    SAFE_ALLOCATE(denominator(1:size))
    numerator   = M_ZERO
    denominator = M_ZERO

    do j = 0, gorder-1
      numerator = numerator + (gbase_coeff(j+1)-gcoeff(j+1))*z**j
    end do

    do j = 1, gorder
      denominator = denominator + j*( gbase_coeff(j+1)-(gbase_coeff(j+1)-gcoeff(j+1))*t )*z**(j-1)
    end do

    res = numerator/denominator

    SAFE_DEALLOCATE_A(numerator)
    SAFE_DEALLOCATE_A(denominator)

    POP_SUB(func_ws)

  end subroutine func_ws


  ! ---------------------------------------------------------
  !> Newton-Raphson scheme can only be used in the real case.
  subroutine droot_newton(rs, func, root, startval, success)
    type(root_solver_t), intent(in)  :: rs
    FLOAT,               intent(out) :: root(:)        !< root we are searching
    FLOAT,               intent(in)  :: startval(:)    !< start value for the search
    logical,             intent(out) :: success
    interface
      subroutine func(z, f, jf)
        FLOAT, intent(in)  :: z(:)
        FLOAT, intent(out) :: f(:), jf(:, :)
      end subroutine func
    end interface

    integer :: iter
    FLOAT   :: err
    FLOAT, allocatable :: f(:), jf(:, :), delta(:, :), rhs(:, :)

    ! no push_sub, called too often

    SAFE_ALLOCATE(    f(1:rs%dim))
    SAFE_ALLOCATE(   jf(1:rs%dim, 1:rs%dim))
    SAFE_ALLOCATE(delta(1:rs%dim, 1))
    SAFE_ALLOCATE(  rhs(1:rs%dim, 1))

    root = startval
    call func(root, f, jf)
    err = sum(f(1:rs%dim)**2)

    success = .true.
    iter = 0
    do while(err > rs%abs_tolerance)
      rhs(1:rs%dim, 1) = -f(1:rs%dim)
      call lalg_linsyssolve(rs%dim, 1, jf, rhs, delta)
      root(1:rs%dim) = root(1:rs%dim) + delta(1:rs%dim, 1)
      iter = iter + 1
      if(iter > rs%maxiter) then
        success = .false.
        exit
      end if
      call func(root, f, jf)
      err = sum(f(1:rs%dim)**2)
    end do

    SAFE_DEALLOCATE_A(f)
    SAFE_DEALLOCATE_A(jf)
    SAFE_DEALLOCATE_A(delta)
    SAFE_DEALLOCATE_A(rhs)

  end subroutine droot_newton



#include "undef.F90"
#include "complex.F90"
#include "root_solver_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "root_solver_inc.F90"


end module root_solver_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module root_solver_m
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_parser_m
  use ode_solver_m

  implicit none

  private
  public ::                               &
    root_solver_t,                        &
    droot_solver_init,                    &
!    droot_solver_run,                     &
    droot_solver_end,                     &
    zroot_solver_init,                    &
!    zroot_solver_run,                     &
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
    integer :: solver_type    ! what solver to use (see ROOT_* variables above)_m
    integer :: maxiter        ! maximal number of iterations
    integer :: usediter       ! number of actually performed iterations
    FLOAT   :: abs_tolerance
    FLOAT   :: rel_tolerance
    FLOAT   :: ws_radius      ! radius of circle in complex plane; used for initial values
    logical :: have_polynomial
   integer :: poly_order
  end type root_solver_t

  ! a few variables which we have to define global
  ! for this module
  CMPLX, allocatable :: gbase_coeff(:), gcoeff(:)
  integer            :: gorder

contains

  ! ---------------------------------------------------------
  subroutine droot_bisection
    call push_sub('root_solver.droot_bisection')

    message(1) = 'Not implemented yet.'
    call write_fatal(1)

    call pop_sub()
  end subroutine droot_bisection


  ! ---------------------------------------------------------
!!$subroutine droot_brent(rs, func, root, interval)
!!$  type(root_solver_t), intent(in)  :: rs
!!$  FLOAT,                  intent(out) :: root
!!$  FLOAT,  optional,       intent(in)  :: interval(2)     ! lower and upper boundary of search region
!!$
!!$  interface
!!$     subroutine func(x,s)
!!$       FLOAT :: x,s
!!$     end subroutine func
!!$  end interface
!!$
!!$  call push_sub('root_solver.droot_brent')
!!$
!!$  message(1) = 'Error: Not implemented yet.'
!!$  call write_fatal(1)
!!$
!!$  call pop_sub()
!!$end subroutine droot_brent


  ! ---------------------------------------------------------
  ! Implementation of J. Comp. Phys., 8, (1971), p. 304-308
  subroutine zroot_watterstrom(rs, roots, coeff)
    type(root_solver_t), intent(in)  :: rs
    CMPLX,                  intent(out) :: roots(:)    ! roots we are searching
    CMPLX,                  intent(in)  :: coeff(:)    ! polynomial coefficients

    type(ode_solver_t) :: os
    CMPLX, allocatable    :: base_roots(:)
    FLOAT   :: theta
    integer :: order, j

    call push_sub('root_solver.zroot_watterstrom')

    order  = rs%poly_order
    gorder = order

    ALLOCATE(gbase_coeff(order+1), order+1)
    ALLOCATE(gcoeff     (order+1), order+1)
    ALLOCATE(base_roots (order),   order)

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
    !% The Watterstrom method (cf. J. Comp. Phys., 8, (1971), p. 304-308) transforms
    !% the root finding for n-th order polynomials into the solution of n uncoupled 
    !% ODEs. This variable specifies the ODESolver that should be used for the ODE 
    !% stepping. Valid solver types are the ones that are allowed for ODESolver (cf. 
    !% variable ODESolver).
    !%Option ode_rk4 1
    !% Standard Runge-Kutta 4th order
    !%Option ode_fb78 2
    !% Fehlberg solver
    !%Option ode_vr89 3
    !% Verner solver
    !%Option ode_pd89 4
    !% Prince-Dormand solver
    !%End
    call loct_parse_int(check_inp('WatterstromODESolver'),       ODE_PD89, os%solver_type)

    !%Variable WatterstromODESolverNSteps
    !%Type integer
    !%Default 400
    !%Section Math::General
    !%Description
    !% Number of steps which the chosen ODE solver should perform
    !% in the integration interval [a,b] of the Watterstrom ODE.
    !%End
    call loct_parse_int(check_inp('WatterstromODESolverNSteps'),      400, os%nsteps)

    ! setup ode solver
    os%nsize       = order
    os%tmin        = M_ZERO
    os%tmax        = M_ONE
    call zode_solver_create(os)
    call zode_solver_run(os, func_ws, base_roots, roots)

    deallocate(gbase_coeff, gcoeff, base_roots)

    call pop_sub()

  end subroutine zroot_watterstrom


  ! ---------------------------------------------------------
  subroutine func_ws(size, t, z, res)
    integer, intent(in)  :: size
    FLOAT,   intent(in)  :: t
    CMPLX,   intent(in)  :: z(:)
    CMPLX,   intent(out) :: res(:)

    CMPLX, allocatable   :: numerator(:), denominator(:)
    integer :: j

    ALLOCATE(  numerator(size), size)
    ALLOCATE(denominator(size), size)
    numerator   = M_ZERO
    denominator = M_ZERO

    do j = 0, gorder-1
      numerator = numerator + (gbase_coeff(j+1)-gcoeff(j+1))*z**j
    end do

    do j = 1, gorder
      denominator = denominator + j*( gbase_coeff(j+1)-(gbase_coeff(j+1)-gcoeff(j+1))*t )*z**(j-1)
    end do

    res = numerator/denominator

    deallocate(numerator, denominator)

  end subroutine func_ws


#include "undef.F90"
#include "complex.F90"
#include "root_solver_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "root_solver_inc.F90"


end module root_solver_m

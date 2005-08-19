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

module root_solver
  use global
  use messages
  use syslabels
  use ode_solver

  implicit none

  private
  public :: root_solver_type
  public :: &
       droot_solver_init, droot_solver_run, droot_solver_end, &
       zroot_solver_init, zroot_solver_run, zroot_solver_end, &
       zroot_watterstrom

  integer, public, parameter :: &
       ROOT_BISECTION   =  1,   &
       ROOT_BRENT       =  2,   &
       ROOT_NEWTON      =  3,   &
       ROOT_LAGUERRE    =  4,   &
       ROOT_WATTERSTROM =  5

  type root_solver_type
     integer :: solver_type    ! what solver to use (see ROOT_* variables above)
     integer :: maxiter        ! maximal number of iterations
     integer :: usediter       ! number of actually performed iterations
     FLOAT   :: abs_tolerance
     FLOAT   :: rel_tolerance
     FLOAT   :: ws_radius      ! radius of circle in complex plane; used for initial values
     logical :: have_polynomial
     integer :: poly_order
  end type root_solver_type

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
!!$  type(root_solver_type), intent(in)  :: rs
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
  type(root_solver_type), intent(in)  :: rs
  CMPLX,                  intent(out) :: roots(:)    ! roots we are searching
  CMPLX,                  intent(in)  :: coeff(:)    ! polynomial coefficients

  type(ode_solver_type) :: os
  CMPLX, allocatable    :: base_roots(:)
  FLOAT   :: theta
  integer :: order, j

  call push_sub('root_solver.zroot_watterstrom')

  order  = rs%poly_order
  gorder = order

  allocate(gbase_coeff(order+1), gcoeff(order+1), base_roots(order))

  ! normalize polynomial
  do j = 1, order+1
     gcoeff(j) = coeff(j)/coeff(order+1)
  enddo

  gbase_coeff = M_ZERO
  gbase_coeff(1)       = (rs%ws_radius)**order
  gbase_coeff(order+1) = M_ONE

  do j = 1, order
     theta = (M_TWO*j-M_ONE)*M_PI/order
     base_roots(j) = exp(M_zI*theta)*(rs%ws_radius)
  enddo

  ! setup ode solver
  call loct_parse_int(check_inp('WatterstromODESolver'),       ODE_PD89, os%solver_type)
  call loct_parse_int(check_inp('WatterstromODESolverNSteps'),      400, os%nsteps)
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

  allocate(numerator(size), denominator(size))
  numerator   = M_ZERO
  denominator = M_ZERO

  do j = 0, gorder-1
     numerator = numerator + (gbase_coeff(j+1)-gcoeff(j+1))*z**j
  enddo

  do j = 1, gorder
     denominator = denominator + j*( gbase_coeff(j+1)-(gbase_coeff(j+1)-gcoeff(j+1))*t )*z**(j-1)
  enddo

  res = numerator/denominator

  deallocate(numerator, denominator)

end subroutine func_ws


#include "undef.F90"
#include "complex.F90"
#include "root_solver_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "root_solver_inc.F90"


end module root_solver

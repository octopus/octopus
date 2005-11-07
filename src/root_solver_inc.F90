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

! ---------------------------------------------------------
subroutine X(root_solver_init)(rs)
  type(root_solver_type), intent(out) :: rs

  call push_sub('root_solver_inc.Xroot_solver_init')

  !%Variable RootSolver
  !%Type integer
  !%Default root_newton
  !%Section Math::General
  !%Description
  !% Specifies what kind of root solver will be used.
  !%Option root_bisection 1
  !% Bisection method
  !%Option root_brent 2
  !% Brent method
  !%Option root_newton 3
  !% Newton method
  !%Option root_laguerre 4
  !% Laguerre method
  !%Option root_watterstrom 5
  !% Watterstrom method
  !%End
  call loct_parse_int(check_inp('RootSolver'),        ROOT_NEWTON, rs%solver_type)
  if( rs%solver_type.lt.ROOT_MINVAL.or.rs%solver_type.gt.ROOT_MAXVAL ) then
    call input_error(check_inp('RootSolver'))
  end if

  !%Variable RootSolverMaxIter
  !%Type integer
  !%Default 100
  !%Section Math::General
  !%Description
  !% In case of an interative root solver this variable determines the maximal number
  !% of iteration steps.
  !%End
  call loct_parse_int    (check_inp('RootSolverMaxIter'),               100, rs%maxiter)

  !%Variable RootSolverRelTolerance
  !%Type float
  !%Default 1e-8
  !%Section Math::General
  !%Description
  !% Relative tolerance for the root finding process.
  !%End
  call loct_parse_float  (check_inp('RootSolverRelTolerance'),   CNST(1e-8), rs%rel_tolerance)

  !%Variable RootSolverAbsTolerance
  !%Type float
  !%Default 1e-8
  !%Section Math::General
  !%Description
  !% Relative tolerance for the root finding process.
  !%End
  call loct_parse_float  (check_inp('RootSolverAbsTolerance'),   CNST(1e-8), rs%abs_tolerance)

  !%Variable RootSolverHavePolynomial
  !%Type logical
  !%Default no
  !%Section Math::General
  !%Description
  !%  If set to yes, the coefficients of the polynomial have to be passed to
  !%  the root solver.
  !%End
  call loct_parse_logical(check_inp('RootSolverHavePolynomial'),    .false., rs%have_polynomial)

  !%Variable RootSolverWSRadius
  !%Type float
  !%Default 1.0
  !%Section Math::General
  !%Description
  !% Radius of circle in the complex plane. If RootSolverWSRadius = 1.0
  !% the unit roots of a n-th order polynomial are taken as initial values.
  !%End
  call loct_parse_float  (check_inp('RootSolverWSRadius'),       CNST( 1.0), rs%ws_radius)

  call pop_sub()
end subroutine X(root_solver_init)


!!$! ---------------------------------------------------------
!!$subroutine X(root_solver_run)(rs, func, roots, startval, interval_, coeff)
!!$  type(root_solver_type), intent(inout) :: rs
!!$  R_TYPE,                 intent(out)   :: roots(:)    ! roots we are searchin
!!$  R_TYPE, optional,       intent(in)    :: startval    ! start value for the search
!!$  FLOAT,  optional,       intent(in)    :: interval_(2) ! lower and upper boundary of search region
!!$  R_TYPE, optional,       intent(in)    :: coeff(:)    ! polynomial coefficients
!!$
!!$  FLOAT :: interval(2) 
!!$
!!$  interface
!!$    subroutine func(z,s)
!!$      R_TYPE :: z,s      ! s = f(z)
!!$    end subroutine func
!!$  end interface
!!$
!!$
!!$  call push_sub('root_solver_inc.Xroot_solver_run')
!!$
!!$  ! initialize array
!!$  roots = M_ZERO
!!$
!!$  if(present(interval_)) interval = interval_
!!$
!!$  select case(rs%solver_type)
!!$
!!$  case(ROOT_BISECTION)
!!$#ifdef R_TREAL
!!$    message(1) = 'Info: root_solver: Using bisection.'
!!$    call write_info(1)
!!$    call droot_bisection()
!!$#endif
!!$#ifdef R_TCOMPLEX
!!$    message(1) = 'Error: root_solver: Simple line bisection not defined for complex arithmetic'
!!$    call write_fatal(1)
!!$#endif
!!$
!!$
!!$  case(ROOT_BRENT)
!!$#ifdef R_TREAL
!!$    if(present(interval_)) then
!!$      message(1) = 'Info: root_solver: Using Brent method.'
!!$      call write_info(1)
!!$!!        call droot_brent(rs, func, roots(1), interval)
!!$    else
!!$      message(1) = 'Error: root_solver: search interval required for Brent method.'
!!$      call write_fatal(1)
!!$    end if
!!$#endif
!!$#ifdef R_TCOMPLEX
!!$    message(1) = 'Error: root_solver: Brent method not defined for complex arithmetic'
!!$    call write_fatal(1)
!!$#endif
!!$
!!$
!!$  case(ROOT_NEWTON)
!!$    message(1) = 'Info: root_solver: Using Newton method.'
!!$    call write_info(1)
!!$    ! pass roots(1): only a single root will be returned
!!$!!     call X(root_newton)(rs, func, roots(1), startval)
!!$
!!$
!!$  case(ROOT_LAGUERRE)
!!$    if(present(coeff)) then
!!$      message(1) = 'Info: root_solver: Using Laguerre method.'
!!$      call write_info(1)
!!$      ! pass roots(1): only a single root will be returned
!!$      call X(root_laguerre)(rs, roots(1), startval, coeff)
!!$    else
!!$      message(1) = 'Error: root_solver: Laguerre method only valid for polynomials.'
!!$      call write_fatal(1)
!!$    end if
!!$
!!$  case(ROOT_WATTERSTROM)
!!$#ifdef R_TREAL
!!$    message(1) = 'Error: root_solver: Watterstrom method not defined for pure real arithmetic'
!!$    call write_fatal(1)
!!$#endif
!!$#ifdef R_TCOMPLEX
!!$    if(present(coeff)) then
!!$      message(1) = 'Info: root_solver: Using Watterstrom method.'
!!$      call write_info(1)
!!$      call zroot_watterstrom(rs, roots, coeff)
!!$    else
!!$      message(1) = 'Error: root_solver: Watterstrom method only valid for polynomials.'
!!$      call write_fatal(1)
!!$    end if
!!$#endif
!!$
!!$  case default
!!$    write(message(1), '(a,i4,a)') "Input: '", rs%solver_type, &
!!$      "' is not a valid root solver"
!!$    message(2) = '( root solver =  root_bisection | root_brent | root_newton | root_laguerre | root_watterstrom )'
!!$    call write_fatal(2)
!!$  end select
!!$
!!$
!!$  call pop_sub()
!!$end subroutine X(root_solver_run)


! ---------------------------------------------------------
subroutine X(root_solver_create)()
!!$  type(root_solver_type), intent(in) :: rs

  call push_sub('root_solver_inc.Xroot_solver_create')
  ! do allocation stuff

  call pop_sub()
end subroutine X(root_solver_create)


! ---------------------------------------------------------
subroutine X(root_solver_end)()
  call push_sub('root_solver_inc.Xroot_solver_end')
  ! do deallocation stuff

  call pop_sub()
end subroutine X(root_solver_end)


! ---------------------------------------------------------
!!$subroutine X(root_newton)(rs, func, root, startval)
!!$  type(root_solver_type), intent(in)  :: rs
!!$  R_TYPE,                 intent(out) :: root        ! root we are searching
!!$  R_TYPE,                 intent(in)  :: startval    ! start value for the search
!!$
!!$    interface
!!$       subroutine func(z,s)
!!$         R_TYPE :: z,s
!!$       end subroutine func
!!$    end interface
!!$
!!$  call push_sub('root_solver_inc.Xroot_newton')
!!$
!!$  message(1) = 'Error: root_solver: Not Newton Method not implemented yet.'
!!$  call write_fatal(1)
!!$
!!$  call pop_sub()
!!$end subroutine X(root_newton)


! ---------------------------------------------------------
subroutine X(root_laguerre)(rs, root, startval, coeff)
  type(root_solver_type), intent(inout) :: rs
  R_TYPE,                 intent(out)   :: root        ! root we are searching
  R_TYPE,                 intent(in)    :: startval    ! start value for the search
  R_TYPE,                 intent(in)    :: coeff(:)    ! polynomial coefficients

  R_TYPE,  allocatable  :: b(:), c(:), d(:)
  R_TYPE  :: z, zold, s1, s2, denom1, denom2, lroot
  integer :: order, j, k

  call push_sub('root_solver_inc.Xroot_laguerre')

  order = rs%poly_order

  allocate(b(order), c(order-1), d(order-2))

  z = startval

  do k = 1, rs%maxiter

    ! Horner scheme for function value of polynomial
    b(1) = coeff(1)
    do j = 2, order
      b(j) = coeff(j) + z*b(j-1)
    end do

    ! first derivative
    c(1) = b(1)
    do j = 2, order-1
      c(j) = b(j) + z*c(j-1)
    end do

    s1 = c(order-1)/b(order)

    ! second derivative
    d(1) = c(1)
    do j = 2, order-2
      d(j) = c(j) + z*d(j-1)
    end do

    s2     = M_TWO*d(order-2)/b(order) - s1**2
    lroot  = sqrt( (1-order)*( order*s2 + s1**2 ) )
    denom1 = s1 + lroot
    denom2 = s1 - lroot

    zold = z

    if ( abs(denom1).gt.abs(denom2) ) then
      z = z - order / denom1
    else
      z = z - order / denom2
    end if

    if (abs(zold-z).lt.rs%abs_tolerance) exit

  end do

  rs%usediter = k
  root = z

  deallocate(b, c, d)

  call pop_sub()
end subroutine X(root_laguerre)

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

! ---------------------------------------------------------
subroutine X(root_solver_run)(rs, func, root, success, startval, interval_, coeff)
  type(root_solver_t), intent(inout) :: rs
  R_TYPE,                intent(out)  :: root(:)        ! roots we are searching
  logical,               intent(out)  :: success
  R_TYPE, optional,      intent(in)   :: startval(:)    ! start value for the search
  FLOAT,  optional,      intent(in)   :: interval_(2)   ! lower and upper boundary of search region
  R_TYPE, optional,      intent(in)   :: coeff(:)       ! polynomial coefficients
  interface
    subroutine func(z, f, jf)
      R_TYPE, intent(in)  :: z(:)
      R_TYPE, intent(out) :: f(:), jf(:, :)
    end subroutine func
  end interface

!!$  FLOAT :: interval(2) 

  PUSH_SUB(X(root_solver_run))

  ! Initializations
  root = M_ZERO
  success = .false.
!!$
!!$  if(present(interval_)) interval = interval_
!!$
  select case(rs%solver_type)
!!$
!!$  case(ROOT_BISECTION)
!!$#ifdef R_TREAL
!!$    message(1) = 'Info: root_solver: Using bisection.'
!!$    call messages_info(1)
!!$    call droot_bisection()
!!$#endif
!!$#ifdef R_TCOMPLEX
!!$    message(1) = 'Error: root_solver: Simple line bisection not defined for complex arithmetic'
!!$    call messages_fatal(1)
!!$#endif
!!$
!!$
!!$  case(ROOT_BRENT)
!!$#ifdef R_TREAL
!!$    if(present(interval_)) then
!!$      message(1) = 'Info: root_solver: Using Brent method.'
!!$      call messages_info(1)
!!$!!        call droot_brent(rs, func, root(1), interval)
!!$    else
!!$      message(1) = 'Error: root_solver: search interval required for Brent method.'
!!$      call messages_fatal(1)
!!$    end if
!!$#endif
!!$#ifdef R_TCOMPLEX
!!$    message(1) = 'Error: root_solver: Brent method not defined for complex arithmetic'
!!$    call messages_fatal(1)
!!$#endif
!!$
!!$
#if defined(R_TREAL)
  case(ROOT_NEWTON)
    call droot_newton(rs, func, root, startval, success)
#endif
!!$!!     call X(root_newton)(rs, func, root(1), startval)
!!$
!!$
!!$  case(ROOT_LAGUERRE)
!!$    if(present(coeff)) then
!!$      message(1) = 'Info: root_solver: Using Laguerre method.'
!!$      call messages_info(1)
!!$      ! pass root(1): only a single root will be returned
!!$      call X(root_laguerre)(rs, root(1), startval, coeff)
!!$    else
!!$      message(1) = 'Error: root_solver: Laguerre method only valid for polynomials.'
!!$      call messages_fatal(1)
!!$    end if
!!$
  case(ROOT_WATTERSTROM)
#ifdef R_TREAL
    message(1) = 'Error: root_solver: Watterstrom method not defined for pure real arithmetic'
    call messages_fatal(1)
#endif
#ifdef R_TCOMPLEX
    if(present(coeff)) then
      message(1) = 'Info: root_solver: Using Watterstrom method.'
      call messages_info(1)
      call zroot_watterstrom(rs, root, coeff)
    else
      message(1) = 'Error: root_solver: Watterstrom method only valid for polynomials.'
      call messages_fatal(1)
    end if
#endif

  case default
    write(message(1), '(a,i4,a)') "Input: '", rs%solver_type, &
      "' is not a valid root solver"
    message(2) = '( root solver =  root_bisection | root_brent | root_newton | root_laguerre | root_watterstrom )'
    call messages_fatal(2)
  end select

  POP_SUB(X(root_solver_run))
end subroutine X(root_solver_run)


! ---------------------------------------------------------
subroutine X(root_solver_create)()
!!$  type(root_solver_t), intent(in) :: rs

  PUSH_SUB(X(root_solver_create))
  ! do allocation stuff

  POP_SUB(X(root_solver_create))
end subroutine X(root_solver_create)


! ---------------------------------------------------------
subroutine X(root_solver_end)()
  PUSH_SUB(X(root_solver_end))
  ! do deallocation stuff

  POP_SUB(X(root_solver_end))
end subroutine X(root_solver_end)


! ---------------------------------------------------------
subroutine X(root_laguerre)(rs, root, startval, coeff)
  type(root_solver_t), intent(inout) :: rs
  R_TYPE,                 intent(out)   :: root        ! root we are searching
  R_TYPE,                 intent(in)    :: startval    ! start value for the search
  R_TYPE,                 intent(in)    :: coeff(:)    ! polynomial coefficients

  R_TYPE,  allocatable  :: b(:), c(:), d(:)
  R_TYPE  :: z, zold, s1, s2, denom1, denom2, lroot
  integer :: order, j, k

  PUSH_SUB(X(root_laguerre))

  order = rs%poly_order

  SAFE_ALLOCATE(b(1:order))
  SAFE_ALLOCATE(c(1:order-1))
  SAFE_ALLOCATE(d(1:order-2))

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

  SAFE_DEALLOCATE_A(b)
  SAFE_DEALLOCATE_A(c)
  SAFE_DEALLOCATE_A(d)

  POP_SUB(X(root_laguerre))
end subroutine X(root_laguerre)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

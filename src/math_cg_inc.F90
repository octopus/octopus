!! Copyright (C) 2004-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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



!/* -------------------------------------------------------------------
! The two following subroutines, sym_conjugate_gradients, and bi_conjugate_gradients,
! must be called under a common interface: conjugate_gradients. It provides an
! approximate solution to the linear system problem  Ax = b.
!
!
! $Id$
! Note: the complex valued versions (both CG and BiCG) work only with hermitean
!       operators. This is a property of the CG algorithm. A comment on this can
!       be found in chapter 4.1 of ftp://ftp.netlib.org/templates/templates.ps
!
!
! subroutine conjugate_gradients(np, x, b, op, [opt,] iter [, residue] [, threshold] )
!    integer, intent(in)  :: np    =>  The dimension of the problem.
!    FLOAT, intent(inout) :: x     => On input, an estimate to the solution.
!                                  => On output, the approximate solution.
!    FLOAT, intent(in)    :: b     => The inhomogeneous term of the equation.
!    interface
!      subroutine op(x, y)
!         FLOAT, intent(in)  :: x(:)
!         FLOAT, intent(out) :: y(:)
!      end subroutine op
!    end interface                 => This should be an interface to a procedure that
!                                     computes Ax = y
!    interface
!      subroutine opt(x, y)
!         FLOAT, intent(in)  :: x(:)
!         FLOAT, intent(out) :: y(:)
!      end subroutine opt
!    end interface                 => If present, this should be an interface to a procedure that
!                                     computes A^T x = y. Only useful for non-symmetric
!                                     operators.
!    integer, intent(inout) :: iter => On input, the maximum number of iteratios that
!                                      the procedure is allowed to take.
!                                      On output, the iterations it actually did.
!    FLOAT, intent(out) :: residue  => If present, it measures the final error:
!                                      residue = < Ax - b | Ax - b>
!    FLOAT, intent(in)  :: threshold => If present, it sets the required accuracy
!                                       threshold for the algorithm to stop. If not
!                                       present, this is set to 1.0e-6. [The algorithm
!                                       stops when < Ax - b | Ax - b > <= threshold, or
!                                       iter iterations are reached]
! end subroutine conjugate_gradients
!
! (*) NOTE: The algorithm assumes that the vectors are given in an orthonormal basis.
!
! ------------------------------------------------------------------- */!
subroutine X(sym_conjugate_gradients)(np, x, b, op, dotp, iter, residue, threshold)
  integer, intent(in) :: np
  R_TYPE,                intent(inout) :: x(:)
  R_TYPE,                intent(in)    :: b(:)
  interface
    subroutine op(x, y)
      R_TYPE, intent(inout) :: x(:)
      R_TYPE, intent(out)   :: y(:)
    end subroutine op
    R_TYPE function dotp(x, y) result(res)
      R_TYPE, intent(inout) :: x(:)
      R_TYPE, intent(in)    :: y(:)
    end function dotp
  end interface
  integer, intent(inout)        :: iter
  FLOAT,  optional, intent(in)  :: threshold
  R_TYPE, optional, intent(out) :: residue


  R_TYPE, allocatable :: r(:), ax(:), p(:), ap(:)
  R_TYPE  :: alpha, beta, gamma
  FLOAT   :: threshold_
  integer :: max_iter

  call push_sub('math_cg_inc.Xsym_conjugate_gradients')

  if(present(threshold)) then
    threshold_ = threshold
  else
    threshold_ = CNST(1.0e-6)
  end if

  ALLOCATE( r(np), np)
  ALLOCATE(ax(np), np)
  ALLOCATE( p(np), np)
  ALLOCATE(ap(np), np)

  ! Initial residue
  call op(x, ax)
  r(1:np) = b(1:np) - ax(1:np)

  ! Initial search direction
  p(1:np) = r(1:np)

  max_iter = iter
  iter = 1
  do while(iter < max_iter)
    gamma = dotp(r, r)
    if(abs(gamma) < THRESHOLD_) exit
    call op(p, ap)
    alpha   = gamma/dotp(p, ap)
    r(1:np) = r(1:np) - alpha*ap(1:np)
    x(1:np) = x(1:np) + alpha*p(1:np)
    beta    = dotp(r, r)/gamma
    p(1:np) = r(1:np) + beta*p(1:np)
    iter    = iter + 1
  end do
  if(present(residue)) residue = gamma

  deallocate(r, ax, p, ap)

  call pop_sub()
end subroutine X(sym_conjugate_gradients)


! ---------------------------------------------------------
subroutine X(bi_conjugate_gradients)(np, x, b, op, opt, dotp, iter, residue, threshold)
  integer, intent(in) :: np
  R_TYPE,                intent(inout) :: x(:)
  R_TYPE,                intent(in)    :: b(:)
  interface
    subroutine op(x, y)
      R_TYPE, intent(inout) :: x(:)
      R_TYPE, intent(out)   :: y(:)
    end subroutine op
  end interface
  interface
    subroutine opt(x, y)
      R_TYPE, intent(inout) :: x(:)
      R_TYPE, intent(out)   :: y(:)
    end subroutine opt
    R_TYPE function dotp(x, y) result(res)
      R_TYPE, intent(in) :: x(:)
      R_TYPE, intent(in) :: y(:)
    end function dotp
  end interface
  integer, intent(inout) :: iter
  R_TYPE, optional, intent(out) :: residue
  FLOAT, optional, intent(in) :: threshold

  R_TYPE, allocatable :: r(:), rr(:), ax(:), p(:), pp(:), ap(:), atp(:)
  FLOAT :: alpha, beta, gamma, err, threshold_
  integer :: max_iter

  call push_sub('math_cg_inc.Xbi_conjugate_gradients')

  if(present(threshold)) then
    threshold_ = threshold
  else
    threshold_ = CNST(1.0e-6)
  end if

  ALLOCATE(  r(np), np)
  ALLOCATE( rr(np), np)
  ALLOCATE( ax(np), np)
  ALLOCATE(  p(np), np)
  ALLOCATE( pp(np), np)
  ALLOCATE( ap(np), np)
  ALLOCATE(atp(np), np)

  ! Initial residue
  call op(x, ax)
  r  = b - ax
  rr = r

  ! Initial search direction
  p = r
  pp = p

  max_iter = iter
  iter = 1
  do while(iter < MAX_ITER)
    gamma = dotp(rr, r)
    err = dotp(r, r)
    if(abs(err) < THRESHOLD_) exit
    call op (p,  ap)
    call opt(pp, atp)
    alpha = gamma/dotp(pp, ap)
    r  = r  - alpha*ap
    rr = rr - alpha*atp
    x = x + alpha*p
    beta = dotp(rr, r)/gamma
    p  = r  + beta*p
    pp = rr + beta*pp
    iter = iter + 1
  end do
  if(present(residue)) residue = err

  deallocate(r, rr, ax, p, pp, ap, atp)

  call pop_sub()
end subroutine X(bi_conjugate_gradients)

!! Copyright (C) 2004 M. Marques, A. Castro, A. Rubio, G. Bertsch
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



!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/!
subroutine X(sym_conjugate_gradients)(np, x, b, op, iter, residue, threshold)
  integer, intent(in) :: np
  R_TYPE,                intent(inout) :: x(:)
  R_TYPE,                intent(in)    :: b(:)
  interface
     subroutine op(x, y)
       R_TYPE, intent(in)  :: x(:)
       R_TYPE, intent(out) :: y(:)
     end subroutine op
  end interface
  integer, intent(inout)        :: iter
  FLOAT,  optional, intent(in)  :: threshold
  R_TYPE, optional, intent(out) :: residue
  
  
  R_TYPE, allocatable :: r(:), ax(:), p(:), ap(:)
  R_TYPE  :: alpha, beta, gamma
  FLOAT   :: threshold_
  integer :: max_iter
  
  call push_sub('sym_conjugate_gradients')
  
  if(present(threshold)) then
     threshold_ = threshold
  else
     threshold_ = CNST(1.0e-6)
  endif
  
  allocate(r(np), ax(np), p(np), ap(np))
  
  ! Initial residue
  call op(x, ax)
  r = b - ax
  
  ! Initial search direction
  p = r
  
  max_iter = iter
  iter = 1
  do while(iter < max_iter)
     gamma = dot_product(r, r)
     if(abs(gamma) < THRESHOLD_) exit
     call op(p, ap)
     alpha = gamma/dot_product(p, ap)
     r = r - alpha*ap
     x = x + alpha*p
     beta = dot_product(r, r)/gamma
     p = r + beta*p
     iter  = iter + 1
  enddo
  if(present(residue)) residue = gamma
  
  deallocate(r, ax, p, ap)
  call pop_sub(); return
end subroutine X(sym_conjugate_gradients)

subroutine X(bi_conjugate_gradients)(np, x, b, op, opt, iter, residue, threshold)
  integer, intent(in) :: np
  R_TYPE,                intent(inout) :: x(:)
  R_TYPE,                intent(in)    :: b(:)
  interface
     subroutine op(x, y)
       R_TYPE, intent(in)  :: x(:)
       R_TYPE, intent(out) :: y(:)
     end subroutine op
  end interface
  interface
     subroutine opt(x, y)
       R_TYPE, intent(in)  :: x(:)
       R_TYPE, intent(out) :: y(:)
     end subroutine opt
  end interface
  integer, intent(inout) :: iter
  R_TYPE, optional, intent(out) :: residue
  FLOAT, optional, intent(in) :: threshold
  
  R_TYPE, allocatable :: r(:), rr(:), ax(:), p(:), pp(:), ap(:), atp(:)
  FLOAT :: alpha, beta, gamma, err, threshold_
  integer :: max_iter
  
  call push_sub('bi_conjugate_gradients')
  
  if(present(threshold)) then
     threshold_ = threshold
  else
     threshold_ = CNST(1.0e-6)
  endif
  
  allocate(r(np), rr(np), ax(np), p(np), pp(np), ap(np), atp(np))
  
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
     gamma = dot_product(rr, r)
     err = dot_product(r, r)
     if(abs(err) < THRESHOLD_) exit
     call op (p,  ap)
     call opt(pp, atp)
     alpha = gamma/dot_product(pp, ap)
     r  = r  - alpha*ap
     rr = rr - alpha*atp
     x = x + alpha*p
     beta = dot_product(rr, r)/gamma
     p  = r  + beta*p
     pp = rr + beta*pp
     iter = iter + 1
  enddo
  if(present(residue)) residue = err
  
  deallocate(r, rr, ax, p, pp, ap, atp)
  call pop_sub(); return
end subroutine X(bi_conjugate_gradients)

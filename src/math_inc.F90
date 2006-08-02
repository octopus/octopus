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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

! ----------------------------------------------------------------------
! This subroutine performs polynomial extrapolation (interpolation) of
! orders (1-4). It assumes:
!                      f(t) = sum_{j=0}^{order} a_j * t^k,
! and the inputs are f(0), f(-dt), f(-2*dt), ..., f(-order*dt).
!
! (I know there is a smarter and more general way to set the coefficients,
! but for the moment this works. If someday higher orders are needed, or
! I am bored, I will put a more general formula)
! ----------------------------------------------------------------------
subroutine X(extrapolate)(order, n1, n2, v, vex, dt, t)
  integer, intent(in)  :: order, n1, n2
  R_TYPE,  intent(in)  :: v(:,:,:)  ! (n1, n2, order)
  R_TYPE,  intent(out) :: vex(:,:)  ! (n1, n2)
  FLOAT,   intent(in)  :: dt, t

  integer :: j
  FLOAT :: x
  R_TYPE, allocatable :: c(:)

  x = (t/dt)
  ALLOCATE(c(0:order), order+1)

  ! I got these coefficients from mathematica...
  select case(order)
  case(1)
    c(0) = 1.0 + x
    c(1) =     - x
  case(2)
    c(0) = 1.0 + (3.0/2.0)*x + (1.0/2.0)*x**2
    c(1) =     -      2.0 *x -           x**2
    c(2) =       (1.0/2.0)*x + (1.0/2.0)*x**2
  case(3)
    c(0) = 1.0 + (11.0/6.0)*x +           x**2 + (1.0/6.0)*x**3
    c(1) =            -3.0 *x - (5.0/2.0)*x**2 - (1.0/2.0)*x**3
    c(2) =       ( 3.0/2.0)*x +      2.0 *x**2 + (1.0/2.0)*x**3
    c(3) =     - ( 1.0/3.0)*x - (1.0/2.0)*x**2 - (1.0/6.0)*x**3
  case(4)
    c(0) = 1.0 + (25.0/12.0)*x + (35.0/24.0)*x**2 + ( 5.0/12.0)*x**3 + ( 1.0/24.0)*x**4
    c(1) =     -        4.0 *x - (13.0/ 3.0)*x**2 - ( 3.0/ 2.0)*x**3 - ( 1.0/ 6.0)*x**4
    c(2) =              3.0 *x + (19.0/ 4.0)*x**2 +        2.0 *x**3 + ( 1.0/ 4.0)*x**4
    c(3) =     - ( 4.0/ 3.0)*x - ( 7.0/ 3.0)*x**2 - ( 7.0/ 6.0)*x**3 - ( 1.0/ 6.0)*x**4
    c(4) =       ( 1.0/ 4.0)*x + (11.0/24.0)*x**2 + ( 1.0/ 4.0)*x**3 + ( 1.0/24.0)*x**4
  case default
    message(1) = 'extrapolate: Unsupported order.'
    call write_fatal(1)
  end select

  vex(:,:) = R_TOTYPE(M_ZERO)
  do j = 0, order
    call lalg_axpy(n1, n2, c(j), v(:,:, j+1), vex(:,:))
  end do

  deallocate(c)
end subroutine X(extrapolate)


! ---------------------------------------------------------
subroutine X(shellsort1)(a, x)
  FLOAT, intent(inout) :: a(:)
  R_TYPE, intent(inout) :: x(:, :)

  integer :: i, j, inc, n, m
  FLOAT   :: v
  R_TYPE, allocatable :: b(:)

  n = size(a)
  m = size(x, 1)
  ALLOCATE(b(m), m)

  inc = 1
  do
    inc=3*inc+1
    if (inc > n) exit
  end do

  do
    inc=inc/3
    do i=inc+1,n
      v=a(i)
      b(:) = x(:, i)
      j=i
      do
        if (a(j-inc) <= v) exit
        !if (a(j-inc) >= v) exit
        a(j)=a(j-inc)
        x(:, j) = x(:, j-inc)
        j=j-inc
        if (j <= inc) exit
      end do
      a(j)=v
      x(:, j) = b(:)
    end do
    if (inc <= 1) exit
  end do

end subroutine X(shellsort1)


! ---------------------------------------------------------
subroutine X(shellsort2)(a, x)
  FLOAT, intent(inout) :: a(:)
  R_TYPE, intent(inout) :: x(:, :, :)

  integer :: i, j, inc, n, p, q
  FLOAT   :: v
  R_TYPE, allocatable :: b(:, :)

  n = size(a)
  p = size(x, 1)
  q = size(x, 2)
  ALLOCATE(b(p, q), p*q)

  inc = 1
  do
    inc=3*inc+1
    if (inc > n) exit
  end do

  do
    inc=inc/3
    do i=inc+1,n
      v=a(i)
      b(:, :) = x(:, :, i)
      j=i
      do
        if (a(j-inc) <= v) exit
        !if (a(j-inc) >= v) exit
        a(j)=a(j-inc)
        x(:, :, j) = x(:, :, j-inc)
        j=j-inc
        if (j <= inc) exit
      end do
      a(j)=v
      x(:, :, j) = b(:, :)
    end do
    if (inc <= 1) exit
  end do

end subroutine X(shellsort2)


! ---------------------------------------------------------
! These routines compare numbers or arrays of numbers to
! within a certain threshold (to avoid considering differences
! due to rounding as significant). They are not to be
! accessed directly, but throught the .app. operator.
! ---------------------------------------------------------
logical function X(approximate_equal)(a, b) result(app)
  R_TYPE, intent(in) :: a, b
#if defined(R_TREAL)
  app = abs(a-b) < APP_THRESHOLD
#else
  app = &
    (abs(R_REAL(a)-R_REAL(b))   < APP_THRESHOLD) .and. &
    (abs(R_AIMAG(a)-R_AIMAG(b)) < APP_THRESHOLD)
#endif
end function X(approximate_equal)


! ---------------------------------------------------------
logical function X(approximate_equal_1)(a, b) result(app)
  R_TYPE, intent(in) :: a(:), b(:)
  integer :: i
  app = .false.
  if(size(a).ne.size(b)) return
  do i = 1, size(a)
    app = X(approximate_equal)(a(i), b(i))
    if(.not.app) return
  end do
end function X(approximate_equal_1)


! ---------------------------------------------------------
logical function X(approximate_equal_2)(a, b) result(app)
  R_TYPE, intent(in) :: a(:, :), b(:, :)
  integer :: i
  app = .false.
  if(any(shape(a).ne.shape(b))) return
  do i = 1, size(a, 1)
    app = X(approximate_equal_1)(a(i, :), b(i, :))
    if(.not.app) return
  end do
end function X(approximate_equal_2)


! ---------------------------------------------------------
logical function X(approximate_equal_3)(a, b) result(app)
  R_TYPE, intent(in) :: a(:, :, :), b(:, :, :)
  integer :: i
  app = .false.
  if(any(shape(a).ne.shape(b))) return
  do i = 1, size(a, 1)
    app = X(approximate_equal_2)(a(i, :, :), b(i, :, :))
    if(.not.app) return
  end do
end function X(approximate_equal_3)


! ---------------------------------------------------------
! This is an implementation of the Parker-Traub algorithm
! for the inversion of Vandermonde matrices. 
!   F. Parker, Inverses of Vandermonde matrices, Amer. Math. 
!   Monthly, 71, 1964, p410-411
!   J. Traub, Associated polynomials and uniform methods for 
!   the solution of linear problems
!   SIAM Review, 8, No. 3, 1966, p277-301
! The algorithm inverts a Vandermonde matrix in O(N^2) operations
! 
! 2006-08-02 Heiko Appel
! ---------------------------------------------------------
subroutine X(parker_traub)(nsize, vdm_base, vdm_inverse)
  integer, intent(in)  :: nsize
  R_TYPE,  intent(in)  :: vdm_base(:)       ! Vandermonde base
  R_TYPE,  intent(out) :: vdm_inverse(:,:)  ! Inverse of Vandermonde matrix
  
  integer :: j, k
  R_TYPE, allocatable :: ap(:), an(:), q(:, :), pp(:)
  
  call push_sub('math_inc.Xparker_traub')

  ALLOCATE(an(0:nsize), nsize+1)
  ALLOCATE(ap(0:nsize), nsize+1)
  ALLOCATE(pp(1:nsize), nsize  )
  ALLOCATE( q(0:nsize, 1:nsize), nsize*(nsize+1))
  
  ! compute coefficients a_j
  ap(0) = - vdm_base(1)
  ap(1) =   R_TOTYPE(M_ONE)
  
  do k = 2, nsize
    an(0) = - vdm_base(k)*ap(0)
    an(k) = ap(k-1)
    do j = 1, k-1
      an(j) = ap(j-1) - vdm_base(k)*ap(j)
    end do
    ap = an
  end do
  
  ! compute coefficients q_k,j
  q = R_TOTYPE(M_ZERO)
  
  do j = 1, nsize
    q(0, j) = R_TOTYPE(M_ONE)
    do k = 1, nsize
      q(k, j) = vdm_base(j)*q(k-1, j) + an(nsize-k) 
    end do
  end do
  
  ! compute derivative polynomial at base points
  ! this amounts to Parkers variant of the algorithm
  pp = R_TOTYPE(M_ONE)
  do j = 1, nsize
    do k = 1, nsize
      if(k.ne.j) then
        pp(j) = pp(j) * ( vdm_base(j) - vdm_base(k) )
      end if
    end do
  end do
  
  ! compute inverse 
  do j = 1, nsize
    do k = 0, nsize-1 
      vdm_inverse(k+1, j) = q(nsize-1 - k, j) / pp(j) 
    end do
  end do
  
  deallocate(q, pp, ap, an)
  
  call pop_sub()
end subroutine X(parker_traub)


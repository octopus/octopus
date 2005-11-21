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

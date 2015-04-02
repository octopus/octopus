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
!! $Id$

! ---------------------------------------------------------
subroutine X(shellsort1)(a, x)
  FLOAT,  intent(inout) :: a(:)
  R_TYPE, intent(inout) :: x(:, :)

  integer :: i, j, inc, n, m
  FLOAT   :: v
  R_TYPE, allocatable :: b(:)

  PUSH_SUB(X(shellsort1))

  n = size(a)
  m = size(x, 1)
  SAFE_ALLOCATE(b(1:m))

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

  SAFE_DEALLOCATE_A(b)
  POP_SUB(X(shellsort1))
end subroutine X(shellsort1)


! ---------------------------------------------------------
subroutine X(shellsort2)(a, x)
  FLOAT,  intent(inout) :: a(:)
  R_TYPE, intent(inout) :: x(:, :, :)

  integer :: i, j, inc, n, p, q
  FLOAT   :: v
  R_TYPE, allocatable :: b(:, :)

  PUSH_SUB(X(shellsort2))

  n = size(a)
  p = size(x, 1)
  q = size(x, 2)
  SAFE_ALLOCATE(b(1:p, 1:q))

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

  SAFE_DEALLOCATE_A(b)
  POP_SUB(X(shellsort2))
end subroutine X(shellsort2)

! ---------------------------------------------------------
!> sort the eigenvectors according to eigenvalues
!! with increasing absolute value
subroutine X(matrix_sort)(np, matrix, eigenvals)
  integer, intent(in)    :: np
  R_TYPE,  intent(inout) :: matrix(:, :)
  R_TYPE,  intent(inout) :: eigenvals(:)

  integer              :: i
  R_TYPE, allocatable  :: unsorted_matrix(:, :), unsorted_eigenvals(:)
  FLOAT, allocatable   :: abs_e(:)
  integer, allocatable :: index(:)

  PUSH_SUB(X(matrix_sort))

  SAFE_ALLOCATE( abs_e(1:np) )
  SAFE_ALLOCATE( index(1:np) )
  SAFE_ALLOCATE( unsorted_matrix(1:np, 1:np) )
  SAFE_ALLOCATE( unsorted_eigenvals(1:np) )

  unsorted_matrix(:, :) = matrix(:, :)
  unsorted_eigenvals(:) = eigenvals(:)
  abs_e(:) = abs(unsorted_eigenvals(:))
  call sort(abs_e, index)
  do i=1, np
    eigenvals(i) = unsorted_eigenvals(index(i))
    matrix(:, i) = unsorted_matrix(:, index(i))
  end do
  SAFE_DEALLOCATE_A(abs_e)
  SAFE_DEALLOCATE_A(index)
  SAFE_DEALLOCATE_A(unsorted_matrix)
  SAFE_DEALLOCATE_A(unsorted_eigenvals)

  POP_SUB(X(matrix_sort))
end subroutine X(matrix_sort)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

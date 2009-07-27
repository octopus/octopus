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
!! $Id$

subroutine X(interpolate_2)(xa, ya, x, y)
  FLOAT,  intent(in)  :: xa(:)
  R_TYPE, intent(in)  :: ya(:, :, :)
  FLOAT,  intent(in)  :: x
  R_TYPE, intent(out) :: y(:, :)

  integer :: n, n1, n2, i, k
  FLOAT, allocatable :: c(:)

  n = size(xa)
  n1 = size(ya, 1)
  n2 = size(ya, 2)

  SAFE_ALLOCATE(c(1:n))

  call interpolation_coefficients(n, xa, x, c)

  do k = 1, n2

    y(:, k) = M_ZERO

    do i = 1, n
      call lalg_axpy(n1, c(i), ya(:, k, i), y(:, k)) 
    end do

  end do

  SAFE_DEALLOCATE_A(c)
end subroutine X(interpolate_2)


subroutine X(interpolate_1)(xa, ya, x, y)
  FLOAT,  intent(in)  :: xa(:)
  R_TYPE, intent(in)  :: ya(:, :)
  FLOAT,  intent(in)  :: x
  R_TYPE, intent(out) :: y(:)

  integer :: n, n1, i
  FLOAT, allocatable :: c(:)

  n = size(xa)
  n1 = size(ya, 1)
  SAFE_ALLOCATE(c(1:n))

  call interpolation_coefficients(n, xa, x, c)

  y(:) = M_ZERO
  do i = 1, n
    call lalg_axpy(n1, c(i), ya(:, i), y(:)) 
  end do

  SAFE_DEALLOCATE_A(c)
end subroutine X(interpolate_1)


subroutine X(interpolate_0)(xa, ya, x, y)
  FLOAT,  intent(in)  :: xa(:)
  R_TYPE, intent(in)  :: ya(:)
  FLOAT,  intent(in)  :: x
  R_TYPE, intent(out) :: y

  integer :: n, i
  FLOAT, allocatable :: c(:)

  n = size(xa)
  SAFE_ALLOCATE(c(1:n))

  call interpolation_coefficients(n, xa, x, c)

  y = M_ZERO
  do i = 1, n
    y = y + c(i)*ya(i)
  end do

  SAFE_DEALLOCATE_A(c)
end subroutine X(interpolate_0)

! ---------------------------------------------------------
subroutine X(shellsort1)(a, x)
  FLOAT, intent(inout) :: a(:)
  R_TYPE, intent(inout) :: x(:, :)

  integer :: i, j, inc, n, m
  FLOAT   :: v
  R_TYPE, allocatable :: b(:)

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

  SAFE_ALLOCATE(an(0:nsize))
  SAFE_ALLOCATE(ap(0:nsize))
  SAFE_ALLOCATE(pp(1:nsize))
  SAFE_ALLOCATE( q(0:nsize, 1:nsize))
  
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
  
  SAFE_DEALLOCATE_A(q)
  SAFE_DEALLOCATE_A(pp)
  SAFE_DEALLOCATE_A(ap)
  SAFE_DEALLOCATE_A(an)
  
  call pop_sub()
end subroutine X(parker_traub)


! ---------------------------------------------------------
! Newton-Raphson method for matrices. Starting from a 
! sufficiently good initial guess, this routine can be used 
! to polish inverse matrices.
! 
! 2006-08-02 Heiko Appel
! ---------------------------------------------------------
subroutine X(matrix_newton_raphson)(nsteps, nsize, a, b)
  integer, intent(in)    :: nsteps, nsize
  R_TYPE,  intent(in)    :: a(:,:)  ! matrix whose inverse to improve
  R_TYPE,  intent(inout) :: b(:,:)  ! initial guess for the inverse
  
  integer :: is
  R_TYPE, allocatable :: ab(:,:), bab(:,:)
  
  call push_sub('math_inc.Xmatrix_newton_raphson')

  SAFE_ALLOCATE( ab(1:nsize, 1:nsize))
  SAFE_ALLOCATE(bab(1:nsize, 1:nsize))

  do is = 1, nsteps

    ab  = R_TOTYPE(M_ZERO)
    bab = R_TOTYPE(M_ZERO)
    ab  = matmul(a, b)
    bab = matmul(b, ab)

    ! Newton-Raphson step
    b = R_TOTYPE(M_TWO)*b - bab 

  end do

  SAFE_DEALLOCATE_A(ab)
  SAFE_DEALLOCATE_A(bab)

  call pop_sub()
end subroutine X(matrix_newton_raphson)


! ---------------------------------------------------------
! gives a measure for the quality of the inverse
FLOAT function X(matrix_inv_residual)(nsize, a, b) result(residual)
  integer, intent(in) :: nsize
  R_TYPE,  intent(in) :: a(:,:), b(:,:) ! matrix and approximate inverse

  integer :: i, j
  R_TYPE, allocatable :: ab(:,:)

  SAFE_ALLOCATE(ab(1:nsize, 1:nsize))

  ab = R_TOTYPE(M_ZERO)
  ab = matmul(a, b)

  residual = M_ZERO
  do i = 1, nsize   
    do j = 1, nsize
      if(i.ne.j) then 
        residual = residual + abs(ab(i, j)) 
      end if
    end do
  end do

  ! account for the size of the matrices
  residual = residual / nsize**2

  SAFE_DEALLOCATE_A(ab)

end function X(matrix_inv_residual)

pure function X(cross_product)(a, b) result(c)
  R_TYPE, intent(in) :: a(1:3)
  R_TYPE, intent(in) :: b(1:3)
  R_TYPE :: c(1:3)

  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)

end function X(cross_product)


! ---------------------------------------------------------
! Caclulate infinity-norm of matrix.
FLOAT function X(infinity_norm)(matrix)
  R_TYPE, intent(in) :: matrix(:, :)

  integer :: m_min, m_max, i
  FLOAT   :: norm_old, norm

  call push_sub('math_inc.infinity_norm')

  norm = 0

  m_min = lbound(matrix, 1)
  m_max = ubound(matrix, 1)
  do i = m_min, m_max
    norm_old = norm
    norm     = sum(abs(matrix(i, :)))
    norm     = max(norm, norm_old)
  end do
  
  X(infinity_norm) = norm

  call pop_sub()
end function X(infinity_norm)


! ---------------------------------------------------------
! Takes the average of the matrix and its transposed.
subroutine X(matrix_symmetric_average)(matrix, np)
  R_TYPE,  intent(inout) :: matrix(np, np)
  integer, intent(in)    :: np

  integer             :: j
  R_TYPE, allocatable :: tmp(:, :)

  call push_sub('math_inc.symmetric_average')

  SAFE_ALLOCATE(tmp(1:np, 1:np))

  if(np.gt.1) then
    do j = 1, np
      tmp(:, j) = M_HALF*(matrix(j, :) + matrix(:, j))
    end do
    matrix = tmp
  end if

  SAFE_DEALLOCATE_A(tmp)
  call pop_sub()
end subroutine X(matrix_symmetric_average)


! ---------------------------------------------------------
! Copies the upper triangular matrix into the lower part
subroutine X(matrix_symmetrize)(matrix, np)
  R_TYPE,  intent(inout) :: matrix(np, np)
  integer, intent(in)    :: np

  integer :: j

  call push_sub('math_inc.Xmatrix_symmetrize')

  do j = 1, np-1
    matrix(j+1:np,j) = matrix(j,j+1:np)
  end do

  call pop_sub()
end subroutine X(matrix_symmetrize)


! ---------------------------------------------------------
! sort the eigenvectors accordingly to eigenvalues
! with increasing absolute value
subroutine X(matrix_sort)(np, matrix, eigenvals)
  integer, intent(in)    :: np
  R_TYPE,  intent(inout) :: matrix(1:np, 1:np)
  R_TYPE,  intent(inout) :: eigenvals(1:np)

  integer              :: i
  R_TYPE, allocatable  :: unsorted_matrix(:, :), unsorted_eigenvals(:)
  FLOAT, allocatable   :: abs_e(:)
  integer, allocatable :: index(:)

  call push_sub('math_inc.X(matrix_sort)')

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

  call pop_sub()
end subroutine X(matrix_sort)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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

subroutine X(interpolate_2)(xa, ya, x, y)
  FLOAT,  intent(in)  :: xa(:)
  R_TYPE, intent(in)  :: ya(:, :, :)
  FLOAT,  intent(in)  :: x
  R_TYPE, intent(out) :: y(:, :)

  integer :: n, n1, n2, i, k
  FLOAT, allocatable :: c(:)

  PUSH_SUB(X(interpolate_2))

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
  POP_SUB(X(interpolate_2))
end subroutine X(interpolate_2)


! ---------------------------------------------------------
subroutine X(interpolate_1)(xa, ya, x, y)
  FLOAT,  intent(in)  :: xa(:)
  R_TYPE, intent(in)  :: ya(:, :)
  FLOAT,  intent(in)  :: x
  R_TYPE, intent(out) :: y(:)

  integer :: n, n1, i
  FLOAT, allocatable :: c(:)

  PUSH_SUB(X(interpolate_1))

  n = size(xa)
  n1 = size(ya, 1)
  SAFE_ALLOCATE(c(1:n))

  call interpolation_coefficients(n, xa, x, c)

  y(:) = M_ZERO
  do i = 1, n
    call lalg_axpy(n1, c(i), ya(:, i), y(:))
  end do

  SAFE_DEALLOCATE_A(c)
  POP_SUB(X(interpolate_1))
end subroutine X(interpolate_1)


! ---------------------------------------------------------
subroutine X(interpolate_0)(xa, ya, x, y)
  FLOAT,  intent(in)  :: xa(:)
  R_TYPE, intent(in)  :: ya(:)
  FLOAT,  intent(in)  :: x
  R_TYPE, intent(out) :: y

  integer :: n, i
  FLOAT, allocatable :: c(:)

  ! no push_sub, called too frequently

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
!> These routines compare numbers or arrays of numbers to
!! within a certain threshold (to avoid considering differences
!! due to rounding as significant). They are not to be
!! accessed directly, but through the .app. operator.
logical function X(approximately_equal)(a, b) result(app)
  R_TYPE, intent(in) :: a, b

  PUSH_SUB(X(approximately_equal))
#if defined(R_TREAL)
  app = abs(a-b) < APP_THRESHOLD
#else
  app = &
    (abs(R_REAL(a)-R_REAL(b))   < APP_THRESHOLD) .and. &
    (abs(R_AIMAG(a)-R_AIMAG(b)) < APP_THRESHOLD)
#endif

  POP_SUB(X(approximately_equal))
end function X(approximately_equal)


! ---------------------------------------------------------
logical function X(approximately_equal_1)(a, b) result(app)
  R_TYPE, intent(in) :: a(:), b(:)

  integer :: i

  PUSH_SUB(X(approximately_equal_1))

  app = .false.
  if(size(a) /= size(b)) then
    POP_SUB(X(approximately_equal_1))
    return
  end if
  do i = 1, size(a)
    app = X(approximately_equal)(a(i), b(i))
    if(.not.app) then
      POP_SUB(X(approximately_equal_1))
      return
    end if
  end do

  POP_SUB(X(approximately_equal_1))
end function X(approximately_equal_1)


! ---------------------------------------------------------
logical function X(approximately_equal_2)(a, b) result(app)
  R_TYPE, intent(in) :: a(:, :), b(:, :)

  integer :: i

  PUSH_SUB(X(approximately_equal_2))

  app = .false.
  if(any(shape(a) /= shape(b))) then
    POP_SUB(X(approximately_equal_2))
    return
  end if
  do i = 1, size(a, 1)
    app = X(approximately_equal_1)(a(i, :), b(i, :))
    if(.not.app) then
      POP_SUB(X(approximately_equal_2))
      return
    end if
  end do

  POP_SUB(X(approximately_equal_2))
end function X(approximately_equal_2)


! ---------------------------------------------------------
logical function X(approximately_equal_3)(a, b) result(app)
  R_TYPE, intent(in) :: a(:, :, :), b(:, :, :)

  integer :: i

  PUSH_SUB(X(approximately_equal_3))

  app = .false.
  if(any(shape(a) /= shape(b))) then
    POP_SUB(X(approximately_equal_3))
    return
  end if
  do i = 1, size(a, 1)
    app = X(approximately_equal_2)(a(i, :, :), b(i, :, :))
    if(.not.app) then
      POP_SUB(X(approximately_equal_3))
      return
    end if
  end do

  POP_SUB(X(approximately_equal_3))
end function X(approximately_equal_3)


! ---------------------------------------------------------
pure function X(cross_product)(a, b) result(c)
  R_TYPE, intent(in) :: a(:) !< (3)
  R_TYPE, intent(in) :: b(:) !< (3)

  R_TYPE :: c(1:3)

  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)

end function X(cross_product)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

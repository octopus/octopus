! interchanges two vectors.
subroutine dlalg_swap(n, dx, dy)
  integer, intent(in)    :: n
  FLOAT,   intent(inout) :: dx(*), dy(*)

  integer :: i
  FLOAT :: s

  do i = 1, n
    s = dx(i)
    dx(i) = dy(i)
    dy(i) = s
  end do

end subroutine dlalg_swap

subroutine zlalg_swap(n, zx, zy)
  integer, intent(in)    :: n
  CMPLX,   intent(inout) :: zx(*), zy(*)

  integer :: i
  CMPLX :: s

  do i = 1, n
    s = zx(i)
    zx(i) = zy(i)
    zy(i) = s
  end do

end subroutine zlalg_swap

! scales a vector by a constant
subroutine dlalg_scal(n, da, dx)
  integer, intent(in)    :: n
  FLOAT,   intent(in)    :: da
  FLOAT,   intent(inout) :: dx(*)

  dx(1:n) = dx(1:n)*da

end subroutine dlalg_scal

subroutine zlalg_scal(n, za, zx)
  integer, intent(in)    :: n
  CMPLX,   intent(in)    :: za
  CMPLX,   intent(inout) :: zx(*)

  zx(1:n) = zx(1:n)*za

end subroutine zlalg_scal

! constant times a vector plus a vector
subroutine dlalg_axpy(n, da, dx, dy)
  integer, intent(in)    :: n
  FLOAT,   intent(in)    :: da, dx(*)
  FLOAT,   intent(inout) :: dy(*)

  dy(1:n) = dy(1:n) + da*dx(1:n)

end subroutine dlalg_axpy

subroutine zlalg_axpy(n, za, zx, zy)
  integer, intent(in)    :: n
  CMPLX,   intent(in)    :: za, zx(*)
  CMPLX,   intent(inout) :: zy(*)

  zy(1:n) = zy(1:n) + za*zx(1:n)

end subroutine zlalg_axpy

! forms the dot product of two vectors
FLOAT function dlalg_dot(n, dx, dy)
  integer, intent(in) :: n
  FLOAT,   intent(in) :: dx(*), dy(*)

  dlalg_dot = dot_product(dx(1:n), dy(1:n))

end function dlalg_dot

CMPLX function zlalg_dot(n, zx, zy)
  integer, intent(in) :: n
  CMPLX,   intent(in) :: zx(*), zy(*)

  zlalg_dot = dot_product(zx(1:n), zy(1:n))

end function zlalg_dot

! returns the euclidean norm of a vector
FLOAT function dlalg_nrm2(n, dx)
  integer, intent(in) :: n
  FLOAT,   intent(in) :: dx(*)

  dlalg_nrm2 = sqrt(dot_product(dx(1:n), dx(1:n)))

end function dlalg_nrm2

FLOAT function zlalg_nrm2(n, zx)
  integer, intent(in) :: n
  CMPLX,   intent(in) :: zx(*)

  zlalg_nrm2 = sqrt(real(dot_product(zx(1:n), zx(1:n)), PRECISION))

end function zlalg_nrm2

! copies a vector, x, to a vector, y
subroutine dlalg_copy(n, dx, dy)
  integer, intent(in)  :: n
  FLOAT,   intent(in)  :: dx(*)
  FLOAT,   intent(out) :: dy(*)

  dy(1:n) = dx(1:n)

end subroutine dlalg_copy

subroutine zlalg_copy(n, zx, zy)
  integer, intent(in)  :: n
  CMPLX,   intent(in)  :: zx(*)
  CMPLX,   intent(out) :: zy(*)

  zy(1:n) = zx(1:n)

end subroutine zlalg_copy

! matrix-matrix multiplication plus matrix
subroutine dlalg_gemm(m, n, k, alpha, a, b, beta, c)
  integer,      intent(in)    :: m, n, k
  FLOAT,        intent(in)    :: alpha, beta
  FLOAT,        intent(in)    :: a(*)
  FLOAT,        intent(in)    :: b(*)
  FLOAT,        intent(inout) :: c(*)

  integer :: i, j, l, ls, js1, js2
  real :: temp

  js1 = -m
  js2 = -k
  do j = 1, n
    js1 = js1 + m
    if (beta == 0.0) then
      do i = 1, m
        c(js1 + i) = 0.0
      end do
    else if(beta /= 1.0) then
      do i = 1, m
        c(js1 + i) = beta*c(js1 + i)
      end do
    end if

    js2 = js2 + k
    ls = -m
    do l = 1, k
      ls = ls + m
      if (b(l + js2) /= 0.0) then
        temp = alpha*b(js2 + l)
        do i = 1, m
          c(js1 + i) = c(js1 + i) + temp*a(i + ls)
        end do
      end if
    end do
  end do

end subroutine dlalg_gemm

subroutine zlalg_gemm(m, n, k, alpha, a, b, beta, c)
  integer,      intent(in)    :: m, n, k
  CMPLX,        intent(in)    :: alpha, beta
  CMPLX,        intent(in)    :: a(*)
  CMPLX,        intent(in)    :: b(*)
  CMPLX,        intent(inout) :: c(*)

  integer :: i, j, l, ls, js1, js2
  real :: temp

  js1 = -m
  js2 = -k
  do j = 1, n
    js1 = js1 + m
    if (beta == 0.0) then
      do i = 1, m
        c(js1 + i) = 0.0
      end do
    else if(beta /= 1.0) then
      do i = 1, m
        c(js1 + i) = beta*c(js1 + i)
      end do
    end if

    js2 = js2 + k
    ls = -m
    do l = 1, k
      ls = ls + m
      if (b(l + js2) /= 0.0) then
        temp = alpha*b(js2 + l)
        do i = 1, m
          c(js1 + i) = c(js1 + i) + temp*a(i + ls)
        end do
      end if
    end do
  end do

end subroutine zlalg_gemm
 
! matrix-vector multiplication plus vector
subroutine dlalg_gemv(m, n, alpha, a, x, beta, y)
  integer,      intent(in)    :: m, n
  FLOAT,        intent(in)    :: alpha, beta
  FLOAT,        intent(in)    :: a(*)
  FLOAT,        intent(in)    :: x(*)
  FLOAT,        intent(inout) :: y(*)

  integer :: i, j, is

  do i = 1, m
    y(i) = beta*y(i)
  end do

  is = -m
  do j = 1, n
    is = is + m
    do i = 1, m
      y(i) = y(i) + alpha*a(is + i)*x(j)
    end do
  end do

end subroutine dlalg_gemv

subroutine zlalg_gemv(m, n, alpha, a, x, beta, y)
  integer,      intent(in)    :: m, n
  CMPLX,        intent(in)    :: alpha, beta
  CMPLX,        intent(in)    :: a(*)
  CMPLX,        intent(in)    :: x(*)
  CMPLX,        intent(inout) :: y(*)

  integer :: i, j, is

  do i = 1, m
    y(i) = beta*y(i)
  end do

  is = -m
  do j = 1, n
    is = is + m
    do i = 1, m
      y(i) = y(i) + alpha*a(is + i)*x(j)
    end do
  end do

end subroutine zlalg_gemv

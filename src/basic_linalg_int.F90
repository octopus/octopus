! interchanges two vectors.
subroutine d_swap(n, dx, dy)
  integer, intent(in)    :: n
  FLOAT,   intent(inout) :: dx(*), dy(*)

  integer :: i
  FLOAT :: s

  do i = 1, n
    s = dx(i)
    dx(i) = dy(i)
    dy(i) = s
  end do

end subroutine d_swap

subroutine z_swap(n, zx, zy)
  integer, intent(in)    :: n
  CMPLX,   intent(inout) :: zx(*), zy(*)

  integer :: i
  CMPLX :: s

  do i = 1, n
    s = dx(i)
    zx(i) = zy(i)
    zy(i) = s
  end do

end subroutine z_swap

! scales a vector by a constant
subroutine d_scal(n, da, dx)
  integer, intent(in)    :: n
  FLOAT,   intent(in)    :: da
  FLOAT,   intent(inout) :: dx(*)

  dx(1:n) = dx(1:n)*da

end subroutine d_scal

subroutine z_scal(n, za, zx)
  integer, intent(in)    :: n
  CMPLX,   intent(in)    :: za
  CMPLX,   intent(inout) :: zx(*)

  zx(1:n) = zx(1:n)*za

end subroutine z_scal

! constant times a vector plus a vector
subroutine d_axpy(n, da, dx, dy)
  integer, intent(in)    :: n
  FLOAT,   intent(in)    :: da, dx(*)
  FLOAT,   intent(inout) :: dy(*)

  dy(1:n) = dy(1:n) + da*dx(1:n)

end subroutine d_axpy

subroutine z_axpy(n, za, zx, zy)
  integer, intent(in)    :: n
  CMPLX,   intent(in)    :: za, zx(*)
  CMPLX,   intent(inout) :: zy(*)

  zy(1:n) = zy(1:n) + za*zx(1:n)

end subroutine z_axpy

! forms the dot product of two vectors
FLOAT function d_dot(n, dx, dy)
  integer, intent(in) :: n
  FLOAT,   intent(in) :: dx(*), dy(*)

  d_dot = dot_product(dx(1:n), dy(1:n))

end function d_dot

CMPLX function z_dot(n, zx, zy)
  integer, intent(in) :: n
  CMPLX,   intent(in) :: zx(*), zy(*)

  z_dot = dot_product(zx(1:n), zy(1:n))

end function z_dot

! returns the euclidean norm of a vector
FLOAT function d_nrm2(n, dx)
  integer, intent(in) :: n
  FLOAT,   intent(in) :: dx(*)

  d_nrm2 = sqrt(dot_product(dx(1:n), dx(1:n)))

end function d_nrm2

FLOAT function dz_nrm2(n, zx)
  integer, intent(in) :: n
  CMPLX,   intent(in) :: zx(*)

  dz_nrm2 = sqrt(real(dot_product(zx(1:n), zx(1:n)), PRECISION))

end function dz_nrm2

! copies a vector, x, to a vector, y
subroutine d_copy(n, dx, dy)
  integer, intent(in)  :: n
  FLOAT,   intent(in)  :: dx(*)
  FLOAT,   intent(out) :: dy(*)

  dy(1:n) = dx(1:n)

end subroutine d_copy

subroutine z_copy(n, zx, zy)
  integer, intent(in)  :: n
  CMPLX,   intent(in)  :: zx(*)
  CMPLX,   intent(out) :: zy(*)

  zy(1:n) = zx(1:n)

end subroutine z_copy

! matrix-matrix multiplication plus matrix
subroutine d_gemm(transa, transb, m, n, k, alpha, a, b, beta, c)
  character(1), intent(in)    :: transa, transb
  integer,      intent(in)    :: m, n, k
  FLOAT,        intent(in)    :: alpha, beta
  FLOAT,        intent(in)    :: a(*)
  FLOAT,        intent(in)    :: b(*)
  FLOAT,        intent(inout) :: c(*)

  select case (transa)
  case ('N','n')
    c(1:m*n) = alpha*matmul(a(1:m*k), b(1:k*n)) + beta*c(1:m*n)
  case ('T','t','C','c')
    c(1:m*n) = alpha*matmul(transpose(a(1:m*k)), b(1:k*n)) + beta*c(1:m*n)
  end select

end subroutine d_gemm

subroutine z_gemm(transa, transb, m, n, k, alpha, a, b, beta, c)
  character(1), intent(in)    :: transa, transb
  integer,      intent(in)    :: m, n, k
  CMPLX,        intent(in)    :: alpha, beta
  CMPLX,        intent(in)    :: a(*)
  CMPLX,        intent(in)    :: b(*)
  CMPLX,        intent(inout) :: c(*)

  select case (transa)
  case ('N','n')
    c(1:m*n) = alpha*matmul(a(1:m*k), b(1:k*n)) + beta*c(1:m*n)
  case ('T','t')
    c(1:m*n) = alpha*matmul(transpose(a(1:m*k)), b(1:k*n)) + beta*c(1:m*n)
  case('C','c')
    c(1:m*n) = alpha*matmul(transpose(conjg(a(1:m*k))), b(1:k*n)) + beta*c(1:m*n)
  end select

end subroutine z_gemm
 
! matrix-vector multiplication plus vector
subroutine d_gemv(m, n, alpha, a, x, beta, y)
  integer,      intent(in)    :: m, n
  FLOAT,        intent(in)    :: alpha, beta
  FLOAT,        intent(in)    :: a(*)
  FLOAT,        intent(in)    :: x(*)
  FLOAT,        intent(inout) :: y(*)

  y(1:m) = alpha*matmul(a(1:n*m), x(1:n)) + beta*y(1:m)

end subroutine d_gemv

subroutine z_gemv(m, n, alpha, a, x, beta, y)
  integer,      intent(in)    :: m, n
  CMPLX,        intent(in)    :: alpha, beta
  CMPLX,        intent(in)    :: a(*)
  CMPLX,        intent(in)    :: x(*)
  CMPLX,        intent(inout) :: y(*)

  y(1:m) = alpha*matmul(a(1:n*m), x(1:n)) + beta*y(1:m)

end subroutine z_gemv

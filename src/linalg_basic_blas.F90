#if defined(SINGLE_PRECISION)
#  define DBLAS(x) s ## x
#  define ZBLAS(x) c ## x
#  define ZNRM2    scnrm2
#else
#  define DBLAS(x) d ## x
#  define ZBLAS(x) z ## x
#  define ZNRM2    dznrm2
#endif

! interchanges two vectors.
subroutine dlalg_swap(n, dx, dy)
  integer, intent(in)    :: n
  FLOAT,   intent(inout) :: dx(*), dy(*)

  interface
    subroutine  DBLAS(swap) (n, dx, incx, dy, incy)
      integer, intent(in)    :: n, incx, incy
      FLOAT,   intent(inout) :: dx, dy ! dx(n), dy(n)
    end subroutine DBLAS(swap)
  end interface

  call DBLAS(swap) (n, dx(1), 1, dy(1), 1)

end subroutine dlalg_swap

subroutine zlalg_swap(n, zx, zy)
  integer, intent(in)    :: n
  CMPLX,   intent(inout) :: zx(*), zy(*)

  interface
    subroutine  ZBLAS(swap) (n, zx, incx, zy, incy)
      integer, intent(in)    :: n, incx, incy
      CMPLX,   intent(inout) :: zx, zy ! zx(n), zy(n)
    end subroutine ZBLAS(swap)
  end interface

  call ZBLAS(swap) (n, zx(1), 1, zy(1), 1)

end subroutine zlalg_swap

! scales a vector by a constant
subroutine dlalg_scal(n, da, dx)
  integer, intent(in)    :: n
  FLOAT,   intent(in)    :: da
  FLOAT,   intent(inout) :: dx(*)

  interface
    subroutine DBLAS(scal) (n, da, dx, incx)
      integer, intent(in)    :: n, incx
      FLOAT,   intent(in)    :: da
      FLOAT,   intent(inout) :: dx ! dx(n)
    end subroutine DBLAS(scal)
  end interface

  call DBLAS(scal) (n, da, dx(1), 1)

end subroutine dlalg_scal

subroutine zlalg_scal(n, za, zx)
  integer, intent(in)    :: n
  CMPLX,   intent(in)    :: za
  CMPLX,   intent(inout) :: zx(*)

  interface
    subroutine ZBLAS(scal) (n, za, zx, incx)
      integer, intent(in)    :: n, incx
      CMPLX,   intent(in)    :: za
      CMPLX,   intent(inout) :: zx ! zx(n)
    end subroutine ZBLAS(scal)
  end interface

  call ZBLAS(scal) (n, za, zx(1), 1)

end subroutine zlalg_scal

! constant times a vector plus a vector
subroutine dlalg_axpy(n, da, dx, dy)
  integer, intent(in)    :: n
  FLOAT,   intent(in)    :: da
  FLOAT,   intent(IN)    :: dx(*)
  FLOAT,   intent(inout) :: dy(*)

  interface
    subroutine DBLAS(axpy) (n, da, dx, incx, dy, incy)
      integer, intent(in)    :: n, incx, incy
      FLOAT,   intent(in)    :: da, dx ! dx(n)
      FLOAT,   intent(inout) :: dy     ! dy(n)
    end subroutine DBLAS(axpy)
  end interface

  call DBLAS(axpy) (n, da, dx(1), 1, dy(1), 1)

end subroutine dlalg_axpy

subroutine zlalg_axpy(n, za, zx, zy)
  integer, intent(in)    :: n
  CMPLX,   intent(in)    :: za
  CMPLX,   intent(IN)    :: zx(*)
  CMPLX,   intent(inout) :: zy(*)

  interface
    subroutine ZBLAS(axpy) (n, za, zx, incx, zy, incy)
      integer, intent(in)    :: n, incx, incy
      CMPLX,   intent(in)    :: za, zx ! zx(n)
      CMPLX,   intent(inout) :: zy     ! zy(n)
    end subroutine ZBLAS(axpy)
  end interface

  call ZBLAS(axpy) (n, za, zx(1), 1, zy(1), 1)

end subroutine zlalg_axpy

! forms the dot product of two vectors
FLOAT function dlalg_dot(n, dx, dy)
  integer, intent(in) :: n
  FLOAT,   intent(IN) :: dx(*), dy(*)

  interface
    FLOAT function DBLAS(dot) (n, dx, incx, dy, incy)
      integer, intent(in) :: n, incx, incy
      FLOAT,   intent(in) :: dx, dy ! dx(n), dy(n)
    end function DBLAS(dot)
  end interface

  dlalg_dot = DBLAS(dot) (n, dx(1), 1, dy(1), 1)

end function dlalg_dot

CMPLX function zlalg_dot(n, zx, zy)
  integer, intent(in) :: n
  CMPLX,   intent(IN) :: zx(*), zy(*)

  interface
    CMPLX function ZBLAS(dotc) (n, zx, incx, zy, incy)
      integer, intent(in) :: n, incx, incy
      CMPLX,   intent(in) :: zx, zy ! zx(n), zy(n)
    end function ZBLAS(dotc)
  end interface

  zlalg_dot = ZBLAS(dotc) (n, zx(1), 1, zy(1), 1)

end function zlalg_dot

! returns the euclidean norm of a vector
FLOAT function dlalg_nrm2(n, dx)
  integer, intent(in) :: n
  FLOAT,   intent(IN) :: dx(*)

  interface
    FLOAT function DBLAS(nrm2) (n, dx, incx)
      integer, intent(in) :: n, incx
      FLOAT,   intent(in) :: dx ! dx(n)
    end function DBLAS(nrm2)
  end interface

  dlalg_nrm2 = DBLAS(nrm2) (n, dx(1), 1)

end function dlalg_nrm2

FLOAT function zlalg_nrm2(n, zx)
  integer, intent(in) :: n
  CMPLX,   intent(IN) :: zx(*)

  interface
    FLOAT function ZNRM2 (n, zx, incx)
      integer, intent(in) :: n, incx
      CMPLX,   intent(in) :: zx ! zx(n)
    end function ZNRM2
  end interface

  zlalg_nrm2 = ZNRM2 (n, zx(1), 1)

end function zlalg_nrm2

! copies a vector, x, to a vector, y
subroutine dlalg_copy(n, dx, dy)
  integer, intent(in)  :: n
  FLOAT,   intent(IN)  :: dx(*)
  FLOAT,   intent(out) :: dy(*)

  interface
    subroutine DBLAS(copy) (n, dx, incx, dy, incy)
      integer, intent(in)  :: n, incx, incy
      FLOAT,   intent(in)  :: dx ! dx(n)
      FLOAT,   intent(out) :: dy ! dy(n)
    end subroutine DBLAS(copy)
  end interface

  call DBLAS(copy) (n, dx(1), 1, dy(1), 1)

end subroutine dlalg_copy

subroutine zlalg_copy(n, zx, zy)
  integer, intent(in)  :: n
  CMPLX,   intent(IN)  :: zx(*)
  CMPLX,   intent(out) :: zy(*)

  interface
    subroutine ZBLAS(copy) (n, zx, incx, zy, incy)
      integer,    intent(in)  :: n, incx, incy
      CMPLX,      intent(in)  :: zx ! dz(n)
      CMPLX,      intent(out) :: zy ! zy(n)
    end subroutine ZBLAS(copy)
  end interface

  call ZBLAS(copy) (n, zx(1), 1, zy(1), 1)

end subroutine zlalg_copy

! matrix-matrix multiplication plus matrix
subroutine dlalg_gemm(m, n, k, alpha, a, b, beta, c)
  integer,      intent(in)    :: m, n, k
  FLOAT,        intent(in)    :: alpha, beta
  FLOAT,        intent(IN)    :: a(*)
  FLOAT,        intent(IN)    :: b(*)
  FLOAT,        intent(inout) :: c(*)

  interface
    subroutine DBLAS(gemm) (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      FLOAT,        intent(in)    :: alpha, beta
      FLOAT,        intent(in)    :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      FLOAT,        intent(in)    :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      FLOAT,        intent(inout) :: c ! c(ldc,n) 
    end subroutine DBLAS(gemm)
  end interface

  call DBLAS(gemm) ('N', 'N', m, n, k, alpha, a(1), m, b(1), k, beta, c(1), m)

end subroutine dlalg_gemm

subroutine zlalg_gemm(m, n, k, alpha, a, b, beta, c)
  integer,      intent(in)    :: m, n, k
  CMPLX,        intent(in)    :: alpha, beta
  CMPLX,        intent(IN)    :: a(*)
  CMPLX,        intent(IN)    :: b(*)
  CMPLX,        intent(inout) :: c(*)

  interface
    subroutine ZBLAS(gemm) (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(1), intent(in)    :: transa, transb
      integer,      intent(in)    :: m, n, k, lda, ldb, ldc
      CMPLX,        intent(in)    :: alpha, beta
      CMPLX,        intent(in)    :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
      CMPLX,        intent(in)    :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
      CMPLX,        intent(inout) :: c ! c(ldc,n)
    end subroutine ZBLAS(gemm)
  end interface

  call ZBLAS(gemm) ('N', 'N', m, n, k, alpha, a(1), m, b(1), k, beta, c(1), m)

end subroutine zlalg_gemm
 
! matrix-vector multiplication plus vector
subroutine dlalg_gemv(m, n, alpha, a, x, beta, y)
  integer,      intent(in)    :: m, n
  FLOAT,        intent(in)    :: alpha, beta
  FLOAT,        intent(IN)    :: a(*)
  FLOAT,        intent(IN)    :: x(*)
  FLOAT,        intent(inout) :: y(*)

  interface
    subroutine DBLAS(gemv) (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      FLOAT,        intent(in)    :: alpha, beta
      FLOAT,        intent(in)    :: a ! a(lda,n)
      FLOAT,        intent(in)    :: x ! x(:)
      FLOAT,        intent(inout) :: y ! y(:)
    end subroutine DBLAS(gemv)
  end interface

  call DBLAS(gemv) ('N', m, n, alpha, a(1), m, x(1), 1, beta, y(1), 1)

end subroutine dlalg_gemv

subroutine zlalg_gemv(m, n, alpha, a, x, beta, y)
  integer,      intent(in)    :: m, n
  CMPLX,        intent(in)    :: alpha, beta
  CMPLX,        intent(IN)    :: a(*)
  CMPLX,        intent(IN)    :: x(*)
  CMPLX,        intent(inout) :: y(*)

  interface
    subroutine ZBLAS(gemv) (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: m, n, lda, incx, incy
      CMPLX,        intent(in)    :: alpha, beta
      CMPLX,        intent(in)    :: a ! a(lda,n)
      CMPLX,        intent(in)    :: x ! x(:)
      CMPLX,        intent(inout) :: y ! y(:)
    end subroutine ZBLAS(gemv)
  end interface

  call ZBLAS(gemv) ('N', m, n, alpha, a(1), m, x(1), 1, beta, y(1), 1)

end subroutine zlalg_gemv

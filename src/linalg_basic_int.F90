#if   TYPE == 1
#  define TYPE1 real(4)
#  define TYPE2 real(4)
#elif TYPE == 2
#  define TYPE1 real(8)
#  define TYPE2 real(8)
#elif TYPE == 3
#  define TYPE1 complex(4)
#  define TYPE2 real(4)
#elif TYPE == 4
 #  define TYPE1 complex(8)
 #  define TYPE2 real(8)
#endif

#define FNAME(x) xFNAME(x, TYPE)
#define xFNAME(x,y) yFNAME(x,y)
#define yFNAME(x,y) x ## _ ## y

! ------------------------------------------------------------------
! BLAS level I
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! swap two vectors
! ------------------------------------------------------------------

subroutine FNAME(swap_1)(n1, dx, dy)
  integer, intent(in)    :: n1
  TYPE1,   intent(inout) :: dx(:), dy(:)

  TYPE1, allocatable :: dz(:)
  allocate(dz(n1))
  dz = dx; dx = dy; dy = dz
  deallocate(dz) 
  
end subroutine FNAME(swap_1)

subroutine FNAME(swap_2)(n1, n2, dx, dy)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(inout) :: dx(:,:), dy(:,:)

  TYPE1, allocatable :: dz(:,:)
  allocate(dz(n1,n2))
  dz = dx; dx = dy; dy = dz
  deallocate(dz) 
  
end subroutine FNAME(swap_2)

subroutine FNAME(swap_3)(n1, n2, n3, dx, dy)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(inout) :: dx(:,:,:), dy(:,:,:)

  TYPE1, allocatable :: dz(:,:,:)
  allocate(dz(n1,n2,n3))
  dz = dx; dx = dy; dy = dz
  deallocate(dz)
  
end subroutine FNAME(swap_3)

subroutine FNAME(swap_4)(n1, n2, n3, n4, dx, dy)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(inout) :: dx(:,:,:,:), dy(:,:,:,:)

  TYPE1, allocatable :: dz(:,:,:,:)
  allocate(dz(n1,n2,n3,n4))
  dz = dx; dx = dy; dy = dz
  deallocate(dz)
  
end subroutine FNAME(swap_4)

! ------------------------------------------------------------------
! scales a vector by a constant
! ------------------------------------------------------------------

subroutine FNAME(scal_1)(n1, da, dx)
  integer, intent(in)    :: n1
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:)

  dx(1:n1) = da*dx(1:n1)

end subroutine FNAME(scal_1)

subroutine FNAME(scal_2)(n1, n2, da, dx)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:,:)

  dx(1:n1, 1:n2) = da*dx(1:n1,1:n2)

end subroutine FNAME(scal_2)

subroutine FNAME(scal_3)(n1, n2, n3, da, dx)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:,:,:)

  dx(1:n1,1:n2,1:n3) = da*dx(1:n1,1:n2,1:n3)

end subroutine FNAME(scal_3)

subroutine FNAME(scal_4)(n1, n2, n3, n4, da, dx)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:,:,:,:)

  dx(1:n1,1:n2,1:n3,1:n4) = da*dx(1:n1,1:n2,1:n3,1:n4)

end subroutine FNAME(scal_4)

! ------------------------------------------------------------------
! constant times a vector plus a vector
! ------------------------------------------------------------------

subroutine FNAME(axpy_1)(n1, da, dx, dy)
  integer, intent(in)    :: n1
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(IN)    :: dx(:)
  TYPE1,   intent(inout) :: dy(:)

  dy(1:n1) = dy(1:n1) + da*dx(1:n1)

end subroutine FNAME(axpy_1)

subroutine FNAME(axpy_2)(n1, n2, da, dx, dy)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(IN)    :: dx(:,:)
  TYPE1,   intent(inout) :: dy(:,:)

  dy(1:n1,1:n2) = dy(1:n1,1:n2) + da*dx(1:n1,1:n2)

end subroutine FNAME(axpy_2)

subroutine FNAME(axpy_3)(n1, n2, n3, da, dx, dy)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(IN)    :: dx(:,:,:)
  TYPE1,   intent(inout) :: dy(:,:,:)

  dy(1:n1,1:n2,1:n3) = dy(1:n1,1:n2,1:n3) + da*dx(1:n1,1:n2,1:n3)

end subroutine FNAME(axpy_3)

subroutine FNAME(axpy_4)(n1, n2, n3, n4, da, dx, dy)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(IN)    :: dx(:,:,:,:)
  TYPE1,   intent(inout) :: dy(:,:,:,:)

  dy(1:n1,1:n2,1:n3,1:n4) = dy(1:n1,1:n2,1:n3,1:n4) + da*dx(1:n1,1:n2,1:n3,1:n4)

end subroutine FNAME(axpy_4)

! ------------------------------------------------------------------
! Copies a vector x, to a vector y
! ------------------------------------------------------------------

subroutine FNAME(copy_1)(n1, dx, dy)
  integer, intent(in)  :: n1
  TYPE1,   intent(IN)  :: dx(:)
  TYPE1,   intent(out) :: dy(:)

  dy(1:n1) = dx(1:n1)
end subroutine FNAME(copy_1)

subroutine FNAME(copy_2)(n1, n2, dx, dy)
  integer, intent(in)  :: n1, n2
  TYPE1,   intent(IN)  :: dx(:,:)
  TYPE1,   intent(out) :: dy(:,:)

  dy(1:n1,1:n2) = dx(1:n1,1:n2)
end subroutine FNAME(copy_2)

subroutine FNAME(copy_3)(n1, n2, n3, dx, dy)
  integer, intent(in)  :: n1, n2, n3
  TYPE1,   intent(IN)  :: dx(:,:,:)
  TYPE1,   intent(out) :: dy(:,:,:)

  dy(1:n1,1:n2,1:n3) = dx(1:n1,1:n2,1:n3)
end subroutine FNAME(copy_3)

subroutine FNAME(copy_4)(n1, n2, n3, n4, dx, dy)
  integer, intent(in)  :: n1, n2, n3, n4
  TYPE1,   intent(IN)  :: dx(:,:,:,:)
  TYPE1,   intent(out) :: dy(:,:,:,:)

  dy(1:n1,1:n2,1:n3,1:n4) = dx(1:n1,1:n2,1:n3,1:n4)
end subroutine FNAME(copy_4)

! ------------------------------------------------------------------
! Forms the dot product of two vectors
! ------------------------------------------------------------------

TYPE1 function FNAME(dot) (n, dx, dy) result(dot)
  integer, intent(in) :: n
  TYPE1,   intent(IN) :: dx(:), dy(:)

  dot = dot_product(dx(1:n), dy(1:n))

end function FNAME(dot)

! ------------------------------------------------------------------
! Returns the euclidean norm of a vector
! ------------------------------------------------------------------

TYPE2 function FNAME(nrm2)(n, dx) result(nrm2)
  integer, intent(in) :: n
  TYPE1,   intent(IN) :: dx(:)

  nrm2 = sqrt(real(dot_product(dx(1:n), dx(1:n)), PRECISION))

end function FNAME(nrm2)

! ------------------------------------------------------------------
! BLAS level II
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! matrix-matrix multiplication plus matrix
! ------------------------------------------------------------------

subroutine FNAME(gemm)(m, n, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m, n, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(IN)    :: a(:,:)  ! a(m, k)
  TYPE1,   intent(IN)    :: b(:,:)  ! b(k, n)
  TYPE1,   intent(inout) :: c(:,:)  ! c(m, n)

  c(1:m,1:n) = alpha*matmul(a(1:m,1:k), b(1:m,1:k)) + beta*c(1:m,1:n)

end subroutine FNAME(gemm)

! ------------------------------------------------------------------
! matrix-vector multiplication plus vector
! ------------------------------------------------------------------
 
subroutine FNAME(gemv)(m, n, alpha, a, x, beta, y)
  integer, intent(in)    :: m, n
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(IN)    :: a(:,:)
  TYPE1,   intent(IN)    :: x(:)
  TYPE1,   intent(inout) :: y(:)

  y(1:n) = alpha*matmul(a(1:m, 1:n), x(1:n)) + beta*y(1:n)

end subroutine FNAME(gemv)

! ------------------------------------------------------------------
! Clean up preprocessor directives
! ------------------------------------------------------------------

#undef ARG_LIST
#undef ARG_CALL
#undef TYPE1
#undef TYPE2
#undef FNAME
#undef xFNAME
#undef yFNAME

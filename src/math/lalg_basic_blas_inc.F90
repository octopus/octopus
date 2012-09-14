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

! ------------------------------------------------------------------
! Preprocessor directives
! ------------------------------------------------------------------

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

  if (n1 < 1) return

  ASSERT(ubound(dx, dim = 1) >= n1)
  ASSERT(ubound(dy, dim = 1) >= n1)

  call blas_swap(n1, dx(1), 1, dy(1), 1)

end subroutine FNAME(swap_1)

subroutine FNAME(swap_2)(n1, n2, dx, dy)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(inout) :: dx(:,:), dy(:,:)

  if (n1*n2 < 1) return

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) >= n2)
  ASSERT(ubound(dy, dim = 2) >= n2)

  call blas_swap(n1*n2, dx(1,1), 1, dy(1,1), 1)

end subroutine FNAME(swap_2)

subroutine FNAME(swap_3)(n1, n2, n3, dx, dy)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(inout) :: dx(:,:,:), dy(:,:,:)

  if (n1*n2*n3 < 1) return

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dy, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) >= n3)
  ASSERT(ubound(dy, dim = 3) >= n3)

  call blas_swap(n1*n2*n3, dx(1,1,1), 1, dy(1,1,1), 1)

end subroutine FNAME(swap_3)

subroutine FNAME(swap_4)(n1, n2, n3, n4, dx, dy)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(inout) :: dx(:,:,:,:), dy(:,:,:,:)

  if (n1*n2*n3*n4 < 1) return

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dy, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) == n3)
  ASSERT(ubound(dy, dim = 3) == n3)
  ASSERT(ubound(dx, dim = 4) >= n4)
  ASSERT(ubound(dy, dim = 4) >= n4)

  call blas_swap(n1*n2*n3*n4, dx(1,1,1,1), 1, dy(1,1,1,1), 1)

end subroutine FNAME(swap_4)

! ------------------------------------------------------------------
! scales a vector by a constant
! ------------------------------------------------------------------

subroutine FNAME(scal_1)(n1, da, dx)
  integer, intent(in)    :: n1
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:)

  ASSERT(ubound(dx, dim = 1) >= n1)

  if (n1 < 1) return

  call blas_scal(n1, da, dx(1), 1)

end subroutine FNAME(scal_1)

subroutine FNAME(scal_2)(n1, n2, da, dx)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:, :)

  if (n1*n2 < 1) return

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) >= n2)

  call blas_scal(n1*n2, da, dx(1,1), 1)

end subroutine FNAME(scal_2)

subroutine FNAME(scal_3)(n1, n2, n3, da, dx)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:, :, :)

  if (n1*n2*n3 < 1) return

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) >= n3)

  call blas_scal(n1*n2*n3, da, dx(1,1,1), 1)

end subroutine FNAME(scal_3)

subroutine FNAME(scal_4)(n1, n2, n3, n4, da, dx)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:, :, :, :)

  if (n1*n2*n3*n4 < 1) return

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) == n3)
  ASSERT(ubound(dx, dim = 4) >= n4)

  call blas_scal(n1*n2*n3*n4, da, dx(1,1,1,1), 1)

end subroutine FNAME(scal_4)

#if TYPE == 3 || TYPE == 4
subroutine FNAME(scal_5)(n1, da, dx)
  integer, intent(in)    :: n1
  TYPE2,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:)

  ASSERT(ubound(dx, dim = 1) >= n1)

  if (n1 < 1) return

  call blas_scal(n1, da, dx(1))

end subroutine FNAME(scal_5)
#endif

! ------------------------------------------------------------------
! constant times a vector plus a vector
! ------------------------------------------------------------------

subroutine FNAME(axpy_1)(n1, da, dx, dy)
  integer, intent(in)    :: n1
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:)
  TYPE1,   intent(inout) :: dy(:)

  ASSERT(ubound(dx, dim = 1) >= n1)
  ASSERT(ubound(dy, dim = 1) >= n1)

  if (n1 < 1) return
  
  call profiling_in(axpy_profile, "BLAS_AXPY")

  call blas_axpy(n1, da, dx(1), 1, dy(1), 1)
  
#if TYPE == 1 || TYPE == 2
  call profiling_count_operations(n1*2)
#else
  call profiling_count_operations(n1*8)
#endif

  call profiling_out(axpy_profile)

end subroutine FNAME(axpy_1)

subroutine FNAME(axpy_2)(n1, n2, da, dx, dy)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:, :)
  TYPE1,   intent(inout) :: dy(:, :)

  if (n1*n2 < 1) return

  call profiling_in(axpy_profile, "BLAS_AXPY")

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) >= n2)
  ASSERT(ubound(dy, dim = 2) >= n2)

  call blas_axpy(n1*n2, da, dx(1,1), 1, dy(1,1), 1)

#if TYPE == 1 || TYPE == 2
  call profiling_count_operations(n1*n2*2)
#else
  call profiling_count_operations(n1*n2*8)
#endif

  call profiling_out(axpy_profile)
end subroutine FNAME(axpy_2)

subroutine FNAME(axpy_3)(n1, n2, n3, da, dx, dy)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:, :, :)
  TYPE1,   intent(inout) :: dy(:, :, :)

  if (n1*n2*n3 < 1) return

  call profiling_in(axpy_profile, "BLAS_AXPY")

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dy, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) >= n3)
  ASSERT(ubound(dy, dim = 3) >= n3)

  call blas_axpy(n1*n2*n3, da, dx(1,1,1), 1, dy(1,1,1), 1)

#if TYPE == 1 || TYPE == 2
  call profiling_count_operations(n1*n2*n3*2)
#else
  call profiling_count_operations(n1*n2*n3*8)
#endif

  call profiling_out(axpy_profile)
end subroutine FNAME(axpy_3)

subroutine FNAME(axpy_4)(n1, n2, n3, n4, da, dx, dy)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:, :, :, :)
  TYPE1,   intent(inout) :: dy(:, :, :, :)

  if (n1*n2*n3*n4 < 1) return

  call profiling_in(axpy_profile, "BLAS_AXPY")

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dy, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) == n3)
  ASSERT(ubound(dy, dim = 3) == n3)
  ASSERT(ubound(dx, dim = 4) >= n4)
  ASSERT(ubound(dy, dim = 4) >= n4)

  call blas_axpy(n1*n2*n3*n4, da, dx(1,1,1,1), 1, dy(1,1,1,1), 1)

#if TYPE == 1 || TYPE == 2
  call profiling_count_operations(n1*n2*n3*n4*2)
#else
  call profiling_count_operations(n1*n2*n3*n4*8)
#endif

  call profiling_out(axpy_profile)

end subroutine FNAME(axpy_4)

#if TYPE == 3 || TYPE == 4
subroutine FNAME(axpy_5)(n1, da, dx, dy)
  integer, intent(in)    :: n1
  TYPE2,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:)
  TYPE1,   intent(inout) :: dy(:)

  ASSERT(ubound(dx, dim = 1) >= n1)
  ASSERT(ubound(dy, dim = 1) >= n1)

  if (n1 < 1) return
  
  call profiling_in(axpy_profile, "BLAS_AXPY")

  call blas_axpy(n1, da, dx(1), dy(1))

  call profiling_count_operations(n1*4)

  call profiling_out(axpy_profile)

end subroutine FNAME(axpy_5)

subroutine FNAME(axpy_6)(n1, n2, da, dx, dy)
  integer, intent(in)    :: n1
  integer, intent(in)    :: n2
  TYPE2,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:, :)
  TYPE1,   intent(inout) :: dy(:, :)

  if (n1 < 1) return
  
  call profiling_in(axpy_profile, "BLAS_AXPY")

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) >= n2)
  ASSERT(ubound(dy, dim = 2) >= n2)

  call blas_axpy(n1*n2, da, dx(1, 1), dy(1, 1))

  call profiling_count_operations(n1*n2*4)

  call profiling_out(axpy_profile)

end subroutine FNAME(axpy_6)
#endif

! ------------------------------------------------------------------
! Copies a vector x, to a vector y
! ------------------------------------------------------------------

subroutine FNAME(copy_1)(n1, dx, dy)
  integer, intent(in)  :: n1
  TYPE1,   intent(in)  :: dx(:)
  TYPE1,   intent(out) :: dy(:)

#ifdef HAVE_OPENMP
  integer :: ini, nn_loc
#endif

  ASSERT(ubound(dx, dim = 1) >= n1)
  ASSERT(ubound(dy, dim = 1) >= n1)

  if (n1 < 1) return

  call profiling_in(copy_profile, "BLAS_COPY")

  call blas_copy(n1, dx(1), 1, dy(1), 1)

  call profiling_count_transfers(n1, dx(1))

  call profiling_out(copy_profile)

end subroutine FNAME(copy_1)

subroutine FNAME(copy_2)(n1, n2, dx, dy)
  integer, intent(in)  :: n1, n2
  TYPE1,   intent(in)  :: dx(:,:)
  TYPE1,   intent(out) :: dy(:,:)

  if (n1*n2 < 1) return

  call profiling_in(copy_profile, "BLAS_COPY")

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) >= n2)
  ASSERT(ubound(dy, dim = 2) >= n2)

  call blas_copy(n1*n2, dx(1,1), 1, dy(1,1), 1)

  call profiling_out(copy_profile)

end subroutine FNAME(copy_2)

subroutine FNAME(copy_3)(n1, n2, n3, dx, dy)
  integer, intent(in)  :: n1, n2, n3
  TYPE1,   intent(in)  :: dx(:,:,:)
  TYPE1,   intent(out) :: dy(:,:,:)

  if(n1*n2*n3 < 1) return

  call profiling_in(copy_profile, "BLAS_COPY")

  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dy, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) >= n3)
  ASSERT(ubound(dy, dim = 3) >= n3)

  call blas_copy (n1*n2*n3, dx(1,1,1), 1, dy(1,1,1), 1)

  call profiling_out(copy_profile)

end subroutine FNAME(copy_3)

subroutine FNAME(copy_4)(n1, n2, n3, n4, dx, dy)
  integer, intent(in)  :: n1, n2, n3, n4
  TYPE1,   intent(in)  :: dx(:,:,:,:)
  TYPE1,   intent(out) :: dy(:,:,:,:)

  if (n1*n2*n3*n4 < 1) return

  call profiling_in(copy_profile, "BLAS_COPY")
 
  ASSERT(ubound(dx, dim = 1) == n1)
  ASSERT(ubound(dy, dim = 1) == n1)
  ASSERT(ubound(dx, dim = 2) == n2)
  ASSERT(ubound(dy, dim = 2) == n2)
  ASSERT(ubound(dx, dim = 3) == n3)
  ASSERT(ubound(dy, dim = 3) == n3)
  ASSERT(ubound(dx, dim = 4) >= n4)
  ASSERT(ubound(dy, dim = 4) >= n4)

  call blas_copy (n1*n2*n3*n4, dx(1,1,1,1), 1, dy(1,1,1,1), 1)

  call profiling_out(copy_profile)

end subroutine FNAME(copy_4)

! ------------------------------------------------------------------
! Returns the euclidean norm of a vector
! ------------------------------------------------------------------

TYPE2 function FNAME(nrm2)(n, dx) result(nrm2)
  integer, intent(in) :: n
  TYPE1,   intent(in) :: dx(:)

  ASSERT(ubound(dx, dim = 1) >= n)

  nrm2 = CNST(0.0)
  if (n < 1) return

  nrm2 = blas_nrm2(n, dx(1), 1)

end function FNAME(nrm2)

! ------------------------------------------------------------------
! BLAS level II
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Matrix-vector multiplication plus vector.
! ------------------------------------------------------------------

subroutine FNAME(symv_1)(n, alpha, a, x, beta, y)
  integer, intent(in)    :: n
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:, :)
  TYPE1,   intent(in)    :: x(:)
  TYPE1,   intent(inout) :: y(:)

  call profiling_in(symv_profile, 'BLAS_SYMV')
  call blas_symv('U', n, alpha, a(1, 1), n, x(1), 1, beta, y(1), 1)
  call profiling_out(symv_profile)

end subroutine FNAME(symv_1)

subroutine FNAME(symv_2)(n1, n2, alpha, a, x, beta, y)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:, :, :)
  TYPE1,   intent(in)    :: x(:)
  TYPE1,   intent(inout) :: y(:, :)

  call profiling_in(symv_profile, 'BLAS_SYMV')
  call blas_symv('U', n1*n2, alpha, a(1, 1, 1), n1*2, x(1), 1, beta, y(1, 1), 1)
  call profiling_out(symv_profile)

end subroutine FNAME(symv_2)

subroutine FNAME(gemv_1)(m, n, alpha, a, x, beta, y)
  integer, intent(in)    :: m, n
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:,:)
  TYPE1,   intent(in)    :: x(:)
  TYPE1,   intent(inout) :: y(:)

  call profiling_in(gemv_profile, "BLAS_GEMV")
  call blas_gemv('N', m, n, alpha, a(1,1), m, x(1), 1, beta, y(1), 1)
  call profiling_out(gemv_profile)

end subroutine FNAME(gemv_1)

subroutine FNAME(gemv_2)(m1, m2, n, alpha, a, x, beta, y)
  integer, intent(in)    :: m1, m2, n
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:,:,:)
  TYPE1,   intent(in)    :: x(:)
  TYPE1,   intent(inout) :: y(:,:)

  call profiling_in(gemv_profile, "BLAS_GEMV")
  call blas_gemv('N', m1*m2, n, alpha, a(1,1,1), m1*m2, x(1), 1, beta, y(1,1), 1)
  call profiling_out(gemv_profile)

end subroutine FNAME(gemv_2)


! ------------------------------------------------------------------
!> BLAS level III
! ------------------------------------------------------------------

! ------------------------------------------------------------------
!> Matrix-matrix multiplication plus matrix.
! ------------------------------------------------------------------

subroutine FNAME(gemm_1)(m, n, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m, n, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:,:)  !< a(m, k)
  TYPE1,   intent(in)    :: b(:,:)  !< b(k, n)
  TYPE1,   intent(inout) :: c(:,:)  !< c(m, n)

  PUSH_SUB(FNAME(gemm_1))

  call blas_gemm('N', 'N', m, n, k, alpha, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), beta, c(1, 1), lead_dim(c))

  POP_SUB(FNAME(gemm_1))
end subroutine FNAME(gemm_1)

subroutine FNAME(gemm_2)(m, n, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m, n, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:, :, :)  !< a(m, k)
  TYPE1,   intent(in)    :: b(:, :)     !< b(k, n)
  TYPE1,   intent(inout) :: c(:, :, :)  !< c(m, n)

  PUSH_SUB(FNAME(gemm_2))

  call blas_gemm('N', 'N', m, n, k, alpha, a(1, 1, 1), lead_dim(a), &
    b(1, 1), lead_dim(b), beta, c(1, 1, 1), lead_dim(c))

  POP_SUB(FNAME(gemm_2))
end subroutine FNAME(gemm_2)

!> The same as above but with (Hermitian) transpose of a.
subroutine FNAME(gemmt_1)(m, n, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m, n, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:,:)  !< a(k, m)
  TYPE1,   intent(in)    :: b(:,:)  !< b(k, n)
  TYPE1,   intent(inout) :: c(:,:)  !< c(m, n)

  PUSH_SUB(FNAME(gemmt_1))

  call blas_gemm('C', 'N', m, n, k, alpha, a(1, 1), lead_dim(a), b(1, 1), lead_dim(b), beta, c(1, 1), lead_dim(c))

  POP_SUB(FNAME(gemmt_1))
end subroutine FNAME(gemmt_1)

subroutine FNAME(gemmt_2)(m, n, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m, n, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:, :, :)  !< a((k), m)
  TYPE1,   intent(in)    :: b(:, :, :)  !< b((k), n)
  TYPE1,   intent(inout) :: c(:, :)     !< c(m, n)

  PUSH_SUB(FNAME(gemmt_2))

  call blas_gemm('C', 'N', m, n, k, alpha, a(1, 1, 1), lead_dim(a), &
    b(1, 1, 1), lead_dim(b), beta, c(1, 1), lead_dim(c))

  PUSH_SUB(FNAME(gemmt_2))
end subroutine FNAME(gemmt_2)

!> The following matrix multiplications all expect upper triangular matrices for a.
!! For real matrices, a = a^T, for complex matrices a = a^H.
subroutine FNAME(symm_1)(m, n, side, alpha, a, b, beta, c)
  integer,      intent(in)    :: m, n
  character(1), intent(in)    :: side
  TYPE1,        intent(in)    :: alpha, beta, a(:, :), b(:, :)
  TYPE1,        intent(inout) :: c(:, :)

  integer :: lda

  PUSH_SUB(FNAME(symm_1))

  select case(side)
  case('l', 'L')
    lda = max(1, m)
  case('r', 'R')
    lda = max(1, n)
  end select
  
  call blas_symm(side, 'U', m, n, alpha, a(1, 1), lda, b(1, 1), m, beta, c(1, 1), m)

  POP_SUB(FNAME(symm_1))
end subroutine FNAME(symm_1)


subroutine FNAME(symm_2)(m, n, side, alpha, a, b, beta, c)
  integer,      intent(in)    :: m, n
  character(1), intent(in)    :: side
  TYPE1,        intent(in)    :: alpha, beta, a(:, :, :), b(:, :)
  TYPE1,        intent(inout) :: c(:, :, :)

  integer :: lda

  PUSH_SUB(FNAME(symm_2))

  select case(side)
  case('l', 'L')
    lda = max(1, m)
  case('r', 'R')
    lda = max(1, n)
  end select
  
  call blas_symm(side, 'U', m, n, alpha, a(1, 1, 1), lda, b(1, 1), m, beta, c(1, 1, 1), m)

  POP_SUB(FNAME(symm_2))
end subroutine FNAME(symm_2)

! ------------------------------------------------------------------
! Matrix-matrix multiplication.
! ------------------------------------------------------------------

subroutine FNAME(trmm_1)(m, n, uplo, transa, side, alpha, a, b)
  integer,      intent(in)    :: m, n
  character(1), intent(in)    :: side, transa, uplo
  TYPE1,        intent(in)    :: alpha
  TYPE1,        intent(in)    :: a(:, :) ! a(m, m), upper triangular matrix.
  TYPE1,        intent(inout) :: b(:, :) ! b(m, n).

  integer :: lda

  PUSH_SUB(FNAME(trmm_1))

  select case(side)
    case('L', 'l')
      lda = max(1, m)
    case('R', 'r')
      lda = max(1, n)
  end select
      
  call blas_trmm(side, uplo, transa, 'N', m, n, alpha, a(1, 1), lda, b(1, 1), m)

  POP_SUB(FNAME(trmm_1))
end subroutine FNAME(trmm_1)


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


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

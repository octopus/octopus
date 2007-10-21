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

  call blas_swap(n1, dx(1), 1, dy(1), 1)

end subroutine FNAME(swap_1)

subroutine FNAME(swap_2)(n1, n2, dx, dy)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(inout) :: dx(:,:), dy(:,:)

  call blas_swap(n1*n2, dx(1,1), 1, dy(1,1), 1)

end subroutine FNAME(swap_2)

subroutine FNAME(swap_3)(n1, n2, n3, dx, dy)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(inout) :: dx(:,:,:), dy(:,:,:)

  call blas_swap(n1*n2*n3, dx(1,1,1), 1, dy(1,1,1), 1)

end subroutine FNAME(swap_3)

subroutine FNAME(swap_4)(n1, n2, n3, n4, dx, dy)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(inout) :: dx(:,:,:,:), dy(:,:,:,:)

  call blas_swap(n1*n2*n3*n4, dx(1,1,1,1), 1, dy(1,1,1,1), 1)

end subroutine FNAME(swap_4)

! ------------------------------------------------------------------
! scales a vector by a constant
! ------------------------------------------------------------------

subroutine FNAME(scal_1)(n1, da, dx)
  integer, intent(in)    :: n1
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(1:n1)

#ifdef USE_OMP
  integer :: nn_loc, ini

!$omp parallel private(ini, nn_loc)
  call divide_range(n1, omp_get_thread_num(), omp_get_num_threads(), ini, nn_loc)
  call blas_scal(nn_loc, da, dx(ini), 1)
!$omp end parallel
#else
  call blas_scal(n1, da, dx(1), 1)
#endif

end subroutine FNAME(scal_1)

subroutine FNAME(scal_2)(n1, n2, da, dx)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:,:)

  call blas_scal(n1*n2, da, dx(1,1), 1)

end subroutine FNAME(scal_2)

subroutine FNAME(scal_3)(n1, n2, n3, da, dx)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:,:,:)

  call blas_scal(n1*n2*n3, da, dx(1,1,1), 1)

end subroutine FNAME(scal_3)

subroutine FNAME(scal_4)(n1, n2, n3, n4, da, dx)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:,:,:,:)

  call blas_scal(n1*n2*n3*n4, da, dx(1,1,1,1), 1)

end subroutine FNAME(scal_4)

#if TYPE == 3 || TYPE == 4
subroutine FNAME(scal_5)(n1, da, dx)
  integer, intent(in)    :: n1
  TYPE2,   intent(in)    :: da
  TYPE1,   intent(inout) :: dx(:)

  call blas_scal(n1, da, dx(1))

end subroutine FNAME(scal_5)
#endif

! ------------------------------------------------------------------
! constant times a vector plus a vector
! ------------------------------------------------------------------

subroutine FNAME(axpy_1)(n1, da, dx, dy)
  integer, intent(in)    :: n1
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(1:n1)
  TYPE1,   intent(inout) :: dy(1:n1)

#ifdef USE_OMP
  integer :: nn_loc, ini
#endif
  
  call profiling_in(axpy_profile, "BLAS_AXPY")

#ifdef USE_OMP
!$omp parallel private(ini, nn_loc)
  call divide_range(n1, omp_get_thread_num(), omp_get_num_threads(), ini, nn_loc)
  call blas_axpy(nn_loc, da, dx(ini), 1, dy(ini), 1)
!$omp end parallel
#else
  call blas_axpy(n1, da, dx(1), 1, dy(1), 1)
#endif

  call profiling_out(axpy_profile)

end subroutine FNAME(axpy_1)

subroutine FNAME(axpy_2)(n1, n2, da, dx, dy)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(1:n1,1:n2)
  TYPE1,   intent(inout) :: dy(1:n1,1:n2)

  call profiling_in(axpy_profile, "BLAS_AXPY")

  call blas_axpy(n1*n2, da, dx(1,1), 1, dy(1,1), 1)

  call profiling_out(axpy_profile)
end subroutine FNAME(axpy_2)

subroutine FNAME(axpy_3)(n1, n2, n3, da, dx, dy)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(1:n1,1:n2,1:n3)
  TYPE1,   intent(inout) :: dy(1:n1,1:n2,1:n3)

  call profiling_in(axpy_profile, "BLAS_AXPY")

  call blas_axpy(n1*n2*n3, da, dx(1,1,1), 1, dy(1,1,1), 1)

  call profiling_out(axpy_profile)
end subroutine FNAME(axpy_3)

subroutine FNAME(axpy_4)(n1, n2, n3, n4, da, dx, dy)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(1:n1,1:n2,1:n3,1:n4)
  TYPE1,   intent(inout) :: dy(1:n1,1:n2,1:n3,1:n4)

  call profiling_in(axpy_profile, "BLAS_AXPY")

  call blas_axpy(n1*n2*n3*n4, da, dx(1,1,1,1), 1, dy(1,1,1,1), 1)

  call profiling_out(axpy_profile)

end subroutine FNAME(axpy_4)

#if TYPE == 3 || TYPE == 4
subroutine FNAME(axpy_5)(n1, da, dx, dy)
  integer, intent(in)    :: n1
  TYPE2,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(1:n1)
  TYPE1,   intent(inout) :: dy(1:n1)

#ifdef USE_OMP
  integer :: ini, nn_loc
#endif

    call profiling_in(axpy_profile, "BLAS_AXPY")

#ifdef USE_OMP  
  !$omp parallel private(ini, nn_loc)
  call divide_range(n1, omp_get_thread_num(), omp_get_num_threads(), ini, nn_loc)
  call blas_axpy(nn_loc, da, dx(ini), dy(ini))
  !$omp end parallel
#else
  call blas_axpy(n1, da, dx(1), dy(1))
#endif

  call profiling_out(axpy_profile)

end subroutine FNAME(axpy_5)
#endif

! ------------------------------------------------------------------
! Copies a vector x, to a vector y
! ------------------------------------------------------------------

subroutine FNAME(copy_1)(n1, dx, dy)
  integer, intent(in)  :: n1
  TYPE1,   intent(in)  :: dx(:)
  TYPE1,   intent(out) :: dy(:)

#ifdef USE_OMP
  integer :: ini, nn_loc
#endif

  call profiling_in(copy_profile, "BLAS_COPY")

#ifdef USE_OMP
  !$omp parallel private(ini, nn_loc)
  call divide_range(n1, omp_get_thread_num(), omp_get_num_threads(), ini, nn_loc)
  call blas_copy(nn_loc, dx(ini), 1, dy(ini), 1)
  !$omp end parallel
#else
  call blas_copy(n1, dx(1), 1, dy(1), 1)
#endif

  call profiling_out(copy_profile)

end subroutine FNAME(copy_1)

subroutine FNAME(copy_2)(n1, n2, dx, dy)
  integer, intent(in)  :: n1, n2
  TYPE1,   intent(in)  :: dx(:,:)
  TYPE1,   intent(out) :: dy(:,:)

  call profiling_in(copy_profile, "BLAS_COPY")

  call blas_copy(n1*n2, dx(1,1), 1, dy(1,1), 1)

  call profiling_out(copy_profile)

end subroutine FNAME(copy_2)

subroutine FNAME(copy_3)(n1, n2, n3, dx, dy)
  integer, intent(in)  :: n1, n2, n3
  TYPE1,   intent(in)  :: dx(:,:,:)
  TYPE1,   intent(out) :: dy(:,:,:)

  call profiling_in(copy_profile, "BLAS_COPY")

  call blas_copy (n1*n2*n3, dx(1,1,1), 1, dy(1,1,1), 1)

  call profiling_out(copy_profile)

end subroutine FNAME(copy_3)

subroutine FNAME(copy_4)(n1, n2, n3, n4, dx, dy)
  integer, intent(in)  :: n1, n2, n3, n4
  TYPE1,   intent(in)  :: dx(:,:,:,:)
  TYPE1,   intent(out) :: dy(:,:,:,:)

  call profiling_in(copy_profile, "BLAS_COPY")

  call blas_copy (n1*n2*n3*n4, dx(1,1,1,1), 1, dy(1,1,1,1), 1)

  call profiling_out(copy_profile)

end subroutine FNAME(copy_4)

! ------------------------------------------------------------------
! Forms the dot product of two vectors
! ------------------------------------------------------------------

TYPE1 function FNAME(dot) (n, dx, dy) result(dot)
  integer, intent(in) :: n
  TYPE1,   intent(in) :: dx(:), dy(:)

#ifdef USE_OMP
  TYPE1   :: dot_loc
  integer :: nn_loc, ini
  
  dot = CNST(0.0)

  !$omp parallel private(ini, nn_loc, dot_loc)
  call divide_range(n, omp_get_thread_num(), omp_get_num_threads(), ini, nn_loc)
  dot_loc = blas_dot(nn_loc, dx(ini), 1, dy(ini), 1)

  !$omp atomic
  dot = dot + dot_loc

  !$omp end parallel

#else
  dot = blas_dot(n, dx(1), 1, dy(1), 1)
#endif

end function FNAME(dot)

! ------------------------------------------------------------------
! Returns the euclidean norm of a vector
! ------------------------------------------------------------------

TYPE2 function FNAME(nrm2)(n, dx) result(nrm2)
  integer, intent(in) :: n
  TYPE1,   intent(in) :: dx(:)

#ifdef USE_OMP
  TYPE2   :: nrm2_loc
  integer :: nn_loc, ini
  
  nrm2 = CNST(0.0)

  !$omp parallel private(ini, nn_loc, nrm2_loc)
  call divide_range(n, omp_get_thread_num(), omp_get_num_threads(), ini, nn_loc)
  nrm2_loc = blas_nrm2(nn_loc, dx(ini), 1)

  !$omp critical
  nrm2 = hypot(nrm2, nrm2_loc)
  !$omp end critical

  !$omp end parallel
  
#else
  nrm2 = blas_nrm2(n, dx(1), 1)
#endif

end function FNAME(nrm2)

! ------------------------------------------------------------------
! BLAS level II
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Matrix-vector multiplication plus vector.
! ------------------------------------------------------------------

subroutine FNAME(gemv_1)(m, n, alpha, a, x, beta, y)
  integer, intent(in)    :: m, n
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:,:)
  TYPE1,   intent(in)    :: x(:)
  TYPE1,   intent(inout) :: y(:)

  call blas_gemv('N', m, n, alpha, a(1,1), m, x(1), 1, beta, y(1), 1)

end subroutine FNAME(gemv_1)

subroutine FNAME(gemv_2)(m1, m2, n, alpha, a, x, beta, y)
  integer, intent(in)    :: m1, m2, n
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:,:,:)
  TYPE1,   intent(in)    :: x(:)
  TYPE1,   intent(inout) :: y(:,:)

  call blas_gemv('N', m1*m2, n, alpha, a(1,1,1), m1*m2, x(1), 1, beta, y(1,1), 1)

end subroutine FNAME(gemv_2)


! ------------------------------------------------------------------
! BLAS level III
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Matrix-matrix multiplication plus matrix.
! ------------------------------------------------------------------

subroutine FNAME(gemm_1)(m, n, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m, n, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:,:)  ! a(m, k)
  TYPE1,   intent(in)    :: b(:,:)  ! b(k, n)
  TYPE1,   intent(inout) :: c(:,:)  ! c(m, n)

  call blas_gemm('N', 'N', m, n, k, alpha, a(1, 1), m, b(1, 1), k, beta, c(1, 1), m)

end subroutine FNAME(gemm_1)

subroutine FNAME(gemm_2)(m, n, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m, n, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:, :, :)  ! a(m, k)
  TYPE1,   intent(in)    :: b(:, :)     ! b(k, n)
  TYPE1,   intent(inout) :: c(:, :, :)  ! c(m, n)

  call blas_gemm('N', 'N', m, n, k, alpha, a(1, 1, 1), m, b(1, 1), k, beta, c(1, 1, 1), m)

end subroutine FNAME(gemm_2)

! The same as above but with (Hermitian) transposed of a.

subroutine FNAME(gemmt_1)(m, n, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m, n, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:,:)  ! a(m, k)
  TYPE1,   intent(in)    :: b(:,:)  ! b(k, n)
  TYPE1,   intent(inout) :: c(:,:)  ! c(m, n)

  call blas_gemm('C', 'N', m, n, k, alpha, a(1, 1), k, b(1, 1), k, beta, c(1, 1), m)
end subroutine FNAME(gemmt_1)

subroutine FNAME(gemmt_2)(m, n, k, alpha, a, b, beta, c)
  integer, intent(in)    :: m, n, k
  TYPE1,   intent(in)    :: alpha, beta
  TYPE1,   intent(in)    :: a(:, :, :)  ! a(m, k)
  TYPE1,   intent(in)    :: b(:, :)     ! b(k, n)
  TYPE1,   intent(inout) :: c(:, :, :)  ! c(m, n)

  call blas_gemm('C', 'N', m, n, k, alpha, a(1, 1, 1), k, b(1, 1), k, beta, c(1, 1, 1), m)
end subroutine FNAME(gemmt_2)

! The following matrix multiplications all expect upper triangular matrices for a.
! For real matrices, a = a^T, for complex matrices a = a^H.

subroutine FNAME(hemm_1)(m, n, side, alpha, a, b, beta, c)
  integer,      intent(in)    :: m, n
  character(1), intent(in)    :: side
  TYPE1,        intent(in)    :: alpha, beta, a(:, :), b(:, :)
  TYPE1,        intent(inout) :: c(:, :)

  integer :: lda

  select case(side)
    case('L', 'l')
      lda = max(1, m)
    case('R', 'r')
      lda = max(1, n)
  end select

  call blas_hemm(side, 'U', m, n, alpha, a(1, 1), lda, b(1, 1), m, beta, c(1, 1), m)
end subroutine FNAME(hemm_1)

subroutine FNAME(hemm_2)(m, n, side, alpha, a, b, beta, c)
  integer,      intent(in)    :: m, n
  character(1), intent(in)    :: side
  TYPE1,        intent(in)    :: alpha, beta, a(:, :, :), b(:, :)
  TYPE1,        intent(inout) :: c(:, :, :)

  integer :: lda

  select case(side)
    case('L', 'l')
      lda = max(1, m)
    case('R', 'r')
      lda = max(1, n)
  end select

  call blas_hemm(side, 'U', m, n, alpha, a(1, 1, 1), lda, b(1, 1), m, beta, c(1, 1, 1), m)
end subroutine FNAME(hemm_2)

! Expects upper triangular matrix for c.
! trans = 'N' => C <- alpha*A*A^H + beta*C
! trans = 'C' => C <- alpha*A^H*A + beta*C
subroutine FNAME(herk_1)(n, k, trans, alpha, a, beta, c)
  integer,      intent(in)    :: n, k
  character(1), intent(in)    :: trans                ! 'N', 'C'
  TYPE1,        intent(in)    :: alpha, a(:, :), beta
  TYPE1,        intent(inout) :: c(:, :)              ! c(n, n)

  integer :: lda, i, j

  select case(trans)
    case('N', 'n')
      lda = max(1, n)
    case('C', 'c', 'T', 't')
      lda = max(1, k)
  end select

  call blas_herk('U', trans, n, k, alpha, a(1, 1), lda, beta, c(1, 1), n)

  ! Fill lower triangular matrix.
  do i = 2, n
    do j = 1, i-1
#if TYPE == 3 || TYPE == 4
      c(i, j) = conjg(c(j, i))
#else
      c(i, j) = c(j, i)
#endif
    end do
  end do
end subroutine FNAME(herk_1)


! ------------------------------------------------------------------
! Matrix-matrix multiplication.
! ------------------------------------------------------------------

subroutine FNAME(trmm_1)(m, n, side, alpha, a, b)
  integer,      intent(in)    :: m, n
  character(1), intent(in)    :: side
  TYPE1,        intent(in)    :: alpha
  TYPE1,        intent(in)    :: a(:, :) ! a(m, m), upper triangular matrix.
  TYPE1,        intent(inout) :: b(:, :) ! b(m, n).

  integer :: lda

  select case(side)
    case('L', 'l')
      lda = max(1, m)
    case('R', 'r')
      lda = max(1, n)
  end select
      
  call blas_trmm(side, 'U', 'N', 'N', m, n, alpha, a(1, 1), lda, b(1, 1), m)
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

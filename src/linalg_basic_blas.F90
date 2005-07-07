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
  TYPE1,   intent(inout) :: dx(:)

  call blas_scal(n1, da, dx(1), 1)

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

! ------------------------------------------------------------------
! constant times a vector plus a vector
! ------------------------------------------------------------------

subroutine FNAME(axpy_1)(n1, da, dx, dy)
  integer, intent(in)    :: n1
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:)
  TYPE1,   intent(inout) :: dy(:)

  call blas_axpy(n1, da, dx(1), 1, dy(1), 1)

end subroutine FNAME(axpy_1)

subroutine FNAME(axpy_2)(n1, n2, da, dx, dy)
  integer, intent(in)    :: n1, n2
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:,:)
  TYPE1,   intent(inout) :: dy(:,:)

  call blas_axpy(n1*n2, da, dx(1,1), 1, dy(1,1), 1)

end subroutine FNAME(axpy_2)

subroutine FNAME(axpy_3)(n1, n2, n3, da, dx, dy)
  integer, intent(in)    :: n1, n2, n3
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:,:,:)
  TYPE1,   intent(inout) :: dy(:,:,:)

  call blas_axpy(n1*n2*n3, da, dx(1,1,1), 1, dy(1,1,1), 1)

end subroutine FNAME(axpy_3)

subroutine FNAME(axpy_4)(n1, n2, n3, n4, da, dx, dy)
  integer, intent(in)    :: n1, n2, n3, n4
  TYPE1,   intent(in)    :: da
  TYPE1,   intent(in)    :: dx(:,:,:,:)
  TYPE1,   intent(inout) :: dy(:,:,:,:)

  call blas_axpy(n1*n2*n3*n4, da, dx(1,1,1,1), 1, dy(1,1,1,1), 1)

end subroutine FNAME(axpy_4)

! ------------------------------------------------------------------
! Copies a vector x, to a vector y
! ------------------------------------------------------------------

subroutine FNAME(copy_1)(n1, dx, dy)
  integer, intent(in)  :: n1
  TYPE1,   intent(in)  :: dx(:)
  TYPE1,   intent(out) :: dy(:)

  call blas_copy(n1, dx(1), 1, dy(1), 1)

end subroutine FNAME(copy_1)

subroutine FNAME(copy_2)(n1, n2, dx, dy)
  integer, intent(in)  :: n1, n2
  TYPE1,   intent(in)  :: dx(:,:)
  TYPE1,   intent(out) :: dy(:,:)

  call blas_copy(n1*n2, dx(1,1), 1, dy(1,1), 1)

end subroutine FNAME(copy_2)

subroutine FNAME(copy_3)(n1, n2, n3, dx, dy)
  integer, intent(in)  :: n1, n2, n3
  TYPE1,   intent(in)  :: dx(:,:,:)
  TYPE1,   intent(out) :: dy(:,:,:)

  call blas_copy (n1*n2*n3, dx(1,1,1), 1, dy(1,1,1), 1)

end subroutine FNAME(copy_3)

subroutine FNAME(copy_4)(n1, n2, n3, n4, dx, dy)
  integer, intent(in)  :: n1, n2, n3, n4
  TYPE1,   intent(in)  :: dx(:,:,:,:)
  TYPE1,   intent(out) :: dy(:,:,:,:)

  call blas_copy (n1*n2*n3*n4, dx(1,1,1,1), 1, dy(1,1,1,1), 1)

end subroutine FNAME(copy_4)

! ------------------------------------------------------------------
! Forms the dot product of two vectors
! ------------------------------------------------------------------

TYPE1 function FNAME(dot) (n, dx, dy) result(dot)
  integer, intent(in) :: n
  TYPE1,   intent(in) :: dx(:), dy(:)

  dot = blas_dot(n, dx(1), 1, dy(1), 1)

end function FNAME(dot)

! ------------------------------------------------------------------
! Returns the euclidean norm of a vector
! ------------------------------------------------------------------

TYPE2 function FNAME(nrm2)(n, dx) result(nrm2)
  integer, intent(in) :: n
  TYPE1,   intent(in) :: dx(:)

  nrm2 = blas_nrm2(n, dx(1), 1)

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
  TYPE1,   intent(in)    :: a(:,:)  ! a(m, k)
  TYPE1,   intent(in)    :: b(:,:)  ! b(k, n)
  TYPE1,   intent(inout) :: c(:,:)  ! c(m, n)

  call blas_gemm('N', 'N', m, n, k, alpha, a(1,1), m, b(1,1), k, beta, c(1,1), m)

end subroutine FNAME(gemm)

! ------------------------------------------------------------------
! matrix-vector multiplication plus vector
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
! Clean up preprocessor directives
! ------------------------------------------------------------------

#undef ARG_LIST
#undef ARG_CALL
#undef TYPE1
#undef TYPE2
#undef FNAME
#undef xFNAME
#undef yFNAME

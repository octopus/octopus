/*
 Copyright (C) 2006 X. Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

 $Id: dablas.c 2146 2006-05-23 17:36:00Z xavier $
*/

#include <config.h>

#include <stdio.h>

#include <fortran_types.h>

/* declare blas functions */
void FC_FUNC(sscal, SSCAL)(const fint * n, const float  * a, const float  * x, const fint * incx);
void FC_FUNC(dscal, DSCAL)(const fint * n, const double * a, const double * x, const fint * incx);
void FC_FUNC(saxpy, SAXPY)(const fint * n, const float  * a, const float  * x, const fint * incx, float  * y, const fint * incy);
void FC_FUNC(daxpy, DAXPY)(const fint * n, const double * a, const double * x, const fint * incx, double * y, const fint * incy);

/* interface to apply blas with a real constant over complex vectors */

void FC_FUNC(sazscal, SAZSCAL)(const fint * n, 
			       const float * restrict a,
			       float * restrict x){
  
  const fint twon = 2*n[0];
  const fint one = 1;
  
  FC_FUNC(sscal, SSCAL)(&twon, a, x, &one);
}

void FC_FUNC(dazscal, DAZSCAL)(const fint * n, 
			       const double * restrict a,
			       double * restrict x){
  
  const fint twon = 2*n[0];
  const fint one = 1;
  
  FC_FUNC(dscal, DSCAL)(&twon, a, x, &one);
}

void FC_FUNC(dazaxpy, DAZAXPY)(const fint * n, 
			       const double * restrict a,
			       const double * restrict x,
			       double * restrict y){

  const fint twon = 2*n[0];
  const fint one = 1;

  FC_FUNC(daxpy, DAXPY)(&twon, a, x, &one, y, &one);

}

void FC_FUNC(sazaxpy, SAZAXPY)(const fint * n, 
			       const float * restrict a,
			       const float * restrict x,
			       float * restrict y){

  const fint twon = 2*n[0];
  const fint one = 1;

  FC_FUNC(saxpy, SAXPY)(&twon, a, x, &one, y, &one);

}


/* declare blas functions */
void FC_FUNC(sgemm, SGEMM)(const char * transa, const char * transb,
			   const fint * m, const fint * n, const fint * k,
			   const float * alpha, const float * a, const fint * lda,
			   const float * b, const fint * ldb, const float * beta, 
			   float * c, const fint * ldc);

void FC_FUNC(dgemm, DGEMM)(const char * transa, const char * transb,
			   const fint * m, const fint * n, const fint * k,
			   const double * alpha, const double * a, const fint * lda,
			   const double * b, const fint * ldb, const double * beta, 
			   double * c, const fint * ldc);

/* interface to apply dgemm passing complex vectors
   the same as dgemm, but allows giving each an appropriate Fortan interface
   in which alpha, beta, a, b, c are actually complex in Fortran 
   Could be inline, but in that case pgcc will not put it in the symbol table. */
void FC_FUNC(zdgemm, ZDGEMM)(const char * transa, const char * transb,
			     const fint * m, const fint * n, const fint * k,
			     const double * alpha, const double * restrict a, const fint * lda,
			     const double * restrict b, const fint * ldb, const double * beta, 
			     double * restrict c, const fint * ldc) {
  FC_FUNC(dgemm, DGEMM)(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

void FC_FUNC(csgemm, CSGEMM)(const char * transa, const char * transb,
			     const fint * m, const fint * n, const fint * k,
			     const float * alpha, const float * restrict a, const fint * lda,
			     const float * restrict b, const fint * ldb, const float * beta, 
			     float * restrict c, const fint * ldc) {

  FC_FUNC(sgemm, SGEMM)(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

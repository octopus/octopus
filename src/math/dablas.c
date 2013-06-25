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
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: dablas.c 2146 2006-05-23 17:36:00Z xavier $
*/

#include <config.h>

#include <stdio.h>


/* declare blas functions */
void FC_FUNC(sscal, SSCAL)(const int * n, const float  * a, const float  * x, const int * incx);
void FC_FUNC(dscal, DSCAL)(const int * n, const double * a, const double * x, const int * incx);
void FC_FUNC(saxpy, SAXPY)(const int * n, const float  * a, const float  * x, const int * incx, float  * y, const int * incy);
void FC_FUNC(daxpy, DAXPY)(const int * n, const double * a, const double * x, const int * incx, double * y, const int * incy);

/* interface to apply blas with a real constant over complex vectors */

void FC_FUNC(sazscal, SAZSCAL)(const int * n, 
			       const float * restrict a,
			       float * restrict x){
  
  const int twon = 2*n[0];
  const int one = 1;
  
  FC_FUNC(sscal, SSCAL)(&twon, a, x, &one);
}

void FC_FUNC(dazscal, DAZSCAL)(const int * n, 
			       const double * restrict a,
			       double * restrict x){
  
  const int twon = 2*n[0];
  const int one = 1;
  
  FC_FUNC(dscal, DSCAL)(&twon, a, x, &one);
}

void FC_FUNC(dazaxpy, DAZAXPY)(const int * n, 
			       const double * restrict a,
			       const double * restrict x,
			       double * restrict y){

  const int twon = 2*n[0];
  const int one = 1;

  FC_FUNC(daxpy, DAXPY)(&twon, a, x, &one, y, &one);

}

void FC_FUNC(sazaxpy, SAZAXPY)(const int * n, 
			       const float * restrict a,
			       const float * restrict x,
			       float * restrict y){

  const int twon = 2*n[0];
  const int one = 1;

  FC_FUNC(saxpy, SAXPY)(&twon, a, x, &one, y, &one);

}


/* declare blas functions */
void FC_FUNC(sgemm, SGEMM)(const char * transa, const char * transb,
			   const int * m, const int * n, const int * k,
			   const float * alpha, const float * a, const int * lda,
			   const float * b, const int * ldb, const float * beta, 
			   float * c, const int * ldc);

void FC_FUNC(dgemm, DGEMM)(const char * transa, const char * transb,
			   const int * m, const int * n, const int * k,
			   const double * alpha, const double * a, const int * lda,
			   const double * b, const int * ldb, const double * beta, 
			   double * c, const int * ldc);

/* interface to apply dgemm passing complex vectors
   the same as dgemm, but allows giving each an appropriate Fortan interface
   in which alpha, beta, a, b, c are actually complex in Fortran 
   Could be inline, but in that case pgcc will not put it in the symbol table. */
void FC_FUNC(zdgemm, ZDGEMM)(const char * transa, const char * transb,
			     const int * m, const int * n, const int * k,
			     const double * alpha, const double * restrict a, const int * lda,
			     const double * restrict b, const int * ldb, const double * beta, 
			     double * restrict c, const int * ldc) {
  FC_FUNC(dgemm, DGEMM)(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

void FC_FUNC(csgemm, CSGEMM)(const char * transa, const char * transb,
			     const int * m, const int * n, const int * k,
			     const float * alpha, const float * restrict a, const int * lda,
			     const float * restrict b, const int * ldb, const float * beta, 
			     float * restrict c, const int * ldc) {

  FC_FUNC(sgemm, SGEMM)(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

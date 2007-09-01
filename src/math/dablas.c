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

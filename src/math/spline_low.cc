/*
 Copyright (C) 2016 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira, X. Andrade

// Spline routines produced at the Lawrence Livermore National
// Laboratory.  Written by Xavier Andrade (xavier@llnl.gov), Erik
// Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).

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

*/

#include <config.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_spline.h>

#include "string_f.h"

#include <fortran_types.h>

#include <assert.h>
#include <iostream>
#include <vector>

/* Interpolation */
extern "C" void FC_FUNC_(oct_spline_end, OCT_SPLINE_END)(void **spl, void **acc) {
	gsl_spline_free((gsl_spline *)(*spl));
	gsl_interp_accel_free((gsl_interp_accel *)(*acc));
}

extern "C" void FC_FUNC_(oct_spline_fit, OCT_SPLINE_FIT)
	(const fint *nrc, const double *x, const double *y, void **spl, void **acc){

	/* the GSL headers actually specify size_t instead of const int for nrc */
	*acc = (void *)gsl_interp_accel_alloc();
	*spl = (void *)gsl_spline_alloc(gsl_interp_cspline, *nrc);	
	gsl_spline_init((gsl_spline *)(*spl), x, y, *nrc);
	fflush(stdout);
}

extern "C" double FC_FUNC_(oct_spline_eval, OCT_SPLINE_EVAL)
     (const double *x, const void **spl, void **acc){
	/* the GSL headers specify double x instead of const double x */
	return gsl_spline_eval((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}


template <typename Type, int stride>
void oct_spline_eval_array(fint nn, Type * xf, const void **spl, void **acc){
	for(fint ii = 0; ii < nn; ii++){
		xf[stride*ii] = Type(gsl_spline_eval((gsl_spline *)(*spl), xf[stride*ii], (gsl_interp_accel *)(*acc)));
	}
}

extern "C" void FC_FUNC_(oct_spline_eval_array, OCT_SPLINE_EVAL_ARRAY)
     (const fint * nn, double *xf, const void **spl, void **acc){
  oct_spline_eval_array<double, 1>(*nn, xf, spl, acc);
}

/* use a stride of 2 to store into just the real part of a Fortran complex array */
extern "C" void FC_FUNC_(oct_spline_eval_arrayz, OCT_SPLINE_EVAL_ARRAYZ)
  (const fint * nn, double *xf, const void **spl, void **acc){
  oct_spline_eval_array<double, 2>(*nn, xf, spl, acc);
}

/* This function returns the number of points with which a spline
	 was constructed (the size component of the gsl_spline struct). */
extern "C" fint FC_FUNC_(oct_spline_npoints, OCT_SPLINE_NPOINTS)(const void **spl, void **acc){
	return (fint)((gsl_spline *)(*spl))->size;
}

/* This function places in the x array the x values of a given spline spl*/ 
extern "C" void FC_FUNC_(oct_spline_x, OCT_SPLINE_X)(const void **spl, void **acc, double *x){

	int size, i;
	size = (int)((gsl_spline *)(*spl))->size;
	for(i=0; i<size; i++) x[i] = ((gsl_spline *)(*spl))->x[i];
}

/* This function places in the y array the y values of a given spline spl*/ 
extern "C" void FC_FUNC_(oct_spline_y, OCT_SPLINE_Y)(const void **spl, void **acc, double *y){
	int size, i;
	
	size = (int)((gsl_spline *)(*spl))->size;
	for(i=0; i<size; i++) y[i] = ((gsl_spline *)(*spl))->y[i];
}

/* Returns the integral of the spline stored in spl, between a and b */
extern "C" double FC_FUNC_(oct_spline_eval_integ, OCT_SPLINE_EVAL_INTEG)
     (const void **spl, const double *a, const double *b, void **acc){
	/* the GSL headers specify double a, double b */
	return gsl_spline_eval_integ((gsl_spline *)(*spl), *a, *b, (gsl_interp_accel *)(* acc));
}

extern "C" double FC_FUNC_(oct_spline_eval_integ_full, OCT_SPLINE_EVAL_INTEG_FULL)
     (const void **spl, void **acc) {
	/* the GSL headers specify double a, double b */
	const int size = (int)((gsl_spline *)(*spl))->size;
	const double a = ((gsl_spline *)(*spl))->x[0];
	const double b = ((gsl_spline *)(*spl))->x[size - 1];
	return gsl_spline_eval_integ((gsl_spline *)(*spl), a, b, (gsl_interp_accel *)(* acc));
}

/* Performs the derivative of a spline */
extern "C" double FC_FUNC_(oct_spline_eval_der, OCT_SPLINE_EVAL_DER)
     (const double *x, const void **spl, void **acc){
	/* the GSL headers specify double x */
	return gsl_spline_eval_deriv((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}

/* Performs the second derivative of a spline */
extern "C" double FC_FUNC_(oct_spline_eval_der2, OCT_SPLINE_EVAL_DER2)
	(const double *x, const void **spl, void **acc){
  /* the GSL headers specify double x */
  return gsl_spline_eval_deriv2((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}

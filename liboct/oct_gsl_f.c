/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

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
*/

#include "config.h"

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_chebyshev.h>


/* ---------------------- Interface to GSL functions ------------------------ */


/* Mathematical Functions */
double FC_FUNC_(oct_asinh, OCT_ASINH)
		 (double *x)
{
  return gsl_asinh(*x);
}

/* Special Functions */
double FC_FUNC_(oct_gamma, OCT_GAMMA)
		 (double *x)
{
  return gsl_sf_gamma(*x);
}

double FC_FUNC_(oct_sph_bessel, OCT_SPH_BESSEL)
     (int *l, double*x)
{
  return gsl_sf_bessel_jl(*l, *x);
}

double FC_FUNC_(oct_bessel, OCT_BESSEL)
     (int *n, double *x)
{
  return gsl_sf_bessel_Jn(*n, *x);
}

double FC_FUNC_(oct_bessel_j0, OCT_BESSEL_J0)
     (double *x)
{
  return gsl_sf_bessel_J0(*x);
}

double FC_FUNC_(oct_bessel_j1, OCT_BESSEL_J1)
     (double *x)
{
  return gsl_sf_bessel_J1(*x);
}

double FC_FUNC_(oct_bessel_k0, OCT_BESSEL_K0)
     (double *x)
{
  return gsl_sf_bessel_K0(*x);
}

double FC_FUNC_(oct_bessel_k1, OCT_BESSEL_K1)
     (double *x)
{
  return gsl_sf_bessel_K1(*x);
}

double FC_FUNC_(oct_erfc, OCT_ERFC)
		 (double *x)
{
	return gsl_sf_erfc(*x);
}

double FC_FUNC_(oct_erf, OCT_ERF)
		 (double *x)
{
	return gsl_sf_erf(*x);
}


/* Vectors and Matrices */


/* Permutations */


/* Linear Algebra */


/* Random Number Generation */
void FC_FUNC_(oct_ran_init, OCT_RAN_INIT)
     (gsl_rng **r)
{
  gsl_rng_env_setup();
  *r = gsl_rng_alloc(gsl_rng_default);
}

void FC_FUNC_(oct_ran_end, OCT_RAN_END)
     (gsl_rng **r)
{
  gsl_rng_free(*r);
}


/* Random Number Distributions */ 
double FC_FUNC_(oct_ran_gaussian, OCT_RAN_GAUSSIAN)
		(gsl_rng **r, double *sigma)
{
  return gsl_ran_gaussian(*r, *sigma);
}


/* Interpolation */
void FC_FUNC_(oct_spline_end, OCT_SPLINE_END)
		 (void **spl, void **acc)
{
	gsl_spline_free((gsl_spline *)(*spl));
	gsl_interp_accel_free((gsl_interp_accel *)(*acc));
}

void FC_FUNC_(oct_spline_fit, OCT_SPLINE_FIT)
		 (int *nrc, double *x, double *y, void **spl, void **acc)
{
	*acc = (void *)gsl_interp_accel_alloc();
	*spl = (void *)gsl_spline_alloc(gsl_interp_cspline, *nrc);	
	gsl_spline_init((gsl_spline *)(*spl), x, y, *nrc);
	fflush(stdout);
}

double FC_FUNC_(oct_spline_eval, OCT_SPLINE_EVAL)
		 (double *x, void **spl, void **acc)
{
	return gsl_spline_eval((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}


/* This function returns the number of points with which a spline
	 was constructed (the size component of the gsl_spline struct). */
int FC_FUNC_(oct_spline_npoints, OCT_SPLINE_NPOINTS)
		 (void **spl)
{
	return (int)((gsl_spline *)(*spl))->size;
}

/* This function places in the x array the x values of a given spline spl*/ 
void FC_FUNC_(oct_spline_x, OCT_SPLINE_X)
     (void **spl, double *x)
{
  int size, i;
	
  size = (int)((gsl_spline *)(*spl))->size;
  for(i=0; i<size; i++)
		x[i] = ((gsl_spline *)(*spl))->x[i];
}

/* This function places in the y array the y values of a given spline spl*/ 
void FC_FUNC_(oct_spline_y, OCT_SPLINE_Y)
     (void **spl, double *y)
{
  int size, i;
	
  size = (int)((gsl_spline *)(*spl))->size;
  for(i=0; i<size; i++)
		y[i] = ((gsl_spline *)(*spl))->y[i];
}

/* Returns the integral of the spline stored in spl, between a and b */
double FC_FUNC_(oct_spline_eval_integ, OCT_SPLINE_EVAL_INTEG)
     (void **spl, double *a, double *b, void **acc)
{
  return gsl_spline_eval_integ((gsl_spline *)(*spl), *a, *b, (gsl_interp_accel *)(* acc));
}

/* Performs the derivative of a spline */
double FC_FUNC_(oct_spline_eval_der, OCT_SPLINE_EVAL_DER)
     (double *x, void **spl, void **acc)
{
  return gsl_spline_eval_deriv((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}


/* Chebyshev Approximations */
/*
void FC_FUNC_(oct_chebyshev_coeffs, OCT_CHEBYSHEV_COEFFS)
     (gsl_complex *coeffs, int *order)
{
  int i;
  double f (double x, void *p){return cos(x);}
  double g (double x, void *p){return -sin(x);}
  gsl_cheb_series *cs = gsl_cheb_alloc (*order);
  gsl_function F;
  F.function = f;
  F.params = 0;
  gsl_cheb_init (cs, &F, -1.0, 1.0);
  for(i=0; i<=12; i++){GSL_SET_REAL(&coeffs[i], (*cs).c[i]);}
  F.function = g;
  F.params = 0;
  gsl_cheb_init (cs, &F, -1.0, 1.0);
  for(i=0; i<=12; i++){GSL_SET_IMAG(&coeffs[i], (*cs).c[i]);}    
}
*/

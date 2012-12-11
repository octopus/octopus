/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira

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

 $Id$
*/

#include <config.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_deriv.h>

#include "string_f.h"

/* Numerical threshold for oct_bessel_k0 and oct_bessel_k1 */
#define  BESSEL_K_THRES  1.0e2

/* ---------------------- Interface to GSL functions ------------------------ */


/* Mathematical Functions */
double FC_FUNC_(oct_asinh, OCT_ASINH)
     (const double *x)
{
  return gsl_asinh(*x);
}

/* Special Functions */
double FC_FUNC_(oct_gamma, OCT_GAMMA)
     (const double *x)
{
  return gsl_sf_gamma(*x);
}

double FC_FUNC_(oct_hypergeometric, OCT_HYPERGEOMETRIC)
     (const double *a, const double*b, const double *x)
{
  return gsl_sf_hyperg_U(*a, *b, *x);
}

double FC_FUNC_(oct_incomplete_gamma, OCT_INCOMPLETE_GAMMA)
     (const double *a, const double *x)
{
  return gsl_sf_gamma_inc_Q(*a, *x);
}

double FC_FUNC_(oct_sph_bessel, OCT_SPH_BESSEL)
     (const int *l, const double*x)
{
  return gsl_sf_bessel_jl(*l, *x);
}

double FC_FUNC_(oct_bessel, OCT_BESSEL)
     (const int *n, const double *x)
{
  return gsl_sf_bessel_Jn(*n, *x);
}

double FC_FUNC_(oct_bessel_in, OCT_BESSEL_IN)
     (const int *n, const double *x)
{
  return gsl_sf_bessel_In(*n, *x);
}

double FC_FUNC_(oct_bessel_j0, OCT_BESSEL_J0)
     (const double *x)
{
  return gsl_sf_bessel_J0(*x);
}

double FC_FUNC_(oct_bessel_j1, OCT_BESSEL_J1)
     (const double *x)
{
  return gsl_sf_bessel_J1(*x);
}

double FC_FUNC_(oct_bessel_k0, OCT_BESSEL_K0)
     (const double *x)
{
  if( *x > BESSEL_K_THRES )
  {
    return 0.0e0;       
  } else {
    return gsl_sf_bessel_K0(*x);
  }    
}

double FC_FUNC_(oct_bessel_k1, OCT_BESSEL_K1)
     (const double *x)
{
  if( *x > BESSEL_K_THRES )
  {
    return 0.0e0;       
  } else {
    return gsl_sf_bessel_K1(*x);
  }    
}

/* the GSL headers specify double x, not const double x */
double FC_FUNC_(oct_erfc, OCT_ERFC)
     (const double *x)
{
  /* avoid floating invalids in the asymptotic limit */
  if(*x >  20.0) return  0.0;
  if(*x < -20.0) return  2.0;
  /* otherwise call gsl */
  return gsl_sf_erfc(*x);
}

/* the GSL headers specify double x, not const double x */
double FC_FUNC_(oct_erf, OCT_ERF)
     (const double *x)
{
  /* avoid floating invalids in the asymptotic limit */
  if(*x >  20.0) return  1.0;
  if(*x < -20.0) return -1.0;
  /* otherwise call gsl */
  return gsl_sf_erf(*x);
}

double FC_FUNC_(oct_legendre_sphplm, OCT_LEGENDRE_SPHPLM)
     (const int *l, const int *m, const double *x)
{
  return gsl_sf_legendre_sphPlm(*l, *m, *x);
}

double FC_FUNC_(oct_sine_integral, OCT_SINE_INTEGRAL)
     (const double *x)
{
  return gsl_sf_Si(*x);
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
     (const gsl_rng **r, const double *sigma)
{
  return gsl_ran_gaussian(*r, *sigma);
}


double FC_FUNC_(oct_ran_flat, OCT_RAN_FLAT)
     (const gsl_rng **r, const double *a, const double *b)
{
  return gsl_ran_flat(*r, *a, *b);
}


/* Interpolation */
void FC_FUNC_(oct_spline_end, OCT_SPLINE_END)
     (void **spl, void **acc)
{
  gsl_spline_free((gsl_spline *)(*spl));
  gsl_interp_accel_free((gsl_interp_accel *)(*acc));
}

void FC_FUNC_(oct_spline_fit, OCT_SPLINE_FIT)
		 (const int *nrc, const double *x, const double *y, void **spl, void **acc)
{
  /* the GSL headers actually specify size_t instead of const int for nrc */
  *acc = (void *)gsl_interp_accel_alloc();
  *spl = (void *)gsl_spline_alloc(gsl_interp_cspline, *nrc);	
  gsl_spline_init((gsl_spline *)(*spl), x, y, *nrc);
  fflush(stdout);
}

double FC_FUNC_(oct_spline_eval, OCT_SPLINE_EVAL)
     (const double *x, const void **spl, void **acc)
{
  /* the GSL headers specify double x instead of const double x */
  return gsl_spline_eval((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}


void FC_FUNC_(oct_spline_eval_array, OCT_SPLINE_EVAL_ARRAY)
     (const int * nn, double *xf, const void **spl, void **acc)
{
  int ii;
  for(ii = 0; ii < *nn; ii++){
    xf[ii] = gsl_spline_eval((gsl_spline *)(*spl), xf[ii], (gsl_interp_accel *)(*acc));
  }
}

void FC_FUNC_(oct_spline_eval_array4, OCT_SPLINE_EVAL_ARRAY)
     (const int * nn, float *xf, const void **spl, void **acc)
{
  int ii;
  for(ii = 0; ii < *nn; ii++){
    xf[ii] = (float) gsl_spline_eval((gsl_spline *)(*spl), (double) xf[ii], (gsl_interp_accel *)(*acc));
  }
}

/* This function returns the number of points with which a spline
	 was constructed (the size component of the gsl_spline struct). */
int FC_FUNC_(oct_spline_npoints, OCT_SPLINE_NPOINTS)
     (const void **spl)
{
  return (int)((gsl_spline *)(*spl))->size;
}

/* This function places in the x array the x values of a given spline spl*/ 
void FC_FUNC_(oct_spline_x, OCT_SPLINE_X)
     (const void **spl, double *x)
{
  int size, i;
	
  size = (int)((gsl_spline *)(*spl))->size;
  for(i=0; i<size; i++)
    x[i] = ((gsl_spline *)(*spl))->x[i];
}

/* This function places in the y array the y values of a given spline spl*/ 
void FC_FUNC_(oct_spline_y, OCT_SPLINE_Y)
     (const void **spl, double *y)
{
  int size, i;
	
  size = (int)((gsl_spline *)(*spl))->size;
  for(i=0; i<size; i++)
    y[i] = ((gsl_spline *)(*spl))->y[i];
}

/* Returns the integral of the spline stored in spl, between a and b */
double FC_FUNC_(oct_spline_eval_integ, OCT_SPLINE_EVAL_INTEG)
     (const void **spl, const double *a, const double *b, void **acc)
{
  /* the GSL headers specify double a, double b */
  return gsl_spline_eval_integ((gsl_spline *)(*spl), *a, *b, (gsl_interp_accel *)(* acc));
}

/* Performs the derivative of a spline */
double FC_FUNC_(oct_spline_eval_der, OCT_SPLINE_EVAL_DER)
     (const double *x, const void **spl, void **acc)
{
  /* the GSL headers specify double x */
  return gsl_spline_eval_deriv((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}

/* Performs the second derivative of a spline */
double FC_FUNC_(oct_spline_eval_der2, OCT_SPLINE_EVAL_DER2)
     (const double *x, const void **spl, void **acc)
{
  /* the GSL headers specify double x */
  return gsl_spline_eval_deriv2((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}

void FC_FUNC_(oct_strerror, OCT_STRERROR)
     (const int *err, STR_F_TYPE res STR_ARG1)
{
  const char *c;

  c = gsl_strerror(*err);
  TO_F_STR1(c, res);
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

 /* Numerical Derivatives.
    The following is an interface to the GSL support for numerical derivatives,
    in particular the gsl_deriv_central function. It computes a four points
    approximation to the derivative of a function, supplying an error estimation. */

 /* This is a type used to communicate with Fortran; func_nd is the type of the
    interface to a Fortran subroutine that calculates the value of the function.
    The first argument is the function argument, whereas the second is the function
    value. For convenience reasons, it is wrapped around the "param_nd_t" struct. */
typedef void (*func_nd)(double*, double*);
typedef struct{
  func_nd func;
} param_nd_t;

 /* This is the function that is called by the GSL function gsl_deriv_central. It
    receives as first argument the function argument, and as second argument
    a pointer to a params data type (a GSL data type), which in this case should
    be a pointer to a param_nd_t data type, where the address of the Fortran
    subroutine is. */
double function_oct_numerical_derivative (double x, void * params)
{
  double res;
  param_nd_t * p;

  p = (param_nd_t *) params;
  p->func(&x, &res);
  return res;
}
 /* This is the function that should be called by Fortran. The interface is defined
    in loct_math.F90 file. */
void FC_FUNC_(oct_numerical_derivative, OCT_NUMERICAL_DERIVATIVE)
     (const double *x, const double *h, double *result, double *abserr, const func_nd f)
{
  gsl_function F;
  param_nd_t p;

  p.func = f;
  F.function = &function_oct_numerical_derivative;
  F.params = (void *) &p;
  /* the GSL headers specify double x, double h */
  gsl_deriv_central (&F, *x, *h, result, abserr);
  return;
}

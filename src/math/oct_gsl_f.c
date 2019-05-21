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
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <config.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "string_f.h"

#include <fortran_types.h>

/* Numerical threshold for oct_bessel_k0 and oct_bessel_k1 */
#define  BESSEL_K_THRES  1.0e2

/* ---------------------- Interface to GSL functions ------------------------ */


/* Special Functions */
double FC_FUNC_(oct_gamma, OCT_GAMMA)
     (const double *x)
{
  return gsl_sf_gamma(*x);
}

double FC_FUNC_(oct_incomplete_gamma, OCT_INCOMPLETE_GAMMA)
     (const double *a, const double *x)
{
  return gsl_sf_gamma_inc_Q(*a, *x);
}

double FC_FUNC_(oct_sph_bessel, OCT_SPH_BESSEL)
     (const fint *l, const double*x)
{
  return gsl_sf_bessel_jl(*l, *x);
}

double FC_FUNC_(oct_bessel, OCT_BESSEL)
     (const fint *n, const double *x)
{
  return gsl_sf_bessel_Jn(*n, *x);
}

double FC_FUNC_(oct_bessel_in, OCT_BESSEL_IN)
     (const fint *n, const double *x)
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
     (const fint *l, const int *m, const double *x)
{
  return gsl_sf_legendre_sphPlm(*l, *m, *x);
}

/* generalized Laguerre polynomials */
double FC_FUNC_(oct_sf_laguerre_n, OCT_SF_LAGUERRE_N)
     (const int *n, const double *a, const double *x)
{
  return gsl_sf_laguerre_n(*n, *a, *x);
}


/* Combinations */

void FC_FUNC_(oct_combination_init, OCT_COMBINATION_INIT)
     (gsl_combination **c, const fint *n, const fint *k)
{
  *c = gsl_combination_calloc (*n, *k);
}

void FC_FUNC_(oct_combination_next, OCT_COMBINATION_NEXT)
     (gsl_combination **c, fint *next)
{
  *next = gsl_combination_next (((gsl_combination *)(*c)));
}

void FC_FUNC_(oct_get_combination, OCT_GET_COMBINATION)
     (gsl_combination **c, fint *comb)
{
  int i;
  for (i=0;i< ((gsl_combination *)(*c))->k; i++) {
    comb[i]=(fint)((gsl_combination *)(*c))->data[i];  
  }
}

void FC_FUNC_(oct_combination_end, OCT_COMBINATION_END)
     (gsl_combination **c)
{
  gsl_combination_free (((gsl_combination *)(*c)));
}


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

void FC_FUNC_(oct_strerror, OCT_STRERROR)
     (const fint *err, STR_F_TYPE res STR_ARG1)
{
  const char *c;

  c = gsl_strerror(*err);
  TO_F_STR1(c, res);
}

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
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "string_f.h" /* fortran <-> c string compatibility issues */

/* *********************** interface functions ********************** */

void F90_FUNC_(oct_mkdir, OCT_MKDIR)
		 (STR_F_TYPE name STR_ARG1)
{
	struct stat buf;
	char *name_c;

	name_c = TO_C_STR1(name);
	if(!*name_c || stat(name_c, &buf) == 0) return;
	mkdir(name_c, 0775);
	free(name_c);
}

void F90_FUNC_(oct_getcwd, OCT_GETCWD)
		 (STR_F_TYPE name STR_ARG1)
{
	char s[256];
  getcwd(s, 256);
	TO_F_STR1(s, name);
}

/* this function gets a string of the form '1-12, 34' and fills
	 array l with the 1 if the number is in the list, or 0 otherwise */
void F90_FUNC_(oct_wfs_list, OCT_WFS_LIST)
		 (STR_F_TYPE str, int l[32] STR_ARG1)
{
	int i, i1, i2;
	char c[20], *c1, *str_c, *s;

	str_c = TO_C_STR1(str);
	s = str_c;

	/* clear list */
	for(i=0; i<32; i++)
		l[i] = 0;
	
	while(*s){
		/* get integer */
		for(c1 = c; isdigit(*s) || isspace(*s); s++)
			if(isdigit(*s)) *c1++ = *s;
		*c1 = '\0';
		i1 = atoi(c) - 1;

		if(*s == '-'){ /* range */
			s++;
			for(c1 = c; isdigit(*s) || isspace(*s); s++)
				if(isdigit(*s)) *c1++ = *s;
			*c1 = '\0';
			i2 = atoi(c) - 1;
		}else /* single value */
			i2 = i1;

		for(i=i1; i<=i2; i++)
			if(i>=0 && i<1024)
				l[i/32] |= 1 << (i % 32);

		if(*s) s++;
	}

	free(str_c);
}

/* ---------------------- Interface to GSL math functions ------------------------ */
double F90_FUNC_(oct_gamma, OCT_GAMMA)
		 (double *x)
{
  return gsl_sf_gamma(*x);
}

double F90_FUNC_(oct_bessel, OCT_BESSEL)
     (int *n, double *x)
{
  return gsl_sf_bessel_Jn(*n, *x);
}

/* Fortran does not have the asinh intrinsic, */
double F90_FUNC_(oct_asinh, OCT_ASINH)
		 (double *x)
{
  return gsl_asinh(*x);
}

/* complementary error function (we use the one in gsl) */
double F90_FUNC_(oct_erfc, OCT_ERFC)
		 (double *x)
{
	return gsl_sf_erfc(*x);
}

/* error function (we use the one in gsl) */
double F90_FUNC_(oct_erf, OCT_ERF)
		 (double *x)
{
	return gsl_sf_erf(*x);
}

/* from ylm.c */
double ylm(double x, double y, double z, int l, int m);

double F90_FUNC_(oct_ylm, OCT_YLM)
		 (double *x, double *y, double *z, int *l, int *m)
{
	return ylm(*x, *y, *z, *l, *m);
}

/* ----------- Interface to GSL the GSL interpolation functions --------- */
void F90_FUNC_(oct_spline_end, OCT_SPLINE_END)
		 (void **spl, void **acc)
{
	gsl_spline_free((gsl_spline *)(*spl));
	gsl_interp_accel_free((gsl_interp_accel *)(*acc));
}

void F90_FUNC_(oct_spline_fit, OCT_SPLINE_FIT)
		 (int *nrc, double *x, double *y, void **spl, void **acc)
{
	*acc = (void *)gsl_interp_accel_alloc();
	*spl = (void *)gsl_spline_alloc(gsl_interp_cspline, *nrc);	
	gsl_spline_init((gsl_spline *)(*spl), x, y, *nrc);
	fflush(stdout);
}

double F90_FUNC_(oct_spline_eval, OCT_SPLINE_EVAL)
		 (double *x, void **spl, void **acc)
{
	return gsl_spline_eval((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}

/* ---------------------------- Random number generation --------------------- */
double F90_FUNC_(oct_ran_gaussian, OCT_RAN_GAUSSIAN)
		(gsl_rng **r, double *sigma)
{
  return gsl_ran_gaussian(*r, *sigma);
}

void F90_FUNC_(oct_ran_init, OCT_RAN_INIT)
     (gsl_rng **r)
{
  gsl_rng_env_setup();
  *r = gsl_rng_alloc(gsl_rng_default);
}

void F90_FUNC_(oct_ran_end, OCT_RAN_END)
     (gsl_rng **r)
{
  gsl_rng_free(*r);
}

/* ------------------------------ from varia.c ------------------------------- */
#include "varia.h"

void F90_FUNC_(oct_fft_optimize, OCT_FFT_OPTIMIZE)
		 (int *n, int *p, int *par)
{
	fft_optimize(n, *p, *par);
}

void F90_FUNC_(oct_progress_bar, OCT_PROGRESS_BAR)
		 (int *a, int *max)
{
	progress_bar(*a, *max);
}

/* ------------------------------ some stuff  -------------------------------- */
double F90_FUNC_(oct_clock, OCT_CLOCK)
       ()
{
  return (double) clock();
}

/* this function is *not* portable. Should get rid of this! */
int F90_FUNC_(oct_getmem, OCT_GETMEM)
     ()
{
#ifdef linux
	static size_t pagesize = 0;
	FILE *f;
  int pid;
	long mem;
  char s[256];
   
	if(pagesize == 0)
		pagesize = sysconf(_SC_PAGESIZE);

	pid = getpid();
	sprintf(s, "%s%d%s", "/proc/", pid, "/statm");
	if((f = fopen(s, "r")) == NULL) return -1;
	fscanf(f, "%lu", &mem);
	fclose(f);

	return (mem*pagesize) >> 10;
#else
	return -1;
#endif
}


void F90_FUNC_(oct_sysname, OCT_SYSNAME)
		 (STR_F_TYPE name STR_ARG1)
{
	char *name_c;
	
	name_c = TO_C_STR1(name);
	sysname(&name_c);
	TO_F_STR1(name_c, name);
	free(name_c);
}

int F90_FUNC_(number_of_lines, NUMBER_OF_LINES)
     (STR_F_TYPE name STR_ARG1)
{

  FILE *pf;
  int c, i;
  char *name_c;

	name_c = TO_C_STR1(name);
  pf = fopen(name_c, "r");
	free(name_c);

  if (pf != NULL) {
    i = 0;
    while ((c = getc(pf)) != EOF) {
      if (c == '\n') i++;
    }
    fclose(pf); 
    return i;
  }else{
    return -1;
  }
}

/*
void F90_FUNC_(oct_chebyshev_coeffs, OCT_CHEBYSHEV_COEFFS)
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

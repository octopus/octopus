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

#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "symbols.h"
#include "liboct.h"

double F90_FUNC_(oct_gamma, OCT_GAMMA)
		 (double *x)
{
  return gsl_sf_gamma(*x);
}

/* Fortran does not have the asinh intrinsic, 
	 so we use the one from libm.a */
double F90_FUNC_(oct_asinh, OCT_ASINH)
		 (double *x)
{
  return asinh(*x);
}

/* complementary error function (we use the one in gsl) */
double F90_FUNC_(oct_erfc, OCT_ERFC)
		 (double *x)
{
	return gsl_sf_erfc(*x);
}

/* error function (we use the one in gsl) */
double F90_FUNC_(oct_erf, OCT_ERC)
		 (double *x)
{
	return gsl_sf_erf(*x);
}

double F90_FUNC_(oct_ylm, OCT_YLM)
		 (double *x, double *y, double *z, int *l, int *m)
{
	return ylm(*x, *y, *z, *l, *m);
}

/* from varia.c */
void F90_FUNC_(oct_fft_optimize, OCT_FFT_OPTIMIZE)
		 (int *n, int *p, int *par)
{
	fft_optimize(n, *p, *par);
}

/* Interface to the GSL interpolation functions */
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

/* Interface to the parsing routines */
int F90_FUNC_(oct_parse_init, OCT_PARSE_INIT)
		 (char *in, char *out)
{ 
	return parse_init(in, out); 
}

void F90_FUNC_(oct_parse_end, OCT_PARSE_END)
		 ()
{ 
	parse_end(); 
}

int F90_FUNC_(oct_parse_isdef, OCT_PARSE_ISDEF)
		 (char *name)
{ 
	return parse_isdef(name); 
}

void F90_FUNC_(oct_parse_int, OCT_PARSE_INT)
		 (char *name, int *def, int *res)
{ 
	*res = parse_int(name, *def); 
}

void F90_FUNC_(oct_parse_double, OCT_PARSE_DOUBLE)
		 (char *name, double *def, double *res)
{
	*res = parse_double(name, *def); 
}

void F90_FUNC_(oct_parse_complex, OCT_PARSE_COMPLEX)
		 (char *name, gsl_complex *def, gsl_complex *res)
{
	*res = parse_complex(name, *def); 
}

void F90_FUNC_(oct_parse_string, OCT_PARSE_STRING)
		 (char *name, char *def, char *res)
{
	char *c = parse_string(name, def);
	int len = strlen(c);
	strcpy(res, c);
	res[len] = res[len + 1];
}

static void parse_block_error(char *type, char *name, int l, int c){
	fprintf(stderr, "Error: block \"%s\" does not contain a %s in line %d and col %d",
					name, type, l, c);
	exit(1);
}

int F90_FUNC_(oct_parse_block_n, OCT_PARSE_BLOCK_N)
		 (char *name)
{
	return parse_block_n(name);
}

void F90_FUNC_(oct_parse_block_int, OCT_PARSE_BLOCK_INT)
		 (char *name, int *l, int *c, int *res)
{
	if(parse_block_int(name, *l, *c, res) != 0)
		parse_block_error("int", name, *l, *c);
}

void F90_FUNC_(oct_parse_block_double, OCT_PARSE_BLOCK_DOUBLE)
		 (char *name, int *l, int *c, double *res)
{
	if(parse_block_double(name, *l, *c, res) != 0)
		parse_block_error("double", name, *l, *c);
}

void F90_FUNC_(oct_parse_block_complex, OCT_PARSE_BLOCK_COMPLEX)
		 (char *name, int *l, int *c, gsl_complex *res)
{
	if(parse_block_complex(name, *l, *c, res) != 0)
		parse_block_error("complex", name, *l, *c);
}

void F90_FUNC_(oct_parse_block_string, OCT_PARSE_BLOCK_STRING)
		 (char *name, int *l, int *c, char *res)
{
	char *s;
	int len;

	if(parse_block_string(name, *l, *c, &s) != 0)
		parse_block_error("string", name, *l, *c);
	else{
		len = strlen(s);
		strcpy(res, s);
		res[len] = res[len + 1];
	}
}

double  F90_FUNC_(oct_parse_potential, OCT_PARSE_POTENTIAL)
		 (double *x, double *y, double *z, double *r, char *pot)
{
	symrec *rec;
	parse_result c;

	rec = putsym("x", S_CMPLX);
	GSL_SET_COMPLEX(&rec->value.c, *x, 0);
	rec = putsym("y", S_CMPLX);
	GSL_SET_COMPLEX(&rec->value.c, *y, 0);
	rec = putsym("z", S_CMPLX);
	GSL_SET_COMPLEX(&rec->value.c, *z, 0);
	rec = putsym("r", S_CMPLX);
	GSL_SET_COMPLEX(&rec->value.c, *r, 0);

	parse_exp(pot, &c);

	rmsym("x");
	rmsym("y");
	rmsym("z");
	rmsym("r");

	return GSL_REAL(c.value.c);
}

void F90_FUNC_(oct_mkdir, OCT_MKDIR)
		 (char *name)
{
	struct stat buf;
	if(!*name || stat(name, &buf) == 0) return;

	mkdir(name, 0775);
}

void F90_FUNC_(oct_getcwd, OCT_GETCWD)
		 (char *name)
{
  getcwd(name, 256);
}

/* this function gets a string of the form '1-12, 34' and fills
	 array l with the 1 if the number is in the list, or 0 otherwise */
void F90_FUNC_(oct_wfs_list, OCT_WFS_LIST)
		 (char *str, int l[32])
{
	int i, i1, i2;
	char c[20], *c1, *s = str;

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
}

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

double F90_FUNC_(oct_clock, OCT_CLOCK)
       ()
{
  return (double) clock();
}

int F90_FUNC_(oct_getmem, OCT_GETMEM)
     ()
{
  FILE *pf;
  int pid, memory;
  char inst[100];

  mkdir("tmp",0775);
  if ( (pf = fopen("tmp/.tmp","w")) == NULL) return -1;
  pid = getpid();
  sprintf(inst, "ps -p %d -o rsz | tail -n 1 > tmp/.tmp 2> tmp/.tmp", pid );
  if ( system(inst) != 0) return -1;
  fclose(pf);
  if ( (pf = fopen("tmp/.tmp","r")) == NULL) return -1;
  if(fscanf(pf,"%d",&memory) < 1) return -1;
  fclose(pf);
  system("rm -f tmp/.tmp");
  return memory;
}

int F90_FUNC_(print_file, PRINT_FILE)
     (char *name)
{
#define MAXLINEA 256

  FILE *pf;
  char linea[MAXLINEA];

  pf = fopen(name, "r");
  if (pf != NULL) {
    while(fgets(linea, MAXLINEA, pf) != NULL)
      fputs(linea, stdout);
    fclose(pf); 
    fflush(stdout); 
    return 0;
  }else{
    fflush(stdout); 
    return 1;
  }
}

int F90_FUNC_(number_of_lines, NUMBER_OF_LINES)
     (char *name)
{

  FILE *pf;
  int c, i;

  pf = fopen(name, "r");
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

void F90_FUNC_(oct_progress_bar, OCT_PROGRESS_BAR)
		 (int *a, int *max)
{
	progress_bar(*a, *max);
}

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

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "symbols.h"
#include "liboct.h"

/* --------------------- Fortran to C string compatibility ---------------------- */
#if defined(_CRAY)
include <fortran.h>

#define STR_F_TYPE _fcd
#define TO_C_STR1(s) to_c_str(s)
#define TO_C_STR2(s) to_c_str(s)
#define TO_C_STR3(s) to_c_str(s)
#define TO_F_STR1(c, f) to_f_str(c, (_fcd) f)
#define TO_F_STR2(c, f) to_f_str(c, (_fcd) f)
#define TO_F_STR3(c, f) to_f_str(c, (_fcd) f)
#define STR_ARG1
#define STR_ARG2
#define STR_ARG3

char *to_c_str(_fcd f)
{
	char *c;
	int slen;

	slen = _fcdlen(f);
	c = (char *)malloc(slen+1);
	strncpy(c, _fcdtocp(f), slen);
	c[slen] = '\0';
	return c;
}

void to_f_str(char *c, _fcd f)
{
	f = _cptofcd(s, strlen(s));
}
#else

#define STR_F_TYPE char *
#define TO_C_STR1(s) to_c_str(s, l1)
#define TO_C_STR2(s) to_c_str(s, l2)
#define TO_C_STR3(s) to_c_str(s, l3)
#define TO_F_STR1(c, f) to_f_str(c, f, l1)
#define TO_F_STR2(c, f) to_f_str(c, f, l2)
#define TO_F_STR3(c, f) to_f_str(c, f, l3)
#define STR_ARG1     , unsigned long l1
#define STR_ARG2     , unsigned long l1, unsigned long l2
#define STR_ARG3     , unsigned long l1, unsigned long l2, unsigned long l3

char *to_c_str(STR_F_TYPE f, unsigned long l) 
{
	char *c;
	int i;

	for(l--; l>=0; l--)                 // find length of fortran string
		if(f[l] != ' ') break;
	l++;                                // need space for th '\0'
	c = (char *)calloc(l, sizeof(char)); // alloc c string
	for(i=0; i<l; i++) c[i] = f[i];     // copy fortran string onto c string
	c[i] = '\0';                        // add '\0' to the end of the c string
	return c;
}

void to_f_str(char *c, STR_F_TYPE f, unsigned long l)
{
  int i;
	for(i=0; i<l && c[i]!='\0'; i++) // copy string
		f[i] = c[i];
	for(; i<l; i++)                  // fill the rest with whitespace
		f[i] = ' ';
}
#endif

/* --------------------- Interface to the parsing routines ---------------------- */
int F90_FUNC_(oct_parse_init, OCT_PARSE_INIT)
		 (STR_F_TYPE in, STR_F_TYPE out STR_ARG2)
{
	int r;

	in  = TO_C_STR1(in);
	out = TO_C_STR2(out);
	r = parse_init(in, out); 
	free(in); free(out);

	return r;
}

void F90_FUNC_(oct_parse_end, OCT_PARSE_END)
		 ()
{
	parse_end(); 
}

int F90_FUNC_(oct_parse_isdef, OCT_PARSE_ISDEF)
		 (STR_F_TYPE name STR_ARG1)
{ 
	int r;

	name = TO_C_STR1(name);
	r = parse_isdef(name); 
	free(name);

	return r;
}

void F90_FUNC_(oct_parse_int, OCT_PARSE_INT)
		 (STR_F_TYPE name, int *def, int *res STR_ARG1)
{ 
	name = TO_C_STR1(name);
	*res = parse_int(name, *def);
	free(name);
}

void F90_FUNC_(oct_parse_double, OCT_PARSE_DOUBLE)
		 (STR_F_TYPE name, double *def, double *res STR_ARG1)
{
	name = TO_C_STR1(name);
	*res = parse_double(name, *def);
	free(name);
}

void F90_FUNC_(oct_parse_complex, OCT_PARSE_COMPLEX)
		 (STR_F_TYPE name, gsl_complex *def, gsl_complex *res STR_ARG1)
{
	name = TO_C_STR1(name);
	*res = parse_complex(name, *def);
	free(name);
}

void F90_FUNC_(oct_parse_string, OCT_PARSE_STRING)
		 (STR_F_TYPE name, STR_F_TYPE def, STR_F_TYPE res STR_ARG3)
{
	char *c;

	name = TO_C_STR1(name);      // convert string to c strings
	def  = TO_C_STR2(def);
	c = parse_string(name, def); 
	TO_F_STR3(c, res);           // convert string to fortran
	free(name); free(def);       // this has to be *after* the to_f_str
                               // or we will have memory problems
}

static void parse_block_error(char *type, char *name, int l, int c){
	fprintf(stderr, "Error: block \"%s\" does not contain a %s in line %d and col %d\n",
					name, type, l, c);
	exit(1);
}

int F90_FUNC_(oct_parse_block_n, OCT_PARSE_BLOCK_N)
		 (STR_F_TYPE name STR_ARG1)
{
	int r;

	name = TO_C_STR1(name);
	r = parse_block_n(name);
	free(name);
	
	return r;
}

void F90_FUNC_(oct_parse_block_int, OCT_PARSE_BLOCK_INT)
		 (STR_F_TYPE name, int *l, int *c, int *res STR_ARG1)
{
	name = TO_C_STR1(name);
	if(parse_block_int(name, *l, *c, res) != 0)
		parse_block_error("int", name, *l, *c);
	free(name);
}

void F90_FUNC_(oct_parse_block_double, OCT_PARSE_BLOCK_DOUBLE)
		 (STR_F_TYPE name, int *l, int *c, double *res STR_ARG1)
{
	name = TO_C_STR1(name);
	if(parse_block_double(name, *l, *c, res) != 0)
		parse_block_error("double", name, *l, *c);
	free(name);
}

void F90_FUNC_(oct_parse_block_complex, OCT_PARSE_BLOCK_COMPLEX)
		 (STR_F_TYPE name, int *l, int *c, gsl_complex *res STR_ARG1)
{
	name = TO_C_STR1(name);
	if(parse_block_complex(name, *l, *c, res) != 0)
		parse_block_error("complex", name, *l, *c);
	free(name);
}

void F90_FUNC_(oct_parse_block_string, OCT_PARSE_BLOCK_STRING)
		 (STR_F_TYPE name, int *l, int *c, STR_F_TYPE res STR_ARG2)
{
	char *s;
	int len;

	name = TO_C_STR1(name);
	if(parse_block_string(name, *l, *c, &s) != 0)
		parse_block_error("string", name, *l, *c);
	else{
		TO_F_STR2(s, res);
	}
	free(name);
}

double  F90_FUNC_(oct_parse_potential, OCT_PARSE_POTENTIAL)
		 (double *x, double *y, double *z, double *r, STR_F_TYPE pot STR_ARG1)
{
	symrec *rec;
	parse_result c;

	pot = TO_C_STR1(pot);

	rec = putsym("x", S_CMPLX);
	GSL_SET_COMPLEX(&rec->value.c, *x, 0);
	rec = putsym("y", S_CMPLX);
	GSL_SET_COMPLEX(&rec->value.c, *y, 0);
	rec = putsym("z", S_CMPLX);
	GSL_SET_COMPLEX(&rec->value.c, *z, 0);
	rec = putsym("r", S_CMPLX);
	GSL_SET_COMPLEX(&rec->value.c, *r, 0);

	parse_exp(pot, &c);

	free(pot);  // clean up
	rmsym("x"); rmsym("y");	rmsym("z");	rmsym("r");

	return GSL_REAL(c.value.c);
}

void F90_FUNC_(oct_mkdir, OCT_MKDIR)
		 (STR_F_TYPE name STR_ARG1)
{
	struct stat buf;

	name = TO_C_STR1(name);
	if(!*name || stat(name, &buf) == 0) return;
	mkdir(name, 0775);
	free(name);
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
	char c[20], *c1, *s = str;

	str = TO_C_STR1(str);

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

	free(str);
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
void F90_FUNC_(oct_fft_optimize, OCT_FFT_OPTIMIZE)
		 (int *n, int *p, int *par)
{
	fft_optimize(n, *p, *par);
}

double F90_FUNC_(oct_clock, OCT_CLOCK)
       ()
{
  return (double) clock();
}

/* this function is *not* portable. Should get rid of this! */
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

int F90_FUNC_(number_of_lines, NUMBER_OF_LINES)
     (STR_F_TYPE name STR_ARG1)
{

  FILE *pf;
  int c, i;

	name = TO_C_STR1(name);
  pf = fopen(name, "r");
	free(name);

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

#include "config.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_spline.h>

#include "f77_func.h"
#include "liboct.h"

/* Fortran does not have the asinh intrinsic, 
	 so we use the one from libm.a */
double F77_FUNC_(asinh, ASINH)
		 (double *x)
{
  return asinh(*x);
}

/* complementary error function (we use the one in gsl) */
double F77_FUNC_(erfc, ERFC)
		 (double *x)
{
	return gsl_sf_erfc(*x);
}

/* error function (we use the one in gsl) */
double F77_FUNC_(erf, ERC)
		 (double *x)
{
	return gsl_sf_erf(*x);
}

double F77_FUNC_(ylm, YLM)
		 (double *x, double *y, double *z, int *l, int *m)
{
	return ylm(*x, *y, *z, *l, *m);
}

/* from varia.c */
void F77_FUNC_(fft_optimize, FFT_OPTIMIZE)
		 (int *n, int *p, int *par)
{
	fft_optimize(n, *p, *par);
}

/* Interface to the GSL interpolation functions */
void F77_FUNC_(spline_end, SPLINE_END)
		 (void **spl, void **acc)
{
	gsl_spline_free((gsl_spline *)(*spl));
	gsl_interp_accel_free((gsl_interp_accel *)(*acc));
}

void F77_FUNC_(spline_fit, SPLINE_FIT)
		 (int *nrc, double *x, double *y, void **spl, void **acc)
{
	*acc = (void *)gsl_interp_accel_alloc();
	*spl = (void *)gsl_spline_alloc(gsl_interp_cspline, *nrc);	
	gsl_spline_init((gsl_spline *)(*spl), x, y, *nrc);
	fflush(stdout);
}

double F77_FUNC_(spline_eval, SPLINE_EVAL)
		 (double *x, void **spl, void **acc)
{
	return gsl_spline_eval((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}

/* Interface to the parsing routines */
int F77_FUNC_(parse_init, PARSE_INIT)
		 (char *in, char *out)
{ 
	return parse_init(in, out); 
}

void F77_FUNC_(parse_end, PARSE_END)
		 ()
{ 
	parse_end(); 
}

int F77_FUNC_(parse_isdef, PARSE_ISDEF)
		 (char *name)
{ 
	return parse_isdef(name); 
}

void F77_FUNC_(parse_int, PARSE_INT)
		 (char *name, int *def, int *res)
{ 
	*res = parse_int(name, *def); 
}

void F77_FUNC_(parse_double, PARSE_DOUBLE)
		 (char *name, double *def, double *res)
{
	*res = parse_double(name, *def); 
}

void F77_FUNC_(parse_complex, PARSE_COMPLEX)
		 (char *name, gsl_complex *def, gsl_complex *res)
{
	*res = parse_complex(name, *def); 
}

void F77_FUNC_(parse_string, PARSE_STRING)
		 (char *name, char *def, char *res)
{
	char *c = parse_string(name, def);
	strcpy(res, c);
	res[strlen(res)] = ' ';
}

static void parse_block_error(char *type, char *name, int l, int c){
	fprintf(stderr, "Error: block \"%s\" does not contain a %s in line %d and col %d",
					name, type, l, c);
	exit(1);
}

int F77_FUNC_(parse_block_n, PARSE_BLOCK_N)
		 (char *name)
{
	return parse_block_n(name);
}

void F77_FUNC_(parse_block_int, PARSE_BLOCK_INT)
		 (char *name, int *l, int *c, int *res)
{
	if(parse_block_int(name, *l, *c, res) != 0)
		parse_block_error("int", name, *l, *c);
}

void F77_FUNC_(parse_block_double, PARSE_BLOCK_DOUBLE)
		 (char *name, int *l, int *c, double *res)
{
	if(parse_block_double(name, *l, *c, res) != 0)
		parse_block_error("double", name, *l, *c);
}

void F77_FUNC_(parse_block_complex, PARSE_BLOCK_COMPLEX)
		 (char *name, int *l, int *c, gsl_complex *res)
{
	if(parse_block_complex(name, *l, *c, res) != 0)
		parse_block_error("complex", name, *l, *c);
}

void F77_FUNC_(parse_block_string, PARSE_BLOCK_STRING)
		 (char *name, int *l, int *c, char *res)
{
	char *s;

	if(parse_block_string(name, *l, *c, &s) != 0)
		parse_block_error("string", name, *l, *c);
	else{
		strcpy(res, s);
		res[strlen(res)] = ' ';
	}
}

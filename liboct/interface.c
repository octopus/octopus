#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_spline.h>

#include "config.h"
#include "liboct.h"

// Fortran does not have the asinh intrinsic, 
// so we use the one from libm.a
double PROTO(asinh)(double *x)
{
  return asinh(*x);
}

// complementary error function (we use the one in gsl)
double PROTO(erfc)(double *x)
{
	return gsl_sf_erfc(*x);
}

// error function (we use the one in gsl)
double PROTO(erf)(double *x)
{
	return gsl_sf_erf(*x);
}

double PROTO(ylm)(double *x, double *y, double *z, int *l, int *m)
{
	return ylm(*x, *y, *z, *l, *m);
}

// Interface to the GSL interpolation functions
void PROTO(spline_end)(void **spl, void **acc)
{
	gsl_spline_free((gsl_spline *)(*spl));
	gsl_interp_accel_free((gsl_interp_accel *)(*acc));
}

void PROTO(spline_fit)(int *nrc, double *x, double *y, void **spl, void **acc)
{
	*acc = (void *)gsl_interp_accel_alloc();
	*spl = (void *)gsl_spline_alloc(gsl_interp_cspline, *nrc);	
	gsl_spline_init((gsl_spline *)(*spl), x, y, *nrc);
	fflush(stdout);
}

double PROTO(spline_eval)(double *x, void **spl, void **acc)
{
	return gsl_spline_eval((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}

// Interface to the parsing routines
int  PROTO(parse_init) (char *in, char *out) { return parse_init(in, out); }
void PROTO(parse_end)  ()                    { parse_end(); }
int  PROTO(parse_isdef)(char *name)          { return parse_isdef(name); }

void PROTO(parse_int)   (char *name, int *def, int *res) { 
	*res = parse_int(name, *def); 
}

void PROTO(parse_double)(char *name, double *def, double *res){
	*res = parse_double(name, *def); 
}

void PROTO(parse_complex)(char *name, gsl_complex *def, gsl_complex *res) { 
	*res = parse_complex(name, *def); 
}

void PROTO(parse_string)(char *name, char *def, char *res){
	char *c = parse_string(name, def);
	strcpy(res, c);
	res[strlen(res)] = ' ';
}

static void parse_block_error(char *type, char *name, int l, int c){
	fprintf(stderr, "Error: block \"%s\" does not contain a %s in line %d and col %d",
					name, type, l, c);
	exit(1);
}

void PROTO(parse_block_n)(char *name){
	return parse_block_n(name);
}

void PROTO(parse_block_int)(char *name, int *l, int *c, int *res){
	if(parse_block_int(name, *l, *c, res) != 0)
		parse_block_error("int", name, *l, *c);
}

void PROTO(parse_block_double)(char *name, int *l, int *c, double *res){
	if(parse_block_double(name, *l, *c, res) != 0)
		parse_block_error("double", name, *l, *c);
}

void PROTO(parse_block_complex)(char *name, int *l, int *c, gsl_complex *res){
	if(parse_block_complex(name, *l, *c, res) != 0)
		parse_block_error("complex", name, *l, *c);
}

void PROTO(parse_block_string)(char *name, int *l, int *c, char *res){
	char *s;

	if(parse_block_string(name, *l, *c, &s) != 0)
		parse_block_error("string", name, *l, *c);
	else{
		strcpy(res, s);
		res[strlen(res)] = ' ';
	}
}

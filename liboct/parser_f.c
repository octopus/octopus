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
#include "parser.h"
#include "string_f.h"

/* --------------------- Interface to the parsing routines ---------------------- */
int F90_FUNC_(oct_parse_init, OCT_PARSE_INIT)
		 (STR_F_TYPE in, STR_F_TYPE out STR_ARG2)
{
	int r;
	char *in_c, *out_c;

	in_c  = TO_C_STR1(in);
	out_c = TO_C_STR2(out);
	r = parse_init(in_c, out_c); 
	free(in_c); free(out_c);

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
	char *name_c;

	name_c = TO_C_STR1(name);
	r = parse_isdef(name_c); 
	free(name_c);

	return r;
}

void F90_FUNC_(oct_parse_int, OCT_PARSE_INT)
		 (STR_F_TYPE name, int *def, int *res STR_ARG1)
{ 
	char *name_c;
	name_c = TO_C_STR1(name);
	*res = parse_int(name_c, *def);
	free(name_c);
}

void F90_FUNC_(oct_parse_double, OCT_PARSE_DOUBLE)
		 (STR_F_TYPE name, double *def, double *res STR_ARG1)
{
	char *name_c;
	name_c = TO_C_STR1(name);
	*res = parse_double(name_c, *def);
	free(name_c);
}

void F90_FUNC_(oct_parse_complex, OCT_PARSE_COMPLEX)
		 (STR_F_TYPE name, gsl_complex *def, gsl_complex *res STR_ARG1)
{
	char *name_c;
	name_c = TO_C_STR1(name);
	*res = parse_complex(name_c, *def);
	free(name_c);
}

void F90_FUNC_(oct_parse_string, OCT_PARSE_STRING)
		 (STR_F_TYPE name, STR_F_TYPE def, STR_F_TYPE res STR_ARG3)
{
	char *c, *name_c, *def_c;

	name_c = TO_C_STR1(name);      /* convert string to c strings */
	def_c  = TO_C_STR2(def);
	c = parse_string(name_c, def_c); 
	TO_F_STR3(c, res);             /* convert string to fortran */
	free(name_c); free(def_c);     /* this has to be *after* the to_f_str
                                          or we will have memory problems */
}

static void parse_block_error(char *type, char *name, int l, int c){
  fprintf(stderr, "Error: block \"%s\" does not contain a %s in line %d and col %d\n",
	  name, type, l, c);
  /* exit(1); */
}

int F90_FUNC_(oct_parse_block_n, OCT_PARSE_BLOCK_N)
		 (STR_F_TYPE name STR_ARG1)
{
	int r;
	char *name_c;

	name_c = TO_C_STR1(name);
	r = parse_block_n(name_c);
	free(name_c);
	
	return r;
}

int F90_FUNC_(oct_parse_block_cols, OCT_PARSE_BLOCK_COLS)
		 (STR_F_TYPE name, int *l STR_ARG1)
{
	int r;
	char *name_c;

	name_c = TO_C_STR1(name);
	r = parse_block_cols(name_c, *l);
	free(name_c);
	
	return r;
}

int F90_FUNC_(oct_parse_block_int, OCT_PARSE_BLOCK_INT)
		 (STR_F_TYPE name, int *l, int *c, int *res STR_ARG1)
{
	char *name_c;
	int r;

	name_c = TO_C_STR1(name);
	r = parse_block_int(name_c, *l, *c, res);
	if(r!= 0) parse_block_error("int", name_c, *l, *c);
	free(name_c);
	return r;
		
}

void F90_FUNC_(oct_parse_block_double, OCT_PARSE_BLOCK_DOUBLE)
		 (STR_F_TYPE name, int *l, int *c, double *res STR_ARG1)
{
	char *name_c;

	name_c = TO_C_STR1(name);
	if(parse_block_double(name_c, *l, *c, res) != 0)
		parse_block_error("double", name_c, *l, *c);
	free(name_c);
}

void F90_FUNC_(oct_parse_block_complex, OCT_PARSE_BLOCK_COMPLEX)
		 (STR_F_TYPE name, int *l, int *c, gsl_complex *res STR_ARG1)
{
	char *name_c;

	name_c = TO_C_STR1(name);
	if(parse_block_complex(name_c, *l, *c, res) != 0)
		parse_block_error("complex", name_c, *l, *c);
	free(name_c);
}

void F90_FUNC_(oct_parse_block_string, OCT_PARSE_BLOCK_STRING)
		 (STR_F_TYPE name, int *l, int *c, STR_F_TYPE res STR_ARG2)
{
	char *s, *name_c;
	int len;

	name_c = TO_C_STR1(name);
	if(parse_block_string(name_c, *l, *c, &s) != 0)
		parse_block_error("string", name_c, *l, *c);
	else{
		TO_F_STR2(s, res);
	}
	free(name_c);
}

double  F90_FUNC_(oct_parse_potential, OCT_PARSE_POTENTIAL)
		 (double *x, double *y, double *z, double *r, STR_F_TYPE pot STR_ARG1)
{
	symrec *rec;
	parse_result c;
	char *pot_c;

	pot_c = TO_C_STR1(pot);

	rec = putsym("x", S_CMPLX);
	GSL_SET_COMPLEX(&rec->value.c, *x, 0);
	rec = putsym("y", S_CMPLX);
	GSL_SET_COMPLEX(&rec->value.c, *y, 0);
	rec = putsym("z", S_CMPLX);
	GSL_SET_COMPLEX(&rec->value.c, *z, 0);
	rec = putsym("r", S_CMPLX);
	GSL_SET_COMPLEX(&rec->value.c, *r, 0);

	parse_exp(pot_c, &c);

	free(pot_c);  /* clean up */
	rmsym("x"); rmsym("y");	rmsym("z");	rmsym("r");

	return GSL_REAL(c.value.c);
}

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

 $Id$
*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "liboct_parser.h"
#include "string_f.h"

/* --------------------- Interface to the parsing routines ---------------------- */

/* --------------------------------------------------------- */
/* Initialization of the library                             */
/* --------------------------------------------------------- */
int FC_FUNC_(oct_parse_init, OCT_PARSE_INIT)
  (STR_F_TYPE s, int *dont_write STR_ARG1)
{
  int r;
  char *s_c;
  
  s_c  = TO_C_STR1(s);
  r = parse_init(s_c, dont_write); 
  free(s_c);
  
  return r;
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_putsym_int, OCT_PARSE_PUTSYM_INT)
     (STR_F_TYPE s, int *i  STR_ARG1)
{
  char *s_c = TO_C_STR1(s);
  parse_putsym_int(s_c, *i);
  free(s_c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_putsym_double, OCT_PARSE_PUTSYM_DOUBLE)
     (STR_F_TYPE s, double *d STR_ARG1)
{
  char *s_c = TO_C_STR1(s);
  parse_putsym_double(s_c, *d);
  free(s_c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_putsym_complex, OCT_PARSE_PUTSYM_COMPLEX)
     (STR_F_TYPE s, gsl_complex *c STR_ARG1)
{
  char *s_c = TO_C_STR1(s);
  parse_putsym_complex(s_c, *c);
  free(s_c);
}


/* --------------------------------------------------------- */
int FC_FUNC_(oct_parse_input, OCT_PARSE_INIT)
     (STR_F_TYPE s STR_ARG1)
{
  int r;
  char *s_c;
  
  s_c  = TO_C_STR1(s);
  r = parse_input(s_c); 
  free(s_c);
  
  return r;
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_end, OCT_PARSE_END)
		 ()
{
  parse_end(); 
}


/* --------------------------------------------------------- */
/* Parser functions                                          */
/* --------------------------------------------------------- */
int FC_FUNC_(oct_parse_isdef, OCT_PARSE_ISDEF)
     (STR_F_TYPE name STR_ARG1)
{ 
  int r;
  char *name_c;
  
  name_c = TO_C_STR1(name);
  r = parse_isdef(name_c); 
  free(name_c);
  
  return r;
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_int, OCT_PARSE_INT)
     (STR_F_TYPE name, int *def, int *res STR_ARG1)
{ 
  char *name_c;
  name_c = TO_C_STR1(name);
  *res = parse_int(name_c, *def);
  free(name_c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_double, OCT_PARSE_DOUBLE)
     (STR_F_TYPE name, double *def, double *res STR_ARG1)
{
  char *name_c;
  name_c = TO_C_STR1(name);
  *res = parse_double(name_c, *def);
  free(name_c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_complex, OCT_PARSE_COMPLEX)
     (STR_F_TYPE name, gsl_complex *def, gsl_complex *res STR_ARG1)
{
  char *name_c;
  name_c = TO_C_STR1(name);
  *res = parse_complex(name_c, *def);
  free(name_c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_string, OCT_PARSE_STRING)
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


/* --------------------------------------------------------- */
static void parse_block_error(char *type, char *name, int l, int c){
  fprintf(stderr, "Error: block \"%s\" does not contain a %s in line %d and col %d\n",
	  name, type, l, c);
  /* exit(1); */
}
static char *block_name = NULL;


/* --------------------------------------------------------- */
int FC_FUNC_(oct_parse_block, OCT_PARSE_BLOCK)
     (STR_F_TYPE name, sym_block **blk STR_ARG1)
{
	int r;

  if(block_name != NULL){
          free(block_name);
          block_name = NULL;
  }

  block_name = TO_C_STR1(name);
  r = parse_block(block_name, blk);

	if(r != 0){
		free(block_name);
		block_name = NULL;
	}
	return r;
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_block_end, OCT_PARSE_BLOCK_END)
     (sym_block **blk)
{
  assert(block_name!=NULL);

  parse_block_end(blk);
  free(block_name);
	block_name = NULL;
}


/* --------------------------------------------------------- */
int FC_FUNC_(oct_parse_block_n, OCT_PARSE_BLOCK_N)
     (sym_block **blk)
{
  return parse_block_n(*blk);
}


/* --------------------------------------------------------- */
int FC_FUNC_(oct_parse_block_cols, OCT_PARSE_BLOCK_COLS)
     (sym_block **blk, int *l)
{
  return parse_block_cols(*blk, *l);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_block_int, OCT_PARSE_BLOCK_INT)
     (sym_block **blk, int *l, int *c, int *res)
{
  if(parse_block_int(*blk, *l, *c, res) != 0)
    parse_block_error("int", block_name, *l, *c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_block_double, OCT_PARSE_BLOCK_DOUBLE)
     (sym_block **blk, int *l, int *c, double *res)
{
  if(parse_block_double(*blk, *l, *c, res) != 0)
    parse_block_error("double", block_name, *l, *c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_block_complex, OCT_PARSE_BLOCK_COMPLEX)
     (sym_block **blk, int *l, int *c, gsl_complex *res)
{
  if(parse_block_complex(*blk, *l, *c, res) != 0)
    parse_block_error("complex", block_name, *l, *c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(oct_parse_block_string, OCT_PARSE_BLOCK_STRING)
     (sym_block **blk, int *l, int *c, STR_F_TYPE res STR_ARG1)
{
  char *s;
  
  if(parse_block_string(*blk, *l, *c, &s) != 0)
    parse_block_error("string", block_name, *l, *c);
  else{
    TO_F_STR1(s, res);
  }
}


/* --------------------------------------------------------- */
double FC_FUNC_(oct_parse_potential, OCT_PARSE_POTENTIAL)
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

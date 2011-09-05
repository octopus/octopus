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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>

#include "symbols.h"
#include "gsl_userdef.h"



/* The symbol table: a chain of `struct symrec'.  */
symrec *sym_table = (symrec *)0;

char *str_tolower(char *in)
{
  char *s;
  for(s=in; *s; s++)
    *s = tolower(*s);
  return in;
}

symrec *putsym (char *sym_name, symrec_type sym_type)
{
  symrec *ptr;
  ptr = (symrec *)malloc(sizeof(symrec));

  /* names are always lowercase */
  ptr->name = strdup(sym_name);
  str_tolower(ptr->name);
  
  ptr->def  = 0;
  ptr->type = sym_type;
  GSL_SET_COMPLEX(&ptr->value.c, 0, 0); /* set value to 0 even if fctn.  */
  ptr->next = (struct symrec *)sym_table;
  sym_table = ptr;
  return ptr;
}

symrec *getsym (char *sym_name)
{
  symrec *ptr;
  for (ptr = sym_table; ptr != (symrec *) 0;
       ptr = (symrec *)ptr->next)
    if (strcasecmp(ptr->name,sym_name) == 0)
      return ptr;
  return (symrec *) 0;
}

int rmsym (char *sym_name)
{
  symrec *ptr, *prev;
  for (prev = (symrec *) 0, ptr = sym_table; ptr != (symrec *) 0;
       prev = ptr, ptr = ptr->next)
    if (strcasecmp(ptr->name,sym_name) == 0){
      if(prev == (symrec *) 0)
	sym_table = ptr->next;
      else
	prev->next = ptr->next;
      free(ptr->name);
      free(ptr);
      
      return 1;
    }
  
  return 0;
}

struct init_fntc{
  char *fname;
  int  nargs;
  gsl_complex (*fnctptr)();
};

void sym_notdef (symrec *sym)
{
  fprintf(stderr, "Input error:\n\tsymbol '%s' used before being defined\n", sym->name);
  exit(1);
}

void sym_wrong_arg (symrec *sym)
{
  fprintf(stderr, "Input error:\n\tfunction '%s' accepts %d argument\n", sym->name, sym->nargs);
  exit(1);
}

static struct init_fntc arith_fncts[] = {
  {"sqrt",   1, (gsl_complex (*)()) &gsl_complex_sqrt},
  {"exp",    1, (gsl_complex (*)()) &gsl_complex_exp},
  {"ln",     1, (gsl_complex (*)()) &gsl_complex_log},
  {"log",    1, (gsl_complex (*)()) &gsl_complex_log},
  {"log10",  1, (gsl_complex (*)()) &gsl_complex_log10},
  {"logb",   2, (gsl_complex (*)()) &gsl_complex_log_b}, /* takes two arguments logb(z, b) = log_b(z) */

  {"arg",    1, (gsl_complex (*)()) &gsl_complex_carg},
  {"abs",    1, (gsl_complex (*)()) &gsl_complex_cabs},
  {"abs2",   1, (gsl_complex (*)()) &gsl_complex_cabs2},
  {"logabs", 1, (gsl_complex (*)()) &gsl_complex_clogabs},

  {"conjg",  1, (gsl_complex (*)()) &gsl_complex_conjugate},
  {"inv",    1, (gsl_complex (*)()) &gsl_complex_inverse},

  {"sin",    1, (gsl_complex (*)()) &gsl_complex_sin},
  {"cos",    1, (gsl_complex (*)()) &gsl_complex_cos},
  {"tan",    1, (gsl_complex (*)()) &gsl_complex_tan},
  {"sec",    1, (gsl_complex (*)()) &gsl_complex_sec},
  {"csc",    1, (gsl_complex (*)()) &gsl_complex_csc},
  {"cot",    1, (gsl_complex (*)()) &gsl_complex_cot},

  {"asin",   1, (gsl_complex (*)()) &gsl_complex_arcsin},
  {"acos",   1, (gsl_complex (*)()) &gsl_complex_arccos},
  {"atan",   1, (gsl_complex (*)()) &gsl_complex_arctan},
  {"atan2",  2, (gsl_complex (*)()) &gsl_complex_arctan2}, /* takes two arguments atan2(y,x) = atan(y/x) */
  {"asec",   1, (gsl_complex (*)()) &gsl_complex_arcsec},
  {"acsc",   1, (gsl_complex (*)()) &gsl_complex_arccsc},
  {"acot",   1, (gsl_complex (*)()) &gsl_complex_arccot},

  {"sinh",   1, (gsl_complex (*)()) &gsl_complex_sinh},
  {"cosh",   1, (gsl_complex (*)()) &gsl_complex_cosh},
  {"tanh",   1, (gsl_complex (*)()) &gsl_complex_tanh},
  {"sech",   1, (gsl_complex (*)()) &gsl_complex_sech},
  {"csch",   1, (gsl_complex (*)()) &gsl_complex_csch},
  {"coth",   1, (gsl_complex (*)()) &gsl_complex_coth},

  {"asinh",  1, (gsl_complex (*)()) &gsl_complex_arcsinh},
  {"acosh",  1, (gsl_complex (*)()) &gsl_complex_arccosh},
  {"atanh",  1, (gsl_complex (*)()) &gsl_complex_arctanh},
  {"asech",  1, (gsl_complex (*)()) &gsl_complex_arcsech},
  {"acsch",  1, (gsl_complex (*)()) &gsl_complex_arccsch},
  {"acoth",  1, (gsl_complex (*)()) &gsl_complex_arccoth},	
 
/* user defined step function. this is not available in GSL, 
   but we use GSL namespacing and macros here. */
  {"step",   1, (gsl_complex (*)()) &gsl_complex_step_real},

/* Minimum and maximum of two arguments (comparing real parts) */  
  {"min",    2, (gsl_complex (*)()) &gsl_complex_min_real},
  {"max",    2, (gsl_complex (*)()) &gsl_complex_max_real},

  {"erf",    1, (gsl_complex (*)()) &gsl_complex_erf}, 

  {0, 0}
};

struct init_cnst{
	char *fname;
	double re;
	double im;
};

static struct init_cnst arith_cnts[] = {
	{"pi",    M_PI, 0}, 
	{"e",      M_E, 0},
	{"i",        0, 1},
	{"true",     1, 0}, 
	{"t",        1, 0}, 
	{"yes",      1, 0},
	{"false",    0, 0}, 
	{"f",        0, 0}, 
	{"no",       0, 0},
	{0,          0, 0}
};

char *reserved_symbols[] = {
  "x", "y", "z", "r", "w", "t", 0
};

void sym_init_table ()  /* puts arithmetic functions in table. */
{
  int i;
  symrec *ptr;
  for (i = 0; arith_fncts[i].fname != 0; i++){
    ptr = putsym (arith_fncts[i].fname, S_FNCT);
    ptr->def = 1;
    ptr->nargs = arith_fncts[i].nargs;
    ptr->value.fnctptr = arith_fncts[i].fnctptr;
  }

  /* now the constants */
  for (i = 0; arith_cnts[i].fname != 0; i++){
    ptr = putsym(arith_cnts[i].fname, S_CMPLX);
    ptr->def = 1;
    GSL_SET_COMPLEX(&ptr->value.c, arith_cnts[i].re, arith_cnts[i].im);
  }
}

void sym_clear_reserved()
{
  int i;
  for (i = 0; reserved_symbols[i] != 0; i++){
    rmsym(reserved_symbols[i]);
  }
}

void sym_end_table()
{
  symrec *ptr, *ptr2;
  int l, col;

  for (ptr = sym_table; ptr != NULL;){
    free(ptr->name);
    switch(ptr->type){
    case S_STR:
      free(ptr->value.str);
      break;
    case S_BLOCK:
      if(ptr->value.block->n > 0){
	for(l = 0; l < ptr->value.block->n; l++){
	  if(ptr->value.block->lines[l].n > 0){
	    for(col = 0; col < ptr->value.block->lines[l].n; col++)
	      free(ptr->value.block->lines[l].fields[col]);
	    free(ptr->value.block->lines[l].fields);
	  }
	}
	free(ptr->value.block->lines);
      }
      free(ptr->value.block);
      break;
    case S_CMPLX:
    case S_FNCT:
      break;
    }
    ptr2 = ptr->next;
    free(ptr);
    ptr = ptr2;
  }
  
  sym_table = NULL;
}

void sym_output_table()
{
  symrec *ptr;
  for(ptr = sym_table; ptr != NULL; ptr = ptr->next){
    printf("%s", ptr->name);
    switch(ptr->type){
    case S_CMPLX:
      printf(" = (%f,%f)\n", GSL_REAL(ptr->value.c), GSL_IMAG(ptr->value.c));
      break;
    case S_STR:
      printf(" = \"%s\"\n", ptr->value.str);
      break;
    case S_BLOCK:
      printf("%s\n", " <= BLOCK");
      break;
    case S_FNCT:
      printf("%s\n", " <= FUNCTION");
      break;
    }
  }
}

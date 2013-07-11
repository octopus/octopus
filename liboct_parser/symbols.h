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
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

 $Id$
*/

#ifndef _SYMBOLS_H
#define _SYMBOLS_H

#include <gsl/gsl_complex.h>
#include "liboct_parser.h"

typedef enum{
  S_CMPLX, S_STR, S_BLOCK, S_FNCT
} symrec_type;

/* Data type for links in the chain of symbols. */
typedef struct symrec{
  char *name;                  /* name of symbol */
  symrec_type type;            /* type of symbol: either VAR or FNCT */
  int def;                     /* has this symbol been defined */

  int nargs;                   /* if type==FNCT contains the number of arguments of the function */

  union {
    gsl_complex c;             /* value of a VAR */
    char *str;                 /* value of a STRING */
    sym_block *block;          /* to store blocks */
    gsl_complex (*fnctptr)();  /* value of a FNCT */
  } value;

  struct symrec *next;         /* link field */
} symrec;

/* The symbol table: a chain of struct symrec. */
extern symrec *sym_table;
extern char *reserved_symbols[];

symrec *putsym (char *sym_name, symrec_type sym_type);
symrec *getsym (char *sym_name);
int      rmsym (char *sym_name);

void sym_notdef(symrec *sym);
void sym_wrong_arg(symrec *sym);
void sym_init_table();
void sym_clear_reserved();
void sym_end_table();
void sym_output_table();
char *str_tolower(char *in);

#endif

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

#ifndef _LIB_OCT_H
#define _LIB_OCT_H

#include <gsl/gsl_complex.h>

int         parse_init   (char *file_out, int *mpiv_node);
int         parse_input  (char *file_in);
void        parse_end    ();

int         parse_isdef  (char *name);
int         parse_int    (char *name, int def);
double      parse_double (char *name, double def);
gsl_complex parse_complex(char *name, gsl_complex def);
char       *parse_string (char *name, char *def);

/* Now comes stuff for the blocks */
typedef struct sym_block_line{
  int n;
  char **fields;
} sym_block_line;

typedef struct sym_block{
  int n;
  sym_block_line *lines;
} sym_block;

int parse_block        (char *name, sym_block **blk);
int parse_block_end    (sym_block **blk);
int parse_block_n      (sym_block *blk);
int parse_block_cols   (sym_block *blk, int l);
int parse_block_int    (sym_block *blk, int l, int col, int *r);
int parse_block_double (sym_block *blk, int l, int col, double *r);
int parse_block_complex(sym_block *blk, int l, int col, gsl_complex *r);
int parse_block_string (sym_block *blk, int l, int col, char **r);

/* from parse_exp.c */
enum pr_type {PR_NONE,PR_CMPLX, PR_STR};
typedef struct parse_result{
  union {
    gsl_complex c;
    char *s;
  } value;
  enum pr_type type;
} parse_result;

void parse_result_free(parse_result *t);

int parse_exp(char *exp, parse_result *t);

void parse_putsym_int(char *s, int i);
void parse_putsym_double(char *s, double d);
void parse_putsym_complex(char *s, gsl_complex c);

#endif

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

#include <stdio.h>
#include <stdlib.h>

extern char ** environ;

#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "liboct_parser.h"
#include "symbols.h"

static FILE *fout;
static int  disable_write;

#define ROUND(x) ((x)<0 ? (int)(x-0.5) : (int)(x+0.5)) 

void str_trim(char *in)
{
  char *c, *s = in;

  for(c=s; isspace(*c); c++);
  for(; *c != '\0'; *s++=*c++);
  for(s--; s>=in && isspace(*s); s--);
  *(s+1) = '\0';
}

static int parse_get_line(FILE *f, char **s, int *length)
{
  int i, c;

  i = 0;
  do{
    c = getc(f);
    if(c == '#') /* skip comments */
      while(c!=EOF && c!='\n') c = getc(f);
    else if(c != EOF){
      if (i == *length - 1){
	*length *= 2;
	*s = (char *)realloc(*s, *length + 1);
      }
      (*s)[i++] = (char)c;
    }
  }while(c != EOF && c != '\n');
  (*s)[i] = '\0';
	
  str_trim(*s);
  return c;
}

int parse_init(const char *file_out, const int *mpiv_node)
{
  sym_init_table();

  /* only enable writes for node 0 */
  disable_write = *mpiv_node;

  if(disable_write) return 0;

  if(strcmp(file_out, "-") == 0)
    fout = stdout;
  else {
    fout = fopen(file_out, "w");
    if(!fout)
      return -1; /* error opening file */
    setvbuf(fout, NULL, _IONBF, 0);
  }
  fprintf(fout, "# Octopus parser started\n");
  
  return 0;
}

int parse_input(const char *file_in)
{
  FILE *f;
  char *s;
  int c, length = 0;
  
  if(strcmp(file_in, "-") == 0)
    f = stdin;
  else
    f = fopen(file_in, "r");
  
  if(!f)
    return -1; /* error opening file */
  
  /* we now read in the file and parse */
  length = 40;
  s = (char *)malloc(length + 1);
  do{
    c = parse_get_line(f, &s, &length);
    if(*s){
      if(*s == '%'){ /* we have a block */
	*s = ' ';
	str_trim(s);
	if(getsym(s) != NULL){ /* error */
	  fprintf(stderr, "Parser warning: %s \"%s\" %s.\n", "Block", s, "already defined");
	  do{ /* skip block */
	    c = parse_get_line(f, &s, &length);
	  }while(c != EOF && *s != '%');
	}else{ /* parse block */
	  symrec *rec;
	  rec = putsym(s, S_BLOCK);
	  rec->value.block = (sym_block *)malloc(sizeof(sym_block));
	  rec->value.block->n = 0;
	  rec->value.block->lines = NULL;
	  do{
	    c = parse_get_line(f, &s, &length);
	    if(*s && *s != '%'){
	      char *s1, *tok;
	      int l, col;
	      
	      l = rec->value.block->n;
	      rec->value.block->n++;
	      rec->value.block->lines = (sym_block_line *)
		realloc((void *)rec->value.block->lines, sizeof(sym_block_line)*(l+1));
	      rec->value.block->lines[l].n = 0;
	      rec->value.block->lines[l].fields = NULL;
	      
	      /* parse columns */
	      for(s1 = s; (tok = strtok(s1, "|\t")) != NULL; s1 = NULL){
		char *tok2 = strdup(tok);
		str_trim(tok2);
		
		col = rec->value.block->lines[l].n;
		rec->value.block->lines[l].n++;
		rec->value.block->lines[l].fields = (char **)
		  realloc((void *)rec->value.block->lines[l].fields, sizeof(char *)*(col+1));
		rec->value.block->lines[l].fields[col] = tok2;
	      }
	    }
	  }while(c != EOF && *s != '%');
	}
      }else{ /* we can parse it np */
	parse_result pc;
	parse_exp(s, &pc);
      }
    }
  }while(c != EOF);
  
  free(s);
  if(f != stdin)
    fclose(f);
  
#define OCT_ENV_HEADER "OCT_"

  /*now read options from environment variables (by X) */

  if( getenv("OCT_PARSE_ENV")!=NULL ) {
    
    /* environ is an array of C strings with all the environment
       variables, the format of the string is NAME=VALUE, which
       is directly recognized by parse_exp */
    
    char **env = environ;    
    while(*env) {
      /* Only consider variables that begin with OCT_ */
      if( strncmp(OCT_ENV_HEADER, *env, strlen(OCT_ENV_HEADER)) == 0 ){	
	parse_result pc;
	parse_exp( (*env) + strlen(OCT_ENV_HEADER), &pc);
      }
      
      env++;
    }
  }
  
  sym_clear_reserved();

  return 0;
}

void parse_end()
{
  sym_end_table();
  if(!disable_write) {
    fprintf(fout, "# Octopus parser ended\n");
    if(fout != stdout)
      fclose(fout);
  }
}

int parse_isdef(const char *name)
{
  if(getsym(name) == NULL)
    return 0;
  return 1;
}


static void check_is_numerical(const char * name, const symrec * ptr){
  if( ptr->type != S_CMPLX){
    fprintf(stderr, "Parser error: expecting a numerical value for variable '%s' and found a string.\n", name);
    exit(1);
  }
}


int parse_int(const char *name, int def)
{
  symrec *ptr;
  int ret;

  ptr = getsym(name);	
  if(ptr){
    check_is_numerical(name, ptr);
    ret = ROUND(GSL_REAL(ptr->value.c));
    if(!disable_write) {
      fprintf(fout, "%s = %d\n", name, ret);
    }
    if(fabs(GSL_IMAG(ptr->value.c)) > 1e-10) {
       fprintf(stderr, "Parser error: complex value passed for integer variable '%s'.\n", name);
       exit(1);
    }
    if(fabs(ret - GSL_REAL(ptr->value.c)) > 1e-10) {
       fprintf(stderr, "Parser error: non-integer value passed for integer variable '%s'.\n", name);
       exit(1);
    }
  }else{
    ret = def;
    if(!disable_write) {
      fprintf(fout, "%s = %d\t\t# default\n", name, ret);
    }
  }
  return ret;
}

double parse_double(const char *name, double def)
{
  symrec *ptr;
  double ret;

  ptr = getsym(name);	
  if(ptr){
    check_is_numerical(name, ptr);
    ret = GSL_REAL(ptr->value.c);
    if(!disable_write) {
      fprintf(fout, "%s = %g\n", name, ret);
    }
    if(fabs(GSL_IMAG(ptr->value.c)) > 1e-10) {
       fprintf(stderr, "Parser error: complex value passed for real variable '%s'.\n", name);
       exit(1);
    }
  }else{
    ret = def;
    if(!disable_write) {
      fprintf(fout, "%s = %g\t\t# default\n", name, ret);
    }
  }
  return ret;
}

gsl_complex parse_complex(const char *name, gsl_complex def)
{
  symrec *ptr;
  gsl_complex ret;

  ptr = getsym(name);	
  if(ptr){
    check_is_numerical(name, ptr);
    ret = ptr->value.c;
    if(!disable_write) {
      fprintf(fout, "%s = (%g, %g)\n", name, GSL_REAL(ret), GSL_IMAG(ret));
    }
  }else{
    ret = def;
    if(!disable_write) {
      fprintf(fout, "%s = (%g, %g)\t\t# default\n", name, GSL_REAL(ret), GSL_IMAG(ret));
    }
  }
  return ret;
}

/* FIXME: this function should use strcpy rather than assigning char* */
char *parse_string(const char *name, char *def)
{
  symrec *ptr;
  char *ret;
  
  ptr = getsym(name);
  if(ptr){
    if( ptr->type != S_STR){
      fprintf(stderr, "Parser error: expecting a string for variable '%s'.\n", name);
      exit(1);
    }
    ret = ptr->value.str; /* this should use strcpy instead */
    if(!disable_write) {
      fprintf(fout, "%s = \"%s\"\n", name, ret);
    }
  }else{
    ret = def;
    if(!disable_write) {
      fprintf(fout, "%s = \"%s\"\t\t# default\n", name, ret);
    }
  }
  return ret;
}

int parse_block (const char *name, sym_block **blk)
{
  symrec *ptr;

  ptr = getsym(name);
  if(ptr && ptr->type == S_BLOCK){
    *blk = ptr->value.block;
    (*blk)->name = (char *)malloc(strlen(name) + 1);
    strcpy((*blk)->name, name);
    if(!disable_write) {
      fprintf(fout, "Opened block '%s'\n", name);
    }
    return 0;
  }else{
    *blk = NULL;
    return -1;
  }
}

int parse_block_end (sym_block **blk)
{
  if(!disable_write) {
     fprintf(fout, "Closed block '%s'\n", (*blk)->name);
  }
  free((*blk)->name);
  *blk = NULL;
  return 0;
}

int parse_block_n(const sym_block *blk)
{
  assert(blk != NULL);

  return blk->n;
}

int parse_block_cols(const sym_block *blk, int l)
{
  assert(blk!=NULL);
  if(l < 0 || l >= blk->n){
    fprintf(stderr, "Parser error: row %i out of range [0,%i] when parsing block '%s'.\n", l, blk->n-1, blk->name);
    exit(1);
  }
  
  return blk->lines[l].n;
}

static int parse_block_work(const sym_block *blk, int l, int col, parse_result *r)
{
  assert(blk!=NULL);
  assert(l>=0 && l<blk->n);

  if(col < 0 || col >= blk->lines[l].n){
    fprintf(stderr, "Parser error: column %i out of range [0,%i] when parsing block '%s'.\n", col, blk->lines[l].n-1, blk->name);
    exit(1);
  }
  
  return parse_exp(blk->lines[l].fields[col], r);
}


int parse_block_int(const sym_block *blk, int l, int col, int *r)
{
  int o;
  parse_result pr;

  o = parse_block_work(blk, l, col, &pr);

  if(o == 0) {
    if(pr.type == PR_CMPLX){
      *r = ROUND(GSL_REAL(pr.value.c));
      if(!disable_write) {
	fprintf(fout, "  %s (%d, %d) = %d\n", blk->name, l, col, *r);
      }
      if(fabs(GSL_IMAG(pr.value.c)) > 1e-10) {
        fprintf(stderr, "Parser error: complex value passed for integer field in block '%s'.\n", blk->name);
	exit(1);
      }
      if(fabs(*r - GSL_REAL(pr.value.c)) > 1e-10) {
	fprintf(stderr, "Parser error: non-integer value passed for integer field in block '%s'.\n", blk->name);
        exit(1);
      }
    } else {
      o = -1;
    }
  }

  parse_result_free(&pr);
  return o;
}

int parse_block_double(const sym_block *blk, int l, int col, double *r)
{
  int o;
  parse_result pr;
  
  o = parse_block_work(blk, l, col, &pr);
  
  if(o == 0) {
    if(pr.type == PR_CMPLX){
      *r = GSL_REAL(pr.value.c);
      if(!disable_write) {
        fprintf(fout, "  %s (%d, %d) = %g\n", blk->name, l, col, *r);
      }
      if(fabs(GSL_IMAG(pr.value.c)) > 1e-10) {
        fprintf(stderr, "Parser error: complex value passed for real field in block '%s'.\n", blk->name);
	exit(1);
      }
    } else {
      o = -1;
    }
  }
  
  parse_result_free(&pr);
  return o;
}

int parse_block_complex(const sym_block *blk, int l, int col, gsl_complex *r)
{
  int o;
  parse_result pr;

  o = parse_block_work(blk, l, col, &pr);

  if(o == 0) {
    if(pr.type == PR_CMPLX){
      *r = pr.value.c;
      if(!disable_write) {
        fprintf(fout, "  %s (%d, %d) = (%g,%g)\n", blk->name, l, col, GSL_REAL(*r), GSL_IMAG(*r));
      }
    } else {
      o = -1;
    }
  }

  parse_result_free(&pr);
  return o;
}

int parse_block_string(const sym_block *blk, int l, int col, char **r)
{
  int o;
  parse_result pr;

  o = parse_block_work(blk, l, col, &pr);

  if(o == 0) {
    if(pr.type == PR_STR){
      *r = strdup(pr.value.s);
      if(!disable_write) {
        fprintf(fout, "  %s (%d, %d) = \"%s\"\n", blk->name, l, col, *r);
      }
    } else {
      o = -1;
    }
  }

  parse_result_free(&pr);
  return o;
}


void parse_result_free(parse_result *t)
{
  if(t->type == PR_STR)
    free(t->value.s);
  t->type = PR_NONE;
}


void parse_putsym_int(const char *s, int i)
{
  symrec *rec = putsym(s, S_CMPLX);
  GSL_SET_COMPLEX(&rec->value.c, (double)i, 0);
  rec->def = 1;
  rec->used = 1;
}

void parse_putsym_double(const char *s, double d)
{
  symrec *rec =  putsym(s, S_CMPLX);
  GSL_SET_COMPLEX(&rec->value.c, d, 0);
  rec->def = 1;
  rec->used = 1;
}

void parse_putsym_complex(const char *s, gsl_complex c)
{
  symrec *rec =  putsym(s, S_CMPLX);
  rec->value.c = c;
  rec->def = 1;
  rec->used = 1;
}

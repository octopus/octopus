#ifndef _LIB_OCT_H
#define _LIB_OCT_H

#include <gsl/gsl_complex.h>

// from ylm.c
double ylm(double x, double y, double z, int l, int m);

// from varia.c
void fft_optimize(int *n, int p, int par);

// from parse.c
int parse_init(char *file_in, char *file_out);
void parse_end();

int parse_isdef(char *name);

int parse_int(char *name, int def);
double parse_double(char *name, double def);
gsl_complex parse_complex(char *name, gsl_complex def);
char *parse_string(char *name, char *def);

int parse_block_n(char *name);
int parse_block_int(char *name, int l, int col, int *r);
int parse_block_double(char *name, int l, int col, double *r);
int parse_block_complex(char *name, int l, int col, gsl_complex *r);
int parse_block_string(char *name, int l, int col, char **r);

// from parse_exp.c
typedef struct parse_result{
  union {
	  gsl_complex c;
		char *s;
	} value;
	enum {PR_CMPLX, PR_STR} type;
} parse_result;

int parse_exp(char *exp, parse_result *t);

#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include "symbols.h"

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

	// names are always lowercase
	ptr->name = strdup(sym_name);
	str_tolower(ptr->name);

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
	return NULL;
}

struct init_fntc{
	char *fname;
	gsl_complex (*fnct)();
};

static struct init_fntc arith_fncts[] = {
	{"sqrt",  gsl_complex_sqrt},
	{"exp",   gsl_complex_exp},
	{"ln",    gsl_complex_log},
	{"log",   gsl_complex_log},
	{"log10", gsl_complex_log10},

	{"sin",  gsl_complex_sin},
	{"cos",  gsl_complex_cos},
	{"tan",  gsl_complex_tan},
	{"sec",  gsl_complex_sec},
	{"csc",  gsl_complex_csc},
	{"cot",  gsl_complex_cot},

	{"asin", gsl_complex_arcsin},
	{"acos", gsl_complex_arccos},
	{"atan", gsl_complex_arctan},
	{"asec", gsl_complex_arcsec},
	{"acsc", gsl_complex_arccsc},
	{"acot", gsl_complex_arccot},

	{"sinh", gsl_complex_sinh},
	{"cosh", gsl_complex_cosh},
	{"tanh", gsl_complex_tanh},
	{"sech", gsl_complex_sech},
	{"csch", gsl_complex_csch},
	{"coth", gsl_complex_coth},
	
	{"asinh", gsl_complex_arcsinh},
	{"acosh", gsl_complex_arccosh},
	{"atanh", gsl_complex_arctanh},
	{"asech", gsl_complex_arcsech},
	{"acsch", gsl_complex_arccsch},
	{"acoth", gsl_complex_arccoth},	

	{0, 0}
};

struct init_cnst{
	char *fname;
	double c;
};

static struct init_cnst arith_cnts[] = {
	{"pi",      M_PI},
	{"e",       M_E},
  {".true.",  1},
  {".t.",     1},
	{"yes", 1},
	{".false.", 0},
	{".f.",     0},
	{"no",      0}
};

void sym_init_table ()  /* puts arithmetic functions in table. */
{
	int i;
	symrec *ptr;
	for (i = 0; arith_fncts[i].fname != 0; i++){
		ptr = putsym (arith_fncts[i].fname, S_FNCT);
		ptr->value.fnctptr = arith_fncts[i].fnct;
	}

	/* now the constants */
	for (i = 0; arith_cnts[i].fname != 0; i++){
		ptr = putsym(arith_cnts[i].fname, S_CMPLX);
		GSL_SET_COMPLEX(&ptr->value.c, arith_cnts[i].c, 0);
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
			printf(" = (%lf,%lf)\n", GSL_REAL(ptr->value.c), GSL_IMAG(ptr->value.c));
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

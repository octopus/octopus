#ifndef _SYMBOLS_H
#define _SYMBOLS_H

#include <gsl/gsl_complex.h>

typedef struct sym_block_line{
	int n;
	char **fields;
}sym_block_line;

typedef struct sym_block{
	int n;
	sym_block_line *lines;
}sym_block;

typedef enum{
	S_CMPLX, S_STR, S_BLOCK, S_FNCT
}symrec_type;

/* Data type for links in the chain of symbols. */
typedef struct symrec{
	char *name;                  /* name of symbol */
	symrec_type type;            /* type of symbol: either VAR or FNCT */

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

symrec *putsym (char *sym_name, symrec_type sym_type);
symrec *getsym (char *sym_name);
void sym_init_table();
void sym_end_table();
void sym_output_table();
char *str_tolower(char *in);

#endif

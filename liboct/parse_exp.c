#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "liboct.h"
#include "symbols.h"

static char *par_string;
static int par_pos;
parse_result par_res;

int yylex();
int yyerror (char *s)  /* Called by yyparse on error */
{
  /* Do nothing */
	/* printf("%s\n", s); */
	return 0;
}

/* include the parser */
#include "grammar.c"

int parse_exp(char *exp, parse_result *r)
{
	int o;

	par_string = exp;
	par_pos = 0;

	o = yyparse();
	if(o == 0){
		r->type = par_res.type;
		if(r->type == PR_CMPLX)
			r->value.c = par_res.value.c;
		else
			r->value.s = par_res.value.s;
	}
	return o;
}

int get_real(char *s, double *d)
{
	int n=0;
	sscanf(s, "%lg", d);
	while(*s && (isdigit(*s) || *s=='.' || *s=='e' || *s=='E')){
		if((*s=='e' || *s=='E') && (*(s+1)=='+' || *(s+1)=='-')) {s++; n++;}
		s++; n++;
	}
	return n;
}

int yylex (){
	int c;
	static char *symbuf = 0;
	static int length = 0;
     
	/* Ignore whitespace, get first nonwhite character.  */
	while ((c = par_string[par_pos++]) == ' ' || c == '\t');
	
	if (c == '\0')
		return '\n';
	
	/* Char starts a number => parse the number.         */
	if (c == '.' || isdigit (c)){
		par_pos--;
		par_pos += get_real(&par_string[par_pos], &GSL_REAL(yylval.val));
		return NUM;
	}
     
	/* Char starts an identifier => read the name.       */
	if (isalpha (c) || c == '\'' || c == '\"'){
		symrec *s;
		char startc = c;
		int i;
		
		/* Initially make the buffer long enough
			 for a 40-character symbol name.  */
		if (length == 0)
			length = 40, symbuf = (char *)malloc (length + 1);
		
		if(startc == '\'' || startc == '\"')
			c = par_string[par_pos++];
		else
			startc = 0; /* false */

		i = 0;
		do{
			/* If buffer is full, make it bigger.        */
			if (i == length){
				length *= 2;
				symbuf = (char *)realloc (symbuf, length + 1);
			}
			/* Add this character to the buffer.         */
			symbuf[i++] = c;
			/* Get another character.                    */
			c = par_string[par_pos++];
		}while (c != '\0' && ((startc && c!=startc) || 
				      (!startc && (isalnum(c) || c == '_' ))));
		
		if(!startc) par_pos--;
		symbuf[i] = '\0';
		
		if(!startc){
			s = getsym (symbuf);
			if (s == 0)
				s = putsym (symbuf, S_CMPLX);
			yylval.tptr = s;
			if(s->type == S_CMPLX)
				return VAR;
			else
				return FNCT;
		}else{
			yylval.str = strdup(symbuf);
			return STR;
		}
	}
	
	/* Any other character is a token by itself.        */
	return c;
}

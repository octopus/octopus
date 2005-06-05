# $Id$

noinst_LIBRARIES = liboct_parser.a

liboct_parser_a_SOURCES = symbols.c parse.c parse_exp.c string_f.c parser_f.c gsl_userdef.c
## bison is not always installed, so it is better to avoid the possibility of regenerating grammar.c
## EXTRA_liboct_parser_a_SOURCES = grammar.y

noinst_HEADERS = grammar.c liboct_parser.h symbols.h string_f.h gsl_userdef.h

CLEANFILES = *~ *.bak
DISTCLEANFILES = Makefile.in


# $Id$

lib_LIBRARIES = liboct_parser.a

liboct_parser_a_SOURCES = symbols.c parse.c parse_exp.c parser_f.c gsl_userdef.c
## bison is not always installed, so it is better to avoid the possibility of regenerating grammar.c
## EXTRA_liboct_parser_a_SOURCES = grammar.y

noinst_HEADERS = grammar.c symbols.h gsl_userdef.h
include_HEADERS = liboct_parser.h

CLEANFILES = *~ *.bak
DISTCLEANFILES = Makefile.in

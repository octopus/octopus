## Process this file with automake to produce Makefile.in

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id$

lib_LIBRARIES = liboct_parser.a

liboct_parser_a_SOURCES = symbols.c parse.c parse_exp.c parser_f.c gsl_userdef.c
## bison is not always installed, so it is better to avoid the possibility of regenerating grammar.c
## EXTRA_liboct_parser_a_SOURCES = grammar.y

noinst_HEADERS = grammar.c symbols.h gsl_userdef.h string_f.h
include_HEADERS = liboct_parser.h string_f.h

AM_CPPFLAGS = $(GSL_CFLAGS)

CLEANFILES = *~ *.bak

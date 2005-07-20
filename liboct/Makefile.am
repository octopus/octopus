# $Id$

noinst_LIBRARIES = liboct.a

liboct_a_SOURCES = recipes.c ylm.c varia.c oct_f.c oct_gsl_f.c write_iter.c varinfo.c
AM_CPPFLAGS = -I../liboct_parser

noinst_HEADERS = varia.h

CLEANFILES = *~ *.bak
DISTCLEANFILES = Makefile.in

config.h: ../config.h
	ln -sf ../config.h .

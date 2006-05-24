# $Id$

noinst_LIBRARIES = liboct.a

liboct_a_SOURCES = recipes.c ylm.c varia.c oct_f.c oct_gsl_f.c write_iter.c varinfo.c
AM_CPPFLAGS = -I$(top_srcdir)/liboct_parser -I$(top_srcdir)/external_libs/metis-4.0/

noinst_HEADERS = varia.h

CLEANFILES = *~ *.bak
DISTCLEANFILES = Makefile.in

config.h: ../config.h
	ln -sf ../config.h .

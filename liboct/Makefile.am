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

noinst_LIBRARIES = liboct.a

liboct_a_SOURCES = recipes.c ylm.c varia.c \
	oct_f.c oct_gsl_f.c oct_gdlib_f.c  \
	write_iter.c varinfo.c \
	signals.c

AM_CPPFLAGS = -I$(top_srcdir)/liboct_parser -I$(top_srcdir)/external_libs/metis-4.0 \
	-I$(top_srcdir)/libstring_f
AM_CFLAGS = @GD_CFLAGS@

noinst_HEADERS = varia.h

CLEANFILES = *~ *.bak

config.h: ../config.h
	ln -sf ../config.h .


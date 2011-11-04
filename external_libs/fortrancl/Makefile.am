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
## $Id: Makefile.am 2995 2007-06-13 17:49:22Z xavier $

noinst_LIBRARIES = libfortrancl.a

libfortrancl_a_SOURCES =     \
	cl.F90

if COMPILE_OPENCL
 libfortrancl_a_SOURCES += cl_low.c
endif

noinst_HEADERS =               \
	cl_constants_inc.F90   \
	string_f.h

cl.o : cl.F90 cl_constants_inc.F90 config_F90.h string_f.h

config_F90.h: $(top_builddir)/config.h
	grep '^#define' $^ > config_F90.h

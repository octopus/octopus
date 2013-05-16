## Process this file with automake to produce Makefile.in

## Copyright (C) 2013 J. Alberdi-Rodriguez
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

#if COMPILE_METIS_5
  noinst_LIBRARIES = libmetis5.a
#endif

libmetis5_a_SOURCES = \
	libgk.a libmetis_tmp.a

noinst_HEADERS = \
	include/metis.h

SUBDIRS =        \
	 	GKlib	libmetis

CLEANFILES = *~ *.bak

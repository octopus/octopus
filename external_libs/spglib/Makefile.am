## Process this file with automake to produce Makefile.in

## Copyright (C) 2008 X. Andrade
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
## $Id: Makefile.am 2771 2007-03-20 18:29:43Z lorenzen $


# testing is disabled due to problems with gcc-4.3

#TESTS = spg_test

#check_PROGRAMS = spg_test
#spg_test_SOURCES = test.c
#spg_test_LDADD = -L. -lspg

noinst_LIBRARIES = libspg.a

libspg_a_SOURCES =            \
	bravais.c               \
	cell.c                  \
	mathfunc.c              \
	pointgroup.c            \
	primitive.c             \
	spacegroup.c            \
	spacegroup_database.c   \
	spacegroup_data.c       \
	spglib.c                \
	spglib_f.c              \
	symmetry.c

noinst_HEADERS =                \
	bravais.h               \
	cell.h			\
	ispglib.h               \
	mathfunc.h              \
	pointgroup.h            \
	primitive.h             \
	spacegroup_database.h   \
	spacegroup_data.h       \
	spacegroup.h            \
	spglib.h                \
	symmetry.h

CLEANFILES = *~ *.bak

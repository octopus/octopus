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
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id$

if COMPILE_METIS
  noinst_LIBRARIES = libgk.a
endif

AM_CPPFLAGS = -DNDEBUG -DNDEBUG2 

libgk_a_SOURCES = \
	b64.c       fkvkselect.c  htable.c    omp.c     seq.c       \
	blas.c      fs.c          io.c        pdb.c     sort.c      \
	csr.c       getopt.c      itemsets.c  pqueue.c  string.c    \
	error.c     gkregex.c     mcore.c     random.c  timers.c    \
	evaluate.c  graph.c       memory.c    rw.c      tokenizer.c \
	util.c

noinst_HEADERS = \
	gk_arch.h     GKlib.h        gk_mkpqueue2.h  gk_mkutils.h  gk_types.h    \
	gk_defs.h     gk_macros.h    gk_mkpqueue.h   gk_proto.h    ms_inttypes.h \
	gk_externs.h  gk_mkblas.h    gk_mkrandom.h   gkregex.h     ms_stat.h     \
	gk_getopt.h   gk_mkmemory.h  gk_mksort.h     gk_struct.h   ms_stdint.h


CLEANFILES = *~ *.bak

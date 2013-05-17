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

noinst_LIBRARIES = libmetis_tmp.a

libmetis_tmp_a_SOURCES = \
	auxapi.c      contig.c   graph.c       mesh.c      options.c    srefine.c \
	balance.c     debug.c    initpart.c    meshpart.c  parmetis.c   stat.c    \
	bucketsort.c  fm.c       kmetis.c      minconn.c   pmetis.c     timing.c  \
	checkgraph.c  fortran.c  kwayfm.c      mincover.c  refine.c     util.c    \
	coarsen.c     frename.c  kwayrefine.c  mmd.c       separator.c  wspace.c  \
	compress.c    gklib.c    mcutil.c      ometis.c    sfm.c

noinst_HEADERS = \
	defs.h        gklib_rename.h  metislib.h  rename.h      struct.h \
	gklib_defs.h  macros.h        proto.h     stdheaders.h
#	../GKlib/GKlib.h

AM_CPPFLAGS = -I../GKlib/ -I../include

CLEANFILES = *~ *.bak

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

noinst_LIBRARIES = libbasic.a

libbasic_a_SOURCES = \
             liboct.F90 string.F90 varinfo.F90 mpi.F90 global.F90 liboct_parser.F90 \
	     datasets.F90 messages.F90 io.F90 units.F90 mpi_debug.F90 profiling.F90 \
             utils.F90 multicomm.F90 periodic_table.F90 getopt_f.c command_line.F90

noinst_HEADERS = \
	mpi_debug_inc.F90

include $(srcdir)/../common-rules.make

# how to compile global.F90 (have to define some extra variables here)
DATE := $(shell date)
LATEST_SVN := $(shell ../build/svn_release_number.sh 2>/dev/null) 

global.o: global.F90
	@FCCPP@ @CPPFLAGS@ -I. $(AM_CPPFLAGS) \
	        -DSHARE_OCTOPUS='"${pkgdatadir}"' -DLATEST_SVN='"${LATEST_SVN}"' \
		-DOCTOPUS_VERSION='"${VERSION}"' -DBUILD_TIME='"${DATE}"' -DFC='"${FC}"' \
		-DFCFLAGS='"${FCFLAGS} ${FCFLAGS_f90}"' $< > $*_oct.f90
	@FC@ @FCFLAGS@ $(AM_FCFLAGS) -c @FCFLAGS_f90@ -o $@ $*_oct.f90
	@rm -f $*_oct.f90

# This is to generate the dependencies. It makes use of the makedepf90 program of 
# Erik Edelmann (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src directory), substituting
# the dependencies below by the resulting new dependencies. (For normal compilation
# of octopus, there is absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libbasic_a_SOURCES)

mpi_debug.o: mpi_debug_inc.F90

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
	loct.F90 string.F90 varinfo.F90 mpi.F90 global.F90 \
	loct_parser.F90 datasets.F90 messages.F90 io.F90   \
	io_binary.c mpi_lib.F90 c_pointer.F90 \
	units.F90 mpi_debug.F90 profiling.F90 utils.F90      \
	multicomm.F90 periodic_table.F90 getopt_f.c command_line.F90 \
	distributed.F90 \
	hardware.F90 ialist.F90 iihash.F90 calc_mode.F90 \
	recipes.c      \
	varia.c        \
	oct_f.c        \
	oct_gdlib_f.c  \
	write_iter.F90   \
	write_iter_low.c \
	varinfo_low.c    \
	signals.c        \
	subarray.F90     \
	io_binary_f.F90

if COMPILE_PAPI
libbasic_a_SOURCES += papi.c papi_f.F90
endif

noinst_HEADERS =          \
	mpi_debug_inc.F90 \
	mpi_lib_inc.F90   \
	loct_inc.F90      \
	varia.h           \
	io_binary.h

include $(srcdir)/../common-rules.make

# How to compile global.F90 (have to define some extra variables here).
global.o: global.F90
	rev="`$(SHELL) $(top_srcdir)/build/svn_release_number.sh 2> /dev/null`"; \
	truecc="`$(SHELL) $(top_srcdir)/build/true_compiler.sh $(CC) 2> /dev/null`"; \
	truefc="`$(SHELL) $(top_srcdir)/build/true_compiler.sh $(FC) 2> /dev/null`"; \
	@FCCPP@ @CPPFLAGS@ -I. $(AM_CPPFLAGS) \
	        -DSHARE_OCTOPUS='"$(pkgdatadir)"' \
		-DLATEST_SVN="'$$rev'" \
		-DOCTOPUS_VERSION='"$(VERSION)"' \
		-DBUILD_TIME="'`date`'" \
		-DCC="'$(CC) $$truecc'" \
		-DCFLAGS="'$(CFLAGS)'" \
		-DFC="'$(FC) $$truefc'" \
		-DFCFLAGS='"$(FCFLAGS) $(FCFLAGS_f90)"' $< > $*_oct.f90
	@if [ "@F90_ACCEPTS_LINE_NUMBERS@" = "no" ]; then \
		grep -v "^#" $*_oct.f90 > $*_oct.f91; \
		mv -f $*_oct.f91 $*_oct.f90; \
	fi
	@FC@ @FCFLAGS@ $(AM_FCFLAGS) -c @FCFLAGS_f90@ -o $@ $*_oct.f90
	@rm -f $*_oct.f90

# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(libbasic_a_SOURCES)

mpi.o : mpi.F90 
hardware.o : hardware.F90
ialist.o : ialist.F90 messages.o profiling.o
c_pointer.o : c_pointer.F90
command_line.o : command_line.F90 messages.o profiling.o
periodic_table.o : periodic_table.F90
loct.o : loct.F90 loct_inc.F90 c_pointer.o
string.o : string.F90 c_pointer.o
varinfo.o : varinfo.F90 string.o c_pointer.o
global.o : global.F90 mpi.o loct.o hardware.o
loct_parser.o : loct_parser.F90 mpi.o loct.o global.o c_pointer.o
datasets.o : datasets.F90 loct_parser.o c_pointer.o
messages.o : messages.F90 varinfo.o string.o mpi.o loct_parser.o loct.o global.o datasets.o 
io.o : io.F90 mpi.o messages.o loct_parser.o loct.o global.o datasets.o 
mpi_lib.o : mpi_lib.F90 mpi_lib_inc.F90 mpi_debug.o mpi.o messages.o global.o profiling.o
units.o : units.F90 io.o messages.o loct_parser.o global.o datasets.o  varinfo.o
mpi_debug.o : mpi_debug.F90 mpi_debug_inc.F90 mpi.o messages.o loct.o global.o 
profiling.o : profiling.F90 mpi.o messages.o loct.o io.o global.o papi_f.o
utils.o : utils.F90 messages.o global.o 
multicomm.o : multicomm.F90 varinfo.o utils.o mpi.o messages.o loct_parser.o io.o global.o datasets.o profiling.o
iihash.o : iihash.F90 ialist.o global.o
calc_mode.o : calc_mode.F90 multicomm.o
distributed.o: distributed.F90 multicomm.o profiling.o
varia.o : varia.c varia.h
io_binary.o : io_binary.c io_binary.h
getopt_f.o : getopt_f.c
write_iter.o : write_iter.F90 c_pointer.o
subarray.o : subarray.F90 subarray_inc.F90 profiling.o
io_binary_f.o : io_binary_f.F90 io_binary.h
papi.o : papi.c
papi_f.o : papi_f.F90

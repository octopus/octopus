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

# ---------------------------------------------------------------
# Include paths.
# ---------------------------------------------------------------

AM_FCFLAGS = \
	@F90_MODULE_FLAG@$(top_builddir)/src/basic   	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/math    	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/species 	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/grid    	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/poisson 	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/states  	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/xc      	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/h_sys   	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/scf     	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/td          \
	@F90_MODULE_FLAG@$(top_builddir)/src/transport   \
	@F90_MODULE_FLAG@$(top_builddir)/src/opt_control \
	@F90_MODULE_FLAG@$(top_builddir)/src/sternheimer \
	@F90_MODULE_FLAG@$(top_builddir)/libxc/src       \
	@F90_MODULE_FLAG@$(top_builddir)/external_libs/qshep

AM_CPPFLAGS = \
	-I$(top_srcdir)/src/include   \
	-I$(top_builddir)/src/include \
	-I$(top_srcdir)/libstring_f


# ---------------------------------------------------------------
# Define libraries here.
# ---------------------------------------------------------------

octopus_LIBS = \
	$(top_builddir)/src/sternheimer/libsternheimer.a \
	$(top_builddir)/src/opt_control/libopt_control.a \
	$(top_builddir)/src/td/libtd.a                   \
	$(top_builddir)/src/transport/libtransport.a     \
	$(top_builddir)/src/scf/libscf.a                 \
	$(top_builddir)/src/h_sys/libh_sys.a             \
	$(top_builddir)/src/xc/libxc.a                   \
	$(top_builddir)/src/states/libstates.a           \
	$(top_builddir)/src/poisson/libpoisson.a         \
	$(top_builddir)/src/grid/libgrid.a               \
	$(top_builddir)/src/species/libspecies.a         \
	$(top_builddir)/src/math/libmath.a               \
	$(top_builddir)/src/basic/libbasic.a

core_LIBS = \
	$(octopus_LIBS)                               \
	@LIBS_LAPACK@ @LIBS_BLAS@                     \
	$(top_builddir)/liboct/liboct.a               \
	$(top_builddir)/liboct_parser/liboct_parser.a \
	$(top_builddir)/libxc/src/libxc.a             \
	-L$(top_builddir)/libstring_f -lstring_f      \
	@GSL_LIBS@ @GD_LIBS@ @FCEXTRALIBS@

external_LIBS = \
	$(top_builddir)/external_libs/expokit/libexpokit.a \
	$(top_builddir)/external_libs/qshep/libqshep.a     \
	$(top_builddir)/external_libs/poisson_isf/libpoisson_isf.a

if COMPILE_METIS
  external_LIBS += $(top_builddir)/external_libs/metis-4.0/libmetis.a
endif

if COMPILE_LIBNBC
  external_LIBS += $(top_builddir)/external_libs/libnbc/libnbc.a
endif

all_LIBS = $(core_LIBS) @LIBS_FFT@ @LIBS_TRLAN@ @LIBS_ARPACK@ @LIBS_SPARSKIT@ \
  @LIBS_NETCDF@ $(external_LIBS) @LIBS_MPI@ @LIBS_LAPACK@ @LIBS_BLAS@


# ---------------------------------------------------------------
# How to compile F90 files.
# ---------------------------------------------------------------

SUFFIXES = _oct.f90 .F90 .o

# Compilation is a two step process: first we preprocess F90 files
# to generate _oct.f90 files. Then, we compiler this _oct.f90 into
# an object file and delete the intermediate file.
.F90.o:
	@FCCPP@ @CPPFLAGS@ $(AM_CPPFLAGS) -I. $< > $*_oct.f90
	@if [ "@DEBUG@" = "no" ]; then \
		cat $*_oct.f90 | grep -v pop_sub | \
			grep -v push_sub >$*_oct.f91; \
		mv -f $*_oct.f91 $*_oct.f90; \
	fi
	@perl -pi -e 's/\\newline/\n/g; s/\\cardinal/#/g' $*_oct.f90
	@FC@ @FCFLAGS@ @FCFLAGS_NETCDF@ $(AM_FCFLAGS) -c @FCFLAGS_f90@ -o $@ $*_oct.f90
	@rm -f $*_oct.f90

# This rule is basically to create a _oct.f90 file by hand for
# debugging purposes. It is identical to the first part of
# the .F90.o rule.
.F90_oct.f90:
	@FCCPP@ @CPPFLAGS@ $(AM_CPPFLAGS) -I. $< > $*_oct.f90
	@if [ "@DEBUG@" = "no" ]; then \
		cat $*_oct.f90 | grep -v pop_sub | \
			grep -v push_sub >$*_oct.f91; \
		mv -f $*_oct.f91 $*_oct.f90; \
	fi
	@perl -pi -e 's/\\newline/\n/g; s/\\cardinal/#/g' $*_oct.f90


# ---------------------------------------------------------------
# Miscellaneous.
# ---------------------------------------------------------------

# ctags.
CTAGS = ctags-exuberant -e

# Cleaning.
CLEANFILES = *~ *.bak *.mod *.il *.d *.pc* ifc* *_oct.f90 config_F90.h


# Local Variables:
# mode: Makefile
# coding: utf-8
# End:

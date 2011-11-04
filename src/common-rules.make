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
	@F90_MODULE_FLAG@$(top_builddir)/src/ions 	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/grid    	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/poisson 	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/states  	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/xc      	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/system   	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/hamiltonian \
	@F90_MODULE_FLAG@$(top_builddir)/src/scf     	 \
	@F90_MODULE_FLAG@$(top_builddir)/src/td          \
	@F90_MODULE_FLAG@$(top_builddir)/src/opt_control \
	@F90_MODULE_FLAG@$(top_builddir)/src/sternheimer         \
	@F90_MODULE_FLAG@$(top_builddir)/external_libs/qshep     \
	@F90_MODULE_FLAG@$(top_builddir)/external_libs/fortrancl \
	@F90_MODULE_FLAG@$(top_builddir)/external_libs/bpdn

AM_CPPFLAGS = \
	-I$(top_srcdir)/src/include   \
	-I$(top_builddir)/src/include \
        -I$(top_srcdir)/external_libs/spglib/src \
	-I$(top_srcdir)/liboct_parser \
        $(CPPFLAGS_ZOLTAN) \
	-DSHARE_OCTOPUS='"$(pkgdatadir)"'

AM_CCASFLAGS = \
	-I$(top_builddir)/

# ---------------------------------------------------------------
# Define libraries here.
# ---------------------------------------------------------------

octopus_LIBS = \
	$(top_builddir)/src/sternheimer/libsternheimer.a \
	$(top_builddir)/src/opt_control/libopt_control.a \
	$(top_builddir)/src/td/libtd.a                   \
	$(top_builddir)/src/scf/libscf.a                 \
	$(top_builddir)/src/system/libsystem.a           \
	$(top_builddir)/src/hamiltonian/libhamiltonian.a \
	$(top_builddir)/src/xc/libxc.a                   \
	$(top_builddir)/src/states/libstates.a           \
	$(top_builddir)/src/poisson/libpoisson.a         \
	$(top_builddir)/src/grid/libgrid.a               \
	$(top_builddir)/src/ions/libions.a               \
	$(top_builddir)/src/species/libspecies.a         \
	$(top_builddir)/src/math/libmath.a               \
	$(top_builddir)/src/basic/libbasic.a

core_LIBS = \
	$(octopus_LIBS)                               \
	@LIBS_SCALAPACK@ @LIBS_BLACS@                 \
	@LIBS_LAPACK@ @LIBS_BLAS@                     \
	$(top_builddir)/liboct_parser/liboct_parser.a \
	@GSL_LIBS@ @GD_LIBS@ @LIBS_LIBXC@ @FCEXTRALIBS@

external_LIBS = \
	$(top_builddir)/external_libs/qshep/libqshep.a            \
	$(top_builddir)/external_libs/spglib/src/libspglib.a      \
	$(top_builddir)/external_libs/slatec/libslatec.a          \
	$(top_builddir)/external_libs/fortrancl/libfortrancl.a    \
	$(top_builddir)/external_libs/bpdn/libbpdn.a

if COMPILE_METIS
  external_LIBS += $(top_builddir)/external_libs/metis-4.0/libmetis.a
  AM_CPPFLAGS += -I$(top_srcdir)/external_libs/metis-4.0/
endif

if COMPILE_ZOLTAN
  external_LIBS += $(top_builddir)/external_libs/zoltan/libzoltan.a
  AM_CPPFLAGS += -I$(top_srcdir)/external_libs/zoltan/include
endif

if COMPILE_NEWUOA
  external_LIBS += $(top_builddir)/external_libs/newuoa/libnewuoa.a
  AM_FCFLAGS += @F90_MODULE_FLAG@$(top_builddir)/external_libs/newuoa
endif

# Since ETSF_IO depends on netCDF, it must be first in the list
all_LIBS = $(core_LIBS) @LIBS_PFFT@ @LIBS_FFT@ @LIBS_SPARSKIT@ \
  @LIBS_ETSF_IO@ @LIBS_NETCDF@ $(external_LIBS) \
  @LIBS_LIBFM@ @LIBS_MPI@ @LIBS_ZOLTAN@


# ---------------------------------------------------------------
# How to compile F90 files.
# ---------------------------------------------------------------

SUFFIXES = _oct.f90 .F90 .o .S .s

.S.o:
	@CPP@ @CPPFLAGS@ $(INCLUDES) $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $< > $*_oct.s
	@CC@ -c -o $@ $*_oct.s
	@rm -f $*_oct.s

# Compilation is a two-step process: first we preprocess F90 files
# to generate _oct.f90 files. Then, we compile this _oct.f90 into
# an object file and delete the intermediate file.
.F90.o:
	@FCCPP@ @CPPFLAGS@ $(AM_CPPFLAGS) -I. $< > $*_oct.f90
	$(top_srcdir)/build/preprocess.pl $*_oct.f90 \
	  "@DEBUG@" "@F90_ACCEPTS_LINE_NUMBERS@" "@F90_FORALL@"
	@FC@ @FCFLAGS@ @FCFLAGS_NETCDF@ @FCFLAGS_ETSF_IO@ @FCFLAGS_LIBXC@ $(AM_FCFLAGS) -c @FCFLAGS_f90@ -o $@ $*_oct.f90
	@rm -f $*_oct.f90

# This rule is basically to create a _oct.f90 file by hand for
# debugging purposes. It is identical to the first part of
# the .F90.o rule.
.F90_oct.f90:
	@FCCPP@ @CPPFLAGS@ $(AM_CPPFLAGS) -I. $< > $*_oct.f90
	$(top_srcdir)/build/preprocess.pl $*_oct.f90 \
	  "@DEBUG@" "@F90_ACCEPTS_LINE_NUMBERS@" "@F90_FORALL@"


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
